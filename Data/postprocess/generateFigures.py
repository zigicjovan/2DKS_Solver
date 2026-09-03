#!/usr/bin/env python3
"""Headless, streaming conversion of generateFigures.m. See README.md."""
from __future__ import annotations
import argparse
from contextlib import contextmanager
import filecmp
import json
import logging
from pathlib import Path
import re
import shutil
import tempfile
import time

import numpy as np
from scipy.fft import ifft2
from eigenfunctionValidation import CosineProjector
from expFitDelta import expfit_delta

LOG = logging.getLogger('ksPostprocess')
NUMBER = r'[+-]?(?:\d*\.\d+|\d+\.?\d*)(?:[eE][+-]?\d+)?'


def parse_parameters(name):
    tokens = re.findall(r'(?:^|_)(N1|N2|dt|K|ell1|ell2|T|opt|tol|cont|optT)_('
                        + NUMBER + r')(?=_|$)', name)
    params = {key: float(value) for key, value in tokens}
    missing = {'N1', 'N2', 'ell1', 'ell2'} - params.keys()
    if missing:
        raise ValueError(f'Testcase name is missing parameters: {sorted(missing)}')
    if any(not np.isfinite(v) for v in params.values()):
        raise ValueError('Non-finite testcase parameter.')
    for key in ('N1', 'N2'):
        if params[key] < 2 or params[key] != int(params[key]):
            raise ValueError(f'{key} must be an integer of at least 2.')
        params[key] = int(params[key])
    if min(params['ell1'], params['ell2']) <= 0:
        raise ValueError('Length scales must be positive.')
    return params


def exactly_one(case, directory, pattern):
    paths = sorted((case / directory).glob(pattern))
    if len(paths) != 1:
        raise ValueError(f'Expected one {directory}/{pattern}; found {len(paths)} in {case}. '
                         'Separate runs/ranks into distinct testcase directories.')
    return paths[0]


def sorted_forward_files(case):
    files = list((case / 'ForwardSolution').glob('fwd_*.dat'))
    if not files:
        raise ValueError(f'No ForwardSolution/fwd_*.dat files in {case}')
    def key(path):
        match = re.search(r'_(' + NUMBER + r')\.dat$', path.name)
        if not match:
            raise ValueError(f'Cannot parse final numeric suffix: {path.name}')
        return float(match[1])
    files.sort(key=key)
    if len({key(p) for p in files}) != len(files):
        raise ValueError('Duplicate forward-file time suffixes: ambiguous ordering.')
    return files


def iter_fourier_states(paths, shape, byte_order='<', zero_policy='skip'):
    """Read at most one complex128 state per file read, never an entire file.

    Interleaved real/imaginary float64 is exactly complex128 storage. MATLAB's
    reshape is column-major. All-zero slots are skipped by default like MATLAB;
    --zero-states keep is required when an exact zero is a real saved solution.
    """
    size = int(np.prod(shape))
    frame_bytes = 16 * size
    dtype = np.dtype(byte_order + 'c16')
    for path in paths:
        length = path.stat().st_size
        if length % frame_bytes:
            raise ValueError(f'{path}: {length} bytes is not a multiple of {frame_bytes}. '
                             'Expected full complex128 Fourier grids; file may be partial.')
        skipped = 0
        with path.open('rb') as stream:
            for index in range(length // frame_bytes):
                flat = np.fromfile(stream, dtype=dtype, count=size)
                if flat.size != size:
                    raise ValueError(f'{path} changed/truncated while reading.')
                if not np.all(np.isfinite(flat)):
                    raise ValueError(f'Non-finite Fourier data in {path}, slot {index}.')
                if not np.any(flat):
                    if zero_policy == 'error':
                        raise ValueError(f'All-zero state in {path}, slot {index}.')
                    if zero_policy == 'skip':
                        skipped += 1
                        continue
                yield flat.reshape(shape, order='F')
        if skipped:
            LOG.warning('%s: skipped %d all-zero slots (MATLAB behavior).', path.name, skipped)


def first_field(path, shape, args):
    iterator = iter_fourier_states([path], shape, args.byte_order, args.zero_states)
    try:
        hat = next(iterator)
        return ifft2(hat, workers=args.fft_workers).real
    except StopIteration:
        raise ValueError(f'No usable state in {path}') from None
    finally:
        iterator.close()


def center_shift(u):
    """Initial-field minimum clustering and circular-mean centering.

    Fix the original shortened logical mask being applied to the full index
    vector: the mask here applies to the remaining candidate tail.
    """
    flat = u.ravel(order='F')
    # Only the 20 smallest are needed; stable tie order matches linear indexing.
    count = min(20, flat.size)
    threshold = np.partition(flat, count - 1)[count - 1]
    candidates = np.flatnonzero(flat <= threshold)
    indices = candidates[np.argsort(flat[candidates], kind='stable')[:count]]
    cursor = 0
    while cursor < len(indices) - 1:
        tail = indices[cursor:]
        indices = np.r_[indices[:cursor + 1], tail[np.abs(indices[cursor] - tail) > .05*flat.size]]
        cursor += 1
    rows, cols = np.unravel_index(indices, u.shape, order='F')
    shifts = []
    for coords, n in zip((rows, cols), u.shape):
        angles = 2 * np.pi * coords / n
        circular = np.mod(np.arctan2(np.sin(angles).mean(), np.cos(angles).mean())
                          * n / (2*np.pi), n) + 1
        matlab_index = int(np.floor(circular + .5))
        center = (n + 1) // 2 + 1
        shifts.append(center - matlab_index)
    return tuple(shifts)


def load_diagnostics(case, args):
    energy_path = exactly_one(case, 'EnergyEvolution', 'energy*.dat')
    energy = np.loadtxt(energy_path, ndmin=2)
    if energy.shape[1] != 4 or energy.shape[0] < 1 or not np.isfinite(energy).all():
        raise ValueError('Energy file must have four finite columns: t, L2, H1, H2.')
    if np.any(np.diff(energy[:, 0]) < 0):
        raise ValueError('Energy times are not sorted.')
    spectrum_path = exactly_one(case, 'FourierSpectrumEvolution', 'spectrum*.dat')
    raw = np.loadtxt(spectrum_path, ndmin=2)
    states = raw.shape[1] - 1 - args.extra_initial_steps
    if states < 1 or raw.shape[0] < 3 or not np.isfinite(raw).all():
        raise ValueError('Unexpected spectrum shape or non-finite data.')
    if args.states is not None and args.states != states:
        raise ValueError(f'--states {args.states} disagrees with spectrum: {states} states '
                         f'after removing {args.extra_initial_steps} extra initial spectra.')
    # MATLAB drops the zero-radius row, then four extra initial-time columns.
    nonzero = raw[1:]
    selected = np.r_[1, np.arange(2 + args.extra_initial_steps, raw.shape[1])]
    radii, spectra = nonzero[:, 0], nonzero[:, selected]
    if np.any(radii <= 0) or np.any(np.diff(radii) <= 0) or np.any(spectra < 0):
        raise ValueError('Spectral radii must be increasing/positive and energies nonnegative.')
    if args.state_times:
        times = np.loadtxt(args.state_times, ndmin=1).ravel()
        if times.size != states or not np.isfinite(times).all() or np.any(np.diff(times) < 0):
            raise ValueError('State-time file must contain one finite, ordered time per state.')
        if times[0] < energy[0, 0] or times[-1] > energy[-1, 0]:
            raise ValueError('State times lie outside the energy time range.')
        sampled = np.column_stack((times, *[np.interp(times, energy[:, 0], energy[:, j])
                                         for j in range(1, 4)]))
    else:
        # MATLAB round, not NumPy bankers rounding.
        indices = np.floor(np.linspace(0, energy.shape[0] - 1, states) + .5).astype(int)
        sampled = energy[indices]
    fits = np.full((states, 2), np.nan)
    for i in range(states):
        try:
            # The original additionally excludes the first remaining radial bin.
            fits[i] = expfit_delta(radii[1:], spectra[1:, i], fits[i-1, 1] if i else 1e15)
        except ValueError as exc:
            LOG.warning('State %d: strip fit unavailable: %s', i + 1, exc)
    with np.errstate(over='ignore', invalid='ignore'):
        curves = fits[:, 0][None, :] * np.exp(-2*radii[:, None]*fits[:, 1][None, :])
    return energy, sampled, nonzero, radii, spectra, fits, curves


def radial_groups(modes, length):
    radii = np.hypot(modes[:, 0], modes[:, 1]) / length
    order = np.argsort(radii, kind='stable')
    group = np.empty(len(modes), dtype=int)
    unique, representative = [], []
    for i in order:
        if not unique or abs(radii[i] - unique[-1]) > 1e-10:
            unique.append(float(radii[i]))
            a, b = sorted(np.abs(modes[i]), reverse=True)
            representative.append((0, int(a)) if b == 0 else (int(a), int(b)))
        group[i] = len(unique) - 1
    return np.asarray(unique), group, representative


@contextmanager
def display_cache(states, shape, directory):
    with tempfile.TemporaryDirectory(prefix='ks-display-', dir=directory) as temp:
        cache = np.memmap(Path(temp) / 'fields.dat', mode='w+', dtype='float64',
                          shape=(states, *shape))
        try:
            yield cache
        finally:
            cache.flush()
            cache._mmap.close()


def relocate_legacy_outputs(case, data_output, output):
    """Move only known outputs from the old layout after a successful run.

    Differing existing files are preserved with a .previous-N suffix. Unknown
    files are left alone, and only empty legacy directories are removed.
    """
    stems = {
        'initialData': ('InitialData', ('contourIC', 'surfaceIC')),
        'terminalData': ('TerminalData', ('contourTC', 'surfaceTC')),
        'fourierSpectrumEvolution': ('FourierSpectrumEvolution', ('spectrumIC', 'spectrumTC')),
        'energyEvolution': ('EnergyEvolution', ('diagnostics',)),
        'forwardSolution': ('ForwardSolution', ('movie',)),
    }
    pairs = []
    legacy_data_outputs = (case / 'PythonPostprocessing', case / 'python_postprocessing')
    old_group = output / case.name
    for directory, (legacy_directory, prefixes) in stems.items():
        for prefix in prefixes:
            extensions = ('gif', 'mp4') if prefix == 'movie' else ('pdf', 'png', 'svg')
            for extension in extensions:
                name = f'{prefix}{case.name}.{extension}'
                destination = output / directory / case.name / name
                pairs.append((data_output / directory / name, destination))
                pairs.append((data_output / legacy_directory / name, destination))
                for legacy_data_output in legacy_data_outputs:
                    pairs.append((legacy_data_output / legacy_directory / name, destination))
                pairs.append((old_group / directory / name, destination))
                pairs.append((old_group / legacy_directory / name, destination))
                pairs.append((output / legacy_directory / case.name / name, destination))
        if directory == 'forwardSolution':
            for extension in ('gif', 'mp4'):
                name = f'movie{case.name}.partial.{extension}'
                destination = output / directory / case.name / name
                pairs.append((data_output / directory / name, destination))
                pairs.append((data_output / legacy_directory / name, destination))
                for legacy_data_output in legacy_data_outputs:
                    pairs.append((legacy_data_output / legacy_directory / name, destination))
                pairs.append((old_group / directory / name, destination))
                pairs.append((old_group / legacy_directory / name, destination))
                pairs.append((output / legacy_directory / case.name / name, destination))
    csv_layouts = {
        'radial_weights': ('radialWeights', f'radialWeights{case.name}.csv'),
        'strip_width': ('stripWidth', f'stripWidth{case.name}.csv'),
    }
    for old_stem, (directory, filename) in csv_layouts.items():
        destination = output / directory / filename
        old_filename = f'{old_stem}{case.name}.csv'
        former_directory = directory[0].upper() + directory[1:]
        former_filename = filename[0].upper() + filename[1:]
        for legacy_root in (data_output, *legacy_data_outputs):
            pairs.append((legacy_root / f'{old_stem}.csv', destination))
        pairs.extend((
            (old_group / old_stem / old_filename, destination),
            (old_group / directory / filename, destination),
            (old_group / former_directory / former_filename, destination),
            (output / old_stem / case.name / old_filename, destination),
            (output / directory / case.name / filename, destination),
            (output / old_stem / old_filename, destination),
            (output / former_directory / case.name / former_filename, destination),
            (output / former_directory / former_filename, destination),
        ))
    moved = 0
    for source, destination in pairs:
        if not source.is_file():
            continue
        destination.parent.mkdir(parents=True, exist_ok=True)
        if destination.exists():
            if filecmp.cmp(source, destination, shallow=False):
                source.unlink()  # Exact duplicate of the newly saved output.
                moved += 1
                continue
            index = 1
            previous = destination.with_name(f'{destination.stem}.previous-{index}{destination.suffix}')
            while previous.exists():
                index += 1
                previous = destination.with_name(f'{destination.stem}.previous-{index}{destination.suffix}')
            destination = previous
        shutil.move(str(source), str(destination))
        moved += 1
    # The freshly generated CamelCase copies supersede these known files.
    for legacy_data_output in legacy_data_outputs:
        for filename in ('diagnostics.npz', 'run_metadata.json'):
            old_file = legacy_data_output / filename
            if old_file.is_file() and (data_output / filename).is_file():
                old_file.unlink()
    for legacy_root in (data_output, *legacy_data_outputs):
        for directory, (legacy_directory, _) in stems.items():
            for candidate in (directory, legacy_directory):
                old = legacy_root / candidate
                if old.is_dir() and not any(old.iterdir()):
                    old.rmdir()
    legacy_output_dirs = tuple(name for name, _ in stems.values())
    for directory in (*stems, *legacy_output_dirs, 'radial_weights', 'strip_width',
                      'RadialWeights', 'StripWidth', 'radialWeights', 'stripWidth'):
        old = old_group / directory
        if old.is_dir() and not any(old.iterdir()):
            old.rmdir()
    for directory in ('radial_weights', 'strip_width'):
        old = output / directory / case.name
        if old.is_dir() and not any(old.iterdir()):
            old.rmdir()
        parent = output / directory
        if parent.is_dir() and not any(parent.iterdir()):
            parent.rmdir()
    for directory in (*legacy_output_dirs, 'RadialWeights', 'StripWidth'):
        old = output / directory / case.name
        if old.is_dir() and not any(old.iterdir()):
            old.rmdir()
        parent = output / directory
        if parent.is_dir() and not any(parent.iterdir()):
            parent.rmdir()
    for legacy_data_output in legacy_data_outputs:
        if legacy_data_output.is_dir() and not any(legacy_data_output.iterdir()):
            legacy_data_output.rmdir()
    if old_group.is_dir() and not any(old_group.iterdir()):
        old_group.rmdir()
    if moved:
        LOG.info('Relocated %d outputs from the previous folder layout.', moved)


def generate_figures(testcase, args):
    case = Path(testcase).resolve()
    if not case.is_dir():
        raise ValueError(f'Testcase directory does not exist: {case}')
    started = time.perf_counter()
    params = parse_parameters(case.name)
    # The supplied eigenfunction routine and meshgrid/reshape pairing assume a
    # square isotropic domain. Reject ambiguity rather than silently swap axes.
    if params['N1'] != params['N2'] or not np.isclose(params['ell1'], params['ell2'], rtol=0, atol=1e-12):
        raise ValueError('This compatibility conversion requires N1=N2 and ell1=ell2, '
                         'as eigenfunction_validation.m does. Rectangular cases need the '
                         'C++ writer axis convention and an anisotropic mode definition.')
    shape = params['N1'], params['N2']
    data_output = case / 'pythonPostprocessing'
    output_root = Path(args.output).resolve() if args.output else Path.cwd()
    output = output_root
    if output == data_output or data_output in output.parents:
        raise ValueError('Figure/CSV output must be outside the input pythonPostprocessing directory.')
    data_output.mkdir(parents=True, exist_ok=True)
    output.mkdir(parents=True, exist_ok=True)
    energy, sampled, initial_spectra, radii, spectra, fits, fit_curves = load_diagnostics(case, args)
    states = spectra.shape[1]
    LOG.info('%s: %d states, %d x %d grid', case.name, states, *shape)
    initial = first_field(exactly_one(case, 'InitialData', 'fwdIC*.dat'), shape, args)
    terminal = first_field(exactly_one(case, 'TerminalData', 'fwdTC*.dat'), shape, args)
    shift = center_shift(initial)
    initial, terminal = [np.roll(u, shift, axis=(0, 1)) for u in (initial, terminal)]
    projector = CosineProjector(params['ell1'], params['N1'], workers=args.fft_workers)
    group_radii, mode_groups, labels = radial_groups(projector.modes, params['ell1'])
    mode_weights = np.zeros((len(projector.modes), states))
    radial_weights = np.zeros((len(group_radii), states))
    amplitudes = np.zeros_like(mode_weights)
    frame_limits = np.empty((states, 2))
    # No analytical computation is downsampled; only cached movie rendering.
    stride = max(1, int(np.ceil(shape[0] / args.max_display_grid))) if args.max_display_grid else 1
    display_shape = initial[::stride, ::stride].shape
    need_movie = args.plots == 'all' and args.movie != 'none'
    cache_states = states if need_movie else 1
    with display_cache(cache_states, display_shape, args.cache_dir) as cache:
        files = sorted_forward_files(case)
        seen = 0
        limits = [float(min(initial.min(), terminal.min())), float(max(initial.max(), terminal.max()))]
        for i, hat in enumerate(iter_fourier_states(files, shape, args.byte_order, args.zero_states)):
            if i >= states:
                raise ValueError(f'More than {states} usable forward states; spectrum/state alignment is ambiguous.')
            u = np.roll(ifft2(hat, workers=args.fft_workers).real, shift, axis=(0, 1))
            amp = projector.coefficients(u)
            amplitudes[:, i] = amp
            weights = amp**2
            if weights.sum() > 0:
                weights /= weights.sum()
            mode_weights[:, i] = weights
            radial_weights[:, i] = np.bincount(mode_groups, weights=weights, minlength=len(group_radii))
            if need_movie:
                cache[i] = u[::stride, ::stride]
            frame_limits[i] = float(u.min()), float(u.max())
            limits[0], limits[1] = min(limits[0], float(u.min())), max(limits[1], float(u.max()))
            if i == 0 and not np.allclose(u, initial, rtol=1e-8, atol=1e-10):
                LOG.warning('First forward state differs from fwdIC; check saved-state alignment.')
            if i == states - 1 and not np.allclose(u, terminal, rtol=1e-8, atol=1e-10):
                LOG.warning('Last forward state differs from fwdTC; check saved-state alignment.')
            seen += 1
            if seen == 1 or seen % 10 == 0 or seen == states:
                LOG.info('Analyzed %d/%d states', seen, states)
        if seen != states:
            raise ValueError(f'Found {seen} usable forward states, spectrum expects {states}. '
                             'Check --zero-states and --extra-initial-steps.')
        np.savez_compressed(data_output / 'diagnostics.npz', parameters=json.dumps(params),
            energy=energy, sampled_energy=sampled, radii=radii, spectrum=spectra,
            strip_fits=fits, spectrum_fits=fit_curves, modes=projector.modes,
            cosine_amplitudes=amplitudes, mode_weights=mode_weights,
            radial_radii=group_radii, radial_weights=radial_weights,
            representative_modes=np.asarray(labels, dtype=int).reshape(-1, 2),
            initial_field=initial, terminal_field=terminal, center_shift=shift,
            frame_color_limits=frame_limits)
        for directory in ('stripWidth', 'radialWeights'):
            (output / directory).mkdir(parents=True, exist_ok=True)
        np.savetxt(output / 'stripWidth' / f'stripWidth{case.name}.csv', np.column_stack((sampled[:,0], fits)),
                   delimiter=',', header='time,C,delta', comments='')
        np.savetxt(output / 'radialWeights' / f'radialWeights{case.name}.csv', np.column_stack((sampled[:,0], radial_weights.T)),
                   delimiter=',', header='time' + ''.join(f',r={r:.15g}' for r in group_radii), comments='')
        if args.plots != 'none':
            from plotting import render_static, render_movie
            render_static(output, case.name, params, initial, terminal, radii,
                          initial_spectra, spectra, energy, sampled, fits, radial_weights, labels, args)
            if need_movie:
                render_movie(output, case.name, params, cache, stride, energy, sampled,
                             radii, spectra, fits, fit_curves, radial_weights, labels, limits, frame_limits, args)
        metadata = dict(testcase=str(case), parameters=params, states=states,
                        data_output=str(data_output), figure_csv_output=str(output),
                        options=vars(args), center_shift=shift, display_stride=stride,
                        elapsed_seconds=time.perf_counter() - started,
                        forward_files=[str(p) for p in files],
                        conventions={'binary':'complex128, MATLAB column-major reshape',
                         'ifft_normalization':'1/(N1*N2)',
                         'modal_weights':'normalized squared minimum-aligned cosine correlations',
                         'times':'explicit times with interpolated energy' if args.state_times
                                 else 'MATLAB rounded uniform energy-row sampling'})
        (data_output / 'run_metadata.json').write_text(json.dumps(metadata, indent=2) + '\n')
        relocate_legacy_outputs(case, data_output, output)
    LOG.info('Finished in %.2f s: %s', time.perf_counter() - started, output)
    return output


def make_parser():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('testcases', nargs='+', help='Testcase directories; shell wildcards are supported.')
    parser.add_argument('--states', type=int, help='Expected saved-state count; default inferred from spectrum.')
    parser.add_argument('--extra-initial-steps', type=int, default=4)
    parser.add_argument('--plots', choices=('all', 'figures', 'none'), default='all')
    parser.add_argument('--movie', choices=('gif', 'mp4', 'both', 'none'), default='gif')
    parser.add_argument('--output', help='Parent directory for figure/CSV testcase folders; default is pwd. NPZ/JSON stay beside input data.')
    parser.add_argument('--state-times', help='Text file with actual saved-state times (one testcase only).')
    parser.add_argument('--fft-workers', type=int, default=1)
    parser.add_argument('--max-display-grid', type=int, default=0, help='Optional movie grid cap; default 0 retains every grid point.')
    parser.add_argument('--surface-grid', type=int, default=0, help='Optional 3D surface grid cap; default 0 retains every grid point.')
    parser.add_argument('--fps', type=int, default=20)
    parser.add_argument('--dpi', type=int, default=100)
    parser.add_argument('--formats', nargs='+', choices=('pdf', 'png', 'svg'), default=['pdf', 'png'])
    parser.add_argument('--byte-order', choices=('<', '>', '='), default='<')
    parser.add_argument('--zero-states', choices=('skip', 'keep', 'error'), default='skip')
    parser.add_argument('--cache-dir', help='Directory for temporary disk-backed display cache.')
    parser.add_argument('--movie-field', choices=('surface', 'heatmap'), default='surface')
    parser.add_argument('--color-scale', choices=('frame', 'global'), default='frame',
                        help='Field color limits per frame (default, matching static views) or across the run.')
    return parser


def main():
    parser = make_parser()
    args = parser.parse_args()
    if len(args.testcases) > 1 and args.state_times:
        parser.error('--state-times requires exactly one testcase.')
    seen_cases = {}
    for testcase in args.testcases:
        case = Path(testcase).resolve()
        if case.name in seen_cases and seen_cases[case.name] != case:
            parser.error(f'Different input directories share testcase name {case.name}; their output folders would collide.')
        seen_cases[case.name] = case
    if (args.fft_workers < 1 or args.fps < 1 or args.fps > 100 or args.dpi < 30
        or args.max_display_grid < 0 or args.surface_grid < 0 or args.extra_initial_steps < 0
        or (args.states is not None and args.states < 1)):
        parser.error('Invalid numeric option; require workers>=1, 1<=fps<=100, dpi>=30, nonnegative grid caps.')
    logging.basicConfig(level=logging.INFO, format='%(levelname)s: %(message)s')
    failed = 0
    for case in args.testcases:
        try:
            generate_figures(case, args)
        except Exception:
            LOG.exception('Failed testcase: %s', case)
            failed += 1
    return int(failed > 0)


if __name__ == '__main__':
    raise SystemExit(main())
