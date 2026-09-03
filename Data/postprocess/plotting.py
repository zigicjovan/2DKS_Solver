"""Headless figures and bounded-memory GIF/MP4 output."""
from contextlib import ExitStack
from functools import lru_cache
import logging
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.animation import FFMpegWriter, writers
from matplotlib.colors import Normalize
from matplotlib.collections import Collection
from matplotlib.artist import allow_rasterization
from mpl_toolkits.mplot3d import proj3d
from PIL import Image, GifImagePlugin
from redBlue import redblue_cmap

LOG = logging.getLogger('ksPostprocess')
CMAP = redblue_cmap()
plt.rcParams.update({'font.size': 11, 'axes.titlesize': 12, 'axes.spines.top': False,
                     'axes.spines.right': False, 'savefig.facecolor': 'white',
                     'axes.grid': True, 'axes.axisbelow': True,
                     'grid.color': '.78', 'grid.linewidth': .6, 'grid.alpha': .65})


@lru_cache(maxsize=2)
def surface_triangles(rows, cols):
    """Two triangles per cell; all original grid vertices are retained."""
    indices = np.arange(rows*cols, dtype=np.int32).reshape(rows, cols)
    a, b = indices[:-1,:-1].ravel(), indices[:-1,1:].ravel()
    c, d = indices[1:,1:].ravel(), indices[1:,:-1].ravel()
    return np.concatenate((np.column_stack((a,b,c)), np.column_stack((a,c,d))))


class GouraudSurface(Collection):
    """Interpolated vertex colors for a 3D height field using Matplotlib's renderer.

    mplot3d's plot_surface colors each face uniformly. This collection projects
    and depth-sorts full-resolution triangles, then uses the same Gouraud
    primitive as pcolormesh. Like mplot3d, visibility uses painter ordering,
    not a GPU depth buffer. PDF/SVG use rasterized surface artwork.
    """
    def __init__(self, xx, yy, zz, norm):
        super().__init__(cmap=CMAP, norm=norm, alpha=1.0, rasterized=True)
        self.vertices = np.column_stack((xx.ravel(), yy.ravel(), zz.ravel()))
        self.triangles = surface_triangles(*zz.shape)
        self.set_array(zz.ravel())

    def get_paths(self):
        """Exclude projected 3D triangles from 2D tight-layout bounds."""
        return []

    def do_3d_projection(self):
        px, py, pz = proj3d.proj_transform(*self.vertices.T, self.axes.M)
        self._projected = np.column_stack((px, py))
        depths = pz[self.triangles].mean(axis=1)
        self._depth_order = np.argsort(depths, kind='stable')[::-1]
        return float(np.min(pz))

    @allow_rasterization
    def draw(self, renderer):
        if not self.get_visible():
            return
        colors = self.to_rgba(self.get_array(), alpha=1.0)
        gc = renderer.new_gc()
        self._set_gc_clip(gc)
        transform = self.get_transform().frozen()
        # Bound temporary triangle/color buffers for large grids.
        for start in range(0, len(self._depth_order), 32768):
            triangles = self.triangles[self._depth_order[start:start+32768]]
            renderer.draw_gouraud_triangles(gc, self._projected[triangles],
                                           colors[triangles], transform)
        gc.restore()
        self.stale = False


class StreamingGIF:
    """Write one independently quantized frame at a time using Pillow's plugin API."""
    def __init__(self, path, fps):
        self.path, self.duration = Path(path), max(10, int(round(1000/fps/10))*10)
        self.first = True

    def __enter__(self):
        self.temp = self.path.with_suffix('.partial.gif')
        self.stream = self.temp.open('wb')
        return self

    def append(self, rgb):
        indexed = Image.fromarray(np.asarray(rgb, dtype=np.uint8)).quantize(colors=256)
        if self.first:
            header, _ = GifImagePlugin.getheader( indexed, info={'loop': 0, 'duration': self.duration} )
            headerBytes = b''.join(header)
            if b'NETSCAPE2.0' not in headerBytes:
                headerBytes += b'\x21\xff\x0bNETSCAPE2.0\x03\x01\x00\x00\x00'
            self.stream.write(headerBytes)
        for block in GifImagePlugin.getdata(indexed, duration=self.duration, disposal=2,
                                           include_color_table=not self.first):
            self.stream.write(block)
        self.first = False

    def __exit__(self, kind, value, traceback):
        if kind is None:
            self.stream.write(b';')
        self.stream.close()
        if kind is None:
            self.temp.replace(self.path)
        else:
            self.temp.unlink(missing_ok=True)


def save_figure(fig, output, directory, testcase, stem, formats, dpi):
    folder = output / directory
    folder.mkdir(parents=True, exist_ok=True)
    for extension in formats:
        fig.savefig(folder / f'{stem}.{extension}', dpi=dpi, bbox_inches='tight', transparent=False)
    plt.close(fig)


def nonempty_limits(lower, upper):
    if upper <= lower:
        pad = max(abs(lower)*.01, 1e-12)
        return lower - pad, upper + pad
    return lower, upper


def time_limits(times):
    return nonempty_limits(float(times[0]), float(times[-1]))


def fixed_time_labels(times, precision=6):
    """Fixed-width time labels; use one shared exponent for very small/large times."""
    values = np.asarray(times, dtype=float)
    largest = float(np.max(np.abs(values)))
    exponent = int(np.floor(np.log10(largest))) if largest and (largest < 1e-4 or largest >= 1e4) else 0
    scaled = values / (10.0 ** exponent)
    sign = '+' if np.any(values < 0) else ''
    width = max(len(format(float(v), f'{sign}.{precision}f')) for v in scaled)
    suffix = rf'\times10^{{{exponent}}}' if exponent else ''
    return [r'\mathtt{' + format(float(v), f'{sign}0{width}.{precision}f') + '}' + suffix
            for v in scaled]


def parameter_title(params, time, time_label=None):
    """Mathtext accepts LaTeX-style expressions without a LaTeX installation."""
    dt = params.get('dt', float('nan'))
    if np.isfinite(dt):
        mantissa, exponent = f'{dt:.1e}'.split('e')
        dt_text = rf'{mantissa}\times10^{{{int(exponent)}}}'
    else:
        dt_text = r'\mathrm{n/a}'
    if time_label is None:
        time_label = fixed_time_labels([time])[0]
    return (rf'$N_1={params["N1"]},\;N_2={params["N2"]},\;'
            rf'\Delta t={dt_text},\;'
            rf'K={params.get("K",float("nan")):.3g},\;'
            rf'\ell_1={params["ell1"]:g},\;\ell_2={params["ell2"]:g},\;t={time_label}$')


def spectrum_axis(ax, radii, spectra):
    ax.set(xlabel=r'Radial wavenumber $k$', ylabel=r'$E(k)$', title='Energy spectrum',
           xlim=nonempty_limits(float(radii[0]), float(radii[-1])), yscale='log')
    positive = spectra[spectra > 0]
    ymax = max(1e-19, 1.5 * positive.max()) if positive.size else 1
    ax.set_ylim(1e-20, ymax)


def field_axes(ax, length):
    ax.set(xlabel=r'$x_1/(2\pi)$', ylabel=r'$x_2/(2\pi)$', xlim=(0, length), ylim=(0, length))


def show_field(ax, u, length, norm=None, coordinates=None, **kwargs):
    field_axes(ax, length)
    # Close the periodic mesh with copies of its first row/column; original
    # samples stay at their actual node coordinates (no half-cell offset).
    x = np.arange(u.shape[1])*length/u.shape[1] if coordinates is None else coordinates
    y = np.arange(u.shape[0])*length/u.shape[0] if coordinates is None else coordinates
    ax.grid(False)
    obj = ax.pcolormesh(np.r_[x,length], np.r_[y,length], np.pad(u,((0,1),(0,1)),mode='wrap'),
                        cmap=CMAP, norm=norm, shading='gouraud', rasterized=True,
                        alpha=1.0, **kwargs)
    ax.set_aspect('equal')
    ax.grid(False)
    return obj


def update_field(obj, u):
    obj.set_array(np.pad(u,((0,1),(0,1)),mode='wrap').ravel())


def surface(ax, u, length, cap, norm=None, coordinates=None, zlimits=None):
    n = u.shape[0]
    step = max(1, int(np.ceil(n / cap))) if cap else 1
    x = np.arange(n) * length/n if coordinates is None else coordinates
    xx, yy = np.meshgrid(x[::step], x[::step])
    if norm is None:
        norm = Normalize(*nonempty_limits(float(u.min()), float(u.max())))
    obj = GouraudSurface(xx, yy, u[::step, ::step], norm)
    obj.set_transform(ax.transData)
    ax.add_collection(obj, autolim=False)
    low, high = nonempty_limits(*(zlimits if zlimits is not None else (float(u.min()), float(u.max()))))
    floor = low  # Keep the contour on the lower z boundary without extra padding.
    displayed = u[::step, ::step]
    obj._ks_contour = None
    if min(displayed.shape) >= 2 and np.ptp(displayed) > 0:
        obj._ks_contour = ax.contour(xx, yy, displayed,
            levels=np.linspace(float(displayed.min()), float(displayed.max()), 12)[1:-1],
            zdir='z', offset=floor, colors='0.25', linewidths=.6, alpha=1.0)
    ax.set_zlim(floor, high)
    ax.set_facecolor('white')
    for axis in (ax.xaxis, ax.yaxis, ax.zaxis):
        axis.pane.set_alpha(1.0)
        axis.pane.set_facecolor((1, 1, 1, 1))
        axis.pane.set_edgecolor('white')
    ax.grid(True)
    field_axes(ax, length)
    ax.set_box_aspect((1, 1, .7))
    return obj


def remove_surface(obj):
    """Remove the old surface and its contour before drawing the next frame."""
    contour = getattr(obj, '_ks_contour', None)
    if contour is not None:
        if hasattr(contour, 'remove'):
            contour.remove()
        else:  # Matplotlib 3.7 stores separate contour collections.
            for collection in contour.collections:
                collection.remove()
    obj.remove()


def render_static(output, name, params, initial, terminal, radii, initial_spectra,
                  spectra, energy, sampled, fits, radial_weights, labels, args):
    length = params['ell1']
    for suffix, field, directory, timepoint in [('IC',initial,'initialData',sampled[0,0]),
                                               ('TC',terminal,'terminalData',sampled[-1,0])]:
        fig, ax = plt.subplots(figsize=(8,7), layout='constrained')
        fig.colorbar(show_field(ax, field, length), ax=ax, shrink=.8)
        ax.set_title(suffix + '\n' + parameter_title(params, timepoint), fontsize=10)
        save_figure(fig, output, directory, name, 'contour'+suffix+name, args.formats, args.dpi)
        fig = plt.figure(figsize=(8,7), layout='constrained')
        ax = fig.add_subplot(projection='3d')
        surface(ax, field, length, args.surface_grid)
        ax.set_title(suffix + '\n' + parameter_title(params, timepoint), fontsize=10)
        save_figure(fig, output, directory, name, 'surface'+suffix+name, args.formats, args.dpi)
    fig, ax = plt.subplots(figsize=(8,7), layout='constrained')
    for i in range(min(1 + args.extra_initial_steps, initial_spectra.shape[1] - 1)):
        ax.semilogy(radii, initial_spectra[:,i+1], '--', label=f'State {i+1}')
    spectrum_axis(ax, radii, initial_spectra[:,1:])
    ax.set_title('Initial energy spectra')
    ax.legend(frameon=False)
    save_figure(fig, output, 'fourierSpectrumEvolution', name, 'spectrumIC'+name, args.formats, args.dpi)
    fig, ax = plt.subplots(figsize=(8,7), layout='constrained')
    ax.semilogy(radii, spectra[:,-1], '--', label=r'$E(k)$')
    spectrum_axis(ax, radii, spectra)
    ax.legend(frameon=False)
    save_figure(fig, output, 'fourierSpectrumEvolution', name, 'spectrumTC'+name, args.formats, args.dpi)
    # Useful standalone exports of diagnostics otherwise only visible in MATLAB's GIF.
    fig, axes = plt.subplots(1,3,figsize=(16,4.5), layout='constrained')
    for j, label, color in [(1,r'$L^2$','b'),(2,r'$H^1$','r'),(3,r'$H^2$','g')]:
        axes[0].semilogy(energy[:,0], energy[:,j], color, label=label)
    axes[0].set(title='Energy evolution', xlabel='Time', ylabel='Squared norm')
    axes[0].legend(frameon=False)
    axes[1].plot(sampled[:,0], fits[:,1])
    axes[1].axhline(2*np.pi*length/params['N1'], ls='--', color='gray', label='Grid spacing')
    axes[1].set(title='Analyticity strip width', xlabel='Time', ylabel=r'$\delta$')
    axes[1].legend(frameon=False)
    for i, label in enumerate(labels):
        axes[2].plot(sampled[:,0], radial_weights[i], label=str(label))
    axes[2].set(title='Radial cosine weights', xlabel='Time', ylabel='Normalized weight')
    if labels:
        axes[2].legend(title='Radial groups', loc='upper left', bbox_to_anchor=(1.02,1),
                       borderaxespad=0, frameon=False, fontsize=8,
                       ncols=max(1,int(np.ceil(len(labels)/14))))
    save_figure(fig, output, 'energyEvolution', name, 'diagnostics'+name, args.formats, args.dpi)


def render_movie(output, name, params, cache, stride, energy, sampled, radii, spectra,
                 fits, fit_curves, radial_weights, labels, limits, frame_limits, args):
    need_mp4 = args.movie in ('mp4', 'both')
    if need_mp4 and not writers.is_available('ffmpeg'):
        raise ValueError('MP4 requires FFmpeg on PATH. Install FFmpeg or use --movie gif.')
    folder = output / 'forwardSolution'
    folder.mkdir(parents=True, exist_ok=True)
    fig = plt.figure(figsize=(16,9), dpi=args.dpi)
    fig.subplots_adjust(left=.065,right=.97,bottom=.08,top=.84,wspace=.42,hspace=.55)
    axes = [fig.add_subplot(2,4,1,projection='3d' if args.movie_field=='surface' else None)]
    axes += [fig.add_subplot(2,4,i) for i in range(2,9)]
    length, n = params['ell1'], params['N1']
    color_limits = frame_limits[0] if args.color_scale == 'frame' else limits
    norm = Normalize(*nonempty_limits(*color_limits))
    coords = np.arange(0,n,stride) * length/n
    if args.movie_field == 'surface':
        field_artist = surface(axes[0], cache[0], length, args.surface_grid, norm, coords, frame_limits[0])
    else:
        field_artist = show_field(axes[0], cache[0], length, norm, coords)
    axes[0].set_title('Solution field')
    tiled = show_field(axes[1], np.tile(cache[0],(2,2)), 2*length, norm, np.r_[coords,coords+length])
    axes[1].axvline(length, ls='--', color='black', lw=.8)
    axes[1].axhline(length, ls='--', color='black', lw=.8)
    axes[1].set_title('Tiled solution field')
    energy_markers = []
    for j, label, color in [(3,r'$H^2$','g'),(2,r'$H^1$','r'),(1,r'$L^2$','b')]:
        axes[2].semilogy(energy[:,0], energy[:,j], color, label=label)
        marker, = axes[2].plot([], [], 'ko', ms=4)
        energy_markers.append((marker,j))
    axes[2].set(title='Energy',xlabel='Time',ylabel='Squared norm',xlim=time_limits(sampled[:,0]))
    axes[2].legend(frameon=False, fontsize=9)
    spec, = axes[3].semilogy(radii, spectra[:,0], '--', label=r'$E(k)$')
    fit_line, = axes[3].semilogy(radii, fit_curves[:,0], 'r--', label=r'$Ce^{-2\delta k}$')
    spectrum_axis(axes[3], radii, spectra)
    axes[3].legend(frameon=False, fontsize=9)
    axes[4].plot(sampled[:,0], fits[:,1])
    strip_marker, = axes[4].plot([], [], 'ko', ms=4)
    axes[4].axhline(2*np.pi*length/n, ls='--', color='gray')
    axes[4].set(title='Analyticity strip width', xlabel='Time', ylabel=r'$\delta$',
                xlim=time_limits(sampled[:,0]))
    weights_max = max(1e-12, float(radial_weights.max())) if radial_weights.size else 1
    bars = axes[5].bar(np.arange(len(labels)), radial_weights[:,0])
    axes[5].set(title='Projection coefficient weights',xlabel=r'$(k_1,k_2)$',
                ylabel='Normalized weight',ylim=(0,1.1*weights_max))
    axes[5].set_xticks(np.arange(len(labels)), [str(label) for label in labels], rotation=60)
    axes[5].tick_params(axis='x',labelsize=8)
    for i,label in enumerate(labels):
        axes[6].plot(sampled[:,0],radial_weights[i],label=str(label))
    mode_markers, = axes[6].plot([],[],'ko',ms=4)
    axes[6].set(title='Modal energy',xlabel='Time',ylabel='Normalized weight',
                xlim=time_limits(sampled[:,0]),ylim=(0,1.1*weights_max))
    axes[7].axis('off')
    if labels:
        handles, names = axes[6].get_legend_handles_labels()
        axes[6].legend(handles,names,title='Radial groups',loc='upper left',
                       bbox_to_anchor=(1.02,1),borderaxespad=0,frameon=False,
                       fontsize=8,ncols=max(1,int(np.ceil(len(labels)/14))))
    else:
        axes[5].text(.5,.5,'No unstable modes',transform=axes[5].transAxes,ha='center')
    heading = fig.suptitle('',fontsize=16,y=.97)
    time_labels = fixed_time_labels(sampled[:,0])
    mp4_path = folder / ('movie'+name+'.mp4')
    mp4_temp = mp4_path.with_suffix('.partial.mp4')
    try:
        with ExitStack() as stack:
            gif = stack.enter_context(StreamingGIF(folder / ('movie'+name+'.gif'),args.fps)) \
                  if args.movie in ('gif','both') else None
            if need_mp4:
                writer = FFMpegWriter(fps=args.fps,codec='libx264',
                                     extra_args=['-pix_fmt','yuv420p','-crf','18','-threads','1'])
                stack.enter_context(writer.saving(fig, str(mp4_temp), args.dpi))
            for i, current in enumerate(sampled[:,0]):
                if i:
                    color_limits = frame_limits[i] if args.color_scale == 'frame' else limits
                    norm = Normalize(*nonempty_limits(*color_limits))
                    if args.movie_field == 'surface':
                        remove_surface(field_artist)
                        field_artist = surface(axes[0],cache[i],length,args.surface_grid,norm,coords,frame_limits[i])
                    else:
                        update_field(field_artist, cache[i])
                        field_artist.set_norm(norm)
                    update_field(tiled, np.tile(cache[i],(2,2)))
                    tiled.set_norm(norm)
                for marker,j in energy_markers:
                    marker.set_data([current],[sampled[i,j]])
                spec.set_ydata(spectra[:,i])
                fit_line.set_ydata(fit_curves[:,i])
                strip_marker.set_data([current],[fits[i,1]])
                for bar, weight in zip(bars,radial_weights[:,i]):
                    bar.set_height(weight)
                mode_markers.set_data(np.full(len(labels),current),radial_weights[:,i])
                heading.set_text('Time-dependent solution to the 2D Kuramoto–Sivashinsky equation\n'
                                 + parameter_title(params, current, time_labels[i]))
                if gif:
                    fig.canvas.draw()
                    gif.append(np.asarray(fig.canvas.buffer_rgba())[:,:,:3])
                if need_mp4:
                    writer.grab_frame()
                if i == 0 or (i+1)%10 == 0 or i+1 == len(sampled):
                    LOG.info('Rendered %d/%d movie frames',i+1,len(sampled))
        if need_mp4:
            mp4_temp.replace(mp4_path)
    finally:
        plt.close(fig)
        mp4_temp.unlink(missing_ok=True)
