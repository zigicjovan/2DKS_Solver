# Python post-processing for the 2D KS solver

## Install and run

Extract the archive. Keep all its Python modules together in `ksPostprocess`.

```bash
cd ksPostprocess
python3 -m venv .venv
source .venv/bin/activate
python3 -m pip install numpy>=1.24,<3 scipy>=1.10,<2 matplotlib>=3.7,<4 Pillow>=10,<14

# Replace the path with your actual testcase directory, not its ForwardSolution folder.
python3 generateFigures.py /path/to/Data/_IC_s1_N1_256_N2_256_...
```

The testcase's basename must contain the real parameters, as in your existing
directory names; the `...` above is just a placeholder. A normal run generates
PDFs, PNGs, a GIF, and numerical diagnostics. The default GIF is 20 fps with a
1600-by-900 frame at 100 dpi. Figure windows never open.

```bash
# Specify the saved-state count if you want an explicit consistency check.
python generateFigures.py "$testcase" --states 100

# Static figures and numerical diagnostics, without a movie.
python generateFigures.py "$testcase" --plots figures

# Numerical diagnostics only.
python generateFigures.py "$testcase" --plots none

# MP4, or both formats; FFmpeg must be installed separately.
python generateFigures.py "$testcase" --movie mp4
python generateFigures.py "$testcase" --movie both

# Faster rendering: a heatmap replaces the first panel's 3D surface.
python generateFigures.py "$testcase" --movie-field heatmap

# Optional fixed color limits for comparing amplitudes across time.
# Default colors rescale per frame, matching the individual field figures.
python generateFigures.py "$testcase" --color-scale global

# Full-resolution rendering is now the default; these flags make it explicit.
python generateFigures.py "$testcase" --max-display-grid 0 --surface-grid 0

# Process multiple testcases sequentially.
python generateFigures.py /path/to/Data/_IC_*
```

## Output

The six corresponding static views are `contourIC`, `surfaceIC`, `contourTC`,
`surfaceTC`, `spectrumIC`, and `spectrumTC`, each followed by the testcase name.
The contour-named files are colored field maps, matching the purpose of your
top-down MATLAB surfaces; the rendering is not intended to be pixel-identical.
An additional static diagnostic figure contains energy, strip width and
radial cosine-weight evolution.

The movie preserves the seven analytical panels: field, tiled field, energy,
spectrum and fit, strip width, radial weight bars, and radial-weight evolution.
The radial-mode legend is attached immediately to the right of the modal-energy
axes, using the space available in the eighth grid position.

`diagnostics.npz` contains:

- parameters (JSON string), full energy data, sampled energy/time rows;
- radial wavenumbers, spectra, fitted C and delta, and fitted spectra;
- original mode ordering, signed cosine amplitudes, normalized squared
  amplitudes, exact radial groups and their summed weights;
- centered initial and terminal physical fields, the common centering shift,
  and full-resolution field minimum/maximum values for every frame.
