"""Port of the supplied expfit_delta.m, including its actual fit window."""
import numpy as np


def expfit_delta(v_radii, v_mean, prevdelta=None):
    """Return (C, delta) for E(k) = C exp(-2 delta k).

    prevdelta is accepted for calling compatibility. The MATLAB source tries
    just one window (3:-0.5:3), so it does not meaningfully use prevdelta.
    Invalid fits raise ValueError; the driver records these as NaN plus a warning.
    """
    k = np.asarray(v_radii, dtype=float).ravel()
    energy = np.asarray(v_mean, dtype=float).ravel()
    if k.size != energy.size:
        raise ValueError('Radii and energies must have equal lengths.')
    below = np.flatnonzero(energy < 1e-15)
    vcut = int(below[0]) + 1 if below.size else k.size  # MATLAB one-based
    vstop = min(max(10, vcut), k.size - 2)
    vstart = max(1, int(np.ceil(vstop / 3)))
    k, energy = k[vstart - 1:vstop], energy[vstart - 1:vstop]
    valid = np.isfinite(k) & np.isfinite(energy) & (energy > 0)
    k, energy = k[valid], energy[valid]
    if k.size < 5:
        raise ValueError('Fewer than five valid spectral points in fit window.')
    order = np.argsort(k, kind='stable')
    k, energy = k[order], energy[order]
    valid = energy > 1e-12 * energy.max()
    if valid.sum() >= 5:
        k, energy = k[valid], energy[valid]
    if np.ptp(k) == 0:
        raise ValueError('Fit radii must not all be equal.')
    slope, intercept = np.polyfit(k, np.log(energy), 1)
    with np.errstate(over='ignore'):
        coefficient = float(np.exp(intercept))
    if not np.isfinite(coefficient):
        raise ValueError('Fitted prefactor overflowed.')
    return coefficient, float(-slope / 2)
