"""Cosine projections matching the supplied MATLAB definition.

The production path caches basis phases/norms and uses Fourier coefficients.
It does NOT substitute |u_hat|^2 for the phase-dependent cosine projections.
No unused pairwise orthogonality matrix is constructed.
"""
import math
import numpy as np
from scipy.fft import fft2


def matlab_modes(length, kind='full'):
    """Preserve original candidate generation, ordering, signs and duplicates."""
    if not np.isfinite(length) or length <= 0:
        raise ValueError('Length scale must be positive and finite.')
    count = math.ceil(length) ** 2
    candidates = [None] * count
    k1, k2 = 0, 1
    for i in range(count):
        if candidates[i] is not None:
            continue
        candidates[i] = (k1, k2)
        if k2 > k1:
            if i + 1 == len(candidates):
                candidates.append((k2, k1))
            else:
                candidates[i + 1] = (k2, k1)
            k1 += 1
        else:
            k2 += 1
            k1 = 0
    positive = np.asarray(candidates, dtype=int)
    both = np.concatenate((positive, -positive))
    # MATLAB appends two sign variants per candidate, not two separate blocks.
    extra = [v for a, b in positive if a > 0 and b > 0
             for v in [(-a, b), (a, -b)]]
    if extra:
        both = np.concatenate((both, extra))
    radii = np.hypot(both[:, 0], both[:, 1])
    if kind in ('full', 'active'):
        return both[radii < length]
    if kind != 'dominant':
        raise ValueError('kind must be full, active, or dominant.')
    modes = both[np.abs(length - np.sqrt(2) * radii) < 1e-2]
    special = {3:(2,4), 6:(4,8), 9:(8,10), 13:(10,16),
               17:(16,18), 19:(18,20), 23:(20,26), 29:(26,32)}
    for square, shells in special.items():
        if abs(length - np.sqrt(square)) < 1e-10:
            mask = np.zeros(len(both), dtype=bool)
            for shell in shells:
                mask |= np.abs(np.sqrt(shell) - np.sqrt(2) * radii) < 1e-2
            modes = np.concatenate((modes, both[mask]))
    return modes


def minimum_shift(u):
    """Literal MATLAB minimum alignment: min row goes to last, column to first."""
    row, col = np.unravel_index(np.argmin(u.ravel(order='F')), u.shape, order='F')
    return -(int(row) + 1), -int(col)


class CosineProjector:
    def __init__(self, length, n, kind='full', workers=1):
        self.length, self.n, self.kind, self.workers = length, n, kind, workers
        self.modes = matlab_modes(length, kind)
        if np.any(np.abs(self.modes) >= n / 2):
            raise ValueError('Selected modes reach the Nyquist limit; increase grid resolution.')
        # Keep the original arithmetic order to reduce phase differences at ties.
        points = (2 * np.pi * length) * np.linspace(0, 1 - 1/n, n)
        self.x, self.y = np.meshgrid(points, points)
        self.phases = np.empty(len(self.modes), dtype=complex)
        self.norms = np.empty(len(self.modes))
        for i, mode in enumerate(self.modes):
            basis = self.basis(mode)
            self.norms[i] = np.linalg.norm(basis)
            # Scalar normalization can affect last-bit ties, as in MATLAB.
            scaled = basis / (self.norms[i] * 2 * np.pi * length / n)
            sr, sc = minimum_shift(scaled)
            self.phases[i] = np.exp(2j * np.pi * (mode[0]*sc + mode[1]*sr) / n)

    def basis(self, mode):
        return np.cos(mode[0] * (1/self.length) * self.x
                      + mode[1] * (1/self.length) * self.y)

    def coefficients(self, u):
        u = np.asarray(u, dtype=float)
        if u.shape != (self.n, self.n):
            raise ValueError('Field has the wrong shape for this projector.')
        centered = u - u.mean()
        norm = np.linalg.norm(centered)
        if norm == 0 or not len(self.modes):
            return np.zeros(len(self.modes))
        sr, sc = minimum_shift(u)
        transform = fft2(centered, workers=self.workers)
        m = self.modes
        field_phase = np.exp(-2j * np.pi * (m[:,0]*sc + m[:,1]*sr) / self.n)
        return (transform[m[:,1] % self.n, m[:,0] % self.n]
                * field_phase * self.phases).real / (norm * self.norms)

    def validate(self, u):
        """Return original match_score, amplitudes and modes; match is optional cost."""
        u = np.asarray(u, dtype=float)
        amps = self.coefficients(u)
        modes = self.modes.copy()
        if self.kind == 'active':
            keep = np.abs(amps) >= 1e-3
            amps, modes = amps[keep], modes[keep]
        reconstruction = np.zeros_like(u)
        for amp, mode in zip(amps, modes):
            basis = self.basis(mode)
            basis /= np.linalg.norm(basis) * 2 * np.pi * self.length / self.n
            basis *= amp
            reconstruction += np.roll(basis, minimum_shift(basis), axis=(0, 1))
        shifted = np.roll(u, minimum_shift(u), axis=(0, 1))
        shifted -= shifted.mean()
        denom = np.linalg.norm(shifted) * np.linalg.norm(reconstruction)
        score = float(1 - np.sum(shifted * reconstruction) / denom) if denom else np.nan
        return score, amps, modes


def eigenfunction_validation(u_IC_opt, L_s1, N, type='full'):
    """Drop-in numerical helper; flattened input follows MATLAB column order.

    For many fields construct one CosineProjector and reuse it instead.
    The supplied algorithm is for square grids/domains; no anisotropic
    interpretation is silently introduced.
    """
    u = np.asarray(u_IC_opt, dtype=float).reshape((N, N), order='F')
    return CosineProjector(L_s1, N, type).validate(u)
