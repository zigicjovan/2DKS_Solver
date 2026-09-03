"""NumPy translation of Adam Auton's supplied redblue.m (9 October 2009)."""
import numpy as np
from matplotlib.colors import ListedColormap


def redblue(m=256):
    """Return the original m-by-3 blue-white-red RGB table."""
    if not isinstance(m, (int, np.integer)) or m < 1:
        raise ValueError('Colormap length must be a positive integer.')
    half = m // 2
    if m % 2 == 0:
        ramp = np.arange(half) / max(half - 1, 1)
        r = np.r_[ramp, np.ones(half)]
        g = np.r_[ramp, ramp[::-1]]
    else:
        ramp = np.arange(half) / max(half, 1)
        r = np.r_[ramp, np.ones(half + 1)]
        g = np.r_[ramp, 1, ramp[::-1]]
    return np.column_stack((r, g, r[::-1]))


def redblue_cmap(m=256):
    return ListedColormap(redblue(m), name='redblue')
