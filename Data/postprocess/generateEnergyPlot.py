#!/usr/bin/env python3
"""Plot L2, H1, and H2 energy evolution from 2D-KS .dat files.

Python replacement for ``generateEnergyPlot.m``. Large input files are read
line by line and downsampled to approximately 10,000 points per curve.
"""

from __future__ import annotations

import argparse
import glob
from pathlib import Path
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


ELLS = np.array([1.05, 1.20, 1.35, 1.50, 1.65, 1.80, 1.95])
K_FILE_STRINGS = ["1.0e+04", "3.2e+04", "1.0e+05", "3.2e+05", "1.0e+06"]
FONT_SIZE = 20
MAX_SAMPLES = 100_000
COLORS = np.array(
    [
        [0.0000, 0.4470, 0.7410],
        [0.8500, 0.3250, 0.0980],
        [0.9290, 0.6940, 0.1250],
        [0.4940, 0.1840, 0.5560],
        [0.4660, 0.6740, 0.1880],
        [0.3010, 0.7450, 0.9330],
        [0.6350, 0.0780, 0.1840],
    ]
)


def configure_matplotlib() -> None:
    plt.rcParams.update(
        {
            "font.size": FONT_SIZE,
            "axes.titlesize": FONT_SIZE,
            "axes.labelsize": FONT_SIZE,
            "xtick.labelsize": FONT_SIZE,
            "ytick.labelsize": FONT_SIZE,
            "mathtext.fontset": "cm",
            "font.family": "serif",
            "pdf.fonttype": 42,
        }
    )


def count_rows(path: Path) -> int:
    with path.open("r", encoding="utf-8", errors="replace") as stream:
        return sum(1 for _ in stream)


def read_downsampled(path: Path, max_samples: int = MAX_SAMPLES) -> np.ndarray:
    row_count = count_rows(path)
    # MATLAB round() uses half-away-from-zero; row_count/max_samples >= 0.
    sample_step = max(1, int(row_count / max_samples + 0.5))
    samples: list[list[float]] = []

    with path.open("r", encoding="utf-8", errors="replace") as stream:
        for row_index, line in enumerate(stream):
            if row_index % sample_step == 0:
                line_number = row_index + 1
                try:
                    values = [float(value) for value in line.split()[:4]]
                except ValueError as exc:
                    raise ValueError(f"non-numeric data at {path}:{line_number}") from exc
                if len(values) < 4:
                    raise ValueError(f"fewer than four columns at {path}:{line_number}")
                samples.append(values)

    if not samples:
        raise ValueError(f"no numeric rows found in {path}")
    return np.asarray(samples, dtype=float)


def find_input(data_dir: Path, k_string: str, ell: float) -> Path | None:
    pattern = str(data_dir / f"*K_{k_string}_ell1_{ell:.2f}_ell2_{ell:.2f}_*.dat")
    matches = sorted(Path(name) for name in glob.glob(pattern))
    if not matches:
        print(f"Warning: no file matches {Path(pattern).name}", file=sys.stderr)
        return None
    if len(matches) > 1:
        print(f"Warning: multiple files match {Path(pattern).name}; using {matches[0].name}", file=sys.stderr)
    return matches[0]


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("data_dir", nargs="?", type=Path, default=Path.cwd(), help="directory containing the .dat files (default: current directory)")
    args = parser.parse_args()
    data_dir = args.data_dir.expanduser().resolve()
    if not data_dir.is_dir():
        parser.error(f"not a directory: {data_dir}")

    configure_matplotlib()
    figures_axes = [plt.subplots(figsize=(12, 8)) for _ in range(3)]

    for k_string in reversed(K_FILE_STRINGS):
        for color_index in reversed(range(ELLS.size)):
            ell = ELLS[color_index]
            path = find_input(data_dir, k_string, ell)
            if path is None:
                continue
            try:
                data = read_downsampled(path)
            except (OSError, ValueError) as exc:
                print(f"Warning: skipping {path.name}: {exc}", file=sys.stderr)
                continue
            time = data[:, 0]
            positive_time = time > 0
            color = COLORS[color_index]
            for energy_column, (_, ax) in enumerate(figures_axes, start=1):
                ax.loglog(time[positive_time], data[positive_time, energy_column], color=color, linewidth=1.5)

    subtitle = r"$K=[10^4,10^{4.5},\ldots,10^6]$ and $\ell=[1.05,1.20,\ldots,1.95]$"
    specifications = [
        (r"$L^2$ Energy $\|\phi(t;\widetilde{\varphi})\|^2_{L^2}$", r"$L^2$ Energy Evolution for Energy-Optimized 2D Kuramoto-Sivashinsky", "energyL2evolution.pdf"),
        (r"$H^1$ Energy $\|\phi(t;\widetilde{\varphi})\|^2_{H^1}$", r"$H^1$ Energy Evolution for Energy-Optimized 2D Kuramoto-Sivashinsky", "energyH1evolution.pdf"),
        (r"$H^2$ Energy $\|\phi(t;\widetilde{\varphi})\|^2_{H^2}$", r"$H^2$ Energy Evolution for Energy-Optimized 2D Kuramoto-Sivashinsky", "energyH2evolution.pdf"),
    ]

    for (fig, ax), (ylabel, title, filename) in zip(figures_axes, specifications):
        ax.set_xlim(1e-6, 1e-0)
        #ax.autoscale(axis="x", tight=True)
        ax.autoscale(axis="y", tight=True)
        ax.set_xlabel(r"Time $t$")
        ax.set_ylabel(ylabel)
        ax.set_title(f"{title}\n{subtitle}", pad=14)
        ax.grid(True, which="major", alpha=0.35)
        ax.grid(False, which="minor")
        fig.tight_layout()
        fig.savefig(data_dir / filename, format="pdf", bbox_inches="tight")
        plt.close(fig)


if __name__ == "__main__":
    main()
