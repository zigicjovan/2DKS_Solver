#!/usr/bin/env python3
"""Generate scaling-law and solution-branch plots from 2D-KS .dat files.

Python replacement for ``generatePlots.m``. By default, the script reads data
and writes PDFs in the current directory. An alternative directory may be
given as the sole positional argument.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


ELLS = np.array([1.05, 1.20, 1.35, 1.50, 1.65, 1.80, 1.95])
K_VALUES = np.array([1.0e4, 3.2e4, 1.0e5, 3.2e5, 1.0e6])
K_FILE_STRINGS = ["1.0e+04", "3.2e+04", "1.0e+05", "3.2e+05", "1.0e+06"]
K_LABELS = [r"10^{4.0}", r"10^{4.5}", r"10^{5.0}", r"10^{5.5}", r"10^{6.0}"]
FONT_SIZE = 20
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
            "legend.fontsize": FONT_SIZE - 4,
            "xtick.labelsize": FONT_SIZE,
            "ytick.labelsize": FONT_SIZE,
            "mathtext.fontset": "cm",
            "font.family": "serif",
            "pdf.fonttype": 42,
        }
    )


def load_data(path: Path) -> np.ndarray:
    data = np.genfromtxt(path, comments="#", ndmin=2)
    if data.size == 0 or data.shape[1] < 7:
        raise ValueError(f"{path} must contain at least seven numeric columns")
    return data


def fit_loglog(x: np.ndarray, y: np.ndarray) -> tuple[float, float]:
    valid = np.isfinite(x) & np.isfinite(y) & (x > 0) & (y > 0)
    if np.count_nonzero(valid) < 2:
        raise ValueError("at least two positive finite points are required for a fit")
    slope, intercept = np.polyfit(np.log10(x[valid]), np.log10(y[valid]), 1)
    return float(slope), float(intercept)


def fit_semilog_y(x: np.ndarray, y: np.ndarray) -> tuple[float, float]:
    valid = np.isfinite(x) & np.isfinite(y) & (x > 0) & (y > 0)
    if np.count_nonzero(valid) < 2:
        raise ValueError("at least two positive finite points are required for a fit")
    slope, intercept = np.polyfit(x[valid], np.log10(y[valid]), 1)
    return float(slope), float(intercept)


def new_figure() -> tuple[plt.Figure, plt.Axes]:
    return plt.subplots(figsize=(12, 8))


def finish_figure(
    fig: plt.Figure,
    ax: plt.Axes,
    output: Path,
    *,
    xlabel: str,
    ylabel: str,
    title: str,
    subtitle: str,
    legend_location: str,
    legend_columns: int = 1,
    xlim: tuple[float, float] | None = None,
    box: bool = True,
) -> None:
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(f"{title}\n{subtitle}", pad=14)
    ax.grid(True, which="major", alpha=0.35)
    ax.grid(False, which="minor")
    if not box:
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
    if xlim is not None:
        ax.set_xlim(*xlim)
    if legend_location == "outside below":
        ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.18), ncols=legend_columns)
    else:
        ax.legend(loc=legend_location, ncols=legend_columns)
    fig.tight_layout()
    fig.savefig(output, format="pdf", bbox_inches="tight")
    plt.close(fig)


def print_fit_summary(name: str, slopes: np.ndarray) -> tuple[float, float]:
    mean = float(np.mean(slopes))
    sample_std = float(np.std(slopes, ddof=1))
    print(f"{name}:")
    print(slopes)
    print(f"Mean exponent = {mean:.8f}")
    print(f"Sample std    = {sample_std:.8f}\n")
    return mean, sample_std


def plot_energy_vs_k(data_dir: Path, compensated: bool = False) -> tuple[float, float]:
    fig, ax = new_figure()
    slopes = np.zeros(ELLS.size)
    intercepts = np.zeros(ELLS.size)

    for i, ell in enumerate(ELLS):
        data = load_data(data_dir / f"powerlawK_IC_s1_ell1_{ell:.2f}_ell2_{ell:.2f}.dat")
        k, kstar = data[:, 0], data[:, 6]
        plotted = kstar / k**2 if compensated else kstar
        slopes[i], intercepts[i] = fit_loglog(k, plotted)
        if compensated:
            label = rf"$\widetilde K_{{{ell:.2f}}}/K^2 \approx 10^{{{intercepts[i]:.2f}}}K^{{{slopes[i]:.2f}}}$"
        else:
            label = rf"$\widetilde K_{{{ell:.2f}}} \approx 10^{{{intercepts[i]:.2f}}}K^{{{slopes[i]:.2f}}}$"
        ax.loglog(k, plotted, "-o", color=COLORS[i], linewidth=1.5, markersize=6, label=label)

    mean, std = print_fit_summary("Energy scaling K exponents", slopes)
    if compensated:
        ylabel = r"Compensated $\widetilde{K}_{\ell}/K^2$"
        subtitle = rf"$K=[10^4,\ldots,10^6]$, $\ell=[1.05,\ldots,1.95]$: $\widetilde{{K}}_\ell/K^2\sim K^{{{mean:.2f}\,\pm\,{std:.2f}}}$"
        output = data_dir / "powerlawK_comp.pdf"
        location, columns = "outside below", 3
    else:
        ylabel = r"Peak Transient Energy $\widetilde{K}_{\ell}$"
        subtitle = rf"$K=[10^4,10^{{4.5}},\ldots,10^6]$, $\ell=[1.05,1.20,\ldots,1.95]$: $\widetilde{{K}}_\ell\sim K^{{{mean:.2f}\,\pm\,{std:.2f}}}$"
        output = data_dir / "powerlawK.pdf"
        location, columns = "lower right", 1
    finish_figure(
        fig,
        ax,
        output,
        xlabel=r"Initial Energy $K$",
        ylabel=ylabel,
        title=r"$L^2$ Energy-Scaling Law for 2D Kuramoto-Sivashinsky",
        subtitle=subtitle,
        legend_location=location,
        legend_columns=columns,
        xlim=(1e4, 1e6),
    )
    return mean, std


def plot_energy_vs_ell(
    data_dir: Path, *, power_fit: bool = False, compensated: bool = False
) -> tuple[float, float]:
    fig, ax = new_figure()
    slopes = np.zeros(K_VALUES.size)
    intercepts = np.zeros(K_VALUES.size)

    for i, (k_string, k_label) in enumerate(zip(K_FILE_STRINGS, K_LABELS)):
        data = load_data(data_dir / f"powerlawL_IC_s1_K_{k_string}.dat")
        ell_squared, kstar = data[:, 1] ** 2, data[:, 6]
        plotted = kstar / ell_squared if compensated else kstar
        if power_fit:
            slopes[i], intercepts[i] = fit_loglog(ell_squared, plotted)
            prefix = r"/\ell^2" if compensated else ""
            label = rf"$\widetilde K_{{{k_label}}}{prefix}\approx 10^{{{intercepts[i]:.2f}}}\ell^{{2({slopes[i]:.2f})}}$"
        else:
            slopes[i], intercepts[i] = fit_semilog_y(ell_squared, plotted)
            label = rf"$\widetilde K_{{{k_label}}}\approx 10^{{{intercepts[i]:.2f}}}10^{{{slopes[i]:.2f}\ell^2}}$"
        ax.semilogy(ell_squared, plotted, "-o", color=COLORS[i], linewidth=1.5, markersize=6, label=label)

    mean, std = print_fit_summary("Domain scaling energy exponents", slopes)
    if compensated:
        output = data_dir / "powerlawL_comp.pdf"
        ylabel = r"Compensated $\widetilde{K}_{K}/\ell^2$"
        law = rf"$\widetilde{{K}}_K/\ell^2\sim\ell^{{2({mean:.2f}\,\pm\,{std:.2f})}}$"
    elif power_fit:
        output = data_dir / "powerlawL.pdf"
        ylabel = r"Peak Transient Energy $\widetilde{K}_{K}$"
        law = rf"$\widetilde{{K}}_K\sim\ell^{{2({mean:.2f}\,\pm\,{std:.2f})}}$"
    else:
        output = data_dir / "powerlawL_lin.pdf"
        ylabel = r"Peak Transient Energy $\widetilde{K}_{K}$"
        law = rf"$\widetilde{{K}}_K\sim 10^{{{mean:.2f}\ell^2\,\pm\,{std:.2f}}}$"
    if power_fit:
        ax.set_xscale("log")
    finish_figure(
        fig,
        ax,
        output,
        xlabel=r"Domain Factor $\ell^2$",
        ylabel=ylabel,
        title="Isotropic Domain-Scaling Law for 2D Kuramoto-Sivashinsky",
        subtitle=rf"$K=[10^4,\ldots,10^6]$, $\ell=[1.05,\ldots,1.95]$: {law}",
        legend_location="outside below",
        legend_columns=3,
        box=False,
    )
    return mean, std


def plot_universal_energy(data_dir: Path) -> tuple[float, float]:
    fig, ax = new_figure()
    slopes = np.zeros(ELLS.size)
    intercepts = np.zeros(ELLS.size)
    for i, ell in enumerate(ELLS):
        data = load_data(data_dir / f"powerlawK_IC_s1_ell1_{ell:.2f}_ell2_{ell:.2f}.dat")
        k, kstar = data[:, 0], data[:, 6]
        plotted = kstar / ell**2
        slopes[i], intercepts[i] = fit_loglog(k, plotted)
        label = rf"$\widetilde K_{{{ell:.2f}}}/\ell^2\approx 10^{{{intercepts[i]:.2f}}}K^{{{slopes[i]:.2f}}}$"
        ax.loglog(k, plotted, "-o", color=COLORS[i], linewidth=1.5, markersize=6, label=label)

    mean, std = print_fit_summary("Universal energy scaling K exponents", slopes)
    mean_intercept = float(np.mean(intercepts))
    finish_figure(
        fig,
        ax,
        data_dir / "powerlaw_universal.pdf",
        xlabel=r"Initial Energy $K$",
        ylabel=r"Normalized Energy $\widetilde{K}_{\ell}/\ell^2$",
        title="Normalized Energy-Scaling Law for 2D Kuramoto-Sivashinsky",
        subtitle=rf"$K=[10^4,\ldots,10^6]$, $\ell=[1.05,\ldots,1.95]$: $\widetilde{{K}}_\ell/\ell^2\approx 10^{{{mean_intercept:.2f}}}K^{{{mean:.2f}}}$",
        legend_location="outside below",
        legend_columns=3,
        xlim=(1e4, 1e6),
    )
    return mean, std


def plot_time_vs_k(data_dir: Path) -> tuple[float, float]:
    fig, ax = new_figure()
    slopes = np.zeros(ELLS.size)
    intercepts = np.zeros(ELLS.size)
    for i, ell in enumerate(ELLS):
        data = load_data(data_dir / f"powerlawTK_IC_s1_ell1_{ell:.2f}_ell2_{ell:.2f}.dat")
        k, tstar = data[:, 0], data[:, 3]
        slopes[i], intercepts[i] = fit_loglog(k, tstar)
        label = rf"$\widetilde T_{{{ell:.2f}}}\approx 10^{{{intercepts[i]:.2f}}}K^{{{slopes[i]:.2f}}}$"
        ax.loglog(k, tstar, "-o", color=COLORS[i], linewidth=1.5, markersize=6, label=label)
    mean, std = print_fit_summary("Time-vs-K exponents", slopes)
    finish_figure(
        fig,
        ax,
        data_dir / "powerlawTK.pdf",
        xlabel=r"Initial Energy $K$",
        ylabel=r"Time Window of Peak Transient Energy $\widetilde{T}_{\ell}$",
        title=r"$L^2$ Energy Time-Scaling Law for 2D Kuramoto-Sivashinsky",
        subtitle=rf"$K=[10^4,10^{{4.5}},\ldots,10^6]$, $\ell=[1.05,1.20,\ldots,1.95]$: $\widetilde{{T}}_\ell\sim K^{{{mean:.2f}\,\pm\,{std:.2f}}}$",
        legend_location="lower left",
        xlim=(1e4, 1e6),
    )
    return mean, std


def plot_time_vs_ell(data_dir: Path, power_fit: bool = False) -> tuple[float, float]:
    fig, ax = new_figure()
    slopes = np.zeros(K_VALUES.size)
    intercepts = np.zeros(K_VALUES.size)
    for i, (k_string, k_label) in enumerate(zip(K_FILE_STRINGS, K_LABELS)):
        data = load_data(data_dir / f"powerlawTL_IC_s1_K_{k_string}.dat")
        ell_squared, tstar = data[:, 1] ** 2, data[:, 3]
        if power_fit:
            slopes[i], intercepts[i] = fit_loglog(ell_squared, tstar)
            label = rf"$\widetilde T_{{{k_label}}}\approx 10^{{{intercepts[i]:.2f}}}\ell^{{2({slopes[i]:.2f})}}$"
        else:
            slopes[i], intercepts[i] = fit_semilog_y(ell_squared, tstar)
            label = rf"$\widetilde T_{{{k_label}}}\approx 10^{{{intercepts[i]:.2f}}}10^{{{slopes[i]:.2f}\ell^2}}$"
        ax.semilogy(ell_squared, tstar, "-o", color=COLORS[i], linewidth=1.5, markersize=6, label=label)
    mean, std = print_fit_summary("Time-vs-domain exponents", slopes)
    if power_fit:
        ax.set_xscale("log")
        output = data_dir / "powerlawTL.pdf"
        law = rf"$\widetilde{{T}}_K\sim\ell^{{2({mean:.2f}\,\pm\,{std:.2f})}}$"
    else:
        output = data_dir / "powerlawTL_lin.pdf"
        law = rf"$\widetilde{{T}}_K\sim10^{{{mean:.2f}\ell^2\,\pm\,{std:.2f}}}$"
    finish_figure(
        fig,
        ax,
        output,
        xlabel=r"Domain Factor $\ell^2$",
        ylabel=r"Time Window of Peak Transient Energy $\widetilde{T}_{K}$",
        title="Isotropic Domain Time-Scaling Law for 2D Kuramoto-Sivashinsky",
        subtitle=rf"$K=[10^4,10^{{4.5}},\ldots,10^6]$, $\ell=[1.05,1.20,\ldots,1.95]$: {law}",
        legend_location="outside below",
        legend_columns=3,
        box=False,
    )
    return mean, std


def plot_branches(data_dir: Path) -> None:
    fig, ax = new_figure()
    for i_k, k_string in enumerate(K_FILE_STRINGS):
        for i_ell, ell in enumerate(ELLS):
            data = load_data(data_dir / f"branch_IC_s1_K_{k_string}_ell1_{ell:.2f}_ell2_{ell:.2f}.dat")
            label = f"{ell:.2f}" if i_k == 0 else "_nolegend_"
            ax.loglog(data[:, 3], data[:, 6], "-o", color=COLORS[i_ell], linewidth=1.5, markersize=5, label=label)
    finish_figure(
        fig,
        ax,
        data_dir / "branch.pdf",
        xlabel="Time Window",
        ylabel=r"Maximum $L^2$ Energy Amplification",
        title=r"Solution Branches for $K=[10^4,10^{4.5},\ldots,10^6]$ and $\ell=[1.05,1.20,\ldots,1.95]$",
        subtitle="",
        legend_location="upper right",
        legend_columns=5,
    )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("data_dir", nargs="?", type=Path, default=Path.cwd(), help="directory containing the .dat files (default: current directory)")
    args = parser.parse_args()
    data_dir = args.data_dir.expanduser().resolve()
    if not data_dir.is_dir():
        parser.error(f"not a directory: {data_dir}")

    configure_matplotlib()
    k_mean, k_std = plot_energy_vs_k(data_dir)
    plot_energy_vs_k(data_dir, compensated=True)
    ell_mean, ell_std = plot_energy_vs_ell(data_dir)
    plot_energy_vs_ell(data_dir, power_fit=True)
    plot_energy_vs_ell(data_dir, power_fit=True, compensated=True)
    plot_universal_energy(data_dir)
    tk_mean, tk_std = plot_time_vs_k(data_dir)
    tl_mean, tl_std = plot_time_vs_ell(data_dir)
    plot_time_vs_ell(data_dir, power_fit=True)
    plot_branches(data_dir)

    print("\n============================================================")
    print("SUMMARY OF FITTED SCALINGS")
    print("============================================================")
    print(f"K* vs K:    exponent = {k_mean:.8f} +/- {k_std:.8f}")
    print(f"K* vs ell:  log10 slope = {ell_mean:.8f} +/- {ell_std:.8f}")
    print(f"T* vs K:    exponent = {tk_mean:.8f} +/- {tk_std:.8f}")
    print(f"T* vs ell:  log10 slope = {tl_mean:.8f} +/- {tl_std:.8f}")
    print("============================================================")


if __name__ == "__main__":
    main()
