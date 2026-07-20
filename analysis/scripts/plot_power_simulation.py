"""Plot power simulation results: power curves + null calibration QQ plots.

Reads power_simulation_results.json and produces publication-quality figures.

Usage:
    uv run python experiments/batch6_vector_extension/scripts/plot_power_simulation.py
"""

import json
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


RESULTS_DIR = Path(__file__).resolve().parent.parent / "results"
FIGURES_DIR = Path(__file__).resolve().parent.parent / "figures"


def load_results(path=None):
    if path is None:
        path = RESULTS_DIR / "power_simulation_results.json"
    with open(path) as f:
        return json.load(f)


def plot_power_curves(results, output_dir):
    """Main figure: power of vector Q vs per-trait scalar Q across rotation angle.

    Panel grid: rows = d (trait count), columns = K (strata count).
    Each panel: x = theta, lines colored by rho, solid = vector, dashed = scalar.
    """
    power = results["power_results"]
    grid = results["metadata"]["grid"]

    d_values = sorted(set(r["d"] for r in power))
    K_values = sorted(set(r["K"] for r in power))
    rho_values = sorted(set(r["rho"] for r in power))

    colors = {0.0: "#2c7bb6", 0.3: "#abd9e9", 0.6: "#fdae61", 0.9: "#d7191c"}
    rho_labels = {0.0: r"$\rho = 0$", 0.3: r"$\rho = 0.3$",
                  0.6: r"$\rho = 0.6$", 0.9: r"$\rho = 0.9$"}

    fig, axes = plt.subplots(
        len(d_values), len(K_values),
        figsize=(4 * len(K_values), 3.5 * len(d_values)),
        sharex=True, sharey=True,
    )
    if len(d_values) == 1:
        axes = axes[np.newaxis, :]
    if len(K_values) == 1:
        axes = axes[:, np.newaxis]

    for i, d in enumerate(d_values):
        for j, K in enumerate(K_values):
            ax = axes[i, j]
            for rho in rho_values:
                subset = [r for r in power
                          if r["d"] == d and r["K"] == K and r["rho"] == rho]
                subset.sort(key=lambda x: x["theta"])
                thetas = [r["theta"] for r in subset]
                vec_pow = [r["vector_power"] for r in subset]
                sca_pow = [r["scalar_power"] for r in subset]

                ax.plot(thetas, vec_pow, "-o", color=colors[rho],
                        markersize=3, linewidth=1.5, label=rho_labels[rho])
                ax.plot(thetas, sca_pow, "--s", color=colors[rho],
                        markersize=3, linewidth=1.0, alpha=0.7)

            ax.axhline(0.05, color="gray", linestyle=":", linewidth=0.8, alpha=0.5)
            ax.axhline(0.80, color="gray", linestyle=":", linewidth=0.8, alpha=0.5)
            ax.set_xlim(-2, 92)
            ax.set_ylim(-0.02, 1.02)

            if i == 0:
                ax.set_title(f"$K = {K}$ strata", fontsize=11)
            if j == 0:
                ax.set_ylabel(f"$d = {d}$ traits\nPower", fontsize=10)
            if i == len(d_values) - 1:
                ax.set_xlabel(r"Rotation angle $\theta$ (degrees)", fontsize=10)

    handles = []
    for rho in rho_values:
        h, = axes[0, 0].plot([], [], "-o", color=colors[rho],
                             markersize=3, linewidth=1.5, label=rho_labels[rho])
        handles.append(h)
    h_vec, = axes[0, 0].plot([], [], "k-", linewidth=1.5, label="Vector $Q_V$")
    h_sca, = axes[0, 0].plot([], [], "k--", linewidth=1.0, label=r"Per-trait $Q$ (Bonf.)")
    handles.extend([h_vec, h_sca])

    fig.legend(handles=handles, loc="lower center", ncol=len(rho_values) + 2,
               fontsize=9, frameon=False, bbox_to_anchor=(0.5, -0.02))

    fig.suptitle("Power: vector $Q_V$ (solid) vs per-trait scalar $Q$ (dashed)",
                 fontsize=13, y=1.01)
    fig.tight_layout()

    output_dir.mkdir(parents=True, exist_ok=True)
    path = output_dir / "power_curves_vector_vs_scalar.pdf"
    fig.savefig(path, bbox_inches="tight", dpi=300)
    path_png = output_dir / "power_curves_vector_vs_scalar.png"
    fig.savefig(path_png, bbox_inches="tight", dpi=300)
    plt.close(fig)
    print(f"Saved power curves to {path}")
    return path


def plot_null_calibration(results, output_dir):
    """QQ plot: empirical vs theoretical chi-squared quantiles under the null."""
    null = results["null_calibration"]

    n_cells = len(null)
    ncols = min(3, n_cells)
    nrows = (n_cells + ncols - 1) // ncols

    fig, axes = plt.subplots(nrows, ncols, figsize=(4 * ncols, 4 * nrows))
    if nrows == 1 and ncols == 1:
        axes = np.array([[axes]])
    elif nrows == 1:
        axes = axes[np.newaxis, :]
    elif ncols == 1:
        axes = axes[:, np.newaxis]

    for idx, cell in enumerate(null):
        i, j = divmod(idx, ncols)
        ax = axes[i, j]

        theo = cell["qq_theoretical"]
        emp = cell["qq_empirical"]
        ax.plot(theo, emp, "o", markersize=3, color="#2c7bb6")

        lim = max(max(theo), max(emp)) * 1.05
        ax.plot([0, lim], [0, lim], "k--", linewidth=0.8, alpha=0.5)

        ax.set_xlabel(r"Theoretical $\chi^2_{" + str(cell['df']) + "}$", fontsize=9)
        ax.set_ylabel("Empirical quantiles", fontsize=9)
        ax.set_title(
            f"K={cell['K']}, d={cell['d']} (df={cell['df']})\n"
            f"Type I: {cell['empirical_type1_vector']:.3f} "
            f"(nominal: {cell['nominal_alpha']})",
            fontsize=9,
        )
        ax.set_aspect("equal")

    for idx in range(len(null), nrows * ncols):
        i, j = divmod(idx, ncols)
        axes[i, j].set_visible(False)

    fig.suptitle(r"Null calibration: empirical $Q_V$ vs $\chi^2$ quantiles", fontsize=12)
    fig.tight_layout()

    output_dir.mkdir(parents=True, exist_ok=True)
    path = output_dir / "null_calibration_qq.pdf"
    fig.savefig(path, bbox_inches="tight", dpi=300)
    path_png = output_dir / "null_calibration_qq.png"
    fig.savefig(path_png, bbox_inches="tight", dpi=300)
    plt.close(fig)
    print(f"Saved null calibration to {path}")
    return path


def plot_power_heatmap(results, output_dir):
    """Heatmap: power advantage (vector - scalar) across theta x rho, one panel per (K, d)."""
    power = results["power_results"]

    d_values = sorted(set(r["d"] for r in power))
    K_values = sorted(set(r["K"] for r in power))
    theta_values = sorted(set(r["theta"] for r in power))
    rho_values = sorted(set(r["rho"] for r in power))

    fig, axes = plt.subplots(
        len(d_values), len(K_values),
        figsize=(3.5 * len(K_values), 3 * len(d_values)),
    )
    if len(d_values) == 1:
        axes = axes[np.newaxis, :]
    if len(K_values) == 1:
        axes = axes[:, np.newaxis]

    lookup = {(r["theta"], r["rho"], r["K"], r["d"]): r for r in power}

    for i, d in enumerate(d_values):
        for j, K in enumerate(K_values):
            ax = axes[i, j]
            mat = np.zeros((len(rho_values), len(theta_values)))
            for ri, rho in enumerate(rho_values):
                for ti, theta in enumerate(theta_values):
                    r = lookup.get((theta, rho, K, d))
                    if r:
                        mat[ri, ti] = r["vector_power"] - r["scalar_power"]

            im = ax.imshow(mat, aspect="auto", cmap="RdYlBu_r",
                           vmin=-0.1, vmax=0.8, origin="lower")
            ax.set_xticks(range(len(theta_values)))
            ax.set_xticklabels([str(t) for t in theta_values], fontsize=7)
            ax.set_yticks(range(len(rho_values)))
            ax.set_yticklabels([str(r) for r in rho_values], fontsize=8)

            if i == 0:
                ax.set_title(f"$K = {K}$", fontsize=10)
            if j == 0:
                ax.set_ylabel(f"$d = {d}$\n" + r"$\rho$", fontsize=10)
            if i == len(d_values) - 1:
                ax.set_xlabel(r"$\theta$ (degrees)", fontsize=9)

    fig.colorbar(im, ax=axes, label=r"Power advantage (vector $-$ scalar)",
                 shrink=0.6, pad=0.02)
    fig.suptitle("Power advantage of vector $Q_V$ over per-trait scalar $Q$",
                 fontsize=12, y=1.02)
    fig.tight_layout()

    output_dir.mkdir(parents=True, exist_ok=True)
    path = output_dir / "power_advantage_heatmap.pdf"
    fig.savefig(path, bbox_inches="tight", dpi=300)
    path_png = output_dir / "power_advantage_heatmap.png"
    fig.savefig(path_png, bbox_inches="tight", dpi=300)
    plt.close(fig)
    print(f"Saved heatmap to {path}")
    return path


def main():
    path = None
    if len(sys.argv) > 1:
        path = Path(sys.argv[1])

    results = load_results(path)
    print(f"Loaded: {results['metadata']['n_cells']} power cells, "
          f"{len(results['null_calibration'])} null cells")

    plot_power_curves(results, FIGURES_DIR)
    plot_null_calibration(results, FIGURES_DIR)
    plot_power_heatmap(results, FIGURES_DIR)


if __name__ == "__main__":
    main()
