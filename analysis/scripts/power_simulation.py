"""Power simulation for vector Q vs per-trait scalar Q.

Pure-Python simulation logic, separate from Modal wrapper.
Generates synthetic multi-ancestry multi-trait data with known rotation,
computes both vector Q and per-trait scalar Q, and reports rejection rates.

Also includes null calibration: empirical type-I error check for the
chi-squared approximation at various K, d.

Usage (local, small grid):
    uv run python experiments/batch6_vector_extension/scripts/power_simulation.py --local

Usage (full grid, called from Modal):
    See modal_power_simulation.py
"""

import json
import sys
from datetime import datetime
from pathlib import Path

import numpy as np
from scipy import stats
from tqdm import tqdm


def generate_rotated_betas(K, d, theta_deg, base_effect_size=0.2):
    """Generate (K, d) beta matrix where strata k>0 are rotated by theta from stratum 0.

    Rotation is applied in the (0,1) plane of the d-dimensional space.
    Stratum 0 has a fixed effect vector; strata 1..K-1 are each rotated
    by theta degrees (same direction) to create a clean rotation signal.

    Parameters
    ----------
    K : int
        Number of strata (ancestries).
    d : int
        Number of traits.
    theta_deg : float
        Rotation angle in degrees.
    base_effect_size : float
        Magnitude of the base effect vector.

    Returns
    -------
    beta_true : array (K, d)
    """
    beta_0 = np.zeros(d)
    beta_0[0] = base_effect_size
    if d > 1:
        beta_0[1] = base_effect_size * 0.5

    theta_rad = np.radians(theta_deg)
    cos_t, sin_t = np.cos(theta_rad), np.sin(theta_rad)

    betas = np.zeros((K, d))
    betas[0] = beta_0

    for k in range(1, K):
        b = beta_0.copy()
        b[0] = cos_t * beta_0[0] - sin_t * beta_0[1]
        b[1] = sin_t * beta_0[0] + cos_t * beta_0[1]
        betas[k] = b

    return betas


def generate_correlation_matrix(d, rho):
    """Generate a d x d compound-symmetric correlation matrix.

    All off-diagonal entries are rho. Eigenvalues are (1-rho) with
    multiplicity d-1 and (1 + (d-1)*rho) with multiplicity 1,
    so the matrix is PD when rho > -1/(d-1).
    """
    return (1 - rho) * np.eye(d) + rho * np.ones((d, d))


def simulate_one(K, d, theta_deg, rho, se_scale=0.06, alpha=0.05):
    """Run one simulation replicate.

    Returns
    -------
    dict with keys:
        vector_reject : bool
        any_scalar_reject : bool (Bonferroni-corrected)
        all_scalar_reject : bool
        n_scalar_reject : int
    """
    beta_true = generate_rotated_betas(K, d, theta_deg)

    se_matrix = np.full((K, d), se_scale)

    rho_matrix = generate_correlation_matrix(d, rho)

    beta_obs = np.zeros((K, d))
    for k in range(K):
        cov_k = np.diag(se_matrix[k]) @ rho_matrix @ np.diag(se_matrix[k])
        beta_obs[k] = np.random.multivariate_normal(beta_true[k], cov_k)

    cov_matrices = np.zeros((K, d, d))
    for k in range(K):
        cov_matrices[k] = np.diag(se_matrix[k]) @ rho_matrix @ np.diag(se_matrix[k])

    W = np.array([np.linalg.inv(cov_matrices[k]) for k in range(K)])
    W_sum = np.sum(W, axis=0)
    W_sum_inv = np.linalg.inv(W_sum)
    beta_bar = W_sum_inv @ np.sum([W[k] @ beta_obs[k] for k in range(K)], axis=0)

    Q_V = 0.0
    for k in range(K):
        diff = beta_obs[k] - beta_bar
        Q_V += diff @ W[k] @ diff

    df_vec = d * (K - 1)
    p_vec = 1.0 - stats.chi2.cdf(Q_V, df_vec)
    vector_reject = p_vec < alpha

    alpha_bonf = alpha / d
    n_scalar_reject = 0
    for j in range(d):
        w_j = 1.0 / se_matrix[:, j] ** 2
        beta_bar_j = np.sum(w_j * beta_obs[:, j]) / np.sum(w_j)
        Q_j = np.sum(w_j * (beta_obs[:, j] - beta_bar_j) ** 2)
        p_j = 1.0 - stats.chi2.cdf(Q_j, K - 1)
        if p_j < alpha_bonf:
            n_scalar_reject += 1

    return {
        "vector_reject": bool(vector_reject),
        "any_scalar_reject": n_scalar_reject > 0,
        "n_scalar_reject": n_scalar_reject,
    }


def run_power_grid_chunk(chunk, n_sims=1000, alpha=0.05):
    """Run simulations for a chunk of the parameter grid.

    Parameters
    ----------
    chunk : list of dict
        Each dict has keys: theta, rho, K, d
    n_sims : int
    alpha : float

    Returns
    -------
    list of dict with results per grid cell
    """
    results = []
    for cell in tqdm(chunk, desc="Grid cells"):
        theta = cell["theta"]
        rho = cell["rho"]
        K = cell["K"]
        d = cell["d"]

        vec_rejects = 0
        scalar_rejects = 0
        for _ in range(n_sims):
            r = simulate_one(K, d, theta, rho, alpha=alpha)
            vec_rejects += r["vector_reject"]
            scalar_rejects += r["any_scalar_reject"]

        results.append({
            "theta": theta,
            "rho": rho,
            "K": K,
            "d": d,
            "n_sims": n_sims,
            "vector_power": vec_rejects / n_sims,
            "scalar_power": scalar_rejects / n_sims,
        })

    return results


def run_null_calibration(K_values, d_values, n_sims=10000, alpha=0.05):
    """Check empirical type-I error under the null (theta=0, no heterogeneity).

    Returns
    -------
    list of dict with empirical rejection rates and QQ data
    """
    results = []
    for K in K_values:
        for d in d_values:
            print(f"  Null calibration: K={K}, d={d}, n_sims={n_sims}")
            Q_values = []
            vec_rejects = 0
            scalar_rejects = 0

            for _ in tqdm(range(n_sims), desc=f"Null K={K} d={d}", leave=False):
                se_matrix = np.full((K, d), 0.06)
                beta_obs = np.random.normal(0.2, 0.06, size=(K, d))

                cov_matrices = np.array([np.diag(se_matrix[k] ** 2) for k in range(K)])
                W = np.array([np.linalg.inv(cov_matrices[k]) for k in range(K)])
                W_sum = np.sum(W, axis=0)
                W_sum_inv = np.linalg.inv(W_sum)
                beta_bar = W_sum_inv @ np.sum([W[k] @ beta_obs[k] for k in range(K)], axis=0)

                Q_V = 0.0
                for k in range(K):
                    diff = beta_obs[k] - beta_bar
                    Q_V += diff @ W[k] @ diff

                Q_values.append(Q_V)
                df_vec = d * (K - 1)
                p_vec = 1.0 - stats.chi2.cdf(Q_V, df_vec)
                if p_vec < alpha:
                    vec_rejects += 1

                alpha_bonf = alpha / d
                any_rej = False
                for j in range(d):
                    w_j = 1.0 / se_matrix[:, j] ** 2
                    beta_bar_j = np.sum(w_j * beta_obs[:, j]) / np.sum(w_j)
                    Q_j = np.sum(w_j * (beta_obs[:, j] - beta_bar_j) ** 2)
                    p_j = 1.0 - stats.chi2.cdf(Q_j, K - 1)
                    if p_j < alpha_bonf:
                        any_rej = True
                        break
                if any_rej:
                    scalar_rejects += 1

            Q_arr = np.sort(Q_values)
            df_vec = d * (K - 1)
            theoretical = stats.chi2.ppf(np.linspace(0.01, 0.99, 50), df_vec)
            empirical = np.percentile(Q_arr, np.linspace(1, 99, 50))

            results.append({
                "K": K,
                "d": d,
                "df": df_vec,
                "n_sims": n_sims,
                "empirical_type1_vector": vec_rejects / n_sims,
                "empirical_type1_scalar": scalar_rejects / n_sims,
                "nominal_alpha": alpha,
                "qq_theoretical": theoretical.tolist(),
                "qq_empirical": empirical.tolist(),
            })

    return results


FULL_GRID = {
    "theta": [0, 5, 10, 15, 20, 30, 45, 60, 75, 90],
    "rho": [0.0, 0.3, 0.6, 0.9],
    "K": [4, 6, 8],
    "d": [5, 15, 25],
}

SMALL_GRID = {
    "theta": [0, 15, 30, 60, 90],
    "rho": [0.0, 0.6],
    "K": [4, 6],
    "d": [5, 15],
}


def build_grid(grid_spec):
    cells = []
    for theta in grid_spec["theta"]:
        for rho in grid_spec["rho"]:
            for K in grid_spec["K"]:
                for d in grid_spec["d"]:
                    cells.append({"theta": theta, "rho": rho, "K": K, "d": d})
    return cells


def main():
    print("Vector Q Power Simulation")
    print("=" * 50)
    print(f"Timestamp: {datetime.now().isoformat()}")

    is_local = "--local" in sys.argv
    grid_spec = SMALL_GRID if is_local else FULL_GRID
    n_sims = 200 if is_local else 2000

    cells = build_grid(grid_spec)
    print(f"Grid: {len(cells)} cells x {n_sims} sims = {len(cells) * n_sims:,} total")
    print(f"Parameters: {grid_spec}")
    print()

    results = run_power_grid_chunk(cells, n_sims=n_sims)

    print(f"\nRunning null calibration...")
    null_results = run_null_calibration(
        K_values=grid_spec["K"],
        d_values=grid_spec["d"],
        n_sims=5000 if is_local else 10000,
    )

    print("\n" + "=" * 60)
    print("NULL CALIBRATION (empirical type-I error at alpha=0.05)")
    print("=" * 60)
    for r in null_results:
        status = "OK" if abs(r["empirical_type1_vector"] - 0.05) < 0.015 else "CHECK"
        print(f"  K={r['K']:2d}, d={r['d']:2d} (df={r['df']:3d}): "
              f"vector={r['empirical_type1_vector']:.3f}, "
              f"scalar={r['empirical_type1_scalar']:.3f}  [{status}]")

    print("\n" + "=" * 60)
    print("POWER RESULTS (selection)")
    print("=" * 60)
    for r in sorted(results, key=lambda x: (x["d"], x["K"], x["rho"], x["theta"])):
        if r["theta"] in [0, 30, 60, 90]:
            print(f"  theta={r['theta']:3d}, rho={r['rho']:.1f}, K={r['K']}, d={r['d']:2d}: "
                  f"vector={r['vector_power']:.3f}, scalar={r['scalar_power']:.3f}")

    output_dir = Path(__file__).resolve().parent.parent / "results"
    output_dir.mkdir(exist_ok=True)

    output = {
        "metadata": {
            "timestamp": datetime.now().isoformat(),
            "grid": grid_spec,
            "n_sims": n_sims,
            "n_cells": len(cells),
            "mode": "local" if is_local else "full",
        },
        "power_results": results,
        "null_calibration": null_results,
    }

    output_path = output_dir / "power_simulation_results.json"
    with open(output_path, "w") as f:
        json.dump(output, f, indent=2)
    print(f"\nResults saved to {output_path}")


if __name__ == "__main__":
    main()
