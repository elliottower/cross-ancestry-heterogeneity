"""Vector-valued sheaf Q computation.

Pure-Python implementation of the multivariate sheaf Q statistic for
multi-trait multi-ancestry transportability testing.

Usage:
    from vector_q import scalar_cochran_q, vector_sheaf_q, rotation_angle
"""

import numpy as np
from scipy import stats


def scalar_cochran_q(betas, ses):
    """Standard Cochran's Q for scalar effect estimates.

    Parameters
    ----------
    betas : array (K,)
        Per-stratum effect estimates.
    ses : array (K,)
        Per-stratum standard errors.

    Returns
    -------
    dict with keys: Q, df, p, I2, tau2, beta_pooled
    """
    w = 1.0 / ses**2
    beta_bar = np.sum(w * betas) / np.sum(w)
    Q = np.sum(w * (betas - beta_bar) ** 2)
    df = len(betas) - 1
    p = 1.0 - stats.chi2.cdf(Q, df)
    I2 = max(0, (Q - df) / Q) if Q > 0 else 0.0
    c = np.sum(w) - np.sum(w**2) / np.sum(w)
    tau2 = max(0, (Q - df) / c) if c > 0 else 0.0
    return {"Q": Q, "df": df, "p": p, "I2": I2, "tau2": tau2, "beta_pooled": beta_bar}


def vector_sheaf_q(beta_matrix, cov_matrices):
    """Multivariate sheaf Q for vector-valued stalks.

    Parameters
    ----------
    beta_matrix : array (K, d)
        Per-stratum effect vectors. K strata, d traits.
    cov_matrices : array (K, d, d)
        Per-stratum covariance matrices.

    Returns
    -------
    dict with keys: Q_V, df, p, d, K
    """
    K, d = beta_matrix.shape
    assert cov_matrices.shape == (K, d, d)

    w_matrices = np.array([np.linalg.inv(cov_matrices[k]) for k in range(K)])
    W_sum = np.sum(w_matrices, axis=0)
    W_sum_inv = np.linalg.inv(W_sum)

    beta_bar = W_sum_inv @ np.sum(
        [w_matrices[k] @ beta_matrix[k] for k in range(K)], axis=0
    )

    Q_V = 0.0
    for k in range(K):
        diff = beta_matrix[k] - beta_bar
        Q_V += diff @ w_matrices[k] @ diff

    df = d * (K - 1)
    p = 1.0 - stats.chi2.cdf(Q_V, df)

    return {"Q_V": Q_V, "df": df, "p": p, "d": d, "K": K}


def vector_sheaf_q_diagonal(beta_matrix, se_matrix):
    """Multivariate sheaf Q assuming diagonal per-stratum covariance.

    Simpler interface when cross-trait correlations within each stratum
    are unavailable. Each stratum's covariance is diag(se²).

    Parameters
    ----------
    beta_matrix : array (K, d)
        Per-stratum effect vectors.
    se_matrix : array (K, d)
        Per-stratum standard errors.

    Returns
    -------
    dict with keys: Q_V, df, p, d, K, sum_scalar_Q
    """
    K, d = beta_matrix.shape
    cov_matrices = np.array([np.diag(se_matrix[k] ** 2) for k in range(K)])

    result = vector_sheaf_q(beta_matrix, cov_matrices)

    sum_scalar = sum(
        scalar_cochran_q(beta_matrix[:, j], se_matrix[:, j])["Q"] for j in range(d)
    )
    result["sum_scalar_Q"] = sum_scalar

    return result


def rotation_angle(v1, v2):
    """Angle in degrees between two vectors.

    Returns 0-180 degrees. Returns NaN if either vector is zero.
    """
    n1, n2 = np.linalg.norm(v1), np.linalg.norm(v2)
    if n1 < 1e-12 or n2 < 1e-12:
        return np.nan
    cos_theta = np.clip(np.dot(v1, v2) / (n1 * n2), -1.0, 1.0)
    return np.degrees(np.arccos(cos_theta))


def max_rotation(beta_matrix):
    """Maximum pairwise rotation angle across strata.

    Parameters
    ----------
    beta_matrix : array (K, d)

    Returns
    -------
    dict with keys: theta_max, pair (tuple of stratum indices)
    """
    K = beta_matrix.shape[0]
    theta_max = 0.0
    pair = (0, 0)
    for i in range(K):
        for j in range(i + 1, K):
            theta = rotation_angle(beta_matrix[i], beta_matrix[j])
            if not np.isnan(theta) and theta > theta_max:
                theta_max = theta
                pair = (i, j)
    return {"theta_max": theta_max, "pair": pair}


def classify_locus(beta_matrix, se_matrix, alpha=0.05):
    """Classify a locus as concordant, rotation-only, or component-only.

    Parameters
    ----------
    beta_matrix : array (K, d)
    se_matrix : array (K, d)
    alpha : float
        Significance threshold (before any external correction).

    Returns
    -------
    dict with keys: category, vector_result, scalar_results, rotation
    """
    K, d = beta_matrix.shape

    vector_result = vector_sheaf_q_diagonal(beta_matrix, se_matrix)

    scalar_results = []
    for j in range(d):
        res = scalar_cochran_q(beta_matrix[:, j], se_matrix[:, j])
        scalar_results.append(res)

    alpha_bonf = alpha / d
    any_scalar_sig = any(r["p"] < alpha_bonf for r in scalar_results)
    vector_sig = vector_result["p"] < alpha

    if vector_sig and any_scalar_sig:
        category = "concordant_significant"
    elif not vector_sig and not any_scalar_sig:
        category = "concordant_null"
    elif vector_sig and not any_scalar_sig:
        category = "rotation_only"
    else:
        category = "component_only"

    rot = max_rotation(beta_matrix)

    return {
        "category": category,
        "vector_result": vector_result,
        "scalar_results": scalar_results,
        "rotation": rot,
    }


def format_p(p):
    if p < 1e-300:
        return "< 1e-300"
    if p < 0.001:
        return f"{p:.2e}"
    return f"{p:.4f}"


def vector_sheaf_q_correlated(beta_matrix, se_matrix, rho_matrix):
    """Multivariate sheaf Q with within-stratum trait correlations.

    In biobank data (same individuals, multiple traits), the within-stratum
    covariance of effect estimates at a locus is:
        Sigma_k[j,l] = rho_jl * SE_kj * SE_kl
    where rho_jl is the phenotypic correlation between traits j and l.

    Parameters
    ----------
    beta_matrix : array (K, d)
    se_matrix : array (K, d)
    rho_matrix : array (d, d)
        Phenotypic correlation matrix across traits (shared across strata
        or can be made per-stratum by passing different rho per call).

    Returns
    -------
    dict with Q_V, df, p, d, K, sum_scalar_Q
    """
    K, d = beta_matrix.shape
    cov_matrices = np.zeros((K, d, d))
    for k in range(K):
        for j in range(d):
            for l in range(d):
                cov_matrices[k, j, l] = rho_matrix[j, l] * se_matrix[k, j] * se_matrix[k, l]

    result = vector_sheaf_q(beta_matrix, cov_matrices)

    sum_scalar = sum(
        scalar_cochran_q(beta_matrix[:, j], se_matrix[:, j])["Q"] for j in range(d)
    )
    result["sum_scalar_Q"] = sum_scalar

    return result


if __name__ == "__main__":
    np.random.seed(42)
    print("=== Sanity check: d=1 vector Q should equal scalar Q ===")
    betas_1d = np.array([0.30, 0.25, 0.15])
    ses_1d = np.array([0.05, 0.06, 0.08])
    scalar = scalar_cochran_q(betas_1d, ses_1d)
    vec = vector_sheaf_q_diagonal(betas_1d.reshape(-1, 1), ses_1d.reshape(-1, 1))
    print(f"  Scalar Q = {scalar['Q']:.4f}, p = {format_p(scalar['p'])}")
    print(f"  Vector Q = {vec['Q_V']:.4f}, p = {format_p(vec['p'])}")
    print(f"  Match: {np.isclose(scalar['Q'], vec['Q_V'])}")

    print("\n=== Diagonal covariance: Q_V = sum(Q_scalar) ===")
    beta_rot = np.array([
        [0.30, 0.10],
        [0.10, 0.30],
        [0.20, 0.20],
    ])
    se_rot = np.array([
        [0.04, 0.04],
        [0.04, 0.04],
        [0.04, 0.04],
    ])
    vr_diag = vector_sheaf_q_diagonal(beta_rot, se_rot)
    print(f"  Q_V (diagonal cov) = {vr_diag['Q_V']:.3f}")
    print(f"  sum(Q_scalar)      = {vr_diag['sum_scalar_Q']:.3f}")
    print(f"  Equal: {np.isclose(vr_diag['Q_V'], vr_diag['sum_scalar_Q'])}")
    print("  (Expected: diagonal cov => Q_V = sum of marginals)")

    print("\n=== Correlated traits: Q_V DIFFERS from sum(Q_scalar) ===")
    print("  Same betas, but traits correlated within each ancestry (rho=0.6)")
    rho = np.array([[1.0, 0.6], [0.6, 1.0]])
    vr_corr = vector_sheaf_q_correlated(beta_rot, se_rot, rho)
    print(f"  Q_V (corr cov)     = {vr_corr['Q_V']:.3f}, p = {format_p(vr_corr['p'])}")
    print(f"  sum(Q_scalar)      = {vr_corr['sum_scalar_Q']:.3f}")
    ratio = vr_corr["Q_V"] / vr_corr["sum_scalar_Q"]
    print(f"  Ratio Q_V / sum:   {ratio:.3f}")
    rot = max_rotation(beta_rot)
    print(f"  Max rotation: {rot['theta_max']:.1f}° between strata {rot['pair']}")

    print("\n=== Key demonstration: rotation detected by vector Q, missed by scalar ===")
    print("  5 ancestries, 4 blood traits, effect vector rotates but magnitudes ~stable")
    beta_blood = np.array([
        [0.25, 0.15, 0.10, 0.05],   # EUR
        [0.22, 0.18, 0.08, 0.06],   # EAS
        [0.05, 0.28, 0.20, 0.02],   # AFR — rotated (Duffy-like)
        [0.20, 0.16, 0.12, 0.04],   # AMR
        [0.23, 0.14, 0.09, 0.06],   # SAS
    ])
    se_blood = np.full((5, 4), 0.04)
    rho_blood = np.array([
        [1.0, 0.5, 0.3, 0.1],
        [0.5, 1.0, 0.4, 0.2],
        [0.3, 0.4, 1.0, 0.15],
        [0.1, 0.2, 0.15, 1.0],
    ])
    ancestries = ["EUR", "EAS", "AFR", "AMR", "SAS"]
    traits = ["Neutrophil", "Monocyte", "Basophil", "Platelet"]

    print("  Per-trait scalar Q:")
    for j, t in enumerate(traits):
        r = scalar_cochran_q(beta_blood[:, j], se_blood[:, j])
        sig = "***" if r["p"] < 0.001 else "**" if r["p"] < 0.01 else "*" if r["p"] < 0.05 else ""
        print(f"    {t:12s}: Q={r['Q']:7.2f}, p={format_p(r['p']):>10s} {sig}")

    vr = vector_sheaf_q_correlated(beta_blood, se_blood, rho_blood)
    print(f"\n  Vector Q_V: {vr['Q_V']:.2f}, p = {format_p(vr['p'])} (df={vr['df']})")
    print(f"  sum(Q_scalar): {vr['sum_scalar_Q']:.2f}")
    print(f"  Ratio: {vr['Q_V'] / vr['sum_scalar_Q']:.3f}")

    rot = max_rotation(beta_blood)
    print(f"\n  Max rotation: {rot['theta_max']:.1f}° ({ancestries[rot['pair'][0]]}-{ancestries[rot['pair'][1]]})")
    for i in range(5):
        for j2 in range(i + 1, 5):
            theta = rotation_angle(beta_blood[i], beta_blood[j2])
            print(f"    {ancestries[i]:3s}-{ancestries[j2]:3s}: {theta:.1f}°")
