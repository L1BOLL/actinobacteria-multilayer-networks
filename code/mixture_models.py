#!/usr/bin/env python3
from __future__ import annotations

import os
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/tmp/mpl-cache")
Path(os.environ["MPLCONFIGDIR"]).mkdir(parents=True, exist_ok=True)

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from diptest import diptest
from sklearn.decomposition import PCA
from sklearn.metrics import silhouette_score
from sklearn.mixture import GaussianMixture
from sklearn.model_selection import KFold
from sklearn.preprocessing import StandardScaler

from data_io import save_figure, CATEGORY_MAP, FIG_DIR, REPORT_PATH, SEED, TABLE_DIR, compute_out_degree_matrix, ensure_output_dirs, load_or_recompute_embedding4, load_tensor
from report_utils import dataframe_to_md, replace_section


COV_TYPES = ("full", "diag", "spherical")
K_RANGE = list(range(1, 9))
N_INIT = 50
CV_FOLDS = 10


def gmm_icl(gmm: GaussianMixture, X: np.ndarray) -> float:
    """ICL = BIC + 2 × posterior entropy. Resists full-cov overfit at small N."""
    bic = float(gmm.bic(X))
    if gmm.n_components == 1:
        return bic
    resp = gmm.predict_proba(X)
    eps = 1e-12
    entropy_term = -2.0 * float(np.sum(resp * np.log(resp + eps)))
    return bic + entropy_term


def cv_loglik(X: np.ndarray, k: int, cov: str, seed: int, n_splits: int = CV_FOLDS) -> float:
    """Mean held-out log-lik across n_splits folds. NaN if k > train size or fit fails."""
    if k == 1:
        gmm = GaussianMixture(n_components=1, covariance_type=cov, random_state=seed, n_init=N_INIT)
        try:
            gmm.fit(X)
            return float(gmm.score(X))
        except Exception:
            return float("nan")
    kf = KFold(n_splits=n_splits, shuffle=True, random_state=seed)
    scores = []
    for train_idx, test_idx in kf.split(X):
        Xtr, Xte = X[train_idx], X[test_idx]
        if Xtr.shape[0] < k:
            return float("nan")
        try:
            gmm = GaussianMixture(n_components=k, covariance_type=cov, random_state=seed, n_init=N_INIT, reg_covar=1e-4)
            gmm.fit(Xtr)
            scores.append(gmm.score(Xte))
        except Exception:
            return float("nan")
    return float(np.mean(scores))


def gmm_silhouette(X: np.ndarray, k: int, cov: str, seed: int) -> float:
    if k < 2:
        return float("nan")
    try:
        gmm = GaussianMixture(n_components=k, covariance_type=cov, random_state=seed, n_init=N_INIT, reg_covar=1e-4)
        labels = gmm.fit_predict(X)
        if len(set(labels)) < 2:
            return float("nan")
        return float(silhouette_score(X, labels))
    except Exception:
        return float("nan")


def main() -> None:
    ensure_output_dirs()
    tensor, layer_ids, node_ids = load_tensor()
    emb4 = load_or_recompute_embedding4().loc[node_ids]
    out12 = compute_out_degree_matrix(tensor).loc[node_ids, layer_ids]

    datasets = {
        "4D_codebase_embedding": emb4,
        "12D_outdegree_profile": out12,
    }

    records = []
    dip_rows = []
    best_rows = []
    cluster_rows = []
    membership_rows = []

    fig, axes = plt.subplots(len(datasets), len(COV_TYPES), figsize=(13, 7.4), constrained_layout=True, squeeze=False)

    for di, (name, df) in enumerate(datasets.items()):
        X = StandardScaler().fit_transform(df.values)
        n_samples, n_features = X.shape

        pca = PCA(n_components=min(2, n_features), random_state=SEED)
        Z = pca.fit_transform(X)
        dip1 = diptest(Z[:, 0])
        dip2 = diptest(Z[:, 1]) if Z.shape[1] >= 2 else (float("nan"), float("nan"))
        dip_rows.extend(
            [
                {"embedding": name, "component": "PC1", "dip": dip1[0], "pvalue": dip1[1], "seed": SEED},
                {"embedding": name, "component": "PC2", "dip": dip2[0], "pvalue": dip2[1], "seed": SEED},
            ]
        )

        per_cov_best = {}
        for ci, cov in enumerate(COV_TYPES):
            bic_vals, aic_vals, icl_vals, cv_vals, sil_vals = [], [], [], [], []
            for k in K_RANGE:
                # Fit once; compute BIC/AIC/ICL independently so one failure doesn't zero the others.
                bic = aic = icl = float("nan")
                try:
                    gmm = GaussianMixture(n_components=k, covariance_type=cov, random_state=SEED, n_init=N_INIT, reg_covar=1e-4)
                    gmm.fit(X)
                except Exception:
                    gmm = None
                if gmm is not None:
                    try:
                        bic = float(gmm.bic(X))
                    except Exception:
                        bic = float("nan")
                    try:
                        aic = float(gmm.aic(X))
                    except Exception:
                        aic = float("nan")
                    try:
                        icl = gmm_icl(gmm, X)
                    except Exception:
                        icl = float("nan")
                cv = cv_loglik(X, k, cov, SEED)
                sil = gmm_silhouette(X, k, cov, SEED)
                bic_vals.append(bic)
                aic_vals.append(aic)
                icl_vals.append(icl)
                cv_vals.append(cv)
                sil_vals.append(sil)
                records.append({
                    "embedding": name, "covariance": cov, "k": k,
                    "bic": bic, "aic": aic, "icl": icl,
                    "cv_loglik": cv, "gmm_silhouette": sil,
                    "n_features": n_features, "n_samples": n_samples,
                    "seed": SEED,
                })

            bic_arr = np.array(bic_vals)
            icl_arr = np.array(icl_vals)
            cv_arr = np.array(cv_vals)
            sil_arr = np.array(sil_vals)

            k_bic = K_RANGE[int(np.nanargmin(bic_arr))] if np.any(np.isfinite(bic_arr)) else np.nan
            k_icl = K_RANGE[int(np.nanargmin(icl_arr))] if np.any(np.isfinite(icl_arr)) else np.nan
            k_cv = K_RANGE[int(np.nanargmax(cv_arr))] if np.any(np.isfinite(cv_arr)) else np.nan
            with np.errstate(invalid="ignore"):
                sil_after_k2 = sil_arr[1:]
            if np.any(np.isfinite(sil_after_k2)):
                k_sil = int(np.arange(2, max(K_RANGE) + 1)[int(np.nanargmax(sil_after_k2))])
                sil_best = float(np.nanmax(sil_after_k2))
            else:
                k_sil = np.nan
                sil_best = float("nan")

            per_cov_best[cov] = {
                "k_bic": k_bic, "k_icl": k_icl, "k_cv": k_cv, "k_sil": k_sil, "sil_best": sil_best,
            }

            best_rows.append({
                "embedding": name, "covariance": cov,
                "k_bic": k_bic, "k_icl": k_icl, "k_cv": k_cv,
                "k_gmm_silhouette": k_sil, "best_gmm_silhouette": sil_best,
                "dip_pc1_p": float(dip1[1]), "dip_pc2_p": float(dip2[1]),
            })

            ax = axes[di][ci]
            ax.plot(K_RANGE, bic_vals, marker="o", label="BIC", color="#1f77b4")
            ax.plot(K_RANGE, icl_vals, marker="D", label="ICL", color="#8e44ad")
            ax.set_title(f"{name.replace('_', ' ')}\ncov={cov}", fontsize=9.5)
            ax.set_xlabel("k")
            ax.set_ylabel("BIC / ICL")
            ax.grid(axis="y", color="#e5e7eb", linewidth=0.5)
            ax2 = ax.twinx()
            ax2.plot(K_RANGE, cv_vals, marker="s", linestyle=":", color="#d62728", label="CV log-lik")
            ax2.set_ylabel("CV log-lik")
            if di == 0 and ci == len(COV_TYPES) - 1:
                lines = ax.get_lines() + ax2.get_lines()
                ax.legend(lines, [l.get_label() for l in lines], loc="upper right", frameon=False, fontsize=8)

        # Save centroids at ICL-preferred k only.
        for cov in COV_TYPES:
            k_icl = per_cov_best[cov]["k_icl"]
            if not isinstance(k_icl, int) and not (isinstance(k_icl, np.integer)):
                continue
            if k_icl <= 1:
                continue
            try:
                gmm = GaussianMixture(n_components=int(k_icl), covariance_type=cov, random_state=SEED, n_init=N_INIT, reg_covar=1e-4)
                labels = gmm.fit_predict(X)
            except Exception:
                continue
            centers_orig = pd.DataFrame(
                StandardScaler().fit(df.values).inverse_transform(gmm.means_),
                columns=df.columns,
            )
            for cidx in range(int(k_icl)):
                row = {"embedding": name, "covariance": cov, "k_used": int(k_icl), "cluster": cidx, "n_strains": int(np.sum(labels == cidx))}
                for col in df.columns:
                    row[f"{col}_centroid"] = float(centers_orig.loc[cidx, col])
                cluster_rows.append(row)
                for strain in pd.Index(df.index)[labels == cidx]:
                    membership_rows.append({"embedding": name, "covariance": cov, "k_used": int(k_icl), "cluster": cidx, "strain": strain})

    save_figure(fig, "figure_mixture_selection", FIG_DIR)
    plt.close(fig)

    bic_df = pd.DataFrame(records)
    dip_df = pd.DataFrame(dip_rows)
    best_df = pd.DataFrame(best_rows)
    cluster_df = pd.DataFrame(cluster_rows)
    membership_df = pd.DataFrame(membership_rows)
    bic_df.to_csv(TABLE_DIR / "mixture_model_selection.csv", index=False)
    dip_df.to_csv(TABLE_DIR / "mixture_diptest.csv", index=False)
    best_df.to_csv(TABLE_DIR / "mixture_summary.csv", index=False)
    if not cluster_df.empty:
        cluster_df.to_csv(TABLE_DIR / "mixture_cluster_centroids.csv", index=False)
        membership_df.to_csv(TABLE_DIR / "mixture_cluster_membership.csv", index=False)

    diag_rows = best_df[best_df["covariance"] == "diag"]
    headline = []
    for _, r in diag_rows.iterrows():
        headline.append(
            f"- `{r['embedding']}` (diagonal cov): BIC k={r['k_bic']}, ICL k={r['k_icl']}, CV-LL k={r['k_cv']}, GMM-silhouette k={r['k_gmm_silhouette']} (s={r['best_gmm_silhouette']:.3g})"
        )
    headline_block = "\n".join(headline)

    body = f"""
Categories: DA={{IG, IC, RAC}}, IA={{IAC_RAC, IAC_RDE, IRAC}}, MM={{RP, IRP}}, MC={{CCAM, CCVM, CMP, CS}}. Seed: {SEED}.

GMM k=1..8 fit on z-scored 4D and 12D embeddings, three covariances (full, diag, spherical).
Criteria: BIC, ICL (BIC + 2 × posterior entropy), 10-fold CV log-likelihood, silhouette on GMM hard assignments. Hartigan dip on PC1/PC2.

Diagonal cov (regularized) headline:

{headline_block}

Best-model summary:

{dataframe_to_md(best_df, index=False)}

Rules:
- k=1 rejected only if BIC, ICL, CV agree.
- Full-cov BIC unreliable at N=60. Check `mixture_cluster_centroids.csv` for singletons (≤3 members) → overfit.
- Dip p>0.05 = consistent with unimodality on top axes (necessary, not sufficient).

Dip test:

{dataframe_to_md(dip_df, index=False)}

Outputs:
- mixture_model_selection.csv (k × cov: BIC, ICL, CV, silhouette)
- mixture_diptest.csv
- mixture_summary.csv
- mixture_cluster_centroids.csv (ICL-preferred k only)
- mixture_cluster_membership.csv
- figure_mixture_selection.png
"""
    replace_section(REPORT_PATH, "<!-- mixture_models RESULTS -->", "Discrete groups versus a continuum", body)


if __name__ == "__main__":
    main()
