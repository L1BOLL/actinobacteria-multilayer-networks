#!/usr/bin/env python3
"""Recompute every number the manuscript reports and compare it with the text.

Each entry names the claim, the value the manuscript gives, and a callable that
recomputes it from the regenerated result tables. Run it after run_all.py to
confirm that the text and the outputs still agree.

Quantities that are stated in the manuscript but not derived here - descriptions
of the assay itself, for instance - are listed separately at the end.

Exit status is non-zero if any check fails.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[1]
T = REPO / "results" / "tables"
D = REPO / "data"


def load(name: str, **kw) -> pd.DataFrame:
    return pd.read_csv(T / name, **kw)


CHECKS: list[tuple[str, object, object, float]] = []


def check(label: str, stated, actual, tol: float = 0.0) -> None:
    CHECKS.append((label, stated, actual, tol))


def deviance_split() -> dict:
    """Share of each layer's deviance explained by sender identity alone, and by
    receiver identity alone. Quoted in Methods, Results and Discussion for RAC."""
    import sys as _s
    _s.path.insert(0, str(Path(__file__).resolve().parent))
    import numpy as np
    from data_io import LAYER_FILES, load_layer_df
    from source_paths import find_layer_file

    out = {}
    for layer in ("RAC", "IG"):
        A = (load_layer_df(find_layer_file(LAYER_FILES[layer])).values > 0).astype(int)
        n = A.shape[0]
        off = ~np.eye(n, dtype=bool)
        p = A[off].mean()
        ll0 = (A[off] * np.log(p) + (1 - A[off]) * np.log(1 - p)).sum()
        def marginal(by_row: bool) -> float:
            tot = 0.0
            for i in range(n):
                x = A[i][off[i]] if by_row else A[:, i][off[:, i]]
                q = float(np.clip(x.mean(), 1e-9, 1 - 1e-9))
                tot += (x * np.log(q) + (1 - x) * np.log(1 - q)).sum()
            return (1 - tot / ll0) * 100
        # rows are the strain the phenotype is scored in, columns the partner
        out[layer] = (round(marginal(True), 1), round(marginal(False), 1))
    return out


def main() -> None:
    # ---- §1 / Table 1: layer architecture -------------------------------- #
    p7 = load("layer_descriptors.csv").set_index("layer")
    check("RAC density = 0.474", 0.474, p7.loc["RAC", "edge_density"], 5e-4)
    check("RAC R_p = 0.297", 0.297, p7.loc["RAC", "reciprocity_pairs"], 5e-4)
    check("IG density = 0.217", 0.217, p7.loc["IG", "edge_density"], 5e-4)
    check("IG R_p = 0.150", 0.150, p7.loc["IG", "reciprocity_pairs"], 5e-4)
    check("IAC_RAC density = 0.108", 0.108, p7.loc["IAC_RAC", "edge_density"], 5e-4)
    check("IRAC density = 0.091", 0.091, p7.loc["IRAC", "edge_density"], 5e-4)

    for layer, e in [("RAC", 1677), ("IG", 767), ("IAC_RAC", 383), ("IRAC", 321),
                     ("CS", 265), ("CMP", 260), ("CCVM", 252), ("RP", 145),
                     ("IC", 117), ("CCAM", 102), ("IAC_RDE", 67), ("IRP", 19)]:
        check(f"edges {layer} = {e}", e, round(p7.loc[layer, "edge_density"] * 60 * 59))

    # ---- sender versus receiver contribution, quoted in Methods and Results  #
    dev = deviance_split()
    check("RAC deviance from the scored strain = 57.6%", 57.6, dev["RAC"][0], 0.1)
    check("RAC deviance from the partner = 0.8%", 0.8, dev["RAC"][1], 0.1)
    check("IG deviance from the scored strain = 9.3%", 9.3, dev["IG"][0], 0.1)
    check("IG deviance from the partner = 18.1%", 18.1, dev["IG"][1], 0.1)

    # ---- §3 / Table 2: overlap ------------------------------------------- #
    p4 = load("tableS2_layer_overlap_full.csv")
    bio = p4[~p4["methodological_confound"]]
    sig = p4[p4["significant_by"]]
    check("mean J over 66 pairs = 0.047", 0.047, p4["J_obs"].mean(), 5e-4)
    check("mean J excluding confound = 0.043", 0.043, bio["J_obs"].mean(), 5e-4)
    check("n pairs = 66", 66, len(p4))
    check("enriched biological pairs = 8", 8, int(((sig["z"] > 0) & ~sig["methodological_confound"]).sum()))
    check("depleted biological pairs = 1", 1, int(((sig["z"] < 0) & ~sig["methodological_confound"]).sum()))
    conf = p4[p4["methodological_confound"]].iloc[0]
    check("CCAM-CCVM J = 0.336", 0.336, conf["J_obs"], 5e-4)
    check("CCAM-CCVM z = +16.19", 16.19, conf["z"], 5e-3)
    dep = sig[sig["z"] < 0].iloc[0]
    check("IAC_RAC-RAC z = -5.32", -5.32, dep["z"], 5e-3)

    # ---- §4: SCC and the new reciprocity / transitivity nulls ------------ #
    p5 = load("connectivity_null.csv").set_index("layer")
    for layer, val in [("IG", 59), ("RAC", 51), ("IAC_RAC", 46), ("IRAC", 41),
                       ("CMP", 39), ("CS", 26), ("CCVM", 24), ("CCAM", 4)]:
        check(f"largest SCC {layer} = {val}", val, int(p5.loc[layer, "scc_obs"]))
    check("all |z(SCC)| < 1.5", True, bool((p5["z"].abs().dropna() < 1.5).all()))
    check("no SCC significant after BY", 0, int(p5["significant_by"].sum()))

    p14_path = T / "reciprocity_transitivity.csv"
    if p14_path.exists():
        p14 = pd.read_csv(p14_path).set_index("layer")
        rec_sig = p14[p14["reciprocity_significant_by"]]
        tra_sig = p14[p14["transitivity_significant_by"]]
        print("\n  [reciprocity_transitivity] reciprocity significant : "
              + (", ".join(f"{i} ({r.reciprocity_z:+.2f})" for i, r in rec_sig.iterrows()) or "none"))
        print("  [reciprocity_transitivity] transitivity significant: "
              + (", ".join(f"{i} ({r.transitivity_z:+.2f})" for i, r in tra_sig.iterrows()) or "none"))
        for layer, z in [("IG", 6.17), ("IRAC", 4.24), ("CS", 3.10), ("IAC_RAC", 2.82), ("RAC", 2.75)]:
            check(f"reciprocity z {layer} = +{z}", z, p14.loc[layer, "reciprocity_z"], 0.006)
            check(f"reciprocity {layer} significant", True, bool(p14.loc[layer, "reciprocity_significant_by"]))
        for layer, z in [("CS", -4.06), ("RAC", -4.05), ("IG", -3.56), ("CMP", -3.34)]:
            check(f"transitivity z {layer} = {z}", z, p14.loc[layer, "transitivity_z"], 0.006)
            check(f"transitivity {layer} significant", True, bool(p14.loc[layer, "transitivity_significant_by"]))
        check("5 layers significantly more reciprocal", 5, int(p14["reciprocity_significant_by"].sum()))
        check("no layer significantly less reciprocal", 0,
              int((p14["reciprocity_significant_by"] & (p14["reciprocity_z"] < 0)).sum()))
        check("4 layers significantly less transitive", 4, int(p14["transitivity_significant_by"].sum()))
        check("no layer significantly more transitive", 0,
              int((p14["transitivity_significant_by"] & (p14["transitivity_z"] > 0)).sum()))
        check("CCVM transitivity marginal, q = 0.064", 0.064, p14.loc["CCVM", "transitivity_q_by"], 0.0006)
        check("IRP reciprocity null degenerate", False, bool(np.isfinite(p14.loc["IRP", "reciprocity_z"])))
        check("IAC_RDE reciprocity null degenerate", False,
              bool(np.isfinite(p14.loc["IAC_RDE", "reciprocity_z"])))
    else:
        print("\n  [reciprocity_transitivity] NOT YET GENERATED - run run_all.py")

    # ---- §5: strategy space ---------------------------------------------- #
    v12 = load("ordination_variance.csv")
    check("12-D PC1 = 32.5%", 32.5, v12.loc[0, "variance_explained"] * 100, 0.05)
    check("12-D PC2 = 17.3%", 17.3, v12.loc[1, "variance_explained"] * 100, 0.05)
    check("12-D PC1+PC2 = 49.8%", 49.8, v12.loc[1, "cumulative_variance"] * 100, 0.05)
    check("12-D PC3 = 11.8%", 11.8, v12.loc[2, "variance_explained"] * 100, 0.05)
    check("12-D PC4 = 8.2%", 8.2, v12.loc[3, "variance_explained"] * 100, 0.05)
    check("null 95th pct PC1 = 17.0%", 17.0, v12.loc[0, "parallel_null_p95"] * 100, 0.15)
    check("null 95th pct PC2 = 14.4%", 14.4, v12.loc[1, "parallel_null_p95"] * 100, 0.15)
    check("null mean PC3 = 11.6%", 11.6, v12.loc[2, "parallel_null_mean"] * 100, 0.15)
    check("null 95th pct PC3 = 12.6%", 12.6, v12.loc[2, "parallel_null_p95"] * 100, 0.15)
    check("PC3 bootstrap interval 9.8-15.0%",
          (9.8, 15.0), (round(v12.loc[2, "var_ci_lo"] * 100, 1), round(v12.loc[2, "var_ci_hi"] * 100, 1)))
    check("2 components clear the 95th pct", 2, int(v12["significant_vs_parallel"].sum()))
    check("3 components clear the mean", 3, int(v12["significant_vs_parallel_mean_rule"].sum()))
    check("PC4 below the null on both rules", False, bool(v12.loc[3, "significant_vs_parallel"]))

    v4 = load("ordination_variance_4d.csv")
    check("4-D PC1 = 53.0%", 53.0, v4.loc[0, "variance_explained"] * 100, 0.05)
    check("4-D PC1+PC2 = 76.8%", 76.8, v4.loc[1, "cumulative_variance"] * 100, 0.05)

    ld = load("ordination_loadings.csv").set_index("layer")
    for layer, val in [("CMP", 0.41), ("IRP", 0.40), ("CCVM", 0.38), ("CCAM", 0.37), ("IG", 0.33)]:
        check(f"PC1 loading {layer} = {val}", val, ld.loc[layer, "PC1_loading"], 0.005)
    check("PC2 loading RAC = +0.54", 0.54, ld.loc["RAC", "PC2_loading"], 0.005)
    check("PC3 loading IAC_RDE = +0.65", 0.65, ld.loc["IAC_RDE", "PC3_loading"], 0.005)

    sp = load("ordination_pc1_concordance.csv")
    check("12-D/4-D PC1 rho = 0.934", 0.934, sp.loc[0, "rho_pc1_rank"], 5e-4)
    check("that p ~ 1.5e-27", 1.5e-27, sp.loc[0, "pvalue"], 1e-28)

    dip = load("mixture_diptest.csv")
    d4 = dip[(dip.embedding == "4D_codebase_embedding") & (dip.component == "PC1")]["pvalue"].iloc[0]
    d12 = dip[(dip.embedding == "12D_outdegree_profile") & (dip.component == "PC1")]["pvalue"].iloc[0]
    check("dip P (4-D PC1) = 0.53", 0.53, d4, 0.005)
    check("dip P (12-D PC1) = 0.97", 0.97, d12, 0.005)

    # Model order is checked for every covariance parameterization, not just full.
    # Reporting only the full-covariance k = 1 was correction #20; these checks
    # exist so that overstatement cannot come back.
    gmm = load("mixture_model_selection.csv")

    def cvk(emb: str, cov: str) -> int:
        s = gmm[(gmm.embedding == emb) & (gmm.covariance == cov)]
        return int(s.loc[s["cv_loglik"].idxmax(), "k"])

    g12 = gmm[(gmm.embedding == "12D_outdegree_profile") & (gmm.covariance == "full")]
    check("BIC selects k = 5 (12-D, full)", 5, int(g12.loc[g12["bic"].idxmin(), "k"]))
    check("CV k = 1 (12-D, full)", 1, cvk("12D_outdegree_profile", "full"))
    check("CV k = 3 (12-D, diag)", 3, cvk("12D_outdegree_profile", "diag"))
    check("CV k = 2 (12-D, spherical)", 2, cvk("12D_outdegree_profile", "spherical"))
    check("CV k = 1 (4-D, full)", 1, cvk("4D_codebase_embedding", "full"))
    check("CV k = 2 (4-D, diag)", 2, cvk("4D_codebase_embedding", "diag"))
    check("CV k = 2 (4-D, spherical)", 2, cvk("4D_codebase_embedding", "spherical"))
    # The overall CV optimum is NOT k = 1; the paper must not imply otherwise.
    d12 = gmm[gmm.embedding == "12D_outdegree_profile"]
    best12 = d12.loc[d12["cv_loglik"].idxmax()]
    check("overall CV best (12-D) is diag k = 3", ("diag", 3), (best12["covariance"], int(best12["k"])))
    sil = d12[(d12.covariance != "full") & (d12.k.isin([2, 3]))]["gmm_silhouette"]
    check("constrained silhouettes within 0.15-0.44", True,
          bool(sil.min() >= 0.14 and sil.max() <= 0.45))

    ident = pd.read_csv(D / "ecological_identity_with_PCs_and_clusters.csv")
    part = ident["participation"]
    check("participation min = 0.60", 0.60, part.min(), 0.005)
    check("participation max = 0.84", 0.84, part.max(), 0.005)
    od = ident[[c for c in ident.columns if c.startswith("outdeg_")]].sum(axis=1)
    check("min total out-degree = 47", 47, int(od.min()))
    check("max total out-degree = 130", 130, int(od.max()))

    # ---- §6: dominance ---------------------------------------------------- #
    c = load("dominance_cross_layer_correlation.csv", index_col=0)
    v = c.values[np.triu_indices_from(c.values, 1)]
    check("David's cross-layer mean = 0.063", 0.063, v.mean(), 5e-4)
    check("David's cross-layer median = 0.043", 0.043, float(np.median(v)), 5e-4)
    check("David's min = -0.66", -0.66, v.min(), 5e-3)
    check("David's max = +0.62", 0.62, v.max(), 5e-3)
    check("pairs with |rho| > 0.5 = 7", 7, int((np.abs(v) > 0.5).sum()))
    check("CCAM-CCVM rho = 0.62", 0.62, c.loc["CCAM", "CCVM"], 5e-3)
    check("IRAC-RAC rho = 0.59", 0.59, c.loc["IRAC", "RAC"], 5e-3)
    check("IG-IC rho = 0.56", 0.56, c.loc["IG", "IC"], 5e-3)
    check("IG-RAC rho = -0.66", -0.66, c.loc["IG", "RAC"], 5e-3)
    check("IAC_RDE-RAC rho = -0.63", -0.63, c.loc["IAC_RDE", "RAC"], 5e-3)
    check("IAC_RDE-IRAC rho = -0.58", -0.58, c.loc["IAC_RDE", "IRAC"], 5e-3)

    # ---- §7: phylogeny ---------------------------------------------------- #
    pg = load("phylogenetic_signal.csv")
    ml = pg[pg["topology"].str.startswith("ML")]
    check("21 traits on the ML tree", 21, len(ml))
    check("63 trait x topology tests", 63, len(pg))
    check("0 of 21 significant on ML tree", 0, int((ml["q_BH"] < 0.05).sum()))
    check("0 of 42 on UPGMA/WPGMA", 0, int((pg[~pg["topology"].str.startswith("ML")]["q_BH"] < 0.05).sum()))
    check("minimum q = 0.39", 0.3948, pg["q_BH"].min(), 0.005)

    vcv = load("phylogenetic_vcv_diagnostics.csv")
    check("60 tips on every topology", True, bool((vcv["n_tips"] == 60).all()))
    check("ML VCV rank 44 of 60", 44, int(vcv[vcv.topology.str.startswith("ML")]["vcv_rank"].iloc[0]))
    check("UPGMA/WPGMA rank 38 (= genotype count)", True,
          bool((vcv[~vcv.topology.str.startswith("ML")]["vcv_rank"] == 38).all()))

    mt = load("mantel.csv")
    check("Mantel rho = -0.041", -0.041, mt.loc[0, "mantel_rho"], 5e-4)
    check("Mantel P = 0.59", 0.5925, mt.loc[0, "p_value"], 0.005)
    per = mt.iloc[1:]
    check("per-category max |rho| <= 0.069", True, bool(per["mantel_rho"].abs().max() <= 0.069))
    check("per-category min P > 0.22", True, bool(per["p_value"].min() > 0.22))

    # ---- §8: metabolism and taxon sensitivity ---------------------------- #
    prof = load("metabolic_niche_profiles.csv")
    check("assimilation capacity mean = 12.1", 12.1, prof["assimilation_capacity"].mean(), 0.05)
    check("54 of 60 strains use >= 7 sugars", 54, int((prof["n_sugars_used"] >= 7).sum()))
    ms = load("metabolic_niche_strain.csv")
    check("18 strain-level tests", 18, len(ms))
    check("0 strain-level tests at q < 0.05", 0, int((ms["q_BH"] < 0.05).sum()))
    md = load("metabolic_niche_dyad.csv").set_index("layer")
    check("dyad IC rho = +0.121", 0.121, md.loc["IC", "mantel_rho"], 5e-4)
    check("dyad IG rho = +0.131", 0.131, md.loc["IG", "mantel_rho"], 5e-4)
    check("dyad IC q = 0.085", 0.085, md.loc["IC", "q_BH"], 5e-4)

    tax = load("tableS4_taxon_sensitivity.csv").set_index("quantity")

    def taxrow(sub: str, col: str):
        hits = [i for i in tax.index if sub in i]
        assert len(hits) == 1, f"ambiguous taxon row for {sub!r}: {hits}"
        return tax.loc[hits[0], col]

    check("n56 density rank correlation = 1.000", 1.0, float(taxrow("Layer density rank", "n56")), 5e-4)
    check("n56 largest-SCC rank correlation = 0.998", 0.998, float(taxrow("Largest-SCC rank", "n56")), 5e-4)
    check("PCs above the 95th-pct null 2 -> 2", ("2", "2"),
          (str(taxrow("PCs above Horn", "n60")), str(taxrow("PCs above Horn", "n56"))))
    check("CV mixture order 1 -> 1", ("1", "1"),
          (str(taxrow("GMM components", "n60")), str(taxrow("GMM components", "n56"))))
    check("draw-matched enriched pairs 8 -> 6", ("8", "6"),
          (str(taxrow("over-overlapping pairs, matched", "n60")),
           str(taxrow("over-overlapping pairs, matched", "n56"))))
    check("draw-matched depleted pairs 1 -> 1", ("1", "1"),
          (str(taxrow("under-overlapping pairs, matched", "n60")),
           str(taxrow("under-overlapping pairs, matched", "n56"))))

    # ---- §9 / Figure 6: spatial configuration ---------------------------- #
    sp = load("spatial_configuration.csv").set_index("layer")
    check("IG close-proximity edges = 785", 785, int(sp.loc["IG", "edges_close"]))
    check("IG separated edges = 261", 261, int(sp.loc["IG", "edges_separated"]))
    check("IG density 0.222 -> 0.074", (0.222, 0.074),
          (round(float(sp.loc["IG", "density_close"]), 3), round(float(sp.loc["IG", "density_separated"]), 3)))
    check("IG delta density = 0.148", 0.148, sp.loc["IG", "density_difference"], 5e-4)
    check("IG delta 95% CI 0.136-0.160", (0.136, 0.160),
          (round(float(sp.loc["IG", "ci_low"]), 3), round(float(sp.loc["IG", "ci_high"]), 3)))
    check("IG persistent/lost/gained 249/536/12", (249, 536, 12),
          tuple(int(sp.loc["IG", k]) for k in ("persistent", "lost", "gained")))
    check("IG loss:gain OR = 44.67", 44.67, sp.loc["IG", "odds_loss_vs_gain"], 5e-3)
    check("IG dyad-randomisation P < 1e-4", True, bool(sp.loc["IG", "p_paired_dyad_permutation"] <= 1e-4))
    check("CS close-proximity edges = 144", 144, int(sp.loc["CS", "edges_close"]))
    check("CS separated edges = 133", 133, int(sp.loc["CS", "edges_separated"]))
    check("CS delta density = 0.0031", 0.0031, sp.loc["CS", "density_difference"], 5e-5)
    check("CS persistent/lost/gained 59/85/74", (59, 85, 74),
          tuple(int(sp.loc["CS", k]) for k in ("persistent", "lost", "gained")))
    check("CS loss:gain OR = 1.15", 1.15, sp.loc["CS", "odds_loss_vs_gain"], 5e-3)
    check("CS dyad-randomisation P = 0.433", 0.433, sp.loc["CS", "p_paired_dyad_permutation"], 5e-3)
    check("spatial IG differs from main IG in 650 entries", 650, int(sp.loc["IG", "entries_differing_from_main_layer"]))
    check("spatial CS differs from main CS in 291 entries", 291, int(sp.loc["CS", "entries_differing_from_main_layer"]))

    # ---- §6: what the IG-RAC dominance correlation means ------------------ #
    from scipy.stats import spearmanr
    from data_io import load_tensor
    tensor, layer_ids, node_ids = load_tensor()
    ds = load("dominance_davids_scores.csv").set_index(load("dominance_davids_scores.csv").columns[0])
    ds = ds.loc[node_ids]
    rac_release = tensor[layer_ids.index("RAC")].sum(axis=1)   # rows = strain scored = releasing strain
    ig_out = tensor[layer_ids.index("IG")].sum(axis=0)         # columns = sender = inhibitor
    check("DS IG vs DS RAC rho = -0.66", -0.66, spearmanr(ds["IG"], ds["RAC"]).statistic, 5e-3)
    check("DS RAC vs RAC release count rho = -0.97", -0.97, spearmanr(ds["RAC"], rac_release).statistic, 5e-3)
    check("IG out-degree vs RAC release count rho = +0.75", 0.75, spearmanr(ig_out, rac_release).statistic, 5e-3)
    check("n56 PC1 = 31.2%", "31.2", str(taxrow("PC1 (% variance)", "n56")))
    check("n56 PC2 = 15.2%", "15.2", str(taxrow("PC2 (% variance)", "n56")))

    # ---- report ----------------------------------------------------------- #
    fails = []
    print(f"\n{'claim':58s} {'paper':>12s} {'computed':>12s}   ")
    print("-" * 92)
    for label, stated, actual, tol in CHECKS:
        if isinstance(stated, (int, float)) and not isinstance(stated, bool) and tol:
            ok = abs(float(actual) - float(stated)) <= tol
            shown = f"{float(actual):.6g}"
        else:
            ok = actual == stated
            shown = str(actual)
        if not ok:
            fails.append(label)
        print(f"{label:58s} {str(stated):>12s} {shown:>12s}   {'ok' if ok else 'FAIL'}")

    print("-" * 92)
    print(f"{len(CHECKS) - len(fails)} of {len(CHECKS)} checks passed")
    if fails:
        print("\nFAILED:")
        for f in fails:
            print("  -", f)

    print("\nStated in the manuscript but not derived from these data:")
    for line in ["Figure 1: what the plate photographs show",
                 "Table 2 'phenotypic interpretation' column: authors' wording",
                 "reference list"]:
        print("  -", line)

    sys.exit(1 if fails else 0)


if __name__ == "__main__":
    main()
