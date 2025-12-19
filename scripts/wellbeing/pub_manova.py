"""
analysis_pipeline_lambda_pubready_perm_cv.py

Pipeline:
• Drug-aware (subject FE) preprocessing and z-scoring
• Permutation MANOVA (within-subject shuffling for validity)
• Univariate follow-ups + FDR
• Canonical analysis (weights, structure coeffs) + main scatter
• Correlation test via full-shuffle permutation
• Repeated-measures correlation (Pearson rmcorr) 
• Cross-validation (row-wise KFold, K=5) with mean r, median r (printed), and permutation p-value

All plots are shown inline and saved to pub_outputs/.
"""

from pathlib import Path
from datetime import datetime
import numpy as np
import pandas as pd
import statsmodels.api as sm
import matplotlib.pyplot as plt
from scipy.stats import zscore, pearsonr
from sklearn.model_selection import KFold
import pingouin as pg  # NEW: for rmcorr
import os
# ------------------
# Config
# ------------------
BASE = Path.home() / "Documents" / "MATLAB" / "project_cb_tom"
OUTDIR = BASE /"data"/ "behaviour" / "pub_outputs"
EXCEL_PATH = BASE / "data" / "behaviour" / "significant_behaviour.xlsx"

EXCLUDE_SUBJECTS = (150, 170)

# Permutation sizes / seeds
N_PERM_MANOVA = 10000     # MANOVA permutations (within-subject)
N_PERM_CORR   = 10000     # Correlation permutations (full shuffle)
N_PERM_CV     = 10000     # Cross-validation permutation p-value
SEED = 11

# CV config
K_ROW = 5                 # row-wise KFold

PSYCHOSOCIAL = [
    "positive_attitudes_about_life",
    "positive_attitudes_about_self",
    "positive_mood_changes",
    "negative_mood_changes",
    "altruistic_positive_social_effect",
    "positive_behavior_change",
    "panas_pos_affect_median_merge",  # scaled  before analysis
]
ACUTE = ["obe", "edi", "aed"]

# ------------------
# Helpers
# ------------------
def residualise_within_subject(df, subj_col, cols, covar):
    """Within each subject, regress out covar (and intercept) from each col; return residuals."""
    out = df.copy()
    for pid, grp in df.groupby(subj_col):
        X = np.column_stack([np.ones(len(grp)), grp[covar].to_numpy()])
        for c in cols:
            y = grp[c].to_numpy()
            if np.unique(X[:, 1]).size > 1:
                beta = np.linalg.lstsq(X, y, rcond=None)[0]
                resid = y - X @ beta
            else:
                resid = y - y.mean()
            out.loc[grp.index, c] = resid
    return out

def fdr_bh(pvals):
    pvals = np.asarray(pvals, dtype=float)
    n = len(pvals)
    order = np.argsort(pvals)
    ranks = np.arange(1, n + 1)
    q = np.empty(n)
    q[order] = pvals[order] * n / ranks
    for i in range(n - 2, -1, -1):
        q[order[i]] = min(q[order[i]], q[order[i + 1]])
    return np.minimum(q, 1.0)

def pillai_trace(Y, X):
    X0 = np.ones((Y.shape[0], 1))
    X1 = np.column_stack([np.ones((Y.shape[0], 1)), X[:, 0]])
    def sscp_resid(Y, Xd):
        B = np.linalg.pinv(Xd.T @ Xd) @ Xd.T @ Y
        R = Y - Xd @ B
        return R.T @ R
    E_full = sscp_resid(Y, X1)
    E_red  = sscp_resid(Y, X0)
    H = E_red - E_full
    E = E_full
    return float(np.trace(H @ np.linalg.pinv(H + E)))

def permute_within_subject(series, subjects, rng):
    """Permute values within each subject (used for MANOVA only)."""
    arr = np.asarray(series)
    out = arr.copy()
    for pid in np.unique(subjects):
        idx = np.where(subjects == pid)[0]
        if len(idx) > 1:
            out[idx] = rng.permutation(arr[idx])
    return out

def permutation_manova(df, y_cols, xcol, subj_col, n_perm=10000, seed=11):
    """Valid repeated-measures MANOVA via within-subject permutation of X."""
    Y = df[y_cols].to_numpy()
    X = df[[xcol]].to_numpy()
    subs = df[subj_col].to_numpy()
    V_obs = pillai_trace(Y, X)
    rng = np.random.default_rng(seed)
    V_perm = np.empty(n_perm, dtype=float)
    for b in range(n_perm):
        Xp = X.copy()
        Xp[:, 0] = permute_within_subject(df[xcol].to_numpy(), subs, rng)
        V_perm[b] = pillai_trace(Y, Xp)
    p = (np.sum(V_perm >= V_obs) + 1) / (n_perm + 1)
    return V_obs, p, V_perm

def standardised_beta(y_series, x_series):
    y = (y_series - y_series.mean()) / y_series.std(ddof=0)
    x = (x_series - x_series.mean()) / x_series.std(ddof=0)
    X = sm.add_constant(x.to_numpy())
    res = sm.OLS(y.to_numpy(), X).fit()
    return float(res.params[1]), float(res.bse[1]), float(res.pvalues[1])

def canonical_analysis(Y, X):
    """First canonical variate for single predictor X (lambda)."""
    beta = np.linalg.pinv(X) @ Y              # 1 x p
    weights = beta[0, :]
    scores = Y @ weights                      # n x 1 (Y-side canonical variate)
    struct = np.array([np.corrcoef(Y[:, j], scores.squeeze())[0, 1] for j in range(Y.shape[1])])
    return weights, struct, scores

def permutation_corr_fullshuffle(x, y, n_perm=10000, seed=11):
    """Full-shuffle permutation test for correlation between x and y."""
    r_obs, _ = pearsonr(x, y)
    rng = np.random.default_rng(seed)
    r_null = np.empty(n_perm)
    for i in range(n_perm):
        r_null[i], _ = pearsonr(rng.permutation(x), y)
    p_perm = (np.sum(np.abs(r_null) >= abs(r_obs)) + 1) / (n_perm + 1)
    return r_obs, p_perm, r_null

def cv_corr_rowwise(x, y, k=5, n_perm=2000, seed=11):
    """
    Row-wise KFold CV for correlation:
    - Fit OLS on train, predict y on test from x.
    - Compute r(y_test, y_pred_test) per fold.
    - Return mean r, median r, permutation p-value (full-shuffle null), fold r's, and null distribution.
    """
    kf = KFold(n_splits=k, shuffle=True, random_state=seed)
    r_folds = []
    for train_idx, test_idx in kf.split(x):
        fit = sm.OLS(y[train_idx], sm.add_constant(x[train_idx])).fit()
        y_pred = fit.predict(sm.add_constant(x[test_idx]))
        r, _ = pearsonr(y[test_idx], y_pred)
        r_folds.append(r)
    mean_r = float(np.mean(r_folds))
    median_r = float(np.median(r_folds))

    # Permutation p-value for the mean r
    rng = np.random.default_rng(seed)
    cv_null = np.empty(n_perm)
    for i in range(n_perm):
        y_perm = rng.permutation(y)
        r_perm_folds = []
        for train_idx, test_idx in kf.split(x):
            fit = sm.OLS(y_perm[train_idx], sm.add_constant(x[train_idx])).fit()
            y_pred = fit.predict(sm.add_constant(x[test_idx]))
            r, _ = pearsonr(y_perm[test_idx], y_pred)
            r_perm_folds.append(r)
        cv_null[i] = np.mean(r_perm_folds)
    p_perm = (np.sum(np.abs(cv_null) >= abs(mean_r)) + 1) / (n_perm + 1)

    return mean_r, median_r, p_perm, r_folds, cv_null

# ------------------
# Main
# ------------------
def main():
    OUTDIR.mkdir(parents=True, exist_ok=True)
    with open(OUTDIR / "SESSION_INFO.txt", "w") as f:
        f.write(f"Run: {datetime.now().isoformat()}\n")

    # Load
    df = pd.read_excel(EXCEL_PATH)
    df.columns = [c.strip() for c in df.columns]
    df = df[~df["anon_participant"].isin(EXCLUDE_SUBJECTS)].copy()
    df["DrugBinary"] = (~df["drug"].astype(str).str.lower().str.contains("placebo")).astype(int)

    # PANAS scaling
    if "panas_pos_affect_median_merge" in df.columns:
        df["panas_pos_affect_median_merge"] = df["panas_pos_affect_median_merge"].astype(float) / 50.0
    else:
        raise ValueError("Missing column 'panas_pos_affect_median_merge'.")

    y_cols = ACUTE + PSYCHOSOCIAL
    needed = ["anon_participant", "lambda_val", "DrugBinary"] + y_cols
    df = df.dropna(subset=needed).copy()

    # Preprocessing: subject-FE residualisation + z-score
    cols = ["lambda_val"] + y_cols
    df_da = residualise_within_subject(df, "anon_participant", cols, "DrugBinary")
    df_da[cols] = df_da[cols].apply(zscore)

    # =======================
    # MANOVA permutation (within-subject)
    # =======================
    V_obs, p_perm, V_perm = permutation_manova(
        df_da, y_cols, "lambda_val", "anon_participant", N_PERM_MANOVA, SEED
    )
    pd.DataFrame({
        "Pillai_trace":[V_obs],
        "Permutation_p":[p_perm],
        "N_outcomes":[len(y_cols)],
        "N_rows":[len(df_da)],
        "N_subjects":[df_da["anon_participant"].nunique()],
        "N_perm":[N_PERM_MANOVA],
        "Seed":[SEED]
    }).to_csv(OUTDIR/"MANOVA_summary.csv", index=False)

    plt.figure(figsize=(7,5))
    plt.hist(V_perm, bins=40)
    plt.axvline(V_obs, linestyle="--", linewidth=2)
    plt.xlabel("Pillai's trace under null (λ permuted within subject)")
    plt.ylabel("Count")
    plt.title("Permutation MANOVA null vs observed")
    plt.tight_layout()
    plt.savefig(OUTDIR/"Fig_MANOVA_PermutationNull.png", dpi=160)
    plt.show(); plt.close()

    # =======================
    # Univariate follow-ups + FDR
    # =======================
    rows = []
    for y in y_cols:
        beta, se, p = standardised_beta(df_da[y], df_da["lambda_val"])
        rows.append({"Outcome": y, "Std_Beta_lambda": beta, "SE": se, "p": p})
    uni = pd.DataFrame(rows)
    uni["p_FDR"] = fdr_bh(uni["p"].values)
    uni = uni.sort_values("Std_Beta_lambda", ascending=False).reset_index(drop=True)
    uni.to_csv(OUTDIR/"Univariate_lambda_FDR.csv", index=False)

    # Lollipop plot
    plt.figure(figsize=(7,5))
    y_order = uni["Outcome"]; x_vals = uni["Std_Beta_lambda"]; y_pos = np.arange(len(y_order))
    plt.hlines(y=y_pos, xmin=0, xmax=x_vals); plt.plot(x_vals, y_pos, "o")
    plt.yticks(y_pos, y_order); plt.axvline(0, linewidth=1)
    for i, (beta, pfdr) in enumerate(zip(uni["Std_Beta_lambda"], uni["p_FDR"])):
        if pfdr < 0.05: plt.plot(beta, i, "s", markersize=6)
    plt.xlabel("Standardised β (λ)  [Drug-aware + Subject-FE]")
    plt.title("Within-subject λ effects per outcome (FDR<.05 marked)")
    plt.tight_layout(); plt.savefig(OUTDIR/"Fig_StdBetas_Lollipop.png", dpi=160)
    plt.show()

    # =======================
    # Canonical analysis (point estimate) + scatter
    # =======================
    Y = df_da[y_cols].to_numpy()
    X = df_da[["lambda_val"]].to_numpy()
    can_w, struct_c, can_scores = canonical_analysis(Y, X)

    pd.DataFrame({"Outcome": y_cols,
                  "CanonicalWeight": can_w,
                  "StructureCoeff": struct_c}).to_csv(OUTDIR/"Canonical_metrics.csv", index=False)

    lam = X.squeeze()
    r_val, p_val = pearsonr(lam, can_scores.squeeze())
    grid = np.linspace(lam.min(), lam.max(), 200)
    Xline = sm.add_constant(lam); fit = sm.OLS(can_scores.squeeze(), Xline).fit()
    yhat = fit.params[0] + fit.params[1]*grid

    plt.figure(figsize=(6,5))
    plt.scatter(lam, can_scores, alpha=0.8)
    plt.plot(grid, yhat, linewidth=2)
    plt.xlabel("λ (self–other mergence), z-scored")
    plt.ylabel("Canonical variate (Y side), z-scored")
    plt.title(f"λ vs canonical variate  (r = {r_val:.2f}, p = {p_val:.4g})")
    plt.tight_layout(); plt.savefig(OUTDIR/"Fig_Scatter_Lambda_vs_Canonical.png", dpi=160)
    plt.show(); plt.close()

    # Structure coefficients (point)
    order = np.argsort(-np.abs(struct_c)); outs_sorted = np.array(y_cols)[order]; vals_sorted = np.array(struct_c)[order]
    plt.figure(figsize=(7,5)); plt.barh(outs_sorted, vals_sorted); plt.axvline(0, linewidth=1); plt.gca().invert_yaxis()
    plt.xlabel("Structure coefficient (r with canonical score)"); plt.title("Correlation of each outcome with canonical variate")
    plt.tight_layout(); plt.savefig(OUTDIR/"Fig_StructureCoeffs_Point.png", dpi=160); plt.show()

    # Canonical weights (point)
    order_w = np.argsort(-np.abs(can_w)); outs_sorted_w = np.array(y_cols)[order_w]
    plt.figure(figsize=(7,5)); plt.barh(outs_sorted_w, np.array(can_w)[order_w]); plt.axvline(0, linewidth=1); plt.gca().invert_yaxis()
    plt.xlabel("Canonical weight"); plt.title("Canonical weights (point estimate)")
    plt.tight_layout(); plt.savefig(OUTDIR/"Fig_CanonicalWeights_Point.png", dpi=160); plt.show()

    # =======================
    # Correlation test via full-shuffle permutation
    # =======================
    r_obs_full, p_full, r_null_full = permutation_corr_fullshuffle(lam, can_scores.squeeze(),
                                                                   n_perm=N_PERM_CORR, seed=SEED)
    pd.DataFrame({
        "r_obs":[r_obs_full], "p_perm":[p_full],
        "N_perm":[N_PERM_CORR], "null":"full_shuffle"
    }).to_csv(OUTDIR/"CorrTest_FullShuffle.csv", index=False)

    plt.figure(figsize=(7,5))
    plt.hist(r_null_full, bins=40)
    plt.axvline(r_obs_full, linestyle="--", linewidth=2)
    plt.title(f"Correlation null (full shuffle) — p = {p_full:.4f}")
    plt.xlabel("r under null"); plt.ylabel("Count")
    plt.tight_layout(); plt.savefig(OUTDIR/"Fig_CorrNull_FullShuffle.png", dpi=160); plt.show()

    # =======================
    # Repeated-measures correlation (Pearson rmcorr)
    # =======================
    # Attach canonical score to df_da (row-aligned with how can_scores was computed)
    df_da = df_da.copy()
    df_da["canonical_score"] = can_scores.squeeze()
    
    # Pandas-native within-subject centering
    s = df_da["anon_participant"]
    x = df_da["lambda_val"].astype(float)
    y = df_da["canonical_score"].astype(float)
    xc = x - x.groupby(s).transform("mean")
    yc = y - y.groupby(s).transform("mean")
    
    # Common slope + intercept (rmcorr line)
    slope = float((xc * yc).sum() / (xc * xc).sum())
    intercept = float(y.mean() - slope * x.mean())
    
    # r and p from pingouin (fallback if unavailable)
    try:
        rmcorr_res = pg.rm_corr(data=df_da, subject="anon_participant",
                                x="lambda_val", y="canonical_score")
        r_rm = float(rmcorr_res["r"].values[0])
        p_rm = float(rmcorr_res["pval"].values[0])
    except Exception as e:
        # Fallback: correlation of centered variables (approximate rmcorr r); p left as NaN
        r_rm = float(np.corrcoef(xc, yc)[0, 1])
        p_rm = float("nan")
        print("WARNING: pingouin.rm_corr failed; using centered r only. Error:", repr(e))
    
    # Fitted line points (for plotting in MATLAB)
    grid2 = np.linspace(float(x.min()), float(x.max()), 200)
    yhat2 = intercept + slope * grid2
    
    # ---- SAVE exports for MATLAB ----
    OUTDIR.mkdir(parents=True, exist_ok=True)
    
    # Scatter with DrugBinary so you can use shapes in MATLAB
    df_da["drug"] = df["drug"].values  # carry the original condition label
    cols_to_save = ["anon_participant","lambda_val","canonical_score","DrugBinary","drug"]
    df_da[cols_to_save].to_csv(OUTDIR / "rmcorr_scatter.csv", index=False)

    # Stats (r, p, slope, intercept)
    pd.DataFrame({"r": [r_rm], "p": [p_rm],
                  "slope": [slope], "intercept": [intercept]}).to_csv(
        OUTDIR / "rmcorr_stats.csv", index=False
    )
    
    # Fitted line points
    pd.DataFrame({"lambda_val": grid2, "yhat": yhat2}).to_csv(
        OUTDIR / "rmcorr_line.csv", index=False
    )
  
    plt.figure(figsize=(6,5))
    for pid, grp in df_da.groupby("anon_participant"):
        plt.plot(grp["lambda_val"], grp["canonical_score"], "o-", alpha=0.35, linewidth=1)
    plt.plot(grid2, yhat2, "k-", linewidth=2)
    title_str = f"Repeated-measures correlation (Pearson)\nr = {r_rm:.2f}"
    if np.isfinite(p_rm):
        title_str += f", p = {p_rm:.4g}"
    plt.title(title_str)
    plt.xlabel("λ (self–other mergence), z-scored")
    plt.ylabel("Canonical variate (Y side), z-scored")
    plt.tight_layout()
    plt.savefig(OUTDIR / "Fig_rmcorr_Pearson.png", dpi=160)
    plt.show(); plt.close()
    


    # =======================
    # Cross-validation (K=5) with permutation p-value
    # =======================
    mean_r, median_r, p_cv, r_folds, r_null_cv = cv_corr_rowwise(lam, can_scores.squeeze(),
                                                                 k=K_ROW, n_perm=N_PERM_CV, seed=SEED)
    # Save CV summaries
    pd.DataFrame({
        "K":[K_ROW], "mean_r":[mean_r], "median_r":[median_r],
        "p_perm":[p_cv], "N_perm":[N_PERM_CV]
    }).to_csv(OUTDIR/"CV_KFold_Summary.csv", index=False)

    pd.DataFrame({"Fold": np.arange(1, K_ROW+1), "r_fold": r_folds}).to_csv(
        OUTDIR/"CV_KFold_Folds.csv", index=False
    )

    plt.figure(figsize=(7,5))
    plt.hist(r_null_cv, bins=40)
    plt.axvline(mean_r, linestyle="--", linewidth=2)
    plt.title(f"Cross-validation null (K={K_ROW}) — p = {p_cv:.4f}")
    plt.xlabel("Mean r under null"); plt.ylabel("Count")
    plt.tight_layout(); plt.savefig(OUTDIR/"Fig_CVNull.png", dpi=160); plt.show()

    # ----------------------
    # Concise console summary
    # ----------------------
    print("\n=== Results ===")
    print(f"MANOVA: Pillai = {V_obs:.4f}, p_perm = {p_perm:.4g} (within-subject null)")
    print(f"Correlation (full shuffle): r_obs = {r_obs_full:.3f}, p_perm = {p_full:.4g}")
    print(f"rmcorr (Pearson): r = {r_rm:.3f}, p = {p_rm:.4g}")
    print(f"Cross-validation (K={K_ROW}): mean r = {mean_r:.3f}, p_perm = {p_cv:.4g}")
    print(f"\nAll outputs saved in: {OUTDIR.resolve()}")

if __name__ == "__main__":
    main()

