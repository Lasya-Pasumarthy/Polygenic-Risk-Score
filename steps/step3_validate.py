#!/usr/bin/env python3
# steps/step3_validate.py
# Statistical validation: AUC, Incremental R², Odds Ratios, Cohen's d, t-test, plots.

import os
import sys
from math import sqrt

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")  # non-interactive backend for server environments
import matplotlib.pyplot as plt
from scipy.stats import ttest_ind
from sklearn.metrics import roc_auc_score, roc_curve
from statsmodels.discrete.discrete_model import Logit
from statsmodels.tools import add_constant


# ── Statistics helpers ──

def cohen_d(x, y):
    nx, ny = len(x), len(y)
    pooled_std = sqrt(
        ((nx - 1) * np.std(x, ddof=1) ** 2 + (ny - 1) * np.std(y, ddof=1) ** 2)
        / (nx + ny - 2)
    )
    return (np.mean(x) - np.mean(y)) / pooled_std


def fit_logit(y, X):
    """Fit a statsmodels Logit model; return model result or None on failure."""
    try:
        return Logit(y, add_constant(X)).fit(disp=False)
    except Exception as e:
        print(f"  [WARN] Model fitting failed: {e}")
        return None


# ── Plotting helpers ──

def plot_roc(y, y_prob_base, y_prob_full, auc_base, auc_full, out_dir):
    fpr_b, tpr_b, _ = roc_curve(y, y_prob_base)
    fpr_f, tpr_f, _ = roc_curve(y, y_prob_full)

    fig, ax = plt.subplots(figsize=(6, 5))
    ax.plot(fpr_b, tpr_b, label=f"Base model  (AUC = {auc_base:.3f})", linestyle="--")
    ax.plot(fpr_f, tpr_f, label=f"+ PRS model (AUC = {auc_full:.3f})")
    ax.plot([0, 1], [0, 1], "k:", linewidth=0.8)
    ax.set_xlabel("False Positive Rate")
    ax.set_ylabel("True Positive Rate")
    ax.set_title("ROC Curve")
    ax.legend(loc="lower right")
    path = os.path.join(out_dir, "roc_curve.png")
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  [PLOT] ROC curve saved at: {path}")


def plot_risk_gradient(data, prs_col, out_dir, n_deciles=10):
    """Risk gradient: prevalence of cases per PRS decile."""
    data = data.copy()
    data["decile"] = pd.qcut(data[prs_col], q=n_deciles, labels=False) + 1
    decile_prev = data.groupby("decile")["Status"].mean()

    fig, ax = plt.subplots(figsize=(7, 4))
    bars = ax.bar(decile_prev.index, decile_prev.values * 100, color="steelblue", edgecolor="white")
    ax.set_xlabel("PRS Decile")
    ax.set_ylabel("Observed Prevalence (%)")
    ax.set_title("Risk Gradient by PRS Decile")
    ax.set_xticks(range(1, n_deciles + 1))
    path = os.path.join(out_dir, "risk_gradient.png")
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  [PLOT] Risk gradient saved at: {path}")


def plot_prs_distribution(data, prs_col, out_dir):
    """Overlapping PRS distribution for cases vs controls."""
    cases    = data[data["Status"] == 1][prs_col]
    controls = data[data["Status"] == 0][prs_col]

    fig, ax = plt.subplots(figsize=(7, 4))
    ax.hist(controls, bins=40, alpha=0.6, label="Controls", color="steelblue", density=True)
    ax.hist(cases,    bins=40, alpha=0.6, label="Cases",    color="tomato",    density=True)
    ax.set_xlabel("PRS")
    ax.set_ylabel("Density")
    ax.set_title("PRS Distribution: Cases vs Controls")
    ax.legend()
    path = os.path.join(out_dir, "prs_distribution.png")
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  [PLOT] PRS distribution saved at: {path}")


# ── Main entry ──

def run(config):
    print("Step 3: Validating PRS\n")

    prs_file = config.get("data_file")
    if not os.path.isfile(prs_file):
        print("Error: File not found.")
        return
    val_output_dir = config.get("val_output_dir","val_output").strip()
    os.makedirs(val_output_dir, exist_ok = True)
    data = pd.read_csv(prs_file, sep="\t" if prs_file.endswith(".txt") else ",")
    if "Sex" in data.columns and data["Sex"].dtype == "object":
         print("Converting 'Sex' from Male/Female to 0/1...")
         data["Sex"] = data["Sex"].map({"Male": 0, "Female": 1})

    print("Detected columns in dataset:")
    print(data.dtypes)
    prs_col = input("Enter the column name for Scores: ").strip()
    if prs_col not in data.columns:
        print("Column not found.")
        return

    covariates = input("\nEnter comma-separated covariate column names from dataset: ").strip()
    cov = [x.strip() for x in covariates.split(",") if x.strip() in data.columns]

    if not cov:
        print("Error: No covariates found.")
        return

    print(f"Using the covariates: {cov}")

    data["PRS"] = data[prs_col]

    X_base = add_constant(data[cov])
    X_full = add_constant(data[["PRS"] + cov])
    y = data["Status"]

    #val_output_dir = config.get("val_output_dir", "val_output").strip()
    #os.makedirs(val_output_dir, exist_ok=True)

    # ── Load data ──
    """data_file = config.get("data_file", "").strip()
    if not data_file:
        data_file = input("Enter absolute path to merged (PRS + phenotype + covariates) file: ").strip()

    if not os.path.isfile(data_file):
        print(f"[ERROR] File not found: {data_file}")
        return

    sep = "\t" if data_file.endswith(".txt") or data_file.endswith(".tsv") else ","
    try:
        data = pd.read_csv(data_file, sep=sep)
    except Exception as e:
        print(f"[ERROR] Could not read file: {e}")
        return

    print(f"\n[INFO] Loaded {len(data)} rows, {len(data.columns)} columns.")
    print("Columns:", ", ".join(data.columns.tolist()))

    # ── Column selection ──
    prs_col = input("\nEnter the column name for PRS scores: ").strip()
    if prs_col not in data.columns:
        print(f"[ERROR] Column '{prs_col}' not found.")
        return

    status_col = input("Enter the column name for case/control status (0=control, 1=case): ").strip()
    if status_col not in data.columns:
        print(f"[ERROR] Column '{status_col}' not found.")
        return

    # Rename to 'Status' internally for clarity
    data = data.rename(columns={status_col: "Status"})
    data["Status"] = pd.to_numeric(data["Status"], errors="coerce")
    data = data.dropna(subset=["Status", prs_col])
    data["Status"] = data["Status"].astype(int)

    covariates_input = input(
        "\nEnter comma-separated covariate column names (e.g. age,sex,PC1,PC2): "
    ).strip()
    cov_cols = [c.strip() for c in covariates_input.split(",") if c.strip() in data.columns]
    missing_covs = [c.strip() for c in covariates_input.split(",") if c.strip() not in data.columns and c.strip()]
    if missing_covs:
        print(f"  [WARN] Covariates not found and will be ignored: {missing_covs}")
    if not cov_cols:
        print("[ERROR] No valid covariates found. At minimum, include age or sex.")
        return

    print(f"\n[INFO] Using covariates: {cov_cols}")"""

    # Drop rows with missing values in any relevant column
    cols_needed = ["Status", prs_col] + cov
    data = data[cols_needed].dropna()
    print(f"[INFO] Analysis sample size after dropping NAs: {len(data)}")

    y = data["Status"]
    X_cov = data[cov]
    X_full = pd.concat([data[[prs_col]], X_cov], axis=1)

    # ── Model fitting ──
    print("\n[INFO] Fitting logistic regression models...")
    model_base = fit_logit(y, X_cov)
    model_full = fit_logit(y, X_full)

    if model_base is None or model_full is None:
        print("[ERROR] Model fitting failed. Check your data for separation or multicollinearity.")
        return

    # ── Incremental R2 (McFadden) ──
    r2_base = 1 - (model_base.llf / model_base.llnull)
    r2_full = 1 - (model_full.llf / model_full.llnull)
    inc_r2  = r2_full - r2_base

    # ── AUC ──
    y_prob_base = model_base.predict(add_constant(X_cov))
    y_prob_full = model_full.predict(add_constant(X_full))
    auc_base = roc_auc_score(y, y_prob_base)
    auc_full = roc_auc_score(y, y_prob_full)

    # ── Odds Ratios ──
    or_vals = np.exp(model_full.params)
    ci      = np.exp(model_full.conf_int())
    or_table = pd.DataFrame({
        "OR":       or_vals,
        "CI Lower": ci.iloc[:, 0],
        "CI Upper": ci.iloc[:, 1],
        "p-value":  model_full.pvalues,
    })

    # ── Cohen's d and t-test ──
    cases    = data[data["Status"] == 1][prs_col].values
    controls = data[data["Status"] == 0][prs_col].values
    d        = cohen_d(cases, controls)
    t_stat, p_val = ttest_ind(cases, controls)

    # ── Print summary ──
    print("\n" + "="*55)
    print("  Validation Results")
    print("="*55)
    print(f"  McFadden R² (base model):        {r2_base:.4f}")
    print(f"  McFadden R² (PRS model):         {r2_full:.4f}")
    print(f"  Incremental R² from PRS:         {inc_r2:.4f}")
    print(f"  AUC (base model):                {auc_base:.4f}")
    print(f"  AUC (+ PRS):                     {auc_full:.4f}")
    print(f"  Cohen's d:                       {d:.4f}")
    print(f"  t-test p-value (cases vs ctrl):  {p_val:.4e}")
    print("\n  Odds Ratios (full model):")
    print(or_table.to_string())
    print("="*55 + "\n")

    # ── Save stats to CSV ──
    stats_path = os.path.join(val_output_dir, "validation_summary.csv")
    summary = pd.DataFrame({
        "metric": [
            "AUC_base", "AUC_full", "McFadden_R2_base",
            "McFadden_R2_full", "Incremental_R2", "Cohens_d", "ttest_pvalue"
        ],
        "value": [
            auc_base, auc_full, r2_base, r2_full, inc_r2, d, p_val
        ]
    })
    summary.to_csv(stats_path, index=False)
    or_table.to_csv(os.path.join(val_output_dir, "odds_ratios.csv"))
    print(f"  [FILE] Summary stats saved at: {stats_path}")

    # ── Plots ──
    print("\n[INFO] Generating plots...")
    plot_roc(y, y_prob_base, y_prob_full, auc_base, auc_full, val_output_dir)
    plot_risk_gradient(data, prs_col, val_output_dir)
    plot_prs_distribution(data, prs_col, val_output_dir)

    print(f"\n[OK] Validation complete. All outputs in: {val_output_dir}")
