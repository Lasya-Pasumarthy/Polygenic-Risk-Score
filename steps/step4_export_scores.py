#!/usr/bin/env python3
# steps/step4_export_scores.py
# Parses pgsc_calc's aggregated scores output and exports a clean, joinable table.

import glob
import gzip
import os
import sys

import pandas as pd


# ── pgsc_calc output detection ──────────────────────────────────────────────

def find_scores_file(output_dir):
    """
    Search for pgsc_calc's aggregated scores file under the output directory.
    Typical location: <outdir>/score/<sampleset>/aggregated_scores.txt.gz
    """
    patterns = [
        os.path.join(output_dir, "**", "aggregated_scores.txt.gz"),
        os.path.join(output_dir, "**", "aggregated_scores.txt"),
    ]
    for pat in patterns:
        matches = glob.glob(pat, recursive=True)
        if matches:
            return matches[0]
    return None


def load_scores(scores_path):
    """Load pgsc_calc aggregated scores file into a DataFrame."""
    if scores_path.endswith(".gz"):
        with gzip.open(scores_path, "rt") as f:
            df = pd.read_csv(f, sep="\t", comment="#")
    else:
        df = pd.read_csv(scores_path, sep="\t", comment="#")
    return df


# ── Column detection ────────────────────────────────────────────────────────

def detect_id_col(df):
    """Return the sample ID column name. pgsc_calc uses 'IID'; fall back to first column."""
    for candidate in ["IID", "sample_id", "SAMPLE", "ID"]:
        if candidate in df.columns:
            return candidate
    return df.columns[0]


def detect_score_cols(df, id_col):
    """Return all non-ID, non-metadata columns (the PRS score columns)."""
    exclude = {id_col, "#IID", "FID", "sampleset", "cohort", "split"}
    return [c for c in df.columns if c not in exclude]


# ── Export ──────────────────────────────────────────────────────────────────

def export(df, export_path):
    """Export to both CSV and XLSX."""
    os.makedirs(os.path.dirname(export_path) or ".", exist_ok=True)

    # CSV
    csv_path = export_path if export_path.endswith(".csv") else export_path.rsplit(".", 1)[0] + ".csv"
    df.to_csv(csv_path, index=False)
    print(f"  [FILE] CSV  → {csv_path}")

    # XLSX
    xlsx_path = export_path if export_path.endswith(".xlsx") else export_path.rsplit(".", 1)[0] + ".xlsx"
    try:
        with pd.ExcelWriter(xlsx_path, engine="openpyxl") as writer:
            df.to_excel(writer, index=False, sheet_name="PRS_Scores")
            # Auto-fit column widths
            ws = writer.sheets["PRS_Scores"]
            for col_cells in ws.columns:
                max_len = max(len(str(cell.value)) if cell.value else 0 for cell in col_cells)
                ws.column_dimensions[col_cells[0].column_letter].width = min(max_len + 2, 40)
        print(f"  [FILE] XLSX → {xlsx_path}")
    except ImportError:
        print("  [WARN] 'openpyxl' not installed — XLSX export skipped. Run: pip install openpyxl")

    return csv_path, xlsx_path


# ── Main entry ──────────────────────────────────────────────────────────────

def run(config):
    print("Step 4: Exporting PRS scores table\n")

    output_dir    = config.get("output_dir", "").strip()
    scores_file   = config.get("scores_file", "").strip()
    export_path   = config.get("export_output", "prs_scores_export.xlsx").strip()

    # ── Locate scores file ────────────────────────────────────────────────────
    if scores_file and os.path.isfile(scores_file):
        scores_path = scores_file
        print(f"[INFO] Using scores file from config: {scores_path}")
    elif output_dir:
        print(f"[INFO] Auto-detecting scores file in: {output_dir}")
        scores_path = find_scores_file(output_dir)
    else:
        scores_path = None

    if not scores_path:
        # Last resort: ask user
        scores_path = input(
            "Could not auto-detect scores file.\n"
            "Enter absolute path to aggregated_scores.txt.gz: "
        ).strip()

    if not os.path.isfile(scores_path):
        print(f"[ERROR] Scores file not found: {scores_path}")
        return

    print(f"[INFO] Loading: {scores_path}")

    # ── Load and inspect ───────────────────────────────────────────────────────
    try:
        df = load_scores(scores_path)
    except Exception as e:
        print(f"[ERROR] Could not read scores file: {e}")
        return

    print(f"[INFO] Raw file: {len(df)} rows, {len(df.columns)} columns")
    print(f"       Columns: {', '.join(df.columns.tolist())}\n")

    id_col     = detect_id_col(df)
    score_cols = detect_score_cols(df, id_col)

    print(f"[INFO] Sample ID column: '{id_col}'")
    print(f"[INFO] Score column(s):  {score_cols}\n")

    if not score_cols:
        print("[ERROR] No score columns detected. Please check the scores file format.")
        return

    # ── Build clean export table ──────────────────────────────────────────────
    export_df = df[[id_col] + score_cols].copy()
    export_df = export_df.rename(columns={id_col: "Patient_ID"})

    # If only one score column, rename it to 'PRS' for clarity
    if len(score_cols) == 1:
        export_df = export_df.rename(columns={score_cols[0]: "PRS"})
        print("[INFO] Single PRS column renamed to 'PRS'.")
    else:
        print(f"[INFO] Multiple score columns found ({len(score_cols)}). "
              "All will be included as-is.")

    print(f"\nPreview (first 5 rows):")
    print(export_df.head().to_string(index=False))

    # ── Confirm and export ────────────────────────────────────────────────────
    print(f"\n[INFO] Export path: {export_path}")
    confirm = input("Proceed with export? [Y/n]: ").strip().lower()
    if confirm == "n":
        print("[ABORT] Export cancelled.")
        return

    csv_path, xlsx_path = export(export_df, export_path)

    print(
        f"\n[OK] Export complete."
        f"\n     Rows exported: {len(export_df)}"
        f"\n\n[HINT] To append PRS to your patient metadata file in Python:"
        f"\n       import pandas as pd"
        f"\n       meta = pd.read_excel('patient_metadata.xlsx')"
        f"\n       prs  = pd.read_csv('{csv_path}')"
        f"\n       merged = meta.merge(prs, left_on='Patient_ID', right_on='Patient_ID', how='left')"
        f"\n       merged.to_excel('patient_metadata_with_PRS.xlsx', index=False)"
    )
