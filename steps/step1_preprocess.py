#!/usr/bin/env python3
# steps/step1_preprocess.py

# Detects genotype files, auto-generates a pgsc_calc samplesheet, and lets user review/edit.

import csv
import glob
import os
import re
import subprocess
import sys
import tempfile


# ── File detection helpers ──────────────────────────────────────────────────

PLINK1_EXTS = (".bed", ".bim", ".fam")
PLINK2_EXTS = (".pgen", ".pvar", ".psam")
VCF_EXTS    = (".vcf.gz", ".vcf")

def detect_format(prefix):
    """Return ('bfile'|'pfile'|'vcf', primary_ext) or (None, None)."""
    if all(os.path.exists(prefix + e) for e in PLINK1_EXTS):
        return "bfile", ".bed"
    if all(os.path.exists(prefix + e) for e in PLINK2_EXTS):
        return "pfile", ".pgen"
    for ext in VCF_EXTS:
        if os.path.exists(prefix + ext):
            return "vcf", ext
    return None, None

def strip_extensions(path):
    """Strip known genotype extensions from a path to get the bare prefix."""
    for ext in [".vcf.gz", ".vcf", ".bed", ".bim", ".fam", ".pgen", ".pvar", ".psam"]:
        if path.endswith(ext):
            return path[: -len(ext)]
    return path

def extract_chrom(prefix):
    """
    Try to extract a chromosome number from the filename.
    Handles patterns like: _c1_, _chr1_, .chr1., c1_b0, etc.
    Returns the chromosome string (e.g. '1', '22', 'X') or None.
    """
    fname = os.path.basename(prefix)
    patterns = [
        r'_c(\d+|X|Y|XY|MT)_',
        r'_chr(\d+|X|Y|XY|MT)[_\.]',
        r'\.chr(\d+|X|Y|XY|MT)\.',
        r'chr(\d+|X|Y|XY|MT)',
        r'_c(\d+|X|Y|XY|MT)$',
    ]
    for pat in patterns:
        m = re.search(pat, fname, re.IGNORECASE)
        if m:
            return m.group(1)
    return None

def find_prefixes(glob_pattern):
    """
    Expand the glob, strip extensions, and return a sorted unique set of prefixes.
    """
    matched = glob.glob(glob_pattern + "*")
    if not matched:
        # Maybe the prefix itself is exact (no glob needed)
        matched = glob.glob(glob_pattern)
    prefixes = set()
    for path in matched:
        prefixes.add(strip_extensions(path))
    return sorted(prefixes)


# ── Samplesheet building ────────────────────────────────────────────────────

def build_samplesheet_rows(prefixes):
    """
    Given a list of detected prefixes, build samplesheet rows.
    pgsc_calc samplesheet columns: sampleset, path_prefix, chrom, format
    VCF files use the full file path (including extension) in path_prefix.
    """
    rows = []
    errors = []

    for prefix in prefixes:
        fmt, _ = detect_format(prefix)
        if fmt is None:
            errors.append(f"  [SKIP] No recognised genotype files at: {prefix}")
            continue

        chrom = extract_chrom(prefix)

        # For VCF, path_prefix should be the full file path
        if fmt == "vcf":
            path_prefix = (prefix + ".vcf.gz") if os.path.exists(prefix + ".vcf.gz") else (prefix + ".vcf")
        else:
            # bfile / pfile: path_prefix is the bare prefix (no extension)
            path_prefix = prefix

        row = {
            "sampleset":   "default",  # overwritten below after user input
            "path_prefix": path_prefix,
            "chrom":       chrom if chrom else "",
            "format":      fmt,
        }
        rows.append(row)

    return rows, errors


# ── User review loop ────────────────────────────────────────────────────────

FIELDNAMES = ["sampleset", "path_prefix", "chrom", "format"]

def write_samplesheet(rows, path):
    os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
    with open(path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=FIELDNAMES)
        writer.writeheader()
        writer.writerows(rows)

def print_samplesheet(path):
    print(f"\n{'─'*60}")
    print(f"  Samplesheet preview  →  {path}")
    print(f"{'─'*60}")
    with open(path) as f:
        for line in f:
            print(" ", line, end="")
    print(f"{'─'*60}\n")

def open_in_editor(path):
    """Open samplesheet in the user's preferred editor."""
    editor = os.environ.get("EDITOR", "")
    if not editor:
        # Try common editors in order
        for ed in ["nano", "vim", "vi", "notepad"]:
            if subprocess.call(["which", ed], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL) == 0:
                editor = ed
                break
    if not editor:
        print("[WARN] No editor found. Please edit the file manually:")
        print(f"       {path}")
        return
    print(f"[INFO] Opening samplesheet in '{editor}'...")
    subprocess.call([editor, path])

def review_loop(rows, samplesheet_path):
    """Write samplesheet, show preview, let user edit, confirm, or abort."""
    write_samplesheet(rows, samplesheet_path)

    while True:
        print_samplesheet(samplesheet_path)
        print("Options:")
        print("  [c] Confirm and continue")
        print("  [e] Edit samplesheet in terminal editor")
        print("  [p] Re-print preview")
        print("  [q] Quit / abort pipeline")
        choice = input("\nYour choice: ").strip().lower()

        if choice == "c":
            print("[OK] Samplesheet confirmed.\n")
            return True
        elif choice == "e":
            open_in_editor(samplesheet_path)
        elif choice == "p":
            print_samplesheet(samplesheet_path)
        elif choice == "q":
            print("[ABORT] Exiting pipeline.")
            sys.exit(0)
        else:
            print("[?] Unrecognised option. Please enter c, e, p, or q.")


# ── Main entry ──────────────────────────────────────────────────────────────

def run(config):
    print("Step 1: Preprocessing check and samplesheet generation\n")

    genotype_prefix = config.get("genotype_prefix", "").strip()
    samplesheet_path = config.get("samplesheet_output", "samplesheet.csv").strip()

    if not genotype_prefix:
        print("[ERROR] 'genotype_prefix' is not set in config.yaml.")
        return

    if not os.path.isabs(genotype_prefix):
        print("[WARN] 'genotype_prefix' is a relative path. An absolute path is recommended.")

    print(f"[INFO] Scanning for genotype files with prefix:\n       {genotype_prefix}\n")

    prefixes = find_prefixes(genotype_prefix)

    if not prefixes:
        print("[ERROR] No files found matching the given prefix. Check config.yaml.")
        return

    print(f"[INFO] Found {len(prefixes)} file prefix(es):")
    for p in prefixes:
        fmt, _ = detect_format(p)
        chrom  = extract_chrom(p)
        status = f"format={fmt}, chrom={chrom or 'not detected'}" if fmt else "❌ unrecognised format"
        print(f"  {p}  →  {status}")

    print()
    rows, errors = build_samplesheet_rows(prefixes)

    if errors:
        print("[WARN] Some prefixes were skipped:")
        for e in errors:
            print(e)
        print()

    if not rows:
        print("[ERROR] No valid rows could be built. Check your input files.")
        return

    # Detect structure: per-chrom vs merged
    has_chrom = any(r["chrom"] for r in rows)
    if has_chrom:
        print(f"[INFO] Detected per-chromosome file structure ({len(rows)} chromosome(s)).")
    else:
        print("[INFO] Detected single merged file structure.")

    # Prompt for sampleset name
    sampleset_name = input(
        "\nEnter a sampleset name for your cohort (default: 'default'): "
    ).strip() or "default"
    for r in rows:
        r["sampleset"] = sampleset_name

    # Build and enter review loop
    confirmed = review_loop(rows, samplesheet_path)

    if confirmed:
        print(f"[OK] Samplesheet saved to: {samplesheet_path}")
        print("     Pass this path as 'samplesheet' when prompted in Step 2.\n")

    print("Step 1 complete.")
