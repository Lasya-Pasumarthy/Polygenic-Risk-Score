#!/usr/bin/env python3
# steps/step2_run_pgsc_calc.py
# Wraps Nextflow + pgsc_calc. pgsc_calc downloads scoring files natively via PGS ID.

import os
import subprocess
import sys


def prompt_optional(flag_name, description, example=""):
    """Prompt user for an optional flag. Returns value or empty string."""
    hint = f" (e.g. {example})" if example else ""
    val = input(f"  --{flag_name}{hint}  [{description}]  (leave blank to skip): ").strip()
    return val


def build_command(samplesheet, pgs_id, target_build, profile, op_dir, optional_flags):
    """Assemble the nextflow run command."""
    cmd = [
        "nextflow", "run", "pgscatalog/pgsc_calc",
        "-profile", profile,
        "--input",        samplesheet,
        "--pgs_id",       pgs_id,
        "--target_build", target_build,
        "--outdir",       op_dir,
    ]
    for flag, value in optional_flags.items():
        if value:
            cmd.extend([f"--{flag}", value])
    return cmd


def run(config):
    print("Step 2: Running pgsc_calc\n")

    # ── Required inputs ──────────────────────────────────────────────────────
    samplesheet = input(
        "Enter absolute path to samplesheet.csv\n"
        "(generated in Step 1, or provide manually): "
    ).strip()

    if not os.path.isfile(samplesheet):
        print(f"[ERROR] Samplesheet not found: {samplesheet}")
        return

    pgs_id = input(
        "\nEnter PGS Catalog ID (e.g. PGS000039) or comma-separated IDs: "
    ).strip()
    if not pgs_id:
        print("[ERROR] PGS ID is required.")
        return

    # ── Config values ────────────────────────────────────────────────────────
    target_build = config.get("target_build", "GRCh37").strip()
    profile      = config.get("profile", "conda").strip()
    op_dir       = config.get("output_dir", "").strip()

    if not op_dir:
        print("[ERROR] 'output_dir' is not set in config.yaml.")
        return

    os.makedirs(op_dir, exist_ok=True)

    # ── Optional flags ───────────────────────────────────────────────────────
    print("\n[Optional flags]  Press Enter to skip any you don't need.\n")

    optional_flags = {}

    optional_flags["min_overlap"] = prompt_optional(
        "min_overlap", "minimum variant overlap fraction", "0.75"
    )
    optional_flags["ancestry_method"] = prompt_optional(
        "ancestry_method", "ancestry inference method", "pca_projection"
    )
    optional_flags["ref"] = prompt_optional(
        "ref", "ancestry reference panel RDS file"
    )
    optional_flags["liftover"] = prompt_optional(
        "liftover", "enable liftover if build mismatch (true/false)", "true"
    )
    optional_flags["keep_ambiguous"] = prompt_optional(
        "keep_ambiguous", "keep ambiguous SNPs (true/false)", "false"
    )

    # ── Build and preview command ────────────────────────────────────────────
    cmd = build_command(samplesheet, pgs_id, target_build, profile, op_dir, optional_flags)

    print("\n[INFO] Command to be run:")
    print("  " + " \\\n    ".join(cmd))

    confirm = input("\nProceed? [y/N]: ").strip().lower()
    if confirm != "y":
        print("[ABORT] Cancelled by user.")
        sys.exit(0)

    # ── Run ──────────────────────────────────────────────────────────────────
    print("\n[INFO] Launching pgsc_calc via Nextflow...\n")
    try:
        subprocess.run(cmd, check=True)
        print("\n[OK] pgsc_calc completed successfully.")
        print(f"     Results are in: {op_dir}")
        print(
            "\n[NEXT] Run Step 4 to export a clean scores table:\n"
            "       python main.py --step export_scores"
        )
    except FileNotFoundError:
        print(
            "[ERROR] 'nextflow' not found. Make sure Nextflow is installed and on your PATH.\n"
            "        Install guide: https://www.nextflow.io/docs/latest/install.html"
        )
    except subprocess.CalledProcessError as e:
        print(f"\n[ERROR] pgsc_calc failed with exit code {e.returncode}.")
        print("        Check the Nextflow log above for details.")
        print(
            "        Common issues:\n"
            "          - Samplesheet format mismatch (re-check Step 1 output)\n"
            "          - PGS ID not found on PGS Catalog\n"
            "          - Genome build mismatch between samplesheet and --target_build\n"
            "          - Insufficient memory/CPU for the selected profile\n"
            "\n        pgsc_calc docs: https://pgsc-calc.readthedocs.io"
        )
