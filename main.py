#!/usr/bin/env python3
# main.py — PRS Pipeline Entry Point

import argparse
import sys
import yaml
from steps import step1_preprocess, step2_run_pgsc_calc, step3_validate, step4_export_scores

STEPS = {
    "preprocess":     step1_preprocess.run,
    "run_pgsc_calc":  step2_run_pgsc_calc.run,
    "validate":       step3_validate.run,
    "export_scores":  step4_export_scores.run,
}

def load_config(path="config.yaml"):
    try:
        with open(path, "r") as f:
            return yaml.safe_load(f)
    except FileNotFoundError:
        print(f"[ERROR] config.yaml not found at '{path}'. Please create one.")
        sys.exit(1)

def main():
    parser = argparse.ArgumentParser(
        description="PRS Pipeline — Polygenic Risk Score calculation workflow",
        formatter_class=argparse.RawTextHelpFormatter
    )
    parser.add_argument(
        "--step", required=True,
        choices=list(STEPS.keys()),
        help=(
            "Pipeline step to run:\n"
            "  preprocess    — validate input files and generate samplesheet\n"
            "  run_pgsc_calc — run pgsc_calc via Nextflow\n"
            "  validate      — compute AUC, OR, R², Cohen's d, and plots\n"
            "  export_scores — export clean PRS scores table (CSV + XLSX)"
        )
    )
    parser.add_argument(
        "--config", default="config.yaml",
        help="Path to config file (default: config.yaml)"
    )
    args = parser.parse_args()

    config = load_config(args.config)
    print(f"\n{'='*55}")
    print(f"  PRS Pipeline  |  Step: {args.step}")
    print(f"{'='*55}\n")

    STEPS[args.step](config)

if __name__ == "__main__":
    main()
