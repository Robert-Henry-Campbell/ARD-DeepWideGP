#!/usr/bin/env python3
"""
script: create_UKB_pheno_casecount (2025)
description: loads the clean UKB ICD10 pheno df and calculates the case counts for each phenotype excluding psuedo_id column.
inputs:
1. cleaned UKB ICD10 pheno df (cleaned dnanexus output)
outputs:
1. UKB ICD10 pheno case count df

example usage:
python scripts/3.1-create_UKB_casecount.py \
  --input data/3-identify-UKB-casecount-ARDs/dna_nexus_wrangling/ukb_icd10_matrix_2025_cleaned.csv \
  --output-dir data/3-identify-UKB-casecount-ARDs/extract_casecount_ukb2025 \
  --output-name UKB_2025_casecount.csv

"""

import argparse
from pathlib import Path

import pandas as pd


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Calculate case counts per phenotype from a cleaned UKB ICD10 "
            "phenotype dataframe."
        )
    )
    parser.add_argument(
        "--input",
        required=True,
        help="Path to the cleaned UKB ICD10 phenotype dataframe (CSV/TSV/Parquet).",
    )
    parser.add_argument(
        "--output-dir",
        required=True,
        help="Directory where the case count dataframe will be written.",
    )
    parser.add_argument(
        "--output-name",
        default="ukb_icd10_case_counts.csv",
        help="Filename for the case count dataframe (default: ukb_icd10_case_counts.csv).",
    )
    parser.add_argument(
        "--id-column",
        default="psuedo_id",
        help="Name of the identifier column to exclude from counting (default: psuedo_id).",
    )
    return parser.parse_args()


def load_dataframe(path: Path) -> pd.DataFrame:
    if path.suffix.lower() in {".parquet", ".pq"}:
        return pd.read_parquet(path)
    # Simple delimiter guess: TSV if tab present in first line, else CSV.
    with path.open("r", encoding="utf-8") as f:
        first_line = f.readline()
    sep = "\t" if "\t" in first_line else ","
    return pd.read_csv(path, sep=sep)


def compute_case_counts(df: pd.DataFrame, id_column: str) -> pd.DataFrame:
    df = df.copy()
    if id_column in df.columns:
        df = df.drop(columns=[id_column])

    date_cols = [c for c in df.columns if c.startswith("Date ")]
    if not date_cols:
        raise ValueError("No phenotype columns found with prefix 'Date '.")

    counts = []
    for col in date_cols:
        # A case is any non-empty, non-NaN cell; actual date value is ignored.
        series = df[col]
        mask = series.notna() & series.astype(str).str.strip().ne("")
        counts.append((col, int(mask.sum())))

    return pd.DataFrame(counts, columns=["phenotype", "case_count"])


def main() -> None:
    args = parse_args()
    input_path = Path(args.input).expanduser().resolve()
    output_dir = Path(args.output_dir).expanduser().resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    output_path = output_dir / args.output_name

    df = load_dataframe(input_path)
    case_counts = compute_case_counts(df, args.id_column)
    case_counts.to_csv(output_path, index=False)
    print(f"Saved case counts to {output_path}")


if __name__ == "__main__":
    main()
