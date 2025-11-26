#!/usr/bin/env python3
"""
script: merge_GBD_ARD_UKB (2025)
description: connects the GBD ARD data with the UKB casecounts via ICD10 codes.
inputs:
1. ICD10 pheno casecount df: UKB_2025_casecount.csv (from 3.1-create_UKB_casecount.py)
2. full df of ARD ICD10 codes from GBD: bothsex_ARD_GBD_ICD10_exploded_20250911_201039.csv (from project 1)
outputs:
1. UKB ICD10 pheno case count df merged with GBD ARD metadata

example usage:
python scripts/3.2-merge_GBD_ARD_UKB.py \
  --casecounts data/3-identify-UKB-casecount-ARDs/extract_casecount_ukb2025/UKB_2025_casecount.csv \
  --ard-icd10 data/3-identify-UKB-casecount-ARDs/bothsex_ARD_GBD_ICD10_exploded_20250911_201039.csv \
  --output-dir data/3-identify-UKB-casecount-ARDs/extract_casecount_ukb2025 \
  --output-name GBD_ARD_UKB_2025_ICD0_casecounts.csv
"""

import argparse
from pathlib import Path
from typing import Optional

import pandas as pd


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Merge GBD ARD ICD10 codes with UKB ICD10 case counts produced by "
            "3.1-create_UKB_casecount.py."
        )
    )
    parser.add_argument(
        "--casecounts",
        required=True,
        help="Path to the UKB ICD10 case count CSV produced by 3.1-create_UKB_casecount.py.",
    )
    parser.add_argument(
        "--ard-icd10",
        required=True,
        help=(
            "Path to the GBD ARD ICD10 exploded file "
            "(e.g., bothsex_ARD_GBD_ICD10_exploded_20250911_201039.csv)."
        ),
    )
    parser.add_argument(
        "--output-dir",
        required=True,
        help="Directory where the merged dataframe will be written.",
    )
    parser.add_argument(
        "--output-name",
        default="UKB_GBD_ARD_casecounts.csv",
        help="Filename for the merged dataframe (default: UKB_GBD_ARD_casecounts.csv).",
    )
    parser.add_argument(
        "--ard-icd10-column",
        default="ICD10_explo",
        help="Column in the ARD dataframe that holds the ICD10 code (default: ICD10_explo).",
    )
    parser.add_argument(
        "--phenotype-column",
        default="phenotype",
        help="Column in the casecount dataframe that holds the UKB phenotype string (default: phenotype).",
    )
    return parser.parse_args()


def extract_icd10_from_phenotype(phenotype: str) -> Optional[str]:
    """
    UKB casecount phenotypes look like: 'Date A00 first reported (cholera)'.
    This helper returns the ICD10 token between 'Date' and 'first reported'.
    """
    if not isinstance(phenotype, str):
        return None

    text = phenotype.strip()
    if text.lower().startswith("date "):
        text = text[5:]

    marker = " first reported"
    if marker in text.lower():
        # Split once, case-insensitive
        lower_text = text.lower()
        idx = lower_text.index(marker)
        code_part = text[:idx]
    else:
        code_part = text

    code = code_part.strip().replace("(", "").replace(")", "")
    return code or None


def load_dataframe(path: Path) -> pd.DataFrame:
    if path.suffix.lower() in {".parquet", ".pq"}:
        return pd.read_parquet(path)
    return pd.read_csv(path)


def prepare_ard_dataframe(df: pd.DataFrame, icd10_col: str) -> pd.DataFrame:
    if icd10_col not in df.columns:
        raise ValueError(f"ICD10 column '{icd10_col}' not found in ARD dataframe.")
    before = len(df)
    df = df.dropna(subset=[icd10_col])
    df = df.drop_duplicates(subset=[icd10_col])
    after = len(df)
    print(f"ARD dataframe: dropped {before - after} duplicate/null ICD10 rows; remaining {after}.")
    return df


def prepare_casecount_dataframe(df: pd.DataFrame, phenotype_col: str) -> pd.DataFrame:
    required_cols = {phenotype_col, "case_count"}
    missing = required_cols - set(df.columns)
    if missing:
        raise ValueError(f"Casecount dataframe missing required columns: {', '.join(sorted(missing))}")

    df = df.rename(columns={phenotype_col: "phenotype"})
    df["ICD10_explo"] = df["phenotype"].apply(extract_icd10_from_phenotype)
    before = len(df)
    df = df.dropna(subset=["ICD10_explo"])
    after = len(df)
    df["case_count"] = df["case_count"].astype(int)
    print(f"Casecount dataframe: extracted ICD10 for {after}/{before} rows.")
    return df[["ICD10_explo", "phenotype", "case_count"]]


def main() -> None:
    args = parse_args()

    casecount_path = Path(args.casecounts).expanduser().resolve()
    ard_path = Path(args.ard_icd10).expanduser().resolve()
    output_dir = Path(args.output_dir).expanduser().resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    output_path = output_dir / args.output_name

    print(f"Loading UKB casecounts from: {casecount_path}")
    casecounts_df = load_dataframe(casecount_path)
    casecounts_df = prepare_casecount_dataframe(casecounts_df, args.phenotype_column)

    print(f"Loading ARD ICD10 dataframe from: {ard_path}")
    ard_df = load_dataframe(ard_path)
    ard_df = prepare_ard_dataframe(ard_df, args.ard_icd10_column)

    merged_df = pd.merge(ard_df, casecounts_df, on="ICD10_explo", how="inner")
    print(f"Merged rows: {len(merged_df)} (inner join on ICD10_explo).")

    merged_df.to_csv(output_path, index=False)
    print(f"Saved merged dataframe to {output_path}")


if __name__ == "__main__":
    main()
