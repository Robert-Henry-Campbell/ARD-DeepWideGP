#!/usr/bin/env python3
"""
CLI: map MVP ARDs to UKB ARDs via ICD10, filter, log, and plot results.


goal: 
1.load up the ARD data for MVP
2. load up the ARD data for UKB
3. map the MVP phenotypes to UKB phenotypes using ICD10 codes
4. filter for rows with existing MVP GWAS + existing UKB GWAS
5. filter for rows where UKB GWAS N > 10000

Additional tasks:
- Add logging of row counts at each step
- Log rows kept vs. dropped
- Compute surviving row counts by cause_level_2
- Save all logs + a bar chart of surviving categories to the SAME directory as the MVP ARD CSV

inputs:
mvp-ard: ARD GBD ICD10 phecode_map MVP_GWAS availability CSV. 
    unique ICD10_y (individual diagnoses)
    duplicates of ICD10_explo (grouped diagnoses)

ukb-ard: GBD ARD UKB 2025 ICD10 case counts CSV.
    - note that the GBD ARD is not strictly speaking needed as you have it in mvp-ard. you doubled your work somwhere. 

Example:
python scripts/one_offs/4.1-map_UKB_MVP.py \
  --mvp-ard data/2-phecode_mapping/2.2-map_and_join_MVP/ARD_GBD_ICD10_phecode_mapped_with_MVP.csv \
  --ukb-ard data/3-identify-UKB-casecount-ARDs/extract_casecount_ukb2025/GBD_ARD_UKB_2025_ICD0_casecounts.csv \
  --output-dir data/4-integrate_MVP_UKB/2025_casecount \
  --ukb-case-threshold 10000
"""

import argparse
from pathlib import Path
from typing import Iterable, List, Tuple

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Merge MVP and UKB ARD data by ICD10, filter on availability and case "
            "counts, and emit filtered CSV + logs + plots."
        )
    )
    parser.add_argument(
        "--mvp-ard",
        default="data/2-phecode_mapping/ARD_GBD_ICD10_phecode_mapped_with_MVP.csv",
        help="Path to MVP ARD ICD10 dataset (default: %(default)s).",
    )
    parser.add_argument(
        "--ukb-ard",
        default="data/3-identify-UKB-casecount-ARDs/task-extract_casecounts_neale/firstocc_ards_clean.csv",
        help="Path to UKB ARD ICD10 dataset (default: %(default)s).",
    )
    parser.add_argument(
        "--output-dir",
        help=(
            "Directory for outputs (CSV, log, counts, plot). "
            "Defaults to the MVP ARD file's directory."
        ),
    )
    parser.add_argument(
        "--output-prefix",
        default="MVP_UKB",
        help="Prefix for output files (default: %(default)s).",
    )
    parser.add_argument(
        "--ukb-case-threshold",
        type=int,
        default=10000,
        help="Keep rows with UKB case counts above this threshold (default: %(default)s).",
    )
    parser.add_argument(
        "--mvp-flag-col",
        default="mvp_match_flag",
        help="Column marking MVP GWAS availability (default: %(default)s).",
    )
    parser.add_argument(
        "--ukb-case-col",
        default="UKB_case_count",
        help="Column holding UKB case counts (default: %(default)s).",
    )
    parser.add_argument(
        "--cause-col",
        default="cause_level_2",
        help="Column used for cause-level counts/plotting (default: %(default)s).",
    )
    parser.add_argument(
        "--no-plot",
        action="store_true",
        help="Skip plotting surviving rows by cause level.",
    )
    return parser.parse_args()


def require_columns(df: pd.DataFrame, cols: Iterable[str], label: str) -> None:
    missing = [c for c in cols if c not in df.columns]
    if missing:
        raise ValueError(f"{label} missing columns: {', '.join(missing)}")


def to_bool_series(series: pd.Series) -> pd.Series:
    return series.apply(lambda x: str(x).strip().lower() in {"1", "true", "t", "yes"})


def prefix_ukb_columns(df: pd.DataFrame) -> pd.DataFrame:
    df_prefixed = df.add_prefix("UKB_")
    return df_prefixed.rename(columns={"UKB_ICD10_explo": "ICD10_explo"})


def filter_rows(
    df: pd.DataFrame, flag_col: str, case_col: str, threshold: int
) -> Tuple[pd.DataFrame, List[str]]:
    logs: List[str] = []
    require_columns(df, [flag_col, case_col], "Merged dataframe")

    df = df.copy()
    df[case_col] = pd.to_numeric(df[case_col], errors="coerce")
    coerced = df[case_col].isna().sum()
    if coerced:
        logs.append(f"{case_col}: {coerced} non-numeric values coerced to NaN before filtering.")

    flags = to_bool_series(df[flag_col])
    before = len(df)
    filtered = df[flags & (df[case_col] > threshold)]
    after = len(filtered)
    logs.append(
        f"Rows after filtering ({flag_col}==True & {case_col}>{threshold}): {after} "
        f"(dropped {before - after} of {before})."
    )
    return filtered, logs


def plot_cause_level(
    df: pd.DataFrame,
    case_col: str,
    cause_col: str,
    plot_path: Path,
    threshold: int,
) -> None:
    df_plot = df.drop_duplicates(subset="ICD10_explo", keep="first")
    df_plot = df_plot.sort_values([cause_col, case_col], ascending=[True, False]).reset_index(drop=True)

    n = len(df_plot)
    x = np.arange(n)
    y = df_plot[case_col].values

    fig, ax = plt.subplots(figsize=(7.2, 4.0), dpi=300)
    ax.bar(x, y, width=0.8, edgecolor="black", linewidth=0.35)
    ax.axhline(threshold, color="red", linewidth=1)

    ax.set_ylabel("# of cases")
    ax.set_title("age-related disease UKB case threshold & MVP discovery GWAS")
    ax.set_xticks([])

    group_labels = []
    group_positions = []
    curr_group = None
    start_index = 0
    causes = df_plot[cause_col].values

    for i, c in enumerate(causes):
        if c != curr_group:
            if curr_group is not None:
                center = (start_index + i - 1) / 2
                group_positions.append(center)
                group_labels.append(curr_group)
            curr_group = c
            start_index = i

    if curr_group is not None:
        center = (start_index + n - 1) / 2
        group_positions.append(center)
        group_labels.append(curr_group)

    ax.set_xticks(group_positions)
    ax.set_xticklabels(group_labels, rotation=45, ha="right", fontsize=8)

    plt.tight_layout()
    plt.savefig(plot_path, dpi=300)
    plt.close()


def main() -> None:
    args = parse_args()

    mvp_path = Path(args.mvp_ard).expanduser().resolve()
    ukb_path = Path(args.ukb_ard).expanduser().resolve()
    output_dir = Path(args.output_dir).expanduser().resolve() if args.output_dir else mvp_path.parent
    output_dir.mkdir(parents=True, exist_ok=True)

    prefix = args.output_prefix
    log_path = output_dir / f"{prefix}_mapping_log.txt"
    plot_path = output_dir / f"{prefix}_surviving_cause_level2_plot.png"
    counts_csv_path = output_dir / f"{prefix}_surviving_cause_level2_counts_uniqueICD10_category.csv"
    filtered_csv_path = output_dir / f"{prefix}_mapped_ards_MVP_GWAS_and_UKB_{args.ukb_case_threshold}case.csv"

    log_lines: List[str] = [
        "=== MVP–UKB ARD Mapping Log ===",
        f"MVP ARD input: {mvp_path}",
        f"UKB ARD input: {ukb_path}",
        f"Output directory: {output_dir}",
        f"UKB case threshold: {args.ukb_case_threshold}",
        "",
    ]

    mvp_df = pd.read_csv(mvp_path)
    ukb_df = pd.read_csv(ukb_path)
    log_lines.append(f"Initial MVP ARD rows: {len(mvp_df)}")
    log_lines.append(f"Initial UKB ARD rows: {len(ukb_df)}\n")

    require_columns(mvp_df, ["ICD10_explo"], "MVP dataframe")
    require_columns(ukb_df, ["ICD10_explo"], "UKB dataframe")

    ukb_prefixed = prefix_ukb_columns(ukb_df)
    log_lines.append(f"After prefixing UKB columns: {len(ukb_prefixed)} rows")
    log_lines.append("Columns now begin with 'UKB_' except for ICD10_explo.\n")

    merged = pd.merge(
        mvp_df, 
        ukb_prefixed, 
        on="ICD10_explo", 
        how="inner",
        validate="many_to_one") #testing the validation not sure it's right
        
    log_lines.append(f"Rows after MVP–UKB ICD10 merge (inner join): {len(merged)}")

    kept_icd10 = set(merged["ICD10_explo"])
    mvp_dropped = len(set(mvp_df["ICD10_explo"]) - kept_icd10)
    ukb_dropped = len(set(ukb_prefixed["ICD10_explo"]) - kept_icd10)
    log_lines.append(f"MVP ICD10 codes retained: {len(kept_icd10)}")
    log_lines.append(f"MVP ICD10 codes dropped: {mvp_dropped}")
    log_lines.append(f"UKB ICD10 codes retained: {len(kept_icd10)}")
    log_lines.append(f"UKB ICD10 codes dropped: {ukb_dropped}\n")

    filtered, filter_logs = filter_rows(merged, args.mvp_flag_col, args.ukb_case_col, args.ukb_case_threshold)
    log_lines.extend(filter_logs)
    log_lines.append("")

    if args.cause_col in filtered.columns:
        cause_counts = (
            filtered.drop_duplicates(subset="ICD10_explo", keep="first")[args.cause_col]
            .value_counts()
            .sort_index()
        )
        cause_counts.to_csv(counts_csv_path, header=["count"])
        log_lines.append("Surviving row counts by cause_level_2:")
        for k, v in cause_counts.items():
            log_lines.append(f"  {k}: {v}")
        log_lines.append(f"Saved counts to: {counts_csv_path}")
    else:
        log_lines.append(f"Column '{args.cause_col}' not found — cannot compute category counts.")
    log_lines.append("")

    filtered.to_csv(filtered_csv_path, index=False)
    log_lines.append(f"Filtered dataset saved to: {filtered_csv_path}")

    if not args.no_plot and len(filtered) > 0 and args.cause_col in filtered.columns:
        plot_cause_level(filtered, args.ukb_case_col, args.cause_col, plot_path, args.ukb_case_threshold)
        log_lines.append(f"Plot saved to: {plot_path}")
    elif args.no_plot:
        log_lines.append("Plot skipped (--no-plot).")
    else:
        log_lines.append("Plot skipped (no rows or cause column missing).")

    with open(log_path, "w", encoding="utf-8") as f:
        f.write("\n".join(log_lines))
    print(f"Done. Log written to {log_path}")


if __name__ == "__main__":
    main()
