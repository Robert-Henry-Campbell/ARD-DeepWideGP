import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path
import numpy as np

"""
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
"""

# -------------------------------------------------------------------
# Setup paths
# -------------------------------------------------------------------
base_dir = Path("data/2-phecode_mapping")
log_path = "data\\4-integrate_MVP_UKB\\2025_casecount\\MVP_UKB_mapping_log.txt"
plot_path = "data\\4-integrate_MVP_UKB\\MVP_UKB_surviving_cause_level2_plot.png"
counts_csv_path = "data\\4-integrate_MVP_UKB\\MVP_UKB_surviving_cause_level2_counts.csv"

# 1. Load ARD data
mvp_ard = pd.read_csv(base_dir / "ARD_GBD_ICD10_phecode_mapped_with_MVP.csv")
UKB_ard = pd.read_csv("data/3-identify UKB case count for ARDs/firstocc_ards_clean.csv")

# Initialize log text
log_lines = []
log_lines.append("=== MVP–UKB ARD Mapping Log ===\n")

# -------------------------------------------------------------------
# Log initial sizes
# -------------------------------------------------------------------
log_lines.append(f"Initial MVP ARD rows: {mvp_ard.shape[0]}")
log_lines.append(f"Initial UKB ARD rows: {UKB_ard.shape[0]}\n")

# -------------------------------------------------------------------
# 2. Prefix UKB columns
# -------------------------------------------------------------------
UKB_ard = UKB_ard.add_prefix("UKB_")
UKB_ard = UKB_ard.rename(columns={"UKB_ICD10_explo": "ICD10_explo"})

log_lines.append(f"After prefixing UKB columns: {UKB_ard.shape[0]} rows")
log_lines.append("Columns now begin with 'UKB_' except for ICD10_explo.\n")

# -------------------------------------------------------------------
# 3. Merge MVP + UKB by ICD10 code
# -------------------------------------------------------------------
merged_ard = pd.merge(mvp_ard, UKB_ard, on="ICD10_explo", how="inner")
log_lines.append(f"Rows after MVP–UKB ICD10 merge (inner join): {merged_ard.shape[0]}")

kept_icd10 = set(merged_ard["ICD10_explo"])
mvp_dropped = len(set(mvp_ard["ICD10_explo"]) - kept_icd10)
ukb_dropped = len(set(UKB_ard["ICD10_explo"]) - kept_icd10)

log_lines.append(f"MVP ICD10 codes retained: {len(kept_icd10)}")
log_lines.append(f"MVP ICD10 codes dropped: {mvp_dropped}")
log_lines.append(f"UKB ICD10 codes retained: {len(kept_icd10)}")
log_lines.append(f"UKB ICD10 codes dropped: {ukb_dropped}\n")

# -------------------------------------------------------------------
# 4. Apply filtering criteria
# -------------------------------------------------------------------
before_filter_rows = merged_ard.shape[0]
filtered_ard = merged_ard[
    (merged_ard["mvp_match_flag"] == True) &
    (merged_ard["UKB_n_cases_EUR"] > 10000)
]

after_filter_rows = filtered_ard.shape[0]

log_lines.append(f"Rows before filtering: {before_filter_rows}")
log_lines.append(f"Rows after filtering (mvp_match_flag==True & UKB_n_cases_EUR>10000): {after_filter_rows}")
log_lines.append(f"Rows dropped by filter: {before_filter_rows - after_filter_rows}\n")

# -------------------------------------------------------------------
# 5. Count surviving rows by level-2 GBD cause
# -------------------------------------------------------------------
if "cause_level_2" in filtered_ard.columns:
    cause_counts = filtered_ard["cause_level_2"].value_counts().sort_index()
    cause_counts.to_csv(counts_csv_path, header=["count"])
    log_lines.append("Surviving row counts by cause_level_2:")
    for k, v in cause_counts.items():
        log_lines.append(f"  {k}: {v}")
    log_lines.append("")
else:
    log_lines.append("Column 'cause_level_2' not found — cannot compute category counts.\n")

# -------------------------------------------------------------------
# Save filtered dataset (unchanged business logic)
# -------------------------------------------------------------------
output_filtered_path = Path("data/4-integrate_MVP_UKB/mapped_ards_MVP_GWAS_and_UKB_10kcase.csv")
filtered_ard.to_csv(output_filtered_path, index=False)

log_lines.append(f"Filtered dataset saved to: {output_filtered_path}\n")

# -------------------------------------------------------------------
# Save log file
# -------------------------------------------------------------------
with open(log_path, "w") as f:
    f.write("\n".join(log_lines))

# -------------------------------------------------------------------
# Plot surviving rows by cause_level_2
# -------------------------------------------------------------------
if after_filter_rows > 0 and "cause_level_2" in filtered_ard.columns:

    df_plot = filtered_ard.copy()

    # ---------------------------------------------------------------
    # REQUIRED NEW STEP: drop duplicates in ICD10_explo
    # ---------------------------------------------------------------
    df_plot = df_plot.drop_duplicates(subset="ICD10_explo", keep="first")

    df_plot.sort_values(["cause_level_2", "UKB_n_cases_EUR"], ascending=[True, False], inplace=True)
    df_plot = df_plot.reset_index(drop=True)

    n = len(df_plot)
    x = np.arange(n)
    y = df_plot["UKB_n_cases_EUR"].values

    fig, ax = plt.subplots(figsize=(7.2, 4.0), dpi=300)

    ax.bar(
        x, y,
        width=0.8,
        edgecolor="black",
        linewidth=0.35
    )

    ax.axhline(10000, color="red", linewidth=1)

    ax.set_ylabel("# of cases")
    ax.set_title("age-related disease UKB 10k case & MVP discovery GWAS")
    ax.set_xticks([])

    group_labels = []
    group_positions = []
    curr_group = None
    start_index = 0

    causes = df_plot["cause_level_2"].values

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

    log_lines.append(f"Plot saved to: {plot_path}")

    with open(log_path, "w") as f:
        f.write("\n".join(log_lines))
