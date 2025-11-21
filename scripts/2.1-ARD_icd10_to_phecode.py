import pandas as pd
import os

"""
This script appends the phecode mapping to the ARD list from the GBD,
adds cleaned PHECODE_ column, and includes detailed logging
of mapping success/failure overall and by GBD cause levels.
"""

# ==============================
# Load Inputs
# ==============================
ARD_icd10 = pd.read_csv('data/1-inputs/ARD_GBD_ICD10_exploded_20250830_214953.csv')
phecode_map = pd.read_csv('data/1-inputs/Phecode_map_v1_2_icd10_beta.csv')
MVP_meta = pd.read_csv('data/1-inputs/MVP_R4.1000G_AGR.DataDictionary_and_Counts_20240719.csv')

# ==============================
# Merge
# ==============================
mapped = pd.merge(
    ARD_icd10,
    phecode_map,
    left_on='ICD10_explo',
    right_on='ICD10',
    how='left'
)

# ==============================
# Add new PHECODE_ column
# ==============================
if 'PHECODE' in mapped.columns:
    mapped['PHECODE_'] = mapped['PHECODE'].astype(str).str.replace('.', '_', regex=False)
else:
    mapped['PHECODE_'] = None

# ==============================
# Logging Framework
# ==============================
log_records = []

def log_event(event, description, df=None):
    """
    Append a log entry with event name, description, and optional counts DataFrame.
    """
    count = df.shape[0] if isinstance(df, pd.DataFrame) else None
    log_records.append({
        'event': event,
        'description': description,
        'count': count
    })

# Basic counts
log_event("initial_ARD_rows", "Number of rows in ARD_icd10 before merge", ARD_icd10)
log_event("after_merge", "Rows after merge with phecode map", mapped)

# Count duplicates
duplicate_rows = mapped[mapped.duplicated()]
log_event("duplicate_rows", "Duplicate rows found", duplicate_rows)

# Drop duplicates
mapped = mapped.drop_duplicates()
log_event("after_dropping_duplicates", "Rows after removing duplicates", mapped)

# Mapping success/failure
mapped['mapped_flag'] = mapped['PHECODE'].notna()

mapped_success = mapped[mapped['mapped_flag']]
mapped_failure = mapped[~mapped['mapped_flag']]

log_event("mapped_success", "Rows successfully mapped to a PHECODE", mapped_success)
log_event("mapped_failure", "Rows not mapped to any PHECODE", mapped_failure)

# ==============================
# GBD Category Logging
# ==============================

def summarize_by_gbd(label, col):
    """
    Summarize mapped/unmapped counts by a GBD category column.
    """
    # total by category
    total = mapped.groupby(col).size().reset_index(name='n_total')

    # mapped counts
    m_succ = mapped_success.groupby(col).size().reset_index(name='n_mapped')

    # unmapped counts
    m_fail = mapped_failure.groupby(col).size().reset_index(name='n_unmapped')

    # merge summaries
    summary = (
        total
        .merge(m_succ, on=col, how='left')
        .merge(m_fail, on=col, how='left')
        .fillna(0)
    )

    summary['pct_mapped'] = summary['n_mapped'] / summary['n_total']
    summary['pct_unmapped'] = summary['n_unmapped'] / summary['n_total']

    log_records.append({
        'event': f"gbd_summary_{label}",
        'description': f"Mapping summary by {col}",
        'count': None,
        'table': summary
    })

    return summary


# Level 1 summary
if 'cause_level_1' in mapped.columns:
    gbd_lvl1_summary = summarize_by_gbd("level_1", "cause_level_1")
else:
    gbd_lvl1_summary = None

# Level 2 summary
if 'cause_level_2' in mapped.columns:
    gbd_lvl2_summary = summarize_by_gbd("level_2", "cause_level_2")
else:
    gbd_lvl2_summary = None

# ==============================
# Produce Log Output
# ==============================
log_df = pd.DataFrame([
    {k: v for k, v in rec.items() if k != 'table'}  # tables stored separately
    for rec in log_records
])

print("\n========== SUMMARY LOG ==========")
print(log_df)

if gbd_lvl1_summary is not None:
    print("\n========== GBD Level 1 Summary ==========")
    print(gbd_lvl1_summary)

if gbd_lvl2_summary is not None:
    print("\n========== GBD Level 2 Summary ==========")
    print(gbd_lvl2_summary)

# ==============================
# Save outputs
# ==============================
mapped.to_csv('data/2-phecode_mapping/ARD_GBD_ICD10_phecode_mapped_nov_19_2025.csv', index=False)
log_df.to_csv('data/2-phecode_mapping/mapping_log_summary.csv', index=False)

if gbd_lvl1_summary is not None:
    gbd_lvl1_summary.to_csv('data/2-phecode_mapping/gbd_level1_summary.csv', index=False)

if gbd_lvl2_summary is not None:
    gbd_lvl2_summary.to_csv('data/2-phecode_mapping/gbd_level2_summary.csv', index=False)

# Store the full log including embedded tables
# (optional, but good for auditability)
import pickle
with open('data/2-phecode_mapping/full_mapping_log.pkl', 'wb') as f:
    pickle.dump(log_records, f)

