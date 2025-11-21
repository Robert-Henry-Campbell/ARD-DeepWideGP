#!/usr/bin/env python3
"""
map_and_join_mvp.py (revised)

Behavior changes from previous version:
 - Do NOT concatenate Trait. Instead, for mapped rows (mapped_flag==True)
   perform a many-to-many join: if MVP has multiple rows for a PHECODE_,
   the mapped row will be duplicated once per matching MVP row.
 - Preserve all original columns from both mapped and MVP (MVP columns are
   prefixed with "MVP_" to avoid name collisions).
 - Unmappable rows (mapped_flag==False) are preserved untouched; they receive
   MVP_* columns filled with NA.
 - Logging records counts of dropped MVP rows without a PHECODE_, number of duplicate keys,
   the expansion factor, and final match counts.

Inputs/outputs same as before.
"""

import os
import pickle
import datetime
import pandas as pd
import numpy as np

# -----------------------
# Config / paths
# -----------------------
MAPPED_IN = "data/2-phecode_mapping/ARD_GBD_ICD10_phecode_mapped_nov_19_2025.csv"
MVP_META_IN = "data/1-inputs/MVP_R4.1000G_AGR.DataDictionary_and_Counts_20240719.csv"
LOG_PICKLE = "data/2-phecode_mapping/full_mapping_log.pkl"

MAPPED_OUT = "data/2-phecode_mapping/ARD_GBD_ICD10_phecode_mapped_with_MVP.csv"
LOG_SUMMARY_CSV = "data/2-phecode_mapping/mapping_log_summary.csv"
GBD_LVL1_OUT = "data/2-phecode_mapping/mvp_join_summary_level1.csv"
GBD_LVL2_OUT = "data/2-phecode_mapping/mvp_join_summary_level2.csv"
FULL_LOG_OUT = LOG_PICKLE  # overwrite the existing pickle with updated logs

# -----------------------
# Utilities: logging helper
# -----------------------
def load_or_init_log(pickle_path):
    if os.path.exists(pickle_path):
        with open(pickle_path, "rb") as f:
            try:
                log = pickle.load(f)
                if not isinstance(log, list):
                    log = [log]
            except Exception:
                log = []
    else:
        log = []
    return log

def append_log_entry(log_list, event, description, count=None, table=None):
    rec = {
        "timestamp": datetime.datetime.utcnow().isoformat() + "Z",
        "event": event,
        "description": description,
        "count": int(count) if (count is not None and not pd.isna(count)) else None
    }
    if table is not None:
        try:
            rec["table"] = table.copy()
        except Exception:
            rec["table"] = table
    log_list.append(rec)
    return rec

# -----------------------
# Step 1: load inputs
# -----------------------
if not os.path.exists(MAPPED_IN):
    raise FileNotFoundError(f"Mapped input file not found: {MAPPED_IN}")
if not os.path.exists(MVP_META_IN):
    raise FileNotFoundError(f"MVP metadata file not found: {MVP_META_IN}")

mapped = pd.read_csv(MAPPED_IN, dtype=str)
mvp = pd.read_csv(MVP_META_IN, dtype=str)

log_records = load_or_init_log(LOG_PICKLE)
append_log_entry(log_records, "load_inputs", f"Loaded inputs. mapped rows={mapped.shape[0]}, MVP rows={mvp.shape[0]}", count=mapped.shape[0])

# -----------------------
# Step 2: ensure PHECODE_ in mapped exists
# -----------------------
if 'PHECODE_' not in mapped.columns:
    if 'PHECODE' in mapped.columns:
        mapped['PHECODE_'] = mapped['PHECODE'].astype(str).str.replace('.', '_', regex=False)
        append_log_entry(log_records, "created_PHECODE_from_PHECODE", "Created PHECODE_ in mapped from PHECODE", count=int(mapped['PHECODE_'].notna().sum()))
    else:
        mapped['PHECODE_'] = None
        append_log_entry(log_records, "created_PHECODE_missing_source", "PHECODE not present; created empty PHECODE_ in mapped", count=0)
else:
    append_log_entry(log_records, "PHECODE__exists", "PHECODE_ already present in mapped", count=int(mapped['PHECODE_'].notna().sum()))

# -----------------------
# Step 3: create PHECODE_ in MVP (split Trait)
# -----------------------
if 'Trait' not in mvp.columns:
    raise KeyError("MVP metadata file does not contain expected column 'Trait'")

def extract_phe_from_trait(trait):
    if pd.isna(trait):
        return None
    s = str(trait)
    parts = s.split('_', 1)
    return parts[1] if len(parts) > 1 else None

mvp['PHECODE_'] = mvp['Trait'].apply(extract_phe_from_trait)
append_log_entry(log_records, "mvp_create_PHECODE_", f"Created PHECODE_ in MVP by splitting Trait; non-null count={int(mvp['PHECODE_'].notna().sum())}", count=int(mvp['PHECODE_'].notna().sum()))

# -----------------------
# Step 4: join only for mapped rows; many-to-many expansion permitted
# -----------------------

# Ensure mapped_flag exists and is boolean
if 'mapped_flag' not in mapped.columns:
    mapped['mapped_flag'] = mapped['PHECODE'].notna()
    append_log_entry(log_records, "mapped_flag_created_from_PHECODE", "Created mapped_flag from presence of PHECODE", count=int(mapped['mapped_flag'].sum()))

_truth_map = {
    'true': True, 't': True, '1': True, 'yes': True, 'y': True,
    'false': False, 'f': False, '0': False, 'no': False, 'n': False,
    'nan': False, 'none': False, '': False
}

def normalize_flag(val):
    if isinstance(val, bool):
        return val
    if pd.isna(val):
        return False
    s = str(val).strip().lower()
    return _truth_map.get(s, False)

mapped['mapped_flag'] = mapped['mapped_flag'].apply(normalize_flag).astype(bool)
append_log_entry(log_records, "mapped_flag_normalized", "Normalized mapped_flag to boolean", count=int(mapped['mapped_flag'].sum()))

# Split mapped
mapped_mappable = mapped[mapped['mapped_flag']].copy()
mapped_unmappable = mapped[~mapped['mapped_flag']].copy()
append_log_entry(log_records, "split_mappable_unmappable", f"Split mapped into mappable={mapped_mappable.shape[0]} and unmappable={mapped_unmappable.shape[0]} rows.", count=int(mapped_mappable.shape[0]))

# Normalize PHECODE_ fields and treat 'nan'/'none'/' ' as missing
def clean_key_series(s):
    s2 = s.astype(str).str.strip()
    s2.loc[s2.str.lower().isin(['nan', 'none', ''])] = np.nan
    return s2

mapped_mappable['PHECODE_'] = clean_key_series(mapped_mappable['PHECODE_'])
mvp['PHECODE_'] = clean_key_series(mvp['PHECODE_'])

# Drop MVP rows with no key (they cannot be matched on PHECODE_)
n_mvp_before = mvp.shape[0]
mvp_valid = mvp[mvp['PHECODE_'].notna()].copy()
n_mvp_after = mvp_valid.shape[0]
append_log_entry(log_records, "mvp_drop_no_key", f"Dropped MVP rows without PHECODE_: before={n_mvp_before}, after={n_mvp_after}", count=(n_mvp_before - n_mvp_after))

# Diagnostic: duplicates on MVP PHECODE_ keys
dupe_counts = mvp_valid['PHECODE_'].value_counts()
n_keys_with_dupes = int((dupe_counts > 1).sum())
n_rows_in_duplicate_keys = int(dupe_counts[dupe_counts > 1].sum()) if (dupe_counts > 1).any() else 0
append_log_entry(log_records, "mvp_dup_stats", f"MVP keys with duplicates={n_keys_with_dupes}, rows_in_those_keys={n_rows_in_duplicate_keys}", count=n_rows_in_duplicate_keys)

# Perform many-to-many join: preserve all mapped_mappable columns and bring in full MVP rows (prefixed)
# Use add_prefix to avoid name collisions; keep all MVP columns
mvp_prefixed = mvp_valid.add_prefix("MVP_")
# The join key on right is "MVP_PHECODE_"
merged_mappable = mapped_mappable.merge(mvp_prefixed, left_on="PHECODE_", right_on="MVP_PHECODE_", how="left")

append_log_entry(log_records, "mvp_join_mappable", f"Performed many-to-many join; pre-join mappable rows={mapped_mappable.shape[0]}, post-join rows={merged_mappable.shape[0]}", count=int(merged_mappable.shape[0]))

# Recombine with unmappable rows: ensure unmappable has the same MVP_* columns (set to NA)
merged_unmappable = mapped_unmappable.copy()
mvp_cols = [c for c in merged_mappable.columns if c.startswith("MVP_")]
for c in mvp_cols:
    if c not in merged_unmappable.columns:
        merged_unmappable[c] = None

# Now concatenate (mapped rows that matched multiple MVP rows are expanded)
merged = pd.concat([merged_mappable, merged_unmappable], axis=0, ignore_index=True)

# Compute expansion factor: ratio of post-join mappable rows to original mappable rows
expansion_factor = merged_mappable.shape[0] / (mapped_mappable.shape[0] if mapped_mappable.shape[0] > 0 else 1)
append_log_entry(log_records, "mvp_join_expansion", f"Expansion factor (post-join mappable rows / pre-join mappable rows) = {expansion_factor:.3f}", count=int(merged_mappable.shape[0]))

# Mark mvp_match_flag by presence of explicit MVP key column
merged['mvp_match_flag'] = merged['MVP_PHECODE_'].notna() if 'MVP_PHECODE_' in merged.columns else False
n_matched = int(merged['mvp_match_flag'].sum())
n_unmatched = int(merged.shape[0] - n_matched)
append_log_entry(log_records, "mvp_match_final_counts", f"MVP matches final={n_matched}, unmatched={n_unmatched}", count=n_matched)

# -----------------------
# Step 5: compute summaries by GBD level 1 and level 2
# -----------------------
def gbd_match_summary(df, gbd_col):
    if gbd_col not in df.columns:
        return None
    total = df.groupby(gbd_col).size().rename("n_total").reset_index()
    matched = df[df['mvp_match_flag']].groupby(gbd_col).size().rename("n_matched").reset_index()
    unm = df[~df['mvp_match_flag']].groupby(gbd_col).size().rename("n_unmatched").reset_index()
    summary = total.merge(matched, on=gbd_col, how='left').merge(unm, on=gbd_col, how='left').fillna(0)
    summary['n_total'] = summary['n_total'].astype(int)
    summary['n_matched'] = summary['n_matched'].astype(int)
    summary['n_unmatched'] = summary['n_unmatched'].astype(int)
    summary['pct_matched'] = (summary['n_matched'] / summary['n_total']).round(4)
    summary['pct_unmatched'] = (summary['n_unmatched'] / summary['n_total']).round(4)
    return summary

lvl1_summary = gbd_match_summary(merged, "cause_level_1")
lvl2_summary = gbd_match_summary(merged, "cause_level_2")

if lvl1_summary is not None:
    append_log_entry(log_records, "mvp_gbd_level1_summary", "GBD level 1 summary for MVP join", count=int(lvl1_summary.shape[0]), table=lvl1_summary)
    lvl1_summary.to_csv(GBD_LVL1_OUT, index=False)

if lvl2_summary is not None:
    append_log_entry(log_records, "mvp_gbd_level2_summary", "GBD level 2 summary for MVP join", count=int(lvl2_summary.shape[0]), table=lvl2_summary)
    lvl2_summary.to_csv(GBD_LVL2_OUT, index=False)

# -----------------------
# Step 6: persist merged dataframe and update summary CSVs / pickle
# -----------------------
os.makedirs(os.path.dirname(MAPPED_OUT), exist_ok=True)
merged.to_csv(MAPPED_OUT, index=False)
append_log_entry(log_records, "saved_merged_output", f"Saved merged file to {MAPPED_OUT}", count=int(merged.shape[0]))

# Update simple summary CSV
summary_rows = []
for rec in log_records:
    summary_rows.append({
        "timestamp": rec.get("timestamp"),
        "event": rec.get("event"),
        "description": rec.get("description"),
        "count": rec.get("count")
    })
summary_df = pd.DataFrame(summary_rows)
summary_df.to_csv(LOG_SUMMARY_CSV, index=False)

# Save full log
with open(FULL_LOG_OUT, "wb") as f:
    pickle.dump(log_records, f)

# Export embedded log tables
for rec in log_records:
    if rec.get("table") is not None:
        ev = rec.get("event", "table")
        fname = f"data/2-phecode_mapping/log_table_{ev}.csv"
        try:
            rec['table'].to_csv(fname, index=False)
        except Exception:
            try:
                pd.DataFrame(rec['table']).to_csv(fname, index=False)
            except Exception:
                pass

append_log_entry(log_records, "complete", "Completed MVP integration and persisted outputs", count=int(merged.shape[0]))

with open(FULL_LOG_OUT, "wb") as f:
    pickle.dump(log_records, f)

# final prints
print("MVP join completed.")
print(f"Mapped rows (input): {mapped.shape[0]}")
print(f"Mappable rows (input): {mapped_mappable.shape[0] if 'mapped_mappable' in locals() else 'NA'}")
print(f"Merged rows (output): {merged.shape[0]}, MVP matches: {n_matched}, unmatched: {n_unmatched}")
if lvl1_summary is not None:
    print(f"GBD level1 groups: {lvl1_summary.shape[0]} (CSV: {GBD_LVL1_OUT})")
if lvl2_summary is not None:
    print(f"GBD level2 groups: {lvl2_summary.shape[0]} (CSV: {GBD_LVL2_OUT})")
print(f"Full audit log saved to: {FULL_LOG_OUT} and summary saved to: {LOG_SUMMARY_CSV}")
