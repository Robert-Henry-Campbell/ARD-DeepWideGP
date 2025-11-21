# make_ard_gbd_plot_styled.py
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path
import numpy as np

# ---- CONFIG ----
csv_path = Path("data/2-phecode_mapping/ARD_GBD_ICD10_phecode_mapped_with_MVP.csv")
out_png = Path("data/2-phecode_mapping/MVP_ARD_casecount.png")
out_points_csv = Path("data/2-phecode_mapping/MVP_ARD_casecount_points.csv")
min_cases = 5000
hline_y = 10000
figsize = (7.2, 4.0)   # Nature double-column approx 7.2 in width
dpi = 300
bar_width = 0.8        # <1.0 provides slight separation between bars
bar_edgewidth = 0.35   # thin outline to make bars distinct
# ----------------

out_png.parent.mkdir(parents=True, exist_ok=True)

# Read CSV
df = pd.read_csv(csv_path)

# Validate required columns
required_cols = {"MVP_num_cases.EUR", "cause_level_2"}
missing = required_cols - set(df.columns)
if missing:
    raise ValueError(f"Missing required columns in CSV: {missing}")

# Filter and sort
df_filtered = df[df["MVP_num_cases.EUR"] >= min_cases].copy()
df_filtered["original_index"] = df_filtered.index
df_filtered.sort_values(["cause_level_2", "MVP_num_cases.EUR"], ascending=[True, False], inplace=True)
df_filtered.reset_index(drop=True, inplace=True)

# Plot data
n = len(df_filtered)
x = np.arange(n)
y = df_filtered["MVP_num_cases.EUR"].values

fig, ax = plt.subplots(figsize=figsize, dpi=dpi)

# Bars with slight separation (width < 1) and thin black outline for clarity
bars = ax.bar(x, y, width=bar_width, align="center",
              edgecolor="black", linewidth=bar_edgewidth)

# Remove per-row x ticks / labels
ax.set_xticks([])

# Horizontal red line at y = hline_y
ax.axhline(hline_y, color="red", linewidth=1.2)

# Axis labels / title (updated as requested)
ax.set_ylabel("# of cases")
ax.set_xlabel("")   # no per-row x label
ax.set_title("MVP age-related disease EUR case count")

# Compute group centers and set them as xticks/xticklabels
group_positions = []
group_labels = []
current_label = None
start_idx = 0

labels_array = df_filtered["cause_level_2"].values
for idx, label in enumerate(labels_array):
    if label != current_label:
        if current_label is not None:
            center = (start_idx + idx - 1) / 2.0
            group_positions.append(center)
            group_labels.append(current_label)
        current_label = label
        start_idx = idx
# finalize last group
if current_label is not None:
    center = (start_idx + n - 1) / 2.0
    group_positions.append(center)
    group_labels.append(current_label)

# Place group names as xtick labels at their centers
ax.set_xticks(group_positions)
ax.set_xticklabels(group_labels, rotation=45, ha="right", fontsize=8)

# Expand x limits slightly for aesthetics
ax.set_xlim(-0.5, n - 0.5)

# Improve spacing so labels are not clipped
plt.tight_layout()
plt.savefig(out_png, bbox_inches="tight", dpi=dpi)
plt.close(fig)

# Save plotted points CSV (ordered)
points_df = pd.DataFrame({
    "plot_index": x,
    "MVP_num_cases.EUR": y,
    "cause_level_2": df_filtered["cause_level_2"].values,
    "original_index": df_filtered["original_index"].values
})
points_df.to_csv(out_points_csv, index=False)

print(f"Saved figure to: {out_png}")
print(f"Saved plot-points CSV to: {out_points_csv}")
print(f"Plotted {n} rows across {len(group_labels)} groups.")
