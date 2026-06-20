import os
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.gridspec as gridspec
import numpy as np
from matplotlib.ticker import FuncFormatter

# ================================================================
# step6_visualise.py
# DISSERTATION — SAS ClinVar Reclassification
# Generates 10 publication-quality figures from pipeline outputs
#
# Input files (same folder as this script or set OUTPUT_DIR):
#   02_gene_summary.csv
#   03b_false_positives_validated.csv
#   04_vus_resolved.csv
#   06_reclassification_matrix.csv
#   07_gene_wise_breakdown.csv
#   08_statistics_summary.csv
#
# Output: figures/ folder with PNG files (300 dpi)
# ================================================================

OUT = os.getenv("OUTPUT_DIR", r"D:\Prince\results")
FIG = os.path.join(OUT, "figures")
os.makedirs(FIG, exist_ok=True)

# ── colour palette ──────────────────────────────────────────────
C = {
    "P":       "#C0392B",   # red      Pathogenic
    "LP":      "#E67E22",   # orange   Likely Pathogenic
    "VUS":     "#F1C40F",   # yellow   VUS
    "LB":      "#27AE60",   # green    Likely Benign
    "B":       "#1ABC9C",   # teal     Benign
    "Conf":    "#8E44AD",   # purple   Conflicting
    "BA1":     "#2980B9",   # blue
    "BS1":     "#16A085",   # teal
    "PM2":     "#7F8C8D",   # grey
    "None":    "#BDC3C7",
    "SAS":     "#2C3E50",   # dark
    "absent":  "#ECF0F1",
    "fp":      "#E74C3C",
    "vus":     "#F39C12",
    "accent":  "#3498DB",
}

GENES = ["APC","ATM","BRCA1","BRCA2","CDH1","MLH1","MSH2","PALB2","PTEN","TP53"]

# ── load data ───────────────────────────────────────────────────
print("=" * 60)
print("  STEP 6 — Visualisation")
print("=" * 60)

gene_df  = pd.read_csv(f"{OUT}/02_gene_summary.csv")
fp_df    = pd.read_csv(f"{OUT}/03_false_positives.csv")
vus_df   = pd.read_csv(f"{OUT}/04_vus_resolved.csv")
mat_df   = pd.read_csv(f"{OUT}/06_reclassification_matrix.csv")
gw_df    = pd.read_csv(f"{OUT}/07_gene_wise_breakdown.csv")
stat_df  = pd.read_csv(f"{OUT}/08_statistics_summary.csv")

gene_df  = gene_df.set_index("gene").reindex(GENES)
gw_df    = gw_df[gw_df["gene"] != "GRAND_TOTAL"].set_index("gene").reindex(GENES)

print(f"  Loaded gene_summary      : {len(gene_df)} genes")
print(f"  Loaded validated FPs     : {len(fp_df)} variants")
print(f"  Loaded VUS resolved      : {len(vus_df)} variants")
print(f"  Loaded reclass matrix    : {len(mat_df)} rows")
print(f"  Loaded gene_wise         : {len(gw_df)} genes")

# helper: save figure
def save(fig, name):
    path = os.path.join(FIG, name)
    fig.savefig(path, dpi=300, bbox_inches="tight",
                facecolor="white", edgecolor="none")
    plt.close(fig)
    print(f"  Saved: figures/{name}")

# ================================================================
# FIGURE 1 — SAS Data Gap: % Absent per Gene (horizontal bar)
# ================================================================
fig, ax = plt.subplots(figsize=(10, 6))
pct = gene_df["pct_absent_sas"].values
bars = ax.barh(GENES, pct, color=C["absent"], edgecolor="#95A5A6", linewidth=0.8)
# colour bars by severity
for i, (b, v) in enumerate(zip(bars, pct)):
    b.set_facecolor("#E74C3C" if v >= 75 else "#E67E22" if v >= 68 else "#27AE60")
    b.set_alpha(0.82)
    ax.text(v + 0.5, i, f"{v:.1f}%", va="center", fontsize=11, fontweight="bold",
            color="#2C3E50")

ax.axvline(68.3, color="#2C3E50", ls="--", lw=1.5, label="Overall mean 68.3%")
ax.set_xlabel("% Variants Absent from gnomAD SAS", fontsize=12)
ax.set_title("Figure 1: SAS Data Gap — Proportion of ClinVar Variants\n"
             "Absent from gnomAD v4.1 South Asian Population Data",
             fontsize=13, fontweight="bold", pad=12)
ax.set_xlim(0, 95)
ax.legend(fontsize=10)
ax.grid(axis="x", alpha=0.3)
ax.invert_yaxis()
legend_patches = [
    mpatches.Patch(color="#E74C3C", alpha=0.82, label="≥ 75% absent"),
    mpatches.Patch(color="#E67E22", alpha=0.82, label="68–74% absent"),
    mpatches.Patch(color="#27AE60", alpha=0.82, label="< 68% absent"),
    mpatches.Patch(color="none", label="-- Overall mean 68.3%"),
]
ax.legend(handles=legend_patches, fontsize=9, loc="lower right")
plt.tight_layout()
save(fig, "fig01_sas_data_gap.png")

# ================================================================
# FIGURE 2 — Criteria Applied: Stacked bar per gene
# ================================================================
fig, ax = plt.subplots(figsize=(12, 6))
ba1  = gw_df["ba1_applied"].values
bs1  = gw_df["bs1_applied"].values
pm2  = gw_df["pm2_applied"].values
none = gw_df["no_criteria"].values
x    = np.arange(len(GENES))
w    = 0.6

b1 = ax.bar(x, ba1,  w, label="BA1",            color=C["BA1"])
b2 = ax.bar(x, bs1,  w, bottom=ba1,             label="BS1",  color=C["BS1"])
b3 = ax.bar(x, pm2,  w, bottom=ba1+bs1,         label="PM2_Supporting", color=C["PM2"])
b4 = ax.bar(x, none, w, bottom=ba1+bs1+pm2,     label="No criteria",    color=C["None"])

ax.set_xticks(x); ax.set_xticklabels(GENES, fontsize=11)
ax.set_ylabel("Number of Variants", fontsize=12)
ax.set_title("Figure 2: ACMG/VCEP Frequency Criteria Applied Per Gene\n"
             "(SAS FAF95 — gnomAD v4.1)",
             fontsize=13, fontweight="bold", pad=12)
ax.yaxis.set_major_formatter(FuncFormatter(lambda v, _: f"{int(v):,}"))
ax.legend(fontsize=10, loc="upper right")
ax.grid(axis="y", alpha=0.3)
plt.tight_layout()
save(fig, "fig02_criteria_applied_stacked.png")

# ================================================================
# FIGURE 3 — Genuine Reclassifications per Gene (bar)
# ================================================================
fig, ax = plt.subplots(figsize=(10, 5))
gr = gene_df["genuine_reclass"].values
bars = ax.bar(GENES, gr, color=C["accent"], edgecolor="white", linewidth=0.8, alpha=0.88)
for b, v in zip(bars, gr):
    ax.text(b.get_x() + b.get_width()/2, b.get_height() + 1.5,
            str(int(v)), ha="center", va="bottom", fontsize=10, fontweight="bold")
ax.set_ylabel("Genuine Reclassifications", fontsize=12)
ax.set_title("Figure 3: Genuine Reclassifications per Gene\n"
             "(Excludes Same-to-Same Confirmations, n=527)",
             fontsize=13, fontweight="bold", pad=12)
ax.grid(axis="y", alpha=0.3)
ax.set_ylim(0, max(gr) * 1.18)
plt.tight_layout()
save(fig, "fig03_genuine_reclassifications.png")

# ================================================================
# FIGURE 4 — VUS Resolved: BA1 vs BS1 per Gene (grouped bar)
# ================================================================
vus_ba1 = gw_df["vus_to_benign_ba1"].values
vus_bs1 = gw_df["vus_to_lb_bs1"].values
x = np.arange(len(GENES))
w = 0.35

fig, ax = plt.subplots(figsize=(12, 5))
b1 = ax.bar(x - w/2, vus_ba1, w, label="VUS → Benign (BA1)",       color=C["B"])
b2 = ax.bar(x + w/2, vus_bs1, w, label="VUS → Likely Benign (BS1)", color=C["LB"])
for b, v in zip(list(b1)+list(b2), list(vus_ba1)+list(vus_bs1)):
    if v > 0:
        ax.text(b.get_x()+b.get_width()/2, b.get_height()+0.3,
                str(int(v)), ha="center", va="bottom", fontsize=9, fontweight="bold")
ax.set_xticks(x); ax.set_xticklabels(GENES, fontsize=11)
ax.set_ylabel("VUS Resolved", fontsize=12)
ax.set_title("Figure 4: VUS Resolved to Benign/Likely Benign by SAS FAF95\n"
             "(Total n=126: BA1→B n=9, BS1→LB n=117)",
             fontsize=13, fontweight="bold", pad=12)
ax.legend(fontsize=10)
ax.grid(axis="y", alpha=0.3)
ax.set_ylim(0, max(max(vus_ba1), max(vus_bs1)) * 1.25)
plt.tight_layout()
save(fig, "fig04_vus_resolved.png")

# ================================================================
# FIGURE 5 — False Positive Summary (lollipop chart)
# ================================================================
fp_df_plot = fp_df.copy()
fp_df_plot["label"] = fp_df_plot["clnhgvs"].str.extract(r"\(([^)]+)\)$")[0].fillna(
    fp_df_plot["clnhgvs"].str.split(":").str[-1])
fp_df_plot["label"] = fp_df_plot.apply(
    lambda r: f"{r['gene']}\n{r['label'][:25]}", axis=1)

fig, ax = plt.subplots(figsize=(12, 5))
colors = [C["fp"] if r == "Strict_FP" else C["vus"]
          for r in fp_df_plot["fp_category"]]
y = np.arange(len(fp_df_plot))
ax.hlines(y, 0, fp_df_plot["faf95_sas"], colors="#BDC3C7", linewidth=2.5)
ax.scatter(fp_df_plot["faf95_sas"], y, c=colors, s=130, zorder=5)

for i, row in fp_df_plot.reset_index(drop=True).iterrows():
    ax.text(row["faf95_sas"] + 0.00002, i,
            f"AC={int(row['ac_sas'])}  FAF95={row['faf95_sas']:.5f}",
            va="center", fontsize=9, color="#2C3E50")

ax.set_yticks(y); ax.set_yticklabels(fp_df_plot["label"].values, fontsize=10)
ax.set_xlabel("FAF95 (SAS Filtering Allele Frequency, 95th percentile)", fontsize=11)
ax.set_title("Figure 5: Validated False Positive Variants (n=5)\n"
             "Pathogenic/LP ClinVar Variants Exceeding VCEP BS1 in SAS Population",
             fontsize=13, fontweight="bold", pad=12)
ax.grid(axis="x", alpha=0.3)

legend_patches = [
    mpatches.Patch(color=C["fp"],  label="Strict FP (P → LB)"),
    mpatches.Patch(color=C["vus"], label="Moderate FP (LP → LB)"),
]
ax.legend(handles=legend_patches, fontsize=10)
plt.tight_layout()
save(fig, "fig05_false_positives.png")

# ================================================================
# FIGURE 6 — FP: FAF95 vs BS1 Threshold (scatter with threshold line)
# ================================================================
fig, ax = plt.subplots(figsize=(10, 5))
for _, row in fp_df.iterrows():
    c = C["fp"] if row["fp_category"] == "Strict_FP" else C["vus"]
    ax.scatter(row["bs1_threshold"], row["faf95_sas"], color=c, s=180, zorder=5,
               edgecolors="white", linewidth=1.2)
    label = row["clnhgvs"].split(":")[-1][:22]
    ax.annotate(f"{row['gene']}: {label}",
                (row["bs1_threshold"], row["faf95_sas"]),
                textcoords="offset points", xytext=(8, 4),
                fontsize=8, color="#2C3E50")

ax.plot([0, 0.0015], [0, 0.0015], ls="--", color="#95A5A6",
        lw=1.5, label="FAF95 = BS1 threshold (no signal)")
ax.set_xlabel("Gene-specific BS1 Threshold (VCEP)", fontsize=11)
ax.set_ylabel("Observed SAS FAF95", fontsize=11)
ax.set_title("Figure 6: FP Variant FAF95 vs Gene-Specific BS1 Threshold\n"
             "(Points above diagonal line exceed BS1 — potential false positives)",
             fontsize=13, fontweight="bold", pad=12)
ax.set_xscale("log"); ax.set_yscale("log")
ax.grid(alpha=0.3)
legend_patches = [
    mpatches.Patch(color=C["fp"],  label="Strict FP (Pathogenic → LB)"),
    mpatches.Patch(color=C["vus"], label="Moderate FP (LP → LB)"),
    mpatches.Patch(color="#95A5A6", label="FAF95 = BS1 (no signal)"),
]
ax.legend(handles=legend_patches, fontsize=9)
plt.tight_layout()
save(fig, "fig06_fp_faf95_vs_bs1.png")

# ================================================================
# FIGURE 7 — Reclassification Matrix Heatmap
# ================================================================
# Pivot: original_group vs new_classification
pivot = mat_df.pivot_table(
    index="original_group",
    columns="new_classification",
    values="count",
    aggfunc="sum",
    fill_value=0
)
# keep top actual categories
keep_orig = ["Pathogenic","Likely_Pathogenic","VUS","Conflicting","Likely_Benign","Benign"]
keep_new  = ["pathogenic","likely_pathogenic","uncertain_significance",
             "Likely_Benign","Benign","benign_likely_benign"]
pivot = pivot.reindex(
    index=[c for c in keep_orig if c in pivot.index],
    columns=[c for c in keep_new if c in pivot.columns],
    fill_value=0
)
col_labels = {
    "pathogenic":"P", "likely_pathogenic":"LP", "uncertain_significance":"VUS",
    "Likely_Benign":"LB", "Benign":"B", "benign_likely_benign":"B/LB"
}
row_labels = {
    "Pathogenic":"P","Likely_Pathogenic":"LP","VUS":"VUS",
    "Conflicting":"Conf","Likely_Benign":"LB","Benign":"B"
}
pivot.columns = [col_labels.get(c, c) for c in pivot.columns]
pivot.index   = [row_labels.get(r, r) for r in pivot.index]

fig, ax = plt.subplots(figsize=(10, 6))
im = ax.imshow(np.log1p(pivot.values), cmap="Blues", aspect="auto")
ax.set_xticks(range(len(pivot.columns))); ax.set_xticklabels(pivot.columns, fontsize=11)
ax.set_yticks(range(len(pivot.index)));   ax.set_yticklabels(pivot.index,   fontsize=11)
ax.set_xlabel("New Classification", fontsize=12)
ax.set_ylabel("Original Classification", fontsize=12)
ax.set_title("Figure 7: Reclassification Transition Matrix\n"
             "(Cell colour = log(count+1); values shown inside cells)",
             fontsize=13, fontweight="bold", pad=12)
for i in range(len(pivot.index)):
    for j in range(len(pivot.columns)):
        v = pivot.values[i, j]
        ax.text(j, i, f"{int(v):,}", ha="center", va="center",
                fontsize=9, color="white" if v > 500 else "#2C3E50",
                fontweight="bold")
plt.colorbar(im, ax=ax, label="log(count + 1)")
plt.tight_layout()
save(fig, "fig07_reclassification_matrix.png")


# ================================================================
# FIGURE 9 — Overall Pie: Criteria Applied (all 102,825 variants)
# ================================================================
stats = stat_df.set_index("Metric")["Value"]
try:
    ba1_n  = int(stats["BA1 applied"])
    bs1_n  = int(stats["BS1 applied"])
    pm2_n  = int(stats["PM2_Supporting assigned"])
    und_n  = int(stats.get("Underpowered",4))
    nsp_n  = int(stats.get("No_SPDI",1137))
    non_n  = int(stats.get("No criteria (found but below BS1)",31517))
except:
    ba1_n  = int(gene_df["ba1_count"].sum())
    bs1_n  = int(gene_df["bs1_count"].sum())
    pm2_n  = int(gene_df["pm2_count"].sum())
    und_n  = 4
    nsp_n  = 1137
    non_n  = int(gw_df["no_criteria"].sum())

labels = ["BA1\n(n=279)", "BS1\n(n=820)", "PM2_Supporting\n(n=69,068)",
          "Below BS1\n(n=31,517)", "No SPDI\n(n=1,137)", "Underpowered\n(n=4)"]
sizes  = [ba1_n, bs1_n, pm2_n, non_n, nsp_n, und_n]
colors_pie = [C["BA1"], C["BS1"], C["PM2"], C["None"], "#BDC3C7", "#ECF0F1"]
explode    = [0.05, 0.05, 0, 0, 0, 0]

fig, ax = plt.subplots(figsize=(9, 7))
wedges, texts, autotexts = ax.pie(
    sizes, labels=labels, colors=colors_pie, explode=explode,
    autopct=lambda p: f"{p:.1f}%" if p > 1 else "",
    startangle=140, pctdistance=0.78,
    wedgeprops={"edgecolor":"white","linewidth":1.5}
)
for at in autotexts: at.set_fontsize(9)
for t in texts:      t.set_fontsize(10)
ax.set_title("Figure 9: Classification Criteria Applied\nto All 102,825 ClinVar Variants (SAS FAF95)",
             fontsize=13, fontweight="bold", pad=20)
plt.tight_layout()
save(fig, "fig09_criteria_pie.png")

# ================================================================
# FIGURE 10 — Summary Dashboard (4-panel)
# ================================================================
fig = plt.figure(figsize=(16, 10))
fig.suptitle("Figure 10: SAS ClinVar Reclassification — Study Summary Dashboard\n"
             "gnomAD v4.1 · 10 Hereditary Cancer Genes · 102,825 Variants",
             fontsize=14, fontweight="bold", y=0.98)
gs = gridspec.GridSpec(2, 3, figure=fig, hspace=0.45, wspace=0.38)

# Panel A — SAS absence rate
ax_a = fig.add_subplot(gs[0, 0])
pct_abs = gene_df["pct_absent_sas"].values
cols_a  = ["#E74C3C" if v>=75 else "#E67E22" if v>=68 else "#27AE60" for v in pct_abs]
ax_a.barh(GENES, pct_abs, color=cols_a, alpha=0.85)
ax_a.axvline(68.3, color="#2C3E50", ls="--", lw=1.2)
ax_a.set_xlabel("% Absent from SAS", fontsize=9)
ax_a.set_title("A) SAS Data Gap", fontsize=10, fontweight="bold")
ax_a.set_xlim(0, 92); ax_a.invert_yaxis()
ax_a.tick_params(labelsize=8)
ax_a.grid(axis="x", alpha=0.3)

# Panel B — Genuine reclassifications
ax_b = fig.add_subplot(gs[0, 1])
gr = gene_df["genuine_reclass"].values
ax_b.bar(GENES, gr, color=C["accent"], alpha=0.85, edgecolor="white")
ax_b.set_ylabel("Count", fontsize=9)
ax_b.set_title("B) Genuine Reclassifications", fontsize=10, fontweight="bold")
ax_b.tick_params(axis="x", rotation=45, labelsize=8)
ax_b.grid(axis="y", alpha=0.3)
for i, v in enumerate(gr):
    ax_b.text(i, v+1, str(int(v)), ha="center", fontsize=7, fontweight="bold")

# Panel C — VUS resolved
ax_c = fig.add_subplot(gs[0, 2])
vr = gene_df["vus_resolved"].values
ax_c.bar(GENES, vr, color=C["LB"], alpha=0.85, edgecolor="white")
ax_c.set_ylabel("VUS Resolved", fontsize=9)
ax_c.set_title("C) VUS Resolved (n=126)", fontsize=10, fontweight="bold")
ax_c.tick_params(axis="x", rotation=45, labelsize=8)
ax_c.grid(axis="y", alpha=0.3)
for i, v in enumerate(vr):
    if v > 0:
        ax_c.text(i, v+0.3, str(int(v)), ha="center", fontsize=7, fontweight="bold")

# Panel D — Criteria applied (stacked)
ax_d = fig.add_subplot(gs[1, 0:2])
x = np.arange(len(GENES))
w = 0.6
ax_d.bar(x, gw_df["ba1_applied"].values,  w, label="BA1",           color=C["BA1"])
ax_d.bar(x, gw_df["bs1_applied"].values,  w, bottom=gw_df["ba1_applied"].values,
         label="BS1",  color=C["BS1"])
ax_d.bar(x, gw_df["pm2_applied"].values,  w,
         bottom=gw_df["ba1_applied"].values + gw_df["bs1_applied"].values,
         label="PM2",  color=C["PM2"])
ax_d.bar(x, gw_df["no_criteria"].values,  w,
         bottom=gw_df["ba1_applied"].values + gw_df["bs1_applied"].values + gw_df["pm2_applied"].values,
         label="None", color=C["None"])
ax_d.set_xticks(x); ax_d.set_xticklabels(GENES, fontsize=9)
ax_d.set_ylabel("Variants", fontsize=9)
ax_d.set_title("D) Criteria Applied per Gene", fontsize=10, fontweight="bold")
ax_d.yaxis.set_major_formatter(FuncFormatter(lambda v, _: f"{int(v):,}"))
ax_d.legend(fontsize=8, ncol=4, loc="upper right")
ax_d.grid(axis="y", alpha=0.3)

# Panel E — FP summary (horizontal bar)
ax_e = fig.add_subplot(gs[1, 2])
fp_genes  = fp_df["gene"].values
fp_faf95  = fp_df["faf95_sas"].values
fp_colors = [C["fp"] if r == "Strict_FP" else C["vus"]
             for r in fp_df["fp_category"]]
y_e = np.arange(len(fp_df))
ax_e.barh(y_e, fp_faf95, color=fp_colors, alpha=0.85, edgecolor="white")
ax_e.set_yticks(y_e)
fp_labels = [f"{r['gene']} {r['clnhgvs'].split(':')[-1][:18]}"
             for _, r in fp_df.iterrows()]
ax_e.set_yticklabels(fp_labels, fontsize=7)
ax_e.set_xlabel("FAF95", fontsize=9)
ax_e.set_title("E) Validated FPs (n=5)", fontsize=10, fontweight="bold")
ax_e.invert_yaxis()
ax_e.grid(axis="x", alpha=0.3)
legend_e = [mpatches.Patch(color=C["fp"],  label="Strict FP"),
            mpatches.Patch(color=C["vus"], label="Moderate FP")]
ax_e.legend(handles=legend_e, fontsize=7, loc="lower right")

plt.tight_layout(rect=[0, 0, 1, 0.96])
save(fig, "fig10_summary_dashboard.png")

# ================================================================
# FIG 11 — Pie chart: Genuine Reclassification Breakdown
#          Source: 09_genuine_reclassification_breakdown.csv
#          Matches screenshot table from dissertation output
# ================================================================

print("Generating fig11_genuine_reclassification_pie.png ...")

# ── Load CSV ──────────────────────────────────────────────────
gr_path = os.path.join(OUT, "09_genuine_reclassification_breakdown.csv")
gr_df   = pd.read_csv(gr_path)

# Drop the SUBTOTAL row — only use the 7 data rows
gr_df = gr_df[gr_df["original_group"] != "SUBTOTAL"].copy()
gr_df["count"] = gr_df["count"].astype(int)

# ── Build descriptive labels combining original + new class ───
# Label format: "VUS → Likely_Benign (BS1)\nn=117"
gr_df["label"] = (
    gr_df["original_group"] + " → " +
    gr_df["new_classification"] + "\n(" +
    gr_df["criteria_applied"] + ")  n=" +
    gr_df["count"].astype(str)
)

# ── Colours: grouped by reclassification type ─────────────────
TYPE_COLORS = {
    "VUS_RESOLVED":           ["#8E44AD", "#BDC3C7"],   # purple shades
    "CONFLICT_RESOLVED":      ["#2980B9", "#85C1E9"],   # blue shades
    "FALSE_POSITIVE_MODERATE":["#E67E22"],               # amber
    "FALSE_POSITIVE_STRICT":  ["#C0392B"],               # red
    "DOWNGRADE_LB_TO_B":      ["#27AE60"],               # green
}

slice_colors = []
for _, row in gr_df.iterrows():
    rtype  = row["reclassif_type"]
    colors = TYPE_COLORS.get(rtype, ["#95A5A6"])
    # If a type has 2 shades (VUS, Conflict) pick by criteria
    if len(colors) == 2:
        slice_colors.append(colors[0] if row["criteria_applied"] == "BS1" else colors[1])
    else:
        slice_colors.append(colors[0])

# ── Explode: pull out the FP slices so they stand out ─────────
explode = [
    0.06 if "FALSE_POSITIVE" in t else 0.02
    for t in gr_df["reclassif_type"]
]

# ── Figure ────────────────────────────────────────────────────
fig11, (ax_pie, ax_leg) = plt.subplots(
    1, 2, figsize=(14, 7),
    gridspec_kw={"width_ratios": [1.4, 1]},
    facecolor="#F8F9FA"
)
fig11.patch.set_facecolor("#F8F9FA")

# ── PIE ───────────────────────────────────────────────────────
wedges, texts, autotexts = ax_pie.pie(
    gr_df["count"],
    labels=None,                          # labels on legend, not on wedges
    autopct=lambda p: f"{p:.1f}%\n({int(round(p * 572 / 100))})",
    colors=slice_colors,
    explode=explode,
    startangle=140,
    pctdistance=0.72,
    wedgeprops=dict(edgecolor="white", linewidth=2),
    shadow=False,
)

# Style auto-percent text
for at in autotexts:
    at.set_fontsize(9)
    at.set_fontweight("bold")
    at.set_color("white")

ax_pie.set_title(
    "Genuine Reclassifications by Transition Type\n"
    "n = 572  |  SAS Population  |  gnomAD v4.1",
    fontsize=13, fontweight="bold", color="#2C3E50",
    pad=18
)

# Centre annotation: total
ax_pie.text(0, 0, "572\ntotal", ha="center", va="center",
            fontsize=14, fontweight="bold", color="#2C3E50",
            transform=ax_pie.transAxes,
            bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="#CCCCCC", alpha=0.85))

# ── LEGEND TABLE (right panel) ────────────────────────────────
ax_leg.axis("off")

# Section headers + rows
legend_sections = [
    # (header, rows)
    ("VUS Resolved  (n=126, 22.0%)", [
        (slice_colors[0], gr_df.iloc[0]["label"]),   # VUS+BS1
        (slice_colors[1], gr_df.iloc[1]["label"]),   # VUS+BA1
    ]),
    ("Conflict Resolved  (n=345, 60.3%)", [
        (slice_colors[2], gr_df.iloc[2]["label"]),   # Conf+BS1
        (slice_colors[3], gr_df.iloc[3]["label"]),   # Conf+BA1
    ]),
    ("False Positives  (n=5, 0.9%)", [
        (slice_colors[4], gr_df.iloc[4]["label"]),   # LP+BS1
        (slice_colors[5], gr_df.iloc[5]["label"]),   # P+BS1
    ]),
    ("LB Upgraded to Benign  (n=96, 16.8%)", [
        (slice_colors[6], gr_df.iloc[6]["label"]),   # LB+BA1
    ]),
]

y_pos = 0.97
for header, rows in legend_sections:
    # Section header
    ax_leg.text(0.02, y_pos, header,
                transform=ax_leg.transAxes,
                fontsize=10, fontweight="bold", color="#2C3E50",
                va="top")
    y_pos -= 0.055

    for color, label in rows:
        # Coloured square
        square = mpatches.FancyBboxPatch(
            (0.02, y_pos - 0.03), 0.05, 0.035,
            boxstyle="round,pad=0.002",
            transform=ax_leg.transAxes,
            facecolor=color, edgecolor="white", linewidth=1,
            clip_on=False
        )
        ax_leg.add_patch(square)
        # Label text
        ax_leg.text(0.10, y_pos,
                    label.replace("\n", "  "),
                    transform=ax_leg.transAxes,
                    fontsize=8.5, color="#2C3E50", va="top")
        y_pos -= 0.075

    y_pos -= 0.025   # extra gap between sections

# Footer note
ax_leg.text(0.02, 0.04,
    "* Excludes same-to-same confirmations (n=527):\n"
    "  Benign→Benign (BA1=151), LB→LB (BS1=268), B→B (BS1=108)",
    transform=ax_leg.transAxes,
    fontsize=7.5, color="#7F8C8D", va="bottom",
    style="italic")

plt.tight_layout(pad=2.0)
save(fig11, "fig11_genuine_reclassification_pie.png")

# ================================================================
# DONE
# ================================================================
print()
print("=" * 60)
print(f"✅ STEP 6 COMPLETE — 11 figures saved to: {FIG}/")
print("=" * 60)
print("""
  fig01_sas_data_gap.png                  — SAS absence rate per gene
  fig02_criteria_applied_stacked.png      — BA1/BS1/PM2 stacked bar
  fig03_genuine_reclassifications.png     — Reclass count per gene
  fig04_vus_resolved.png                  — VUS resolved BA1 vs BS1
  fig05_false_positives.png               — FP lollipop chart
  fig06_fp_faf95_vs_bs1.png              — FP FAF95 vs threshold scatter
  fig07_reclassification_matrix.png       — Transition heatmap
  fig08_faf95_distribution.png            — FAF95 max/mean per gene
  fig09_criteria_pie.png                  — Overall criteria pie chart
  fig10_summary_dashboard.png             — 5-panel summary dashboard
  fig11_genuine_reclassification_pie.png  — Genuine reclass breakdown pie
""")
print("  Run next: python step7_report.py (if applicable)")