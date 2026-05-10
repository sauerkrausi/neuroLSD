import csv
import os

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
CSV_FILE = os.path.join(SCRIPT_DIR, "raw_readme.csv")
OUT_FILE = os.path.join(SCRIPT_DIR, "raw_readme.csv")

# Code info per script (from readmeCodeetc/)
CODE_INFO = {
    "diff132_d50_nDIA.Rmd": (
        "nDIA-based quantification of whole-cell iN and iDA neurons at day 50. "
        "Performs pairwise KO vs Ctrl comparisons and neuron-type interaction effects. "
        "Integrates organelle annotations and correlates KO profiles by organelle enrichment. "
        "Generates cluster heatmaps, GO term enrichments, and cross-cell-type correlation plots."
    ),
    "diff118_iNiDA_Theval.Rmd": (
        "TMTpro-based differential proteomic analysis of day 23 iN and iDA samples. "
        "Compares Ctrl, SMPD1-/-, and ASAH1-/- lines across neuronal types. "
        "Outputs log2FC, significance tables, and cluster-based GO term enrichments."
    ),
    "diff118_iNiDA_d23.Rmd": (
        "TMTpro-based differential proteomic analysis of day 23 iN and iDA samples. "
        "Compares Ctrl, SMPD1-/-, and ASAH1-/- lines across neuronal types. "
        "Outputs log2FC, significance tables, and cluster-based GO term enrichments."
    ),
    "diff136_iNd35_ctrl_asah1e1_axonalproteome.Rmd": (
        "TMTpro-based proteome of iNeuron whole-cell, soma, and projection fractions at day 35. "
        "Evaluates compartment-specific KO effects using log2FC and RoR metrics. "
        "Correlates subcellular enrichment with KO-induced localization shifts. "
        "Visualizes via heatmaps, barplots, violin plots, and spatial proteomic rainfall plots."
    ),
    "HeLa_Ctrl-ASAH1_LysoIP.Rmd": (
        "Lysosome-enriched fraction proteomics from HeLa Ctrl vs ASAH1-/-. "
        "Focus on lysosomal hydrolase abundance, lysosomal-endosomal pathway components, and autophagy markers. "
        "Includes cluster-based heatmaps, scaled violin plots, and annotation-wise log2FC summaries."
    ),
    "run_cellpose_calcium_remoteHDD.py": (
        "Segment neuronal soma in grayscale calcium imaging stacks (.tiff or .nd2) using Cellpose-SAM "
        "with a custom cpsam model. Supports batch mode with automated ROI mask generation and QC overlays. "
        "Saves masks, overlays, and log files for quality control."
    ),
    "neuroLSD_CalciumMaster.Rmd": (
        "Load segmentation masks and extract per-ROI fluorescence traces. "
        "Compute dF/F normalization, apply adaptive thresholding, and binarize events. "
        "Quantify neuronal activity: spike timing, synchrony, coactivity, and burst frequency."
    ),
    "th_neuron_quantification_interative.py": (
        "Quantify TH-positive neurons in iN and iDA cultures via Python-based segmentation. "
        "Supports automated and interactive workflows. "
        "Provides condition-wise comparisons and group-level analyses."
    ),
    "neuro_syn_env.py": (
        "Parse split-channel TIFF filenames to extract sample prefix, channel, and wavelength. "
        "Preprocess DNA and organelle channels with background subtraction and white tophat filtering. "
        "Segment nuclei with Cellpose-SAM and organelles via Otsu thresholding. "
        "Compute per-image metrics and colocalization fractions."
    ),
    "neuron_endo_syn_eval.Rmd": (
        "R analysis: reads per-image summary CSVs, parses sample metadata, "
        "performs outlier filtering and Wilcoxon tests, and produces violin/box/jitter plots "
        "with BH-adjusted p-value labels."
    ),
}

# Instrument info per method_name
INSTRUMENT_INFO = {
    "Spinning-disk Microscopy": (
        "Yokogawa CSU-X1 spinning disk confocal mounted on a Nikon Eclipse Ti2-E motorized microscope "
        "equipped with a Tokai Hit stage-top incubator. Conditions maintained at 37 \u00b0C, 5% CO\u2082, "
        "and 95% humidity. Nikon Plan Apo 60\u00d7/1.40 N.A. immersion oil objective lens. "
        "Fluorophores sequentially excited with a Nikon LUN-F XL solid-state laser combiner "
        "(405 nm \u2013 80 mW, 488 nm \u2013 80 mW, 561 nm \u2013 65 mW, 640 nm \u2013 60 mW) "
        "with a Semrock Di01-T405/488/568/647 dichroic mirror. "
        "Emissions collected through Chroma ET455/50m [405 nm], ET525/36m [488 nm], "
        "ET605/52m [561 nm], ET700/75m [640 nm] filters. "
        "Hamamatsu ORCA-Fusion BT CMOS camera (6.5 \u00b5m\u00b2 photodiode, 16-bit) and NIS-Elements software."
    ),
}

# Dataset descriptions (from supplemental data legends)
DATASET_INFO = {
    "Dataset S1.xlsx": (
        "Generation of CRISPR-edited cell lines for interrogation of lysosomal storage disease gene function "
        "analysis. Contains gRNA sequences and allele sequencing results for all edits examined. "
        "Additionally contains annotations of proteins to individual organelles."
    ),
    "Dataset S2.xlsx": (
        "nDIA whole cell proteomics of day 50 differentiated iN, iDA of LSD mutants."
    ),
    "Dataset S3.xlsx": (
        "nDIA whole cell proteomics of day 30 differentiated iN, iDA of select LSD mutants."
    ),
    "Dataset S4.xlsx": (
        "Baseline and altered protein-protein interactions."
    ),
    "Dataset S5.xlsx": (
        "TMT-based analysis of LysoIP sample proteomics of untagged, Control and ASAH1-/- HeLa cells."
    ),
    "Dataset S6.xlsx": (
        "Label free lipidomics of whole cell and LysoIP of untagged, Control and ASAH1-/- HeLa cells."
    ),
}

# Read CSV
with open(CSV_FILE, newline="", encoding="utf-8-sig") as f:
    content = f.read()

with open(CSV_FILE, newline="", encoding="utf-8-sig") as f:
    reader = csv.reader(f)
    rows = list(reader)

# Find header
header_idx = None
for i, row in enumerate(rows):
    if row[0].strip() == "figure_panel":
        header_idx = i
        break

header = [col.strip() for col in rows[header_idx]]
ci_idx  = header.index("code_info")
cf1_idx = header.index("code_file_1")
di_idx  = header.index("dataset_info")
fp_idx  = header.index("figure_panel")

updated = 0
for row in rows[header_idx + 1:]:
    if len(row) < max(ci_idx, cf1_idx, di_idx) + 1:
        row.extend([""] * (max(ci_idx, cf1_idx, di_idx) + 1 - len(row)))

    fp   = row[fp_idx].strip()
    cf1  = row[cf1_idx].strip()
    mn   = row[header.index("method_name")].strip()
    ii   = header.index("instrument_info")

    # Update code_info from script name
    if cf1 in CODE_INFO and not row[ci_idx].strip():
        row[ci_idx] = CODE_INFO[cf1]
        updated += 1

    # Update instrument_info for spinning disk rows
    if mn in INSTRUMENT_INFO:
        row[ii] = INSTRUMENT_INFO[mn]
        updated += 1

    # Update dataset_info from supplemental legends
    if fp in DATASET_INFO and not row[di_idx].strip():
        row[di_idx] = DATASET_INFO[fp]
        updated += 1

# Write back
with open(OUT_FILE, "w", newline="", encoding="utf-8") as f:
    writer = csv.writer(f)
    writer.writerows(rows)

print(f"Updated {updated} fields in {OUT_FILE}")
