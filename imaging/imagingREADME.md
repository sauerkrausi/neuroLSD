# Imaging Workflow: Cellpose-SAM, Calcium Analysis & TH Quantification

This folder contains the code and documentation for analyzing live-cell calcium imaging data using Cellpose-SAM for segmentation, R-based pipelines for ΔF/F event analysis, and scripts for TH-positive neuron quantification.

---

## 1. Segmentation and Models (`models/`)
- The `models` folder contains custom trained Cellpose-SAM models for neuronal soma segmentation
- Includes a `MODEL_INDEX.txt` file that lists Zenodo download links for each model

**Requirements:** see `requirements.txt` or `requirements_MPS.txt`
Install via:
```bash
pip install -r requirements.txt
# For Apple Silicon (MPS):
pip install -r requirements_MPS.txt
```


## 2. Calcium Signal Processing (`Calcium/`)
**Python Scripts:** `run_cellpose_calcium_remoteHDD.py`, `run_cellpose_calcium_remoteHDD_MPS.py`
- Segment neuronal soma in grayscale calcium imaging stacks (`.tiff` or `.nd2`) using Cellpose-SAM with the custom `cpsam` model
- Supports batch mode with automated ROI mask generation and QC overlays
- Saves masks, overlays, and log files for quality control
- MPS version optimized for Apple Silicon GPUs

**R Scripts:** `cellposeSAM_CNER_Calcium-liveCell_template.Rmd, neuroLSD_CalciumMaster.Rmd`
- Load segmentation masks and extract per-ROI fluorescence traces
- Compute ΔF/F normalization, apply adaptive thresholding, and binarize events
- Quantify neuronal activity: spike timing, synchrony, coactivity, and burst frequency
- Outputs: raster plots, trace plots, synchrony and correlation matrices, heatmaps, and summary metrics


## 3. TH-Positive Neuron Quantification (`TH/`)

**Scripts:** `diff118_iNiDA_THeval.Rmd, th_neuron_quantification_interative.py`
- Quantify TH-positive neurons in iN and iDA cultures
- Supports both automated and interactive workflows
- Provides condition-wise comparisons and group-level analyses
- Outputs: cell counts, density plots, violin plots, and per-condition summary tables


## 4. Endo-lysosomal Confocal Analysis (`endolyso/`)

Analysis for HeLa endo-lysosomal staining experiments. This folder contains scripts to split channels, run Cellpose‑SAM segmentation, compute EEA1–TFN colocalization, and generate R-based statistics and plots.

**Files:**
- `split_channels_from_composite.py` — Split composite confocal images into channel-specific TIFFs.
- `run_cellpose_EndoLyso_HeLa.py` — Parse split-channel TIFF filenames to extract sample prefix, channel, and wavelength; preprocess DNA and organelle channels with background subtraction, optional denoising, and white tophat filtering; segment nuclei with the custom `cpsam` Cellpose model on Apple MPS/GPU and organelles via Otsu thresholding with size filtering; generate QC overlays for nuclei, organelles, and colocalization pairs; compute per-image metrics (counts, areas, per-nucleus measures) and colocalization fractions; output `per_image_summary.csv` and `pairwise_colocalization.csv` in the `metrics` subfolder.
- `pyEEA1TFNeval.py` — Pair EEA1 and TFN channels, segment, compute per-image metrics and pairwise colocalization; writes `per_image_summary.csv` and `pairwise_colocalization.csv` under `pyeval/metrics/` with QC overlays in `pyeval/qc_overlays/`.
- `HeLa_lyso_eval.Rmd` — End-to-end R analysis: reads the two summary CSVs, parses sample metadata, performs outlier filtering and Wilcoxon tests (Control vs ASAH1 by treatment), and produces violin/box/jitter plots with BH‑adjusted p labels.
- `requirements_MPS.txt` — Python dependencies for Apple Silicon (MPS) execution.

**Typical workflow:**
1. Split channels if needed: `python split_channels_from_composite.py`
2. Segment with Cellpose‑SAM: `python run_cellpose_EndoLyso_HeLa.py`
3. Compute colocalization and per‑image metrics: `python pyEEA1TFNeval.py`
4. Open `HeLa_lyso_eval.Rmd` in RStudio, knit to HTML to generate stats and figures.

**Outputs:**
- Metrics: `pyeval/metrics/per_image_summary.csv`, `pyeval/metrics/pairwise_colocalization.csv`
- QC overlays: `pyeval/qc_overlays/*.png`
- Figures and annotated CSVs: paths configured inside the Rmd (`outdir`); includes plots for `frac_coloc_eea1` and `frac_coloc_tfn` with treatment ordering Control 5 min → Control 45 min → ASAH1 5 min → ASAH1 45 min.