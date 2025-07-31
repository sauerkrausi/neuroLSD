# Imaging Workflow: Cellpose-SAM, Calcium Analysis & TH Quantification

This folder contains the code and documentation for analyzing live-cell calcium imaging data using Cellpose-SAM for segmentation, R-based pipelines for ΔF/F event analysis, and scripts for TH-positive neuron quantification.

---

## 1. Segmentation

**Scripts:** `run_cellpose_calcium_remoteHDD.py`, `run_cellpose_calcium_remoteHDD_MPS.py`  
- Segment grayscale calcium imaging stacks (`.tiff` or `.nd2`) using Cellpose-SAM  
- Supports batch mode with automated ROI mask generation and QC overlays  
- Saves masks, overlays, and log files for quality control  
- MPS version optimized for Apple Silicon GPUs  

**Requirements:** see `requirements.txt` or `requirements_MPS.txt`  
Install via:
```bash
pip install -r requirements.txt
# For Apple Silicon (MPS):
pip install -r requirements_MPS.txt
```



## 2. Calcium Signal Processing

**Scripts:** `cellposeSAM_CNER_Calcium-liveCell_template*.Rmd, neuroLSD_CalciumMaster.Rmd`
- Load segmentation masks and extract per-ROI fluorescence traces
- Compute ΔF/F normalization, apply adaptive thresholding, and binarize events
- Quantify neuronal activity: spike timing, synchrony, coactivity, and burst frequency
- Outputs: raster plots, trace plots, synchrony and correlation matrices, heatmaps, and summary metrics


## 3. TH-Positive Neuron Quantification

**Scripts:** `diff118_iNiDA_THeval.Rmd, th_neuron_quantification_interative.py`
- Quantify TH-positive neurons in iN and iDA cultures
- Supports both automated and interactive workflows
- Provides condition-wise comparisons and group-level analyses
- Outputs: cell counts, density plots, violin plots, and per-condition summary tables