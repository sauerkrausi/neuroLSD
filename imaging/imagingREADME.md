# Imaging Workflow: Cellpose-SAM & Calcium Analysis

This folder contains the code and documentation for analyzing live-cell calcium imaging data using Cellpose-SAM for segmentation and R-based network analysis.

---

## 1. Segmentation

**Script:** `run_cellpose_calcium_remoteHDD.py`  
- Uses Cellpose-SAM to segment grayscale calcium imaging TIFF stacks
- Supports batch processing from a specified input directory
- Saves mask outputs and log files for quality control
- Accepts `.nd2` and `.tiff` formats

**Requirements:** see `requirements.txt`  
Install via:
```bash
pip install -r requirements.txt
```


## 2.  Calcium Signal Processing

**Script:** `cellposeSAM_CNER_Calcium-liveCell.Rmd`
- Loads Cellpose/SAM segmentation masks and raw signal data
- Extracts per-ROI fluorescence over time and computes ΔF/F normalization
- Applies adaptive thresholding and event binarization
- Calculates event timing, synchrony, and coactivity across the network
- Outputs: raster plots, trace plots, synchrony matrices, and summary metrics


## 3. Dataset-Level Meta-Analysis
**Script:** `neuroLSD_CalciumMaster.Rmd`
- Integrates multiple imaging experiments grouped by genotype, condition, or neuron type
- Aggregates per-cell and per-experiment activity metrics
- Computes inter-cell correlation networks and burst frequency distributions
- Produces violin plots, correlation heatmaps, and condition-wise summary plots
