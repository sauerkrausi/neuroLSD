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
- Loads Cellpose/SAM mask outputs and raw calcium signal data
- Computes ΔF/F traces per ROI
- Performs binarization and extracts event synchrony
- Generates per-cell trace plots and community-level summaries
