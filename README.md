![ProjectLogo](/logoNeuroLSD.jpg)
# neuroLSD
This repository contains scripts and evaluation-pipelines related to the **neuroLSD** project.

The **neuroLSD** project investigates how lysosomal storage disorder–associated pathways contribute to neurodegeneration using **stem-cell–derived neuronal models**. We integrate **quantitative proteomics** and **calcium imaging** to profile functional and molecular changes in cortical-like iNeurons (iN) and dopaminergic-like iDA neurons.

🚧 work in progress 🚧

## 🔬 Imaging
Analysis workflows for calcium imaging and TH-positive neuron quantification.  
- **Calcium**: Cellpose-SAM–based segmentation (Python) and R Markdown pipelines for ΔF/F signal extraction, trace binarization, and network-level analysis, including a master pipeline for multi-dataset integration.  
- **TH**: R and Python scripts for quantifying TH-positive neurons in iN/iDA cultures.  

See [`imagingREADME.md`](imaging/imagingREADME.md) for details on workflow steps, input formats, and outputs.



## 🧬 Lipidomics
Lipidomic profiling of HeLa and day‑21 iNeurons (Ctrl vs ASAH1‑/‑) across whole‑cell and organelle‑IP fractions, with outputs including class‑level barplots, volcano plots, and per‑lipid log2FC tables.



## 🧪 proteomics
Contains Rmd scipts for evaluation of stem-cell derived cortical-like iNeurons and stem-cell derived dopamenergic-like iDA neurons at day 50 of in-vitro differentiation.

- `diff118_iNiDA_d23_CtrlSMPD1ASAH1.Rmd`: TMTpro-based proteomic analysis of day 23 iN and iDA cells across control, SMPD1-/-, and ASAH1-/- conditions.
- `diff132_d50_nDIA.Rmd`: nDIA proteomics of whole-cell iNeurons and iDA neurons at day 50, including organelle-level annotation and HeLa cross-comparison.
- `diff136_iNd35_ctrl_asah1e1_axonalproteome.Rmd`: TMTpro analysis of whole-cell, soma, and axonal fractions in iNeurons with axonal enrichment and RoR modeling.
- `HeLa_Ctrl-ASAH1_LysoIP.Rmd`: Proteomic analysis of isolated lysosomes from HeLa cells comparing control and ASAH1-/- lines, with detailed lysosomal and autophagy annotation.


A detailed README for the scripts can be found in [`readmePROTEOME.md`](proteome/readmePROTEOME.md).



## Requirements
- **R ≥ 4.3.0** with: `tidyverse`, `ComplexHeatmap`, `circlize`, `pheatmap`, `ggpubr`, `cowplot`  
- **Python ≥ 3.10** with: `cellpose`, `torch`, `numpy`, `pandas`, `scikit-image`  


## Usage
- **Proteomics:**  
  Run Rmd scripts in RStudio. Input: processed TMTpro/nDIA data.  
- **Calcium Imaging:**  
  Segment with `run_cellpose_calcium_remoteHDD.py`, then analyze with Rmd pipelines. Input: nd2 files.  