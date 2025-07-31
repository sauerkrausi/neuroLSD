![ProjectLogo](/logoNeuroLSD.jpg)
# neuroLSD
The folder contains scripts and evaluations realated to the neuroLSD project.\
🚧 work in progress 🚧

## 🔬 imaging
- `run_cellpose_calcium_remoteHDD.py`: Python script for batch segmentation of time-lapse calcium imaging using Cellpose-SAM models.
- `cellposeSAM_CNER_Calcium-liveCell.Rmd`: R Markdown pipeline for ΔF/F signal extraction, trace binarization, and network-level calcium event analysis.
- `neuroLSD_CalciumMaster.Rmd`: Master analysis pipeline integrating multiple datasets for per-cell event correlation and group-level network dynamics.

See [`imaging/README.md`](imaging/imagingREADME.md) for a detailed explanation of the workflow, input formats, and output structure.



## 🧬 Lipidomics
Lipidomic profiling of HeLa and day‑21 iNeurons (Ctrl vs ASAH1‑/‑) across whole‑cell and organelle‑IP fractions, with outputs including class‑level barplots, volcano plots, and per‑lipid log2FC tables.



## 🧪 proteomics
Contains Rmd scipts for evaluation of stem-cell derived cortical-like iNeurons and stem-cell derived dopamenergic-like iDA neurons at day 50 of in-vitro differentiation.

- `diff118_iNiDA_d23_CtrlSMPD1ASAH1.Rmd`: TMTpro-based proteomic analysis of day 23 iN and iDA cells across control, SMPD1-/-, and ASAH1-/- conditions.
- `diff132_d50_nDIA.Rmd`: nDIA proteomics of whole-cell iNeurons and iDA neurons at day 50, including organelle-level annotation and HeLa cross-comparison.
- `diff136_iNd35_ctrl_asah1e1_axonalproteome.Rmd`: TMTpro analysis of whole-cell, soma, and axonal fractions in iNeurons with axonal enrichment and RoR modeling.
- `HeLa_Ctrl-ASAH1_LysoIP.Rmd`: Proteomic analysis of isolated lysosomes from HeLa cells comparing control and ASAH1-/- lines, with detailed lysosomal and autophagy annotation.


A detailed README for the scripts can be found in [`readmePROTEOME.md`](proteome/readmePROTEOME.md).

