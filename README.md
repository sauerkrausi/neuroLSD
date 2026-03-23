![ProjectLogo](/logoNeuroLSD.jpg)
# neuroLSD
This repository contains scripts and evaluation pipelines related to the **neuroLSD** project.

The **neuroLSD** project investigates how lysosomal storage disorder–associated pathways contribute to neurodegeneration using **stem-cell–derived neuronal models**. We integrate **quantitative proteomics**, **lipidomics**, and **calcium imaging** to profile functional and molecular changes in cortical-like iNeurons (iN) and dopaminergic-like iDA neurons.

Lysosomal Storage Disorders (LSDs) comprise a broad group of diseases affecting lysosomal degradation, ion flux, and lipid catabolism. Within this group, sphingolipidoses genes such as GBA1, SMPD1, and ASAH1 are known or candidate risk factors for Parkinson’s Disease, though the mechanisms remain unclear. Building on prior work in HeLa cells, we generated a human embryonic stem cell library with 23 LSD gene knockouts and differentiated them into cortical and dopaminergic neurons.  

Our proteomic and functional analyses reveal lineage-specific alterations in organellar proteomes, uncovering diverse neuronal vulnerabilities. Notably, GBA1-/- and ASAH1-/- dopaminergic neurons show disruptions in synaptic and mitochondrial compartments, correlating with impaired neuronal firing and presynaptic protein localization.  

Structural analysis of ASAH1-deficient endolysosomes by cryo-electron tomography reveals swollen organelles largely devoid of the dense internal membranes characteristic of wild-type cells, but enriched in intralumenal vesicle compartments.

This LSD mutant toolkit and associated proteomic landscape provide a resource for decoding how lysosomal dysfunction impacts neuronal health and for exploring mechanistic links between lysosomal storage disorders and Parkinson’s Disease.

📖 Current version of Preprint: [BioRxiv]( https://doi.org/10.1101/2025.10.08.681047)


🚧 work in progress 🚧

## 🔬 Imaging
Analysis workflows for calcium imaging and TH-positive neuron quantification.  
- **Calcium:** Segmentation of neuronal soma in calcium imaging using Cellpose-SAM custom models, followed by R pipelines for ΔF/F extraction, trace binarization, network analysis, and multi-dataset integration.
- **TH:** Automated quantification of TH-positive neurons in iN/iDA cultures via Python-based segmentation and R statistical analysis.
- **Models:** Custom Cellpose SAM models for neuronal soma live in imaging/models, with MODEL_INDEX.txt linking to Zenodo downloads.
- **Endo‑lyso HeLa:** Endo‑lysosomal confocal analysis for HeLa in imaging/endolyso, including channel splitting, Cellpose SAM based segmentation, EEA1 and TFN colocalization metrics, and an Rmd that outputs figures and annotated CSVs.

See [`imagingREADME.md`](imaging/imagingREADME.md) for details on workflow steps, input formats, and outputs.



## 🧬 Lipidomics
Lipidomic profiling of HeLa and day‑21 iNeurons (Ctrl vs ASAH1‑/‑) across whole‑cell and organelle‑IP fractions, with outputs including class‑level barplots, volcano plots, and per‑lipid log2FC tables.



## 🧪 Proteomics
Contains Rmd scripts for evaluation of stem-cell derived cortical-like iNeurons and stem-cell derived dopaminergic-like iDA neurons at day 50 of in-vitro differentiation.

- `diff132_d50_nDIA.Rmd:` nDIA proteomics of whole-cell iNeurons and iDA neurons at day 50, including organelle-level annotation and HeLa cross-comparison.
- `HeLa_Ctrl-ASAH1_LysoIP.Rmd:` Proteomic analysis of isolated lysosomes from HeLa cells comparing control and ASAH1-/- lines, with detailed lysosomal and autophagy annotation.


A detailed README for the scripts can be found in [`readmePROTEOME.md`](proteome/readmePROTEOME.md).



## Requirements
- **R 4.5.2** (RStudio Version 2026.01.1+403) with: `bigstatsr`, `broom`, `Cairo`, `CalNetExploreR`, `car`, `circlize`, `ComplexHeatmap`, `ComplexUpset`, `corrplot`, `cowplot`, `data.table`, `devtools`, `dplyr`, `factoextra`, `fmsb`, `forcats`, `furrr`, `ggbiplot`, `ggdendro`, `ggdist`, `gghalves`, `gghighlight`, `ggplot2`, `ggpp`, `ggpmisc`, `ggpubr`, `ggraph`, `ggrepel`, `ggsci`, `ggsignif`, `ggVennDiagram`, `grid`, `gridExtra`, `igraph`, `irlba`, `janitor`, `limma`, `lintr`, `lubridate`, `magick`, `msstats`, `NatParksPalettes`, `org.Hs.eg.db`, `patchwork`, `pheatmap`, `plotly`, `plyr`, `png`, `purrr`, `qs`, `RColorBrewer`, `readr`, `reshape2`, `rlang`, `rstatix`, `scales`, `signal`, `stringr`, `superheat`, `Ternary`, `tibble`, `tidyplots`, `tidyr`, `tidyverse`, `timecourse`, `topGO`, `umap`, `UpSetR`, `uwot`, `viridis`, `zoo`
- **Python ≥ 3.10** with: `cellpose`, `matplotlib`, `nd2`, `numpy`, `pandas`, `pathlib`, `re`, `scipy`, `skimage`, `tifffile`, `torch`, `tqdm`, `typing`

For more detailed package requirements, please check the manuscript.

## Usage
- **Proteomics:**  
  Run Rmd scripts in RStudio. Input: processed TMTpro/nDIA data.  
- **Calcium Imaging:**  
  Segment with `run_cellpose_calcium_remoteHDD.py`, then analyze with Rmd pipelines. Input: nd2 files. 
  
## Versions of preprints
📖 Preprint version 2 (current):[BioRxiv]( https://doi.org/10.1101/2025.10.08.681047) 
📖 Preprint version 1: [BioRxiv](https://www.biorxiv.org/content/10.1101/2025.10.08.681047v1)

Also contains the following datasets (removed from newer versions for clarity):
`diff118_iNiDA_d23.Rmd:` TMTpro-based proteomic analysis of day 23 iN and iDA cells across Control and ASAH1-/- conditions.
`diff136_iNd35_ctrl_asah1e1_axonalproteome.Rmd:` TMTpro analysis of whole-cell, soma, and axonal fractions in iNeurons with axonal enrichment and RoR modeling.


## Acknowledgments & Funding

We thank members of the Harper lab for feedback and support.

This project is supported by:  
- [Aligning Science Across Parkinson’s (ASAP)](https://parkinsonsroadmap.org/) and the [Michael J. Fox Foundation](https://www.michaeljfox.org/)  
- [Bluefield Project](https://www.bluefieldproject.org/)  
- National Institutes of Health (NIH) — (see [nih.gov](https://www.nih.gov/))  
- [Howard Hughes Medical Institute](https://www.hhmi.org/)
- [Warren Alpert Foundation](https://www.warrenalpertfoundation.org/)  
- [Max Planck Society](https://www.mpg.de/en)  
- Fred & Joan Goldberg Post-doctoral Fellowship (see [Goldberg Fellows](https://cellbio.hms.harvard.edu/goldberg-fellows))  
- Boehringer Ingelheim Fonds PhD Fellowship (see [BIFonds PhD Fellowships](https://www.bifonds.de/fellowships-grants/phd-fellowships.html))  

We acknowledge the [Core for Imaging Technology and Education (CITE, HMS)](https://cite.hms.harvard.edu/), the [Harvard Electron Microscopy Core](https://electron-microscopy.hms.harvard.edu/)), and the Central Electron Microscopy Facility at the Max Planck Institute of Biophysics for imaging support.  

## Interested? Learn more
Are you interested to find out more about what we are doing?  
[Harper Lab](https://harper.hms.harvard.edu/everything-protein-and-organelle-quality-control)  
[Gygi Lab](https://gygi.hms.harvard.edu/)  
[Wilfling Lab](https://www.biophys.mpg.de/mechanisms-cellular-quality-control)  
[Farese and Walther Lab](https://www.fwlaboratory.org/)  