![ProjectLogo](/logoNeuroLSD.jpg)
# neuroLSD
The folder contains scripts and evaluations realated to the neuroLSD project.\
🚧 work in progress 🚧


## 🧪 proteome
Contains Rmd scipts for evaluation of stem-cell derived cortical-like iNeurons and stem-cell derived dopamenergic-like iDA neurons at day 50 of in-vitro differentiation.

- `diff132_d50_nDIA.Rmd`: R Markdown pipeline for nDIA proteome analysis of whole cell iNeurons and iDA
- `diff136_iNd35_ctrl_asah1e1_axonalproteome.Rmd`: R Markdown pipeline for TMTpro proteome analysis of whole cell, soma and axonal fractions of iNeurons 


A detailed README for the scripts can be found in [`readmePROTEOME.md`](proteome/readmePROTEOME.md).

## 🔬 imaging
- `imaging/run_cellpose_calcium_remoteHDD.py`: Python script for batch image segmentation using Cellpose-SAM
- `imaging/cellposeSAM_CNER_Calcium-liveCell.Rmd`: R Markdown pipeline for calcium signal binarization and network analysis
- `requirements.txt`: Python environment dependencies

See [`imaging/README.md`](imaging/imagingREADME.md) for a detailed explanation of the workflow, input formats, and output structure.
