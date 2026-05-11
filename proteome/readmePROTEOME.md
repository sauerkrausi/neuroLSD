# Proteome
Contains Rmd scripts for proteomics evaluation of stem-cell derived cortical-like iNeurons and stem-cell derived dopaminergic-like iDA neurons. Includes whole cell analysis of in vitro differentiated neurons at day 23 (TMTpro-DDA), 30 and 50 (LFQ-nDIA), TMTpro proteomics of neuronal projections as well as proteomic analysis of isolated lysosomes. 

## active projects 
🚧 `HeLa_Ctrl-ASAH1_LysoIP.Rmd`: Pipeline for TMTpro proteome analysis of isolated lysosomes from HeLa Ctrl and ASAH1 mutants 
- Lysosome-enriched fraction proteomics from HeLa Ctrl vs ASAH1-/-.
- Focus on lysosomal hydrolase abundance, lysosomal-endosomal pathway components, and autophagy markers.
- Includes cluster-based heatmaps, scaled violin plots, and annotation-wise log2FC summaries.
- Performs GO term enrichment per cluster and compares overlap with iNeuron datasets via Venn diagrams.

🚧 `diff132_d50_nDIA.Rmd`: R Markdown pipeline for nDIA proteome analysis of whole cell iNeurons and iDA
- nDIA-based quantification of whole-cell iN and iDA neurons at day 50.
- Performs pairwise KO vs Ctrl comparisons and neuron-type interaction effects.
- Integrates organelle annotations and correlates KO profiles by organelle enrichment.
- Generates cluster heatmaps, GO term enrichments, and cross-cell-type correlation plots (HeLa vs neuron).

Modules
- Module 0: Setup & Configuration & Experimental Background — package loading, output directories, experimental background for nDIA iN/iDA d30+d50, Circos plot for neuroLSD
- Module 1: Data Import & QC — disease-class annotation, data cleaning, replicate correlation, F-ratio variance analysis, PCA, protein ID and RSD QC plots
- Module 2: Summarized Data & Subcellular Analysis — fold-change import (d30/d50), subcellular annotation, heatmaps, Euler plots, LSD protein abundance, p-value heatmaps, Stouffer correction, ternary plots, half-Circos plots, barplots of significant IDs, lineplots, SV/exocytosis driver analysis, PCA on annotation subsets, iN vs iDA cell-type comparison
- Module 3: Neuro QC & Synaptic Markers — pre/post-synaptic marker plots, neuro annotation by neuron × genotype, differentiation markers, vATPase, SV/SNARE, and SV docking plots
- Module 4: Time Course Analysis (d30 vs d50) — neuromarker abundance, synaptic marker trajectories, NeuroDev protein-list summary, GBA1/ASAH1 time course comparison, statistical tests, barplot function, correlation heatmaps, global d30 vs d50 proteome comparison
- Module 5: Long Diff Analysis ASAH1 (d30–50–70) — timecourse df construction, heatmaps, PCA, violin plots across annotations, pairwise annotation correlations, ChimeraX helper for structure mapping
- Module 6: Violin & Volcano Plots** — violin plot function for LSD organelle proteomics, autophagy machinery heatmap across genotypes × neuron types, volcano plots with annotation highlights
- Module 7: Organelle Correlations Analysis — cross-organelle correlation overview, organelle impact scoring, disease class impact scoring, avg. correlations per annotation/genotype/disease class/nMOST GO cluster, bubble plots iN vs iDA, peroxisomal abundance evaluation
- Module 8: Disease Class Analysis — Sph mutant investigation, organelle correlation analysis, diff118 (d23) vs diff132 (d50) ASAH1 comparison, Euler plots, linear model for GBA1/ASAH1 drivers, GRN d30/d50 DA neuron analysis, disease-class comparison d50 with annotation-based Circos plots
- Module 9: Linear Regression & Modeling — UMAP of neuroLSD dataset, linear regression on gene level across neuron types, stratified modeling of genotype effects, limma modeling, UMAP of linear model results and combined features
- Module 10: nDIA–nMOST Comparative Analysis — nMOST LSD GO signatures, neuron–HeLa proteome correlations, nMOST–nDIA correlation per genotype, HeLa–nMOST vs neuron–nDIA organelle correlation
- Module 11: Master PPI Evidence & Analysis — PPI vulnerability pipeline for iN/iDA LSD mutants; PPI database construction, baseline network, trimer/tetramer detection and validation, complex growing, top-down vs bottom-up comparison, network visualization, stoichiometry analysis, KO vulnerability analysis (significance-filtered), pathway disruption (sphingosine, ceramide, OXPHOS), multi-genotype timecourse summary
- Module 12: PPI Vulnerability Analysis – HeLa nMOST — HeLa data loading/QC, baseline PPI network, trimer/tetramer detection and validation, complex growing, network visualization, KO vulnerability analysis with extended visualizations


## retired projects 
🚧 `Diff118_iNiDA_d23.Rmd`: Pipeline for TMTpro proteome analysis of day 23 iN and iDA whole cell samples
-	TMTpro-based differential proteomic analysis of day 23 iN and iDA samples.
-	Compares Ctrl, SMPD1-/-, and ASAH1-/- lines across neuronal types.
-	Outputs log2FC, significance tables, and cluster-based GO term enrichments.
-	Integrates curated sub-cellular annotations to identify lysosomal and synaptic changes.

🚧 `diff136_iNd35_ctrl_asah1e1_axonalproteome.Rmd`: R Markdown pipeline for TMTpro proteome analysis of neuronal whole cell, soma and projections fractions of iNeurons
-	TMTpro-based proteome of iNeuron whole-cell, soma, and projection fractions at day 35.
-	Evaluates compartment-specific KO effects using log2FC and RoR metrics.
-	Correlates subcellular enrichment (e.g., endo/lyso/synaptic) with KO-induced localization shifts.
-	Visualizes data via heatmaps, barplots, violin plots, and custom “rainfall” plots for spatial proteomic mapping.

Modules
  - Module 0: Setup & Configuration & Experimental Background
  - Module 1: Data Import & add protein organelle annotations
  - Module 2: Heatmaps & Clustering
  - Module 3: Correlation Analysis
  - Module 4: Fraction comparisons (WC/SOMA/PROJECTION)
  - Module 5: Barplots & Comparisons
  - Module 6: Rainbowplots