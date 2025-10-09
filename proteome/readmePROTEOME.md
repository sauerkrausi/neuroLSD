# Proteome
Contains Rmd scripts for proteomics evaluation of stem-cell derived cortical-like iNeurons and stem-cell derived dopaminergic-like iDA neurons. Includes whole cell analysis of in vitro differentiated neurons at day 23 (TMTpro-DDA), 30 and 50 (LFQ-nDIA), TMTpro proteomics of neuronal projections as well as proteomic analysis of isolated lysosomes. 

## active projects 
🚧 `Diff118_iNiDA_d23.Rmd`: Pipeline for TMTpro proteome analysis of day 23 iN and iDA whole cell samples
-	TMTpro-based differential proteomic analysis of day 23 iN and iDA samples.
-	Compares Ctrl, SMPD1-/-, and ASAH1-/- lines across neuronal types.
-	Outputs log2FC, significance tables, and cluster-based GO term enrichments.
-	Integrates curated sub-cellular annotations to identify lysosomal and synaptic changes.

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

🚧 `diff136_iNd35_ctrl_asah1e1_axonalproteome.Rmd`: R Markdown pipeline for TMTpro proteome analysis of neuronal whole cell, soma and projections fractions of iNeurons
-	TMTpro-based proteome of iNeuron whole-cell, soma, and projection fractions at day 35.
-	Evaluates compartment-specific KO effects using log2FC and RoR metrics.
-	Correlates subcellular enrichment (e.g., endo/lyso/synaptic) with KO-induced localization shifts.
-	Visualizes data via heatmaps, barplots, violin plots, and custom “rainfall” plots for spatial proteomic mapping.
