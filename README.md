![ProjectLogo](/logoNeuroLSD.jpg)
# neuroLSD
The folder contains scripts and evaluations realated to the neuroLSD project.

## proteome
Contains Rmd scipts for evaluation of stem-cell derived cortical-like iNeurons and stem-cell derived dopamenergic-like iDA neurons at day 50 of in-vitro differentiation.

> - Diff118_d23_iN-iDA_TMTpro
> - Diff132_d50_nDIA
> - Diff136_axonal proteomics
> - HeLa LysoIP

### Diff132_d50_nDIA

- [ ]	**I. load packackes & set project directory**\
Loads required R libraries and sets the working directory for the proteomics project.

- [ ]	**II. Circos plot for neuroLSD**\
Generates a chord diagram visualizing connections among LSD genes grouped by disease class, highlighting confirmed knockouts.

- [ ]	**III. day50 nDIA whole cell proteomics of iNeurons and iDA**
Performs PCA on day 50 whole-cell proteomics of iN and iDA neurons, with and without QC samples, and generates visualizations by genotype, neuron type, and day.

  - **IIIa. PCA plots of data**\
  Computes and plots PCA projections colored by genotype, neuron type, and day using scaled log2 protein intensities.
  - **IIIb. QC plots for supplements**\
  Generates bar and scatter plots summarizing protein IDs and peptide sequences across neuron types and genotypes for supplemental QC.
  - **IIIc. day50 full dataset for heatmaps**\
  Loads and reshapes full day 30 and 50 proteomics data into matrices for generating annotated heatmaps of log2 fold changes.

- [ ]	IV. Day 30 & Day 50 in-depth evaluation
Performs downstream analyses comparing day 30 and day 50 samples, including subcellular annotations, heatmaps, and expression summaries.

  - **IVa. – Subcell annotation to df for QC and Organelle Violin Plots**\
  Adds subcellular localization annotations to data frames and exports annotated wide-format CSVs for day 30 and day 50.
  
  - **IVb. – LSD KO protein abundnace across genotypes –**\
  Plots raw abundance values of known LSD genes across neuron types and genotypes to confirm knockout efficiency.
  
  - **IVc. – neuro QC plots ––**\
  Plots expression of neuronal QC marker proteins (e.g., synaptic markers) in control iN and iDA samples at day 50.
  
  - **IVd. – Abundance of neuromarkers in Ctrl day 30 & day 50 ––**\
  Compares abundance of general neuronal markers across control iN and iDA neurons between day 30 and day 50.
  
  - **IVe. – Abundance of Pre/Post Synaptic Markers in Ctrl day 30 & day 50 ––**\
  Plots developmental trajectories of pre- and post-synaptic markers in control neurons using log2 abundance values.
  
  - **IVf. – NeuroDev protein trajectory for day 30 & day 50 for select genotypes ––**\
  Analyzes expression changes of pluripotency and neuronal differentiation markers in select genotypes over time.
  
  - **IVg. – Violin Plot Function for LSD Organelle Proteomics ––**\
  Defines a reusable function to visualize log2 fold changes per genotype and neuron type, stratified by subcellular annotation.
  
  - **IVh. – Sph mutant investigation ––**\
  Focuses on sphingolipid disorder mutants, generating heatmaps and correlation analyses of their log2FC profiles and organelle-specific annotations.
      - **IVh1. — Add annotations & heatmap & save csv —**\
    Creates and saves a z-score normalized log2FC heatmap of Sph mutants with organelle annotations.
      - **IVh2. — Correlation analysis on selected annotations —**\
    Computes and visualizes Spearman correlations between selected organelle annotations and Sph mutant log2FC profiles.
      - **IVh3. — Correlation analysis on all organelles —**\
    Generates organelle-specific KO × KO correlation matrices for Sph mutants and arranges them into a multi-panel heatmap PDF.
  
  - **IVi. – GRN d30/50 DA neuron investigation ––**\
  Compares GRN iDA neuron profiles at day 30 and day 50, highlighting unique expression patterns using scatter plots and heatmaps.

- [ ]	**V. Organelle Cross-DF correlations**\
Introduces correlation analysis of log2FC values across organelle annotations and genotypes using annotated day 50 proteomics data.
  - **Va. — User-Select Organelle-annotation correlation**\
  Defines a function to compute and visualize genotype × genotype correlation matrices for user-specified organelle annotations.
  - **Vb. — All Organelle-annotation correlation —**\
  Computes Spearman correlations between genotype effects and each organelle annotation, outputting a heatmap of genotype × organelle associations.  
  - **Vc. — Plot average correlation per annotation across genotypes**\
   Plots average correlation per genotype for selected annotations using a barplot with error bars.  
  - **Vd. — Linear regression on neuro nDIA dataset**\
  Overview section introducing linear modeling strategies for analyzing differential protein expression.
    - **Vd1. — Linear regression on gene level across neurontypes**\
      Fits gene-wise linear models with genotype, neuron type, and interaction terms to identify significant differential effects.
    - **Vd2. — Stratified linear modeling of genotype effects by neuron type —**\
      Runs linear models separately for iN and iDA neurons to assess genotype-specific effects in each neuronal context.
    - **Vd3. — Limma linear modeling —**\
      Uses the limma package to fit an empirical Bayes model testing main and interaction effects of genotype and neuron type across all genes.



- [ ]	**VI. nMOST whole cell plot correlcation matrix per genotype cluster**\
Begins section on correlating nDIA neuronal data with HeLa nMOST proteomics using tidy correlation matrices.
  - **VIa. — nMOST LSD GO signatures of select Clusters for all genotypes —**\
  Plots GO cluster correlation scores across LSD genotypes from nMOST proteome data using boxplots and barplots.
  - **VIb.  — nMOST GO Signatures from select genotype over all GO classes —**\
  Visualizes how a single genotype maps across all GO-derived clusters in nMOST proteomics data.


- [ ]	**VII. iNeuron - Hela Proteome Correlations**\
Introduces global correlation analysis between neuron-derived and HeLa-derived proteomes across matched genotypes.
  - VIIa. **— Neuron–HeLa Proteome Correlation Analysis —**\
  Calculates and visualizes pairwise Pearson correlations across iN, iDA, and HeLa proteomes using Circos plots.
  -	**VIIb. — nMOST-nDIA correlation of select genotype -**\
  Allows user to generate circos plots of neuron–HeLa correlations for a single genotype, separating positive and negative links.
  -	**VIIc1. — Select organelle level —**\
  Computes correlation matrices between neuron and HeLa proteomes for individual organelle annotations and visualizes them.
  -	**VIIc2. — Combined select organelle level —**\
  Combines proteins from multiple selected annotations and calculates a unified neuron–HeLa correlation heatmap.

## imaging
