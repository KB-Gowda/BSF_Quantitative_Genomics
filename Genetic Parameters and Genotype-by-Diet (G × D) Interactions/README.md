# Estimation of Genetic Parameters and Genotype-by-Diet (G × D) Interactions for Growth Traits in Australian Black Soldier Fly (*Hermetia illucens*)

This repository contains the analysis scripts used to evaluate phenotypic variation, population structure, genetic parameters, and genotype × diet (G × D) interactions for growth traits in Australian black soldier fly (*Hermetia illucens*).

## Analyses

### 1. Descriptive Statistics and Fixed-Effect Analysis (GLM) — SAS

Phenotypic descriptive statistics and fixed-effect analyses were performed using SAS.

The analyses evaluated the effects of:

- Diet
- Sampling day
- Tray within diet
- Larval colour

Descriptive statistics, least-squares means, significance tests, and coefficients of determination were obtained for larval body weight (LBW), larval length (LL), larval width (LW), and larval surface area (LSA).

### 2. Population Structure — PCA

Principal component analysis (PCA) was performed using PLINK on the QC-filtered genotype dataset.

The first two principal components were visualised in R using `ggplot2` to assess genetic structure and the distribution of individuals across dietary groups.

The PCA was conducted using the same QC-filtered individuals used for construction of the genomic relationship matrix and subsequent genetic analyses.

### 3. Genetic Parameter Estimation — R

Genetic parameter analyses were performed in R using `ASReml-R`.

Genotype data were subjected to quality-control filtering using PLINK. The QC-filtered genotype dataset was used to construct genomic relationship matrices.

The analyses included:

- Univariate additive genomic models
- Univariate additive–dominance genomic models
- Bivariate genomic models

The models were used to estimate additive genetic variance, dominance variance, residual variance, heritability, dominance ratios, and genetic, dominance, residual, and phenotypic correlations for the growth traits.

### 4. Genotype × Diet Interaction — R

Genotype × diet interactions were estimated using multivariate genomic models implemented in `ASReml-R`.

For each growth trait, measurements from the three dietary environments were treated as correlated traits:

- SYK
- BSG
- FVW

Genetic (co)variances and genetic correlations between diets were estimated to assess genotype re-ranking across dietary environments.

The following pairwise genetic correlations were estimated:

- SYK × BSG
- SYK × FVW
- BSG × FVW

Lower genetic correlations indicate greater genotype re-ranking between diets and stronger evidence of G × D interaction.

Diet-specific heritabilities were estimated from both the multivariate G × D models and the univariate additive and additive–dominance models.

### 5. GEBV Extraction and Reaction-Norm Plots — R

Diet-specific genomic estimated breeding values (GEBVs) were extracted from the multivariate G × D models and stored as CSV files.

Reaction-norm plots were generated using `ggplot2` to visualise individual genetic responses across the three dietary environments.

For each growth trait, individuals were ranked according to their mean GEBV across diets, and the top 5% were visualised to illustrate genetic re-ranking and changes in genetic merit across dietary environments.

## Traits

The primary growth traits analysed were:

- **LBW** — larval body weight (mg)
- **LL** — larval length (mm)
- **LW** — larval width (mm)
- **LSA** — larval surface area (mm²)

## Dietary Treatments

- **SYK** — soy-okara
- **BSG** — brewers' spent grain
- **FVW** — fruit and vegetable waste

## Software

- **SAS** — descriptive statistics and fixed-effect analyses
- **PLINK** — genotype quality control and PCA
- **R** — data processing and visualisation
- **ASReml-R** — genomic mixed models and G × D analyses
- **snpReady** — genomic data processing and relationship matrix construction
- **AGHmatrix** — additive and dominance relationship matrices
- **ggplot2** — visualisation

## Repository Contents

The repository contains the scripts for the analytical workflow, including:

- Descriptive statistics and GLM analysis
- Genotype quality control
- Population structure and PCA
- Genomic relationship matrix construction
- Genetic parameter estimation
- G × D modelling
- GEBV extraction
- GEBV reaction-norm visualisation

Raw genotype and phenotype datasets are not included in the repository.
