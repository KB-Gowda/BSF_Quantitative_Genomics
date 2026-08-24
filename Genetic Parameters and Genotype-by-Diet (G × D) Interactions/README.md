# Estimation of Genetic Parameters and Genotype-by-Diet (G × D) Interactions for Growth Traits in Australian Black soldier fly larvae (*Hermetia illucens*)

This repository contains the two analysis scripts used to evaluate phenotypic variation, genetic parameters, genetic relationships, population structure, and genotype × diet (G × D) interactions for growth traits in Australian black soldier fly larvae (*Hermetia illucens*).

## Analyses

### 1. Descriptive Statistics and Fixed-Effect Analysis - SAS

Phenotypic descriptive statistics and fixed-effect analyses were performed using SAS.

The analyses evaluated the effects of:

* Diet
* Sampling day
* Tray within diet
* Larval colour as covariate

Descriptive statistics, least-squares means, significance tests, and coefficients of determination were obtained for larval body weight (LBW), larval length (LL), larval width (LW), and larval surface area (LSA).

### 2. Genetic Analysis and Genetic Parameter Estimation - R

genetic analyses were performed in R, with genotype quality control and principal component analysis conducted using PLINK and genetic mixed models fitted using `ASReml-R`.

The R script includes:

* Genotype quality-control filtering using PLINK
* Principal component analysis (PCA)
* Additive genetic relationship matrix (GRM) construction
* Dominance relationship matrix construction
* Univariate additive models
* Univariate additive–dominance models
* Bivariate models
* Multivariate genotype × diet models
* GEBV extraction
* Reaction-norm visualisation

The genetic models were used to estimate additive genetic variance, dominance variance, residual variance, heritability, dominance ratio, and genetic, dominance, residual, and phenotypic correlations.

### 3. Population Structure - R/PLINK

PCA was performed using PLINK on the QC-filtered genotype dataset.

The first two principal components were visualised in R using `ggplot2` to assess genetic structure and the distribution of individuals across the dietary groups.

The PCA used the same QC-filtered genotype dataset used for construction of the genetic relationship matrices and subsequent genetic analyses.

### 4. Genotype × Diet Interaction - R

Genotype × diet interactions were estimated using multivariate models implemented in `ASReml-R`.

For each growth trait, measurements from the three dietary environments were treated as correlated traits:

* SYK - soy-okara
* BSG - brewers' spent grain
* FVW - fruit and vegetable waste

Genetic (co)variances and genetic correlations between dietary environments were estimated to assess the consistency of genetic performance and genotype re-ranking across diets.

Pairwise genetic correlations were estimated between:

* SYK × BSG
* SYK × FVW
* BSG × FVW

Lower genetic correlations indicate greater changes in genetic ranking across dietary environments and therefore stronger evidence of genotype × diet interaction.

Diet-specific heritabilities were also estimated from the multivariate G × D models.

### 5. GEBV Extraction and Reaction-Norm Plots - R

Diet-specific genetic estimated breeding values (GEBVs) were extracted from the multivariate G × D models and exported as CSV files.

Reaction-norm plots were generated using `ggplot2` to visualise individual genetic responses across the three dietary environments.

For each growth trait, individuals were ranked according to their mean GEBV across diets, and the top 5% were visualised to illustrate genetic differences and changes in genetic performance across dietary environments.

## Traits

The primary growth traits analysed were:

* **LBW** - larval body weight (mg)
* **LL** - larval length (mm)
* **LW** - larval width (mm)
* **LSA** - larval surface area (mm²)

## Dietary Treatments

* **SYK** - soy-okara
* **BSG** - brewers' spent grain
* **FVW** - fruit and vegetable waste

## Software

* **SAS** - descriptive statistics and fixed-effect analyses
* **PLINK** - genotype quality control and PCA
* **R** - data processing, genetic analyses, and visualisation
* **ASReml-R** - genetic mixed models and G × D analyses
* **snpReady** - genetic data processing and relationship matrix construction
* **AGHmatrix** - additive and dominance relationship matrices
* **ggplot2** - visualisation

## Repository Contents

The repository contains **two analysis scripts** representing the complete analytical workflow:

1. **SAS script** - descriptive statistics and fixed-effect GLM analyses
2. **R script** - genotype quality control, PCA, genetic relationship matrices, genetic parameter estimation, G × D modelling, GEBV extraction, and reaction-norm visualisation

Raw genotype and phenotype datasets are not included in the repository.
