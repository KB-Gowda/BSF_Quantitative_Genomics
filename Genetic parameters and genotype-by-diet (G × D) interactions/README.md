# Genetic parameters and genotype-by-diet interactions for growth traits in Australian black soldier fly larvae: Implications for selective breeding

## Citation

If you use these scripts or find this work helpful, please cite our paper:

> **Gowda, K.B., Septriani, S., Jones, D.B., Jerry, D.R., Tedder, C., & Zenger, K.R.** (2026). Genetic parameters and genotype-by-diet interactions for growth traits in Australian black soldier fly larvae: Implications for selective breeding. *Journal of Insects as Food and Feed* (published online ahead of print 2026). https://doi.org/10.1163/23524588-bja10484

---

This folder contains the analysis scripts used to evaluate phenotypic variation, genetic parameters, genetic relationships, population structure, and genotype × diet (G × D) interactions for growth traits in Australian black soldier fly larvae (*Hermetia illucens*).

## Contents

| File | Description |
|---|---|
| `Descriptive_statistics.sas` | SAS script: descriptive statistics and fixed-effect analyses |
| `Genetic_parameters_GxD.r` | R script: genotype quality control, PCA, relationship matrices, genetic parameter estimation, G × D models, GEBV extraction and reaction-norm plots |
| `data.md` | Description of the input data |
| `Results.md` | Summary of the results |

<!-- CONFIRM: check that the descriptions of data.md and Results.md match their contents. -->

Raw genotype and phenotype datasets are not included in the repository.

---

## Analyses

### 1. Descriptive Statistics and Fixed-Effect Analysis (SAS)

*Script: `Descriptive_statistics.sas`*

Phenotypic descriptive statistics and fixed-effect analyses were performed using SAS.

The analyses evaluated the effects of:

* Diet
* Sampling day
* Tray within diet
* Larval colour (as a covariate)

Descriptive statistics, least-squares means, significance tests, and coefficients of determination were obtained for larval body weight (LBW), larval length (LL), larval width (LW), and larval surface area (LSA).

### 2. Genetic Analysis and Genetic Parameter Estimation (R)

*Script: `Genetic_parameters_GxD.r`*

Genetic analyses were performed in R. Genotype quality control and principal component analysis were conducted using PLINK, and genetic mixed models were fitted using ASReml-R.

The R script includes:

* Genotype quality-control filtering using PLINK
* Principal component analysis (PCA)
* Additive genomic relationship matrix (GRM) construction
* Dominance relationship matrix construction
* Univariate additive models
* Univariate additive–dominance models
* Bivariate models
* Multivariate genotype × diet models
* GEBV extraction
* Reaction-norm visualisation

The genetic models were used to estimate additive genetic variance, dominance variance, residual variance, heritability, dominance ratio, and genetic, dominance, residual and phenotypic correlations.

### 3. Population Structure (R/PLINK)

*Script: `Genetic_parameters_GxD.r`*

PCA was performed using PLINK on the QC-filtered genotype dataset. The first two principal components were visualised in R using ggplot2 to assess genetic structure and the distribution of individuals across the dietary groups.

The PCA used the same QC-filtered genotype dataset that was used to construct the genomic relationship matrices and for the subsequent genetic analyses.

### 4. Genotype × Diet Interaction (R)

*Script: `Genetic_parameters_GxD.r`*

Genotype × diet interactions were estimated using multivariate models implemented in ASReml-R. For each growth trait, measurements from the three dietary environments were treated as correlated traits:

* **SYK** – soy-okara
* **BSG** – brewers' spent grain
* **FVW** – fruit and vegetable waste

Genetic (co)variances and genetic correlations between dietary environments were estimated to assess the consistency of genetic performance and genotype re-ranking across diets. Pairwise genetic correlations were estimated between:

* SYK and BSG
* SYK and FVW
* BSG and FVW

Lower genetic correlations indicate greater changes in genetic ranking across dietary environments and therefore stronger evidence of G × D interaction. Diet-specific heritabilities were also estimated from the multivariate G × D models.

### 5. GEBV Extraction and Reaction-Norm Plots (R)

*Script: `Genetic_parameters_GxD.r`*

Diet-specific genomic estimated breeding values (GEBVs) were extracted from the multivariate G × D models and exported as CSV files.

Reaction-norm plots were generated using ggplot2 to visualise individual genetic responses across the three dietary environments. For each growth trait, individuals were ranked according to their mean GEBV across diets, and the top 5% were visualised to illustrate genetic differences and changes in genetic performance across dietary environments.

---

## Traits

The primary growth traits analysed were:

* **LBW** – larval body weight (mg)
* **LL** – larval length (mm)
* **LW** – larval width (mm)
* **LSA** – larval surface area (mm²)

## Dietary Treatments

* **SYK** – soy-okara
* **BSG** – brewers' spent grain
* **FVW** – fruit and vegetable waste

---

## Software and Packages

* **SAS** – descriptive statistics and fixed-effect analyses (commercial licence required)
* **PLINK** – genotype quality control and PCA
* **R** – data processing, genetic analyses and visualisation
* **ASReml-R** – genetic mixed models and G × D analyses (commercial licence required)
* **snpReady** – genetic data processing and additive relationship matrix construction
* **AGHmatrix** – additive and dominance relationship matrices
* **ggplot2** – visualisation

<!-- CONFIRM: add the PLINK version (1.9 or 2.0), the R version, and the package versions used. -->

## How to Run

1. Run `Descriptive_statistics.sas` for the descriptive statistics and fixed-effect analyses. It is independent of the R script.
2. Run `Genetic_parameters_GxD.r`. PLINK must be installed and available from R, and the file paths at the top of the script must be changed to point to your own copy of the data. The input format is described in `data.md`.

<!-- CONFIRM: check that the run steps match how the scripts work, and state whether PLINK is called from inside R. -->

## Results

A summary of the results is in `Results.md`.
