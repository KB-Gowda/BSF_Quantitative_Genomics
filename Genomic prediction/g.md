# A Genomic Prediction Framework for Growth Traits in Black Soldier Fly (*Hermetia illucens*) Breeding

## Citation

> If you use these scripts or find this work helpful, please cite the associated publication once available:

*Manuscript in preparation. Citation details will be updated upon publication.*

---

## Overview

This repository contains the complete analytical pipeline and reproducible R workflows used to evaluate genomic prediction strategies for larval growth traits in black soldier fly (*Hermetia illucens*). The project investigates predictive performance across:
1. **Dietary environments** (industrial waste streams vs. standard diets)
2. **Genetic architectures** (Additive-only [**Model A**] vs. Additive + Dominance [**Model AD**])
3. **Marker panel densities** (whole-genome vs. density-reduced and trait-targeted panels)
4. **Inter-generational transferability** (forward and backward predictions between temporally separated generations)

The analyses use high-density genotypic and pre-adjusted phenotypic data from a diverse Base breeding population ($n = 2,083$) and an $F_7$ advanced generation ($n = 536$) derived from the same ancestral genetic background.

---

## Key Analyses & Scenarios

### 1. Overall Base Population Prediction (`1.Overall-base/`)
Evaluates 5-fold cross-validation (20 replications) across all four primary growth traits:
- **LBW** – Larval Body Weight (mg)
- **LL** – Larval Length (mm)
- **LW** – Larval Width (mm)
- **LSA** – Larval Surface Area (mm²)

Benchmarked across parametric linear mixed models (GBLUP), variable selection/shrinkage Bayesian alphabets, and non-parametric machine learning.

### 2. Within-Diet Genomic Prediction (`2.Diet-wise-base/`)
Evaluates model performance within individual larval dietary waste streams for **LBW** and **LSA**:
- **SYK** – Soy-okara ($n = 965$)
- **BSG** – Brewers' spent grain ($n = 555$)
- **FVW** – Fruit and vegetable waste ($n = 563$)

### 3. Accounting for Genotype × Diet Interaction (`3.Multivariate-GxD/`)
Evaluates multi-trait models treating performance in SYK, BSG, and FVW as distinct, correlated traits to account for genotype-by-environment ($G \times D$) covariance:
- **GBLUP:** Multi-trait unstructured genetic covariance matrix (`us(trait):vm(Sample_Id, Ginv)`)
- **BayesC:** Multivariate Spike-and-Slab (`BGLR::Multitrait`)
- **Random Forest:** Multi-output joint splitting (`randomForestSRC::rfsrc` with `Multivar() ~ .`)

### 4. Cross-Diet Genomic Prediction (`4.Cross-diet-base/`)
Tests the capacity of models trained on two larval diets to predict breeding/genotypic values in an unobserved third diet:
- **Train:** SYK + BSG $\rightarrow$ **Predict:** FVW
- **Train:** SYK + FVW $\rightarrow$ **Predict:** BSG
- **Train:** BSG + FVW $\rightarrow$ **Predict:** SYK

### 5. Evaluation of Marker Panel Densities (`5.Panels-base/`)
Quantifies predictive decay across reduced-density and targeted SNP panels for **LBW** and **LSA**:
- **Random density reductions:** 4,000, 3,000, 2,000, 1,000, and 500 SNPs (stratified by chromosome length)
- **Targeted panel:** Top 1,000 SNPs selected based on absolute GBLUP-derived allele substitution effects ($\hat{\alpha}$) estimated strictly within training partitions

### 6. Inter-Generational Forward & Backward Prediction (`6.Forward-Backward/`)
Evaluates the cross-generational stability of genomic predictions across temporal and dietary shifts:
- **Forward Prediction (Base $\rightarrow F_7$):**
  - **Scenario 1 (Pooled Baseline):** Models trained on overall Base phenotypes predicting overall $F_7$ performance
  - **Scenario 2 (Diet-Specific & Multi-Diet Averaging):** Models trained on Base diets (SYK, BSG, FVW) predicting $F_7$ performance using individual diet EBVs and composite 3-diet average breeding values
- **Backward Prediction ($F_7 \rightarrow$ Base):**
  - Evaluates models trained on the specialized $F_7$ population ($n = 536$) predicting into the whole Base population ($n = 2,083$) as well as diet-specific Base subsets

---

## Genomic Prediction Models & Software Implementation

### Evaluated Methods
- **GBLUP:** Genomic Best Linear Unbiased Prediction fitted in **ASReml-R** (v4.2)
- **Bayesian Alphabets:** **BayesA**, **BayesB**, **BayesC**, and **Bayesian Lasso (BayesL)** fitted via **BGLR**; **BayesR** fitted via **hibayes**
- **Random Forest (RF):** 
  - Fast univariate implementations via **ranger**
  - Multivariate and joint multi-trait predictions via **randomForestSRC**

### Genetic Effect Partitions
- **Model A (Additive):** Incorporates additive genomic relationships ($\mathbf{G}$) or marker matrices ($\mathbf{X}_a$)
- **Model AD (Additive + Dominance):** Jointly fits orthogonal additive ($\mathbf{G}_a$ / $\mathbf{X}_a$) and dominance relationship/marker matrices ($\mathbf{G}_d$ / $\mathbf{X}_d$) based on the classical Vitezica/Zhu parameterization

### Hyperparameter Tuning
For tree-based models, nested grid searches are evaluated using internal out-of-bag (OOB) error minimization across candidate grids:
- `num.trees` / `ntree`: 150, 200, 300, 400
- `mtry`: candidate subsets evaluated at $\lfloor\sqrt{p}\rfloor$, $\lfloor p/10\rfloor$, $\lfloor p/5\rfloor$, and $\lfloor p/3\rfloor$
- `nodesize` / `min.node.size`: 5, 10

---

## Phenotype Standardization & Validation Metrics

### Phenotypic Pre-Adjustment
Phenotypes are pre-adjusted for non-genetic nuisance effects using linear models fitted in `asreml`, extracting residuals:
- **Base Overall:** `Trait ~ Colour + Sampling_Day + Diet + Tray:Diet`
- **Base Within-Diet:** `Trait ~ Colour + Sampling_Day + Tray + Tray:Sampling_Day`
- **$F_7$ Overall:** `Trait ~ Colour + Tray`

### Validation Schemes
- **5-Fold Cross-Validation (20 Replications):** Applied to Scenarios 1, 2, 3, and 5 using identical cross-validation partitioning folds across all methods
- **Empirical Validation:** Fixed training/testing partitions applied to Cross-Diet (Scenario 4) and Inter-Generational transfers (Scenario 6)

### Evaluation Metrics
1. **Predictive Ability (PA):** Pearson correlation coefficient between observed adjusted phenotypes ($y$) and genomic predictions ($\hat{g}$):
   $$\text{PA} = r(y, \hat{g})$$
2. **Prediction Scale Bias:** Linear regression coefficient ($\beta_1$) of observed phenotypes on predictions ($y = \beta_0 + \beta_1 \hat{g}$). $\beta_1 = 1.0$ indicates unbiased scaling; $\beta_1 < 1.0$ indicates variance inflation.
3. **Mean Squared Error of Prediction (MSEP):** Mean squared difference between centered observed values and centered predictions:
   $$\text{MSEP} = \frac{1}{n}\sum_{i=1}^{n}\Big((y_i - \bar{y}) - (\hat{g}_i - \bar{\hat{g}})\Big)^2$$

---

## Repository Structure & Script Inventory
