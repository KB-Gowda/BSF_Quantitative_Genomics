# A genomic prediction framework for growth traits in black soldier fly (*Hermetia illucens*) breeding

## Citation

> If you use these scripts or find this work helpful, please cite the associated publication once available.

*Manuscript in preparation. Citation details will be added following publication.*

---

This repository contains the analysis scripts used to evaluate genomic prediction methods for growth traits in black soldier fly larvae (*Hermetia illucens*). The analyses investigate prediction performance across dietary environments, genetic effect models, marker panel densities, and populations from different generations to develop a genomic prediction framework for BSF breeding.

The study uses phenotypic and genotypic data from a base population and an $\text{F}_7$ population established from the same original genetic background.

---

## Key Analyses

### 1. Genomic Prediction in the Base Population
Genomic prediction was evaluated in the base population for four growth traits:
* **LBW** – Larval body weight
* **LL** – Larval length
* **LW** – Larval width
* **LSA** – Larval surface area

Different genomic prediction methods and genetic effect models were evaluated to assess their performance in predicting larval growth traits.

### 2. Within-Diet Genomic Prediction
Genomic prediction was evaluated separately within each of the three dietary environments in the base population:
* **SYK** – Soy-okara
* **BSG** – Brewers' spent grain
* **FVW** – Fruit and vegetable waste

This scenario focused on **LBW** and **LSA** to assess prediction performance within individual dietary environments.

### 3. Genomic Prediction Accounting for Genotype × Diet Interaction
Genomic prediction accounting for genotype $\times$ diet ($\text{G} \times \text{D}$) interactions was evaluated using multivariate models.

For LBW and LSA, observations from the three dietary environments were treated as separate but correlated traits. The models incorporated genetic relationships between diets to evaluate genomic prediction while accounting for $\text{G} \times \text{D}$ interactions.

### 4. Cross-Diet Genomic Prediction
Cross-diet genomic prediction was evaluated in the base population by training models using individuals from two dietary environments and predicting individuals from the remaining environment.

The following combinations were considered:
* **Train:** SYK & BSG $\rightarrow$ **Predict:** FVW
* **Train:** SYK & FVW $\rightarrow$ **Predict:** BSG
* **Train:** BSG & FVW $\rightarrow$ **Predict:** SYK

This scenario evaluated the potential to predict growth traits in a target dietary environment using information from other dietary environments.

### 5. Evaluation of SNP Marker Panel Density
The effect of marker panel density on genomic prediction performance was evaluated in the base population.

The full panel of $4,680$ SNPs was compared with reduced panels containing:
* $4,000$ SNPs
* $3,000$ SNPs
* $2,000$ SNPs
* $1,000$ SNPs
* $500$ SNPs

For the reduced panels, markers were randomly selected while accounting for the total number of markers on each chromosome. An additional panel of the top $1,000$ markers was selected based on GBLUP-derived SNP effects using the method described in the study methodology.

### 6. Forward and Backward Genomic Prediction
Forward and backward genomic prediction were evaluated between the base and $\text{F}_7$ populations.

* **Forward Prediction:**
  * **Training Population:** Base population
  * **Prediction Population:** $\text{F}_7$ population
* **Backward Prediction:**
  * **Training Population:** $\text{F}_7$ population
  * **Prediction Population:** Base population

These analyses evaluated the transferability of genomic prediction between populations from different generations.

---

## Genomic Prediction Methods & Implementation

### Evaluated Methods
* Genomic Best Linear Unbiased Prediction (**GBLUP**)
* **BayesA**
* **BayesB$\pi$**
* **BayesC$\pi$**
* **Bayesian Lasso**
* **BayesR**
* Random Forest (**RF**)

Additive (**A**) and additive–dominance (**AD**) genetic effect models were evaluated where applicable.

### Model Implementation Details
* **GBLUP:** Implemented in `ASReml-R` using additive and additive–dominance genomic relationship models.
* **Bayesian Models:** `BayesA`, `BayesB$\pi$`, `BayesC$\pi$`, and `Bayesian Lasso` were implemented using `BGLR`. `BayesR` was implemented using `hibayes`.
* **Random Forest:** Univariate models were implemented using `ranger`, whereas multivariate models were implemented using `randomForestSRC`.
* **Multivariate Prediction:** Multivariate GBLUP, Bayesian, and Random Forest approaches were used to account for genetic relationships between dietary environments.

> **Note on Hyperparameters:** Random forest hyperparameters, including the number of trees ($150, 300$), node size ($5, 10$), and number of candidate split variables ($\text{mtry}$), were optimized within each training dataset using nested $5$-fold cross-validation.

---

## Validation & Evaluation Metrics

Random $5$-fold cross-validation was performed with $20$ replications, except for cross-diet and forward/backward prediction, where the training and validation populations were predefined.

Prediction performance was evaluated using three metrics:
1. **Predictive Ability (PA):** Pearson correlation coefficient between adjusted phenotypes and genomic estimated breeding values (GEBVs) in the validation population.
2. **Mean Squared Error of Prediction (MSEP):** Mean squared difference between adjusted phenotypes and GEBVs.
3. **Prediction Scale Bias:** Regression slope of adjusted phenotypes on GEBVs in the validation population.

---

## Traits & Environments Summary

### Traits Analysed
* **LBW:** Larval body weight ($\text{mg}$)
* **LL:** Larval length ($\text{mm}$)
* **LW:** Larval width ($\text{mm}$)
* **LSA:** Larval surface area ($\text{mm}^2$)

### Dietary Treatments
* **SYK:** Soy-okara
* **BSG:** Brewers' spent grain
* **FVW:** Fruit and vegetable waste
* **CSF:** Commercial chicken starter feed ($15.5\%$ crude protein)

### Study Populations & Genotyping
* **Base Population:** $2,083$ individuals after genotype quality control, reared across three dietary environments (SYK, BSG, and FVW).
* **$\text{F}_7$ Population:** $536$ individuals after genotype quality control, from the seventh generation of the JCU breeding colony and reared on commercial chicken starter feed (CSF; $15.5\%$ crude protein).

Genotyping was performed using a custom-designed $\text{6K}$ Allegro targeted SNP panel. Following genotype quality control, $4,680$ SNP markers were retained across both populations.

---

## Software & Packages

* **R:** Data processing, genomic prediction, and statistical analyses
* **PLINK:** Genotype quality control and population structure analysis
* **ASReml-R:** GBLUP and genomic mixed models
* **BGLR:** Bayesian genomic prediction models
* **hibayes:** BayesR implementation
* **ranger:** Univariate random forest models
* **randomForestSRC:** Multivariate random forest models
* **snpReady:** Additive genomic relationship matrix construction
* **AGHmatrix:** Additive and dominance relationship matrices

---

## Repository Structure

The repository contains scripts organized according to genomic prediction methods and analytical scenarios:

1. Base population genomic prediction
2. Within-diet genomic prediction
3. Multivariate genomic prediction accounting for $\text{G} \times \text{D}$ interactions
4. Cross-diet genomic prediction
5. SNP marker panel density evaluation
6. Forward and backward genomic prediction between the base and $\text{F}_7$ populations

The individual scripts contain the model specifications, validation procedures, and analysis settings for the corresponding scenarios.

> *Note: Raw genotype and phenotype datasets are not included in this repository.*
