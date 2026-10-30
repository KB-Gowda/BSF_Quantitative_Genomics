# A genomic prediction framework for growth traits in black soldier fly (*Hermetia illucens*) breeding

## Citation

> If you use these scripts or find this work helpful, please cite the associated publication once available.

*Manuscript in preparation. Citation details will be added following publication.*

---

This folder contains the analysis scripts used to evaluate genomic prediction methods for growth traits in black soldier fly larvae (*Hermetia illucens*). The analyses investigate prediction performance across dietary environments, genetic effect models, marker panel densities, and populations from different generations to develop a genomic prediction framework for BSF breeding.

The study uses phenotypic and genotypic data from a base population and an F<sub>7</sub> population established from the same original genetic background.

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
Genomic prediction accounting for genotype × diet (G × D) interactions was evaluated using multivariate models.

For LBW and LSA, observations from the three dietary environments were treated as separate but correlated traits. The models incorporated genetic relationships between diets to evaluate genomic prediction while accounting for G × D interactions. Multivariate models were implemented for GBLUP, BayesC and random forest.

### 4. Cross-Diet Genomic Prediction
Cross-diet genomic prediction was evaluated in the base population by training models using individuals from two dietary environments and predicting individuals from the remaining environment.

The following combinations were considered:
* **Train:** SYK & BSG → **Predict:** FVW
* **Train:** SYK & FVW → **Predict:** BSG
* **Train:** BSG & FVW → **Predict:** SYK

This scenario evaluated the potential to predict growth traits in a target dietary environment using information from other dietary environments.

### 5. Evaluation of SNP Marker Panel Density
The effect of marker panel density on genomic prediction performance was evaluated in the base population for **LBW** and **LSA**.

The full panel of 4,680 SNPs was compared with reduced panels containing:
* 4,000 SNPs
* 3,000 SNPs
* 2,000 SNPs
* 1,000 SNPs
* 500 SNPs

For the reduced panels, markers were randomly selected while accounting for the total number of markers on each chromosome. An additional panel of the top 1,000 markers was selected based on GBLUP-derived SNP effects using the method described in the study methodology.

<!-- CONFIRM: say whether the GBLUP SNP effects used to pick the top 1,000 markers were estimated within each training fold only. If they were estimated on all individuals before cross-validation, validation individuals influenced the marker choice and predictive ability is inflated. Add one sentence here stating which was done. -->

### 6. Forward and Backward Genomic Prediction
Forward and backward genomic prediction were evaluated between the base and F<sub>7</sub> populations.

* **Forward Prediction:**
  * **Training Population:** Base population
  * **Prediction Population:** F<sub>7</sub> population
* **Backward Prediction:**
  * **Training Population:** F<sub>7</sub> population
  * **Prediction Population:** Base population

These analyses evaluated the transferability of genomic prediction between populations from different generations.

---

## Genomic Prediction Methods & Implementation

### Evaluated Methods
* Genomic Best Linear Unbiased Prediction (**GBLUP**)
* **BayesA**
* **BayesBπ**
* **BayesCπ**
* **Bayesian Lasso**
* **BayesR**
* Random Forest (**RF**)

Additive (**A**) and additive–dominance (**AD**) genetic effect models were evaluated where applicable. Dominance marker coding was checked with and without imputation of missing genotypes, and the results were nearly identical.

<!-- CONFIRM: list which scenarios use the AD model and which use A only. -->

### Model Implementation Details
* **GBLUP:** Implemented in ASReml-R using additive and additive–dominance genomic relationship models.
* **Bayesian Models:** BayesA, BayesBπ, BayesCπ and Bayesian Lasso were implemented using BGLR. BayesR was implemented using hibayes.
* **Random Forest:** Univariate models were implemented using ranger, whereas multivariate models were implemented using randomForestSRC.
* **Multivariate Prediction:** Multivariate GBLUP, BayesC and random forest approaches were used to account for genetic relationships between dietary environments.

<!-- CONFIRM: BayesBπ and BayesCπ are the right names only if π was estimated (in BGLR this depends on the probIn and counts settings). Otherwise use BayesB and BayesC. -->

> **Note on Hyperparameters:** Random forest hyperparameters, including the number of trees (150, 300), node size (5, 10), and number of candidate split variables (mtry), were optimized within each training dataset using nested 5-fold cross-validation.

<!-- CONFIRM: add the mtry values that were tried. -->

---

## Validation & Evaluation Metrics

Random 5-fold cross-validation was performed with 20 replications, except for cross-diet and forward/backward prediction, where the training and validation populations were predefined. The same cross-validation folds were used for all methods within each scenario.

<!-- CONFIRM: the sentence about shared folds must be true for every scenario. -->

Prediction performance was evaluated using three metrics:
1. **Predictive Ability (PA):** Pearson correlation coefficient between adjusted phenotypes and genomic estimated breeding values (GEBVs) in the validation population.
2. **Mean Squared Error of Prediction (MSEP):** Mean squared difference between adjusted phenotypes and GEBVs.
3. **Prediction Scale Bias:** Regression slope of adjusted phenotypes on GEBVs in the validation population. A slope of 1 indicates no bias.

---

## Traits & Environments Summary

### Traits Analysed
* **LBW:** Larval body weight (mg)
* **LL:** Larval length (mm)
* **LW:** Larval width (mm)
* **LSA:** Larval surface area (mm²)

### Dietary Treatments
* **SYK:** Soy-okara (base population)
* **BSG:** Brewers' spent grain (base population)
* **FVW:** Fruit and vegetable waste (base population)
* **CSF:** Commercial chicken starter feed, 15.5% crude protein (F<sub>7</sub> population only)

### Study Populations & Genotyping
* **Base Population:** 2,083 individuals after genotype quality control, reared across three dietary environments (SYK, BSG, and FVW).
* **F<sub>7</sub> Population:** 536 individuals after genotype quality control, from the seventh generation of the JCU breeding colony and reared on commercial chicken starter feed (CSF; 15.5% crude protein).

Genotyping was performed using a custom-designed 6K Allegro targeted SNP panel. Following genotype quality control, 4,680 SNP markers were retained across both populations.

---

## Software & Packages

* **R:** Data processing, genomic prediction, and statistical analyses
* **PLINK:** Genotype quality control and population structure analysis
* **ASReml-R:** GBLUP and genomic mixed models (requires a commercial licence)
* **BGLR:** Bayesian genomic prediction models
* **hibayes:** BayesR implementation
* **ranger:** Univariate random forest models
* **randomForestSRC:** Multivariate random forest models
* **snpReady:** Additive genomic relationship matrix construction
* **AGHmatrix:** Additive and dominance relationship matrices

<!-- CONFIRM: add the R version and package versions used. -->

---

## Repository Structure

Scripts are grouped by analysis scenario. Each folder has one script per method.

| Folder | Scenario | Scripts |
|---|---|---|
| `Overall-base/` | 1. Genomic prediction in the base population (all diets pooled) | `GBLUP_CV_overall.r`, `BayesA_CV_overall.r`, `BayesB_CV_overall.r`, `BayesC_CV_overall.r`, `BayesL_CV_overall.r`, `BayesR_CV_overall.r`, `RF_CV_overall.r`; data preparation in `BasePop_Prep.md` |
| `Diet-wise-base/` | 2. Within-diet genomic prediction (SYK, BSG and FVW analysed separately) | `GBLUP_CV_dietwise.r`, `BayesA_CV_dietwise.r`, `BayesB_CV_dietwise.r`, `BayesC_CV_dietwise.r`, `BayesL_CV_dietwise.r`, `BayesR_CV_dietwise.r`, `RF_CV_dietwise.r` |
| `Multivariate-GxD/` | 3. Multivariate models accounting for genotype × diet interaction | `GBLUP_CV_GxD.r`, `BayesC_CV_GxD.r`, `RF_CV_GxD.r` |
| `CrossDiet-base/` | 4. Cross-diet prediction | `<Method>_CV_crossdiet.r` for each method |
| `Panels-base/` | 5. SNP marker panel density (4,000 to 500 SNPs and the top 1,000 panel) | `GBLUP_CV_Panels.r`, `BayesA_CV_Panels.r`, `BayesB_CV_Panels.r`, `BayesC_CV_Panels.r`, `BayesL_CV_Panels.r`, `BayesR_CV_Panels.r`, `RF_CV_Panels.r`; panel construction in `Panels_input.md` |
| `Forward-backward/` | 6. Forward and backward prediction between the base and F<sub>7</sub> populations | `<Method>_CV_forwardbackward.r` for each method |

<!-- CONFIRM: use the real folder and file names for scenarios 4 and 6, and list only the methods that have a script. -->

### Script naming

`<Method>_CV_<scenario>.r`, where the method is one of:

| Name in file | Method |
|---|---|
| `GBLUP` | GBLUP (ASReml-R) |
| `BayesA`, `BayesB`, `BayesC`, `BayesL` | BayesA, BayesBπ, BayesCπ and Bayesian Lasso (BGLR) |
| `BayesR` | BayesR (hibayes) |
| `RF` | Random forest (ranger for univariate models, randomForestSRC for multivariate models) |

### Suggested order to run

1. Prepare the base population data (`Overall-base/BasePop_Prep.md`).
2. For the panel scenario, build the marker panels and fold assignments first (`Panels-base/Panels_input.md`).
3. Run the method scripts in each scenario folder. Scripts within a folder are independent of each other.

<!-- CONFIRM: check that this order is right. -->

---

## Data Availability

Raw genotype and phenotype datasets are not included in this repository.

<!-- CONFIRM: add how to request the data, if it can be shared. -->
