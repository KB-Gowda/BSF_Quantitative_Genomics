# Data description: A genomic prediction framework for growth traits in black soldier fly (*Hermetia illucens*) breeding

This file describes the input data used by the scripts in this repository: what each file contains, how the variables are coded, and how many individuals are in each population. The raw data are not included in the repository (see [Data availability](#data-availability)).

Items marked **CONFIRM** could not be taken from the scripts and need to be completed by the data owner.

---

## 1. Overview

| Population | Description | Diets | Phenotyped | After genotype QC |
|---|---|---|---|---|
| **Base** | Base population | SYK, BSG, FVW | 2,097 | 2,083 |
| **F<sub>7</sub>** | Seventh generation of the JCU breeding colony | CSF | 576 | 536 |
| **Combined** | Base + F<sub>7</sub> (used for forward and backward prediction) | SYK, BSG, FVW, CSF | 2,673 | 2,619 |

Each individual has one phenotypic record (one row) and one genotype. Individuals are matched between phenotype and genotype files by `Sample_Id`.

Dietary treatments:

| Code | Diet | Population |
|---|---|---|
| SYK | Soy-okara | Base |
| BSG | Brewers' spent grain | Base |
| FVW | Fruit and vegetable waste | Base |
| CSF | Commercial chicken starter feed (15.5% crude protein) | F<sub>7</sub> |

Records per diet after QC: SYK 965, BSG 555, FVW 563 (Base, total 2,083); CSF 536 (F<sub>7</sub>).

---

## 2. Files

All files are expected in the working folder of the scripts.

### 2.1 Raw input files

| File | Content |
|---|---|
| `Pheno_BSF.txt` | Phenotypes of the base population (2,097 individuals). Tab-delimited text. |
| `Pheno_BSF_F7.txt` | Phenotypes of the F<sub>7</sub> population (576 individuals). Tab-delimited text. |
| `BSF_geno_2097_5562.bed/.bim/.fam` | PLINK binary genotypes of the base population: 2,097 individuals and 5,562 SNPs before QC. |
| `BSF_geno_F7_536.bed/.bim/.fam` | PLINK binary genotypes of the F<sub>7</sub> population (536 individuals). |

### 2.2 Files created by the preparation scripts

| File | Created by | Content |
|---|---|---|
| `BSF_geno_2083_4680.*` | `BasePop_Prep` | Base genotypes after QC: 2,083 individuals and 4,680 SNPs. `.raw` is the additive (0/1/2) coding. |
| `BSF_geno_2619_4680.*` | `FwdBwd_Prep` | Base and F<sub>7</sub> genotypes merged: 2,619 individuals and 4,680 SNPs. `.raw` is the additive (0/1/2) coding. |
| `Pheno_BSF_2097_576.txt` | `FwdBwd_Prep` | Base and F<sub>7</sub> phenotype files combined, with a `Population` column added. |
| `BasePop_prep.RData` | `BasePop_Prep` | Objects for scenarios 1 to 4 (section 5). |
| `BasePop_adjusted_phenotypes.csv` | `BasePop_Prep` | The base-population phenotype table including the adjusted phenotypes. |
| `GP_panel_input.RData` | `Panels_input` | Marker panels, their relationship matrices and cross-validation folds (scenario 5). |
| `snps_panel_<n>.txt`, `snps_top1000_<trait>_<model>.txt` | `Panels_input` | SNP lists of the reduced panels. |
| `BSF_geno_2083_<panel>.raw` | `Panels_input` | Additive genotypes of each reduced panel. |
| `Top1000_SNP_<trait>_<model>.csv` | `Panels_input` | The top-1,000 SNPs with allele frequency, estimated effects and ranking scores. |
| `FwdBwd_prep.RData` | `FwdBwd_Prep` | Objects for forward and backward prediction (scenario 6). |
| `FwdBwd_adjusted_phenotypes.csv` | `FwdBwd_Prep` | The combined phenotype table including the adjusted phenotypes. |

---

## 3. Phenotype files

One row per individual. Columns:

| Column | Type | Description |
|---|---|---|
| `Sl_No` | integer | Row number (serial number) in the phenotype file. |
| `Sample_Id` | text | Unique individual ID, used to link phenotypes to genotypes (for example `BSF_F0_B1_19A01`). |
| `Generation` | text | `Base` or `F7`. Present in the combined data (`FwdBwd_prep.RData`), not in the base-only data. |
| `Tray` | text | Rearing tray. Base trays are coded `B1` to `B9` (9 trays). F<sub>7</sub> trays are coded `T1` to `T3` (3 trays), giving 12 trays in the combined data. |
| `Diet` | text | Diet code: `SYK`, `BSG`, `FVW` (Base) or `CSF` (F<sub>7</sub>). |
| `Sampling_Day` | text | Sampling day of the record. Base values are `13`, `14` and `15`. The combined data has a fourth level, which belongs to the F<sub>7</sub> population. |
| `Colour` | numeric | Larval colour score, continuous. Values are about 120 to 200. Used as a numeric covariate, not as a factor. |
| `Weight` | numeric | Larval body weight (LBW), mg. |
| `Length` | numeric | Larval length (LL), mm. |
| `Width` | numeric | Larval width (LW), mm. |
| `SurfaceArea` | numeric | Larval surface area (LSA), mm². |
| `SYK_<trait>`, `BSG_<trait>`, `FVW_<trait>`, `CSF_<trait>` | numeric | Diet-specific copies of the four traits (`Weight`, `Length`, `Width`, `SurfaceArea`). `CSF_` columns exist only in the combined data. For example `SYK_Weight` equals `Weight` for animals on SYK and is `NA` for all other animals. These columns are used for the diet-wise and multivariate (G × D) analyses. |

Added by the scripts:

| Column | Description |
|---|---|
| `Population` | `Base` or `F7` (factor), added when the two phenotype files are combined. It holds the same information as `Generation`, and the scripts use `Population`. |
| `Adj_Weight`, `Adj_SurfaceArea` | Adjusted phenotypes (residuals), all diets of a population together. |
| `Adj_Length`, `Adj_Width` | Adjusted Length and Width, in the base-only data (`BasePop_prep.RData`) only. |
| `Adj_<diet>_Weight`, `Adj_<diet>_SurfaceArea` | Adjusted phenotypes within one diet (`SYK`, `BSG`, `FVW`). `NA` for animals on other diets and, in the combined file, for all F<sub>7</sub> animals. |

**CONFIRM:** the unit and measurement method of `Colour`, the meaning of `Sampling_Day` (for example days after hatching) and its F<sub>7</sub> value, how and when each trait was measured, and the construction of `Sample_Id`.

### Adjusted phenotypes

Adjusted phenotypes are the residuals from a fixed-effects model with no mean added. They are used as the response in the prediction analyses.

| Data | Models | Fixed effects |
|---|---|---|
| Base population (`BasePop_prep.RData`), all four traits, all diets | Overall | `Colour + Sampling_Day + Diet + Tray:Diet` |
| Base population, Weight and SurfaceArea, each diet separately | Within diet | `Colour + Sampling_Day + Tray + Tray:Sampling_Day` |
| Combined data (`FwdBwd_prep.RData`), Base | Overall and within diet | Same as above, fitted on Base individuals only |
| Combined data (`FwdBwd_prep.RData`), F<sub>7</sub> | Overall | `Colour + Tray`, fitted on F<sub>7</sub> individuals only |

All models were fitted with ASReml-R. In the combined data, Base and F<sub>7</sub> are adjusted separately, and only Weight and SurfaceArea are adjusted. The F<sub>7</sub> population has a single diet, so there is no within-diet adjustment for it.

---

## 4. Genotype data

* **Panel:** custom-designed 6K Allegro targeted SNP panel.
* **Format:** PLINK binary files (`.bed`, `.bim`, `.fam`). Chromosome names are not standard, so PLINK is run with `--allow-extra-chr`.
* **Quality control (base population):** SNPs with more than 10% missing genotypes removed (`--geno 0.1`), individuals with more than 10% missing genotypes removed (`--mind 0.1`), and SNPs with minor allele frequency below 0.005 removed (`--maf 0.005`).
  * Base: 2,097 individuals and 5,562 SNPs before QC, 2,083 individuals and 4,680 SNPs after QC.
  * F<sub>7</sub>: 536 individuals, merged with the QC'd base genotypes, giving 2,619 individuals and 4,680 SNPs.
* **Additive dosage file (`.raw`):** produced by PLINK `--recode A`. The first six columns are `FID`, `IID`, `PAT`, `MAT`, `SEX` and `PHENOTYPE`. The remaining columns hold one SNP each, coded 0, 1 or 2 copies of the counted allele. Column names have the counted allele appended (for example `_A`, `_C`, `_G` or `_T`). `IID` matches `Sample_Id` in the phenotype files.
* **Missing genotypes:** replaced by the marker mean (mean imputation) in the `snpReady` step. The dominance matrix was also checked without imputation, with nearly identical results.
* **Marker coding used in the models:**
  * Additive: 0/1/2 dosage (`Xa_final`, mean-imputed `M_impute`).
  * Dominance: orthogonal coding of Vitezica et al. (`Xd_final`), where genotype 0, 1 and 2 are coded `-2p²`, `2pq` and `-2q²` (*p* and *q* are the allele frequencies).
* **Relationship matrices:** additive `Ga` (VanRaden, via `snpReady`) and dominance `Gd` (Vitezica, via `AGHmatrix`), with their inverses `Ginv` and `Gd.inv` in the format ASReml needs.

**CONFIRM:** the number of chromosomes and the SNP physical positions (from the `.bim` files), and the reference genome version.

### SNP panels (scenario 5)

| Panel | Description |
|---|---|
| 4000, 3000, 2000, 1000, 500 | Random SNPs, sampled in proportion to the number of markers per chromosome. Fixed seed (`set.seed(123)`). |
| `top1000_Weight_A`, `top1000_Weight_AD`, `top1000_SurfaceArea_A`, `top1000_SurfaceArea_AD` | The 1,000 SNPs with the largest variance contribution, from GBLUP (additive, or additive plus dominance) on `Adj_Weight` or `Adj_SurfaceArea`. |

---

## 5. R data objects

### `BasePop_prep.RData` (2,083 individuals)

| Object | Description |
|---|---|
| `pheno_final` | Phenotypes, including all `Adj_` columns. One row per individual, 32 columns. |
| `M` | Raw 0/1/2 genotype matrix (individuals × 4,680 SNPs). |
| `M_impute`, `Xa_final` | Mean-imputed additive marker matrix. |
| `Xd_final` | Dominance marker matrix. |
| `Ga`, `Ginv` | Additive genomic relationship matrix and its inverse. |
| `Gd`, `Gd.inv` | Dominance relationship matrix and its inverse. |
| `trait_info`, `diet_trait_info` | Tables linking each trait to its adjusted column. |
| `diets`, `primary_traits` | Diet codes and trait names. |
| `cv_matrix` | Cross-validation fold labels (1 to 5), 2,083 individuals × 20 replications. Folds are stratified by diet, with `set.seed(421)`. |

### `GP_panel_input.RData` (2,083 individuals)

`pheno_final`, `panel_mats` (a list with, for each of 9 panels, the dosage matrix `Xa`, `Ginv`, `Gd.inv` and the SNP count), `panel_names`, `cv_matrix` (same folds as the full panel), `trait_info`, `diet_trait_info`, `diets`, `primary_traits`.

### `FwdBwd_prep.RData` (2,619 individuals)

| Object | Description |
|---|---|
| `pheno_final` | Combined phenotypes, 36 columns. It has `Generation` and `Population`, the diet-specific columns for all four diets (SYK, BSG, FVW, CSF), `Adj_Weight`, `Adj_SurfaceArea` and the six within-diet `Adj_` columns. It has no `Adj_Length` or `Adj_Width`. |
| `M`, `M_impute`, `Xa_final`, `Xd_final` | Marker matrices for all 2,619 individuals. |
| `Ga`, `Ginv`, `Gd`, `Gd.inv` | Relationship matrices built on all 2,619 individuals together, identical for both scenarios. |
| `idx_base`, `idx_f7` | Row numbers of the Base (2,083) and F<sub>7</sub> (536) individuals. |
| `scenarios` | `forward` (train = Base, valid = F<sub>7</sub>) and `backward` (train = F<sub>7</sub>, valid = Base), as row numbers. |
| `trait_info`, `diet_trait_info`, `diets`, `primary_traits`, `adj_traits` | Trait and diet definitions. |

In every object, rows of the phenotype table and rows of all matrices are in the same order (the scripts check this).

---

## 6. Notes on use

* **Diets are confounded with population.** The Base population has SYK, BSG and FVW; F<sub>7</sub> has CSF only. In forward and backward prediction, differences between the populations include both generation and diet.
* **Trays are population-specific.** Tray codes (`B1` to `B3`, `S1` to `S3`, `V1` to `V3`for Base, `T1` to `T3` for F<sub>7</sub>) do not overlap, and the adjusted phenotypes are computed separately in each population.
* **Combined file name.** The combined phenotype file `Pheno_BSF_2097_576.txt` is named for the 2,097 base and 576 F<sub>7</sub> phenotype records. Only 2,619 of them match a genotype after QC.
* **Dominance and relationship matrices** are computed once on the full set of individuals, not within training folds.

---

## Data availability

Raw genotype and phenotype datasets are not included in this repository.

<!-- CONFIRM: add how to request the data, if it can be shared. -->
