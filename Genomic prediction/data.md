# Phenotypes and SNP Genotypes of Australian Black soldier fly (*Hermetia illucens*) from the base and F<sub>7</sub> populations for genomic prediction

## Description
This dataset contains phenotypic and SNP genotypic data from a Black soldier fly (*Hermetia illucens*) breeding colony: a base population and the seventh generation (F<sub>7</sub>) of the James Cook University (JCU) colony. The F<sub>7</sub> population was established from the same original genetic background as the base population.

The dataset supports genomic prediction analyses of growth traits in Australian Black soldier fly larvae (BSFL): within the base population, within diets, across diets, with different SNP panel densities, and forward and backward between the base and F<sub>7</sub> populations.

This file describes the input data used by the scripts in this repository: what each file contains, how the variables are coded, and how many individuals are in each population. The raw data are not included in the repository (see [Data Availability](#data-availability)).

Items marked **CONFIRM** could not be taken from the scripts and need to be completed by the data owner.

### Populations

| Population | Description | Diets | Phenotyped | After genotype QC |
| :--- | :--- | :--- | ---: | ---: |
| **Base** | Base population | `SYK`, `BSG`, `FVW` | 2,097 | 2,083 |
| **F<sub>7</sub>** | Seventh generation of the JCU breeding colony | `CSF` | 576 | 536 |
| **Combined** | Base + F<sub>7</sub> (used for forward and backward prediction) | `SYK`, `BSG`, `FVW`, `CSF` | 2,673 | 2,619 |

Each individual has one phenotypic record (one row) and one genotype. Individuals are matched between phenotype and genotype files by `Sample_Id`.

### Dietary Treatments

| Code | Diet | Population |
| :--- | :--- | :--- |
| `SYK` | Soy-okara | Base |
| `BSG` | Brewers' spent grain | Base |
| `FVW` | Fruit and vegetable waste | Base |
| `CSF` | Commercial chicken starter feed (15.5% crude protein) | F<sub>7</sub> |

Records per diet after QC: `SYK` 965, `BSG` 555, `FVW` 563 (base, total 2,083); `CSF` 536 (F<sub>7</sub>).

---

## Phenotypic Data

- **Files:**
  - `Pheno_BSF.txt`: base population
  - `Pheno_BSF_F7.txt`: F<sub>7</sub> population
  - `Pheno_BSF_2097_576.txt`: the two files combined (created by `FwdBwd_Prep`)
- **Format:** Tab-delimited text files
- **Number of individuals:** 2,097 (base), 576 (F<sub>7</sub>), 2,673 (combined)
- **Number of variables:** 22 (base), <!-- CONFIRM: number of columns in Pheno_BSF_F7.txt --> (F<sub>7</sub>), 28 (combined)

One row per individual.

### Core Variables

| Variable | Type | Description |
| :--- | :--- | :--- |
| `Sl_No` | integer | Serial number assigned to each individual |
| `Sample_Id` | text | Unique identifier for each individual, used to link phenotypes to genotypes (for example `BSF_F0_B1_19A01`) |
| `Generation` | text | `Base` or `F7`. Present in the combined data (`FwdBwd_prep.RData`), not in the base-only data |
| `Tray` | text | Rearing tray identifier. Base trays are coded `B1` to `B3`, `S1` to `S3` and `V1` to `V3` (9 trays). F<sub>7</sub> trays are coded `T1` to `T3` (3 trays), giving 12 trays in the combined data |
| `Diet` | text | Dietary treatment: `SYK`, `BSG`, `FVW` (base) or `CSF` (F<sub>7</sub>) |
| `Sampling_Day` | text | Day on which the individual was phenotyped. Base values are `13`, `14` and `15`. The combined data has a fourth level, which belongs to the F<sub>7</sub> population |
| `Colour` | numeric | Larval colour measured as a grayscale value, continuous (about 120 to 200). Used as a numeric covariate, not as a factor |
| `Weight` | numeric | Larval body weight (LBW), measured in mg |
| `Length` | numeric | Larval length (LL), measured in mm |
| `Width` | numeric | Larval width (LW), measured in mm |
| `SurfaceArea` | numeric | Larval surface area (LSA), measured in mm² |

### Diet-Specific Phenotypes

Phenotypic values are present for the diet to which the individual was assigned and are `NA` for the other diets. For example, `SYK_Weight` equals `Weight` for animals on `SYK`. These columns are used for the diet-wise and multivariate (genotype × diet) analyses.

| Variable | Description |
| :--- | :--- |
| `SYK_Weight`, `BSG_Weight`, `FVW_Weight`, `CSF_Weight` | Larval body weight measured under each diet (mg) |
| `SYK_Length`, `BSG_Length`, `FVW_Length`, `CSF_Length` | Larval length measured under each diet (mm) |
| `SYK_Width`, `BSG_Width`, `FVW_Width`, `CSF_Width` | Larval width measured under each diet (mm) |
| `SYK_SurfaceArea`, `BSG_SurfaceArea`, `FVW_SurfaceArea`, `CSF_SurfaceArea` | Larval surface area measured under each diet (mm²) |

The `CSF_` columns exist only in the combined data.

### Summary of Growth Traits

- **LBW:** Larval body weight (mg)
- **LL:** Larval length (mm)
- **LW:** Larval width (mm)
- **LSA:** Larval surface area (mm²)

### Variables Added by the Scripts

| Variable | Description |
| :--- | :--- |
| `Population` | `Base` or `F7` (factor), added when the two phenotype files are combined. It holds the same information as `Generation`, and the scripts use `Population` |
| `Adj_Weight`, `Adj_SurfaceArea` | Adjusted phenotypes (residuals), all diets of a population together |
| `Adj_Length`, `Adj_Width` | Adjusted Length and Width, in the base-only data (`BasePop_prep.RData`) only |
| `Adj_<diet>_Weight`, `Adj_<diet>_SurfaceArea` | Adjusted phenotypes within one diet (`SYK`, `BSG`, `FVW`). `NA` for animals on other diets and, in the combined file, for all F<sub>7</sub> animals |

**CONFIRM:** the unit and measurement method of `Colour`, the meaning of `Sampling_Day` (for example days after hatching) and its F<sub>7</sub> value, how and when each trait was measured, and the construction of `Sample_Id`.

### Adjusted Phenotypes

Adjusted phenotypes are the residuals from a fixed-effects model with no mean added. They are used as the response in the prediction analyses.

| Data | Model | Fixed effects |
| :--- | :--- | :--- |
| Base population (`BasePop_prep.RData`), all four traits, all diets | Overall | `Colour + Sampling_Day + Diet + Tray:Diet` |
| Base population, Weight and SurfaceArea, each diet separately | Within diet | `Colour + Sampling_Day + Tray + Tray:Sampling_Day` |
| Combined data (`FwdBwd_prep.RData`), base | Overall and within diet | Same as above, fitted on base individuals only |
| Combined data (`FwdBwd_prep.RData`), F<sub>7</sub> | Overall | `Colour + Tray`, fitted on F<sub>7</sub> individuals only |

All models were fitted with ASReml-R. In the combined data, base and F<sub>7</sub> are adjusted separately, and only Weight and SurfaceArea are adjusted. The F<sub>7</sub> population has a single diet, so there is no within-diet adjustment for it.

---

## SNP Genotypic Data

The SNP genotype data are provided in PLINK binary format. They were generated with a custom-designed 6K Allegro targeted SNP panel. Chromosome names are not standard, so PLINK is run with `--allow-extra-chr`.

### Raw Genotype Files

- **`BSF_geno_2097_5562.bed/.bim/.fam`:** base population before QC (2,097 individuals, 5,562 SNPs, genotyping rate 93.28%).
- **`BSF_geno_F7_536.bed/.bim/.fam`:** F<sub>7</sub> population (536 individuals).

Each `.bed` file holds the binary genotypes, each `.bim` file the marker information (chromosome, marker ID, position and alleles), and each `.fam` file the individual information.

### Quality Control and Merging

- **Base population:** SNPs with more than 10% missing genotypes removed (`--geno 0.1`), individuals with more than 10% missing genotypes removed (`--mind 0.1`), and SNPs with minor allele frequency below 0.005 removed (`--maf 0.005`).
  - Base: 2,097 individuals and 5,562 SNPs before QC; 2,083 individuals and 4,680 SNPs after QC.
- **Merging:** the F<sub>7</sub> genotypes (536 individuals) were merged with the quality-controlled base genotypes (`--bmerge`), giving 2,619 individuals and 4,680 SNPs.
- **Processing tool:** PLINK v1.9

### Additive Dosage File (`.raw`)

Produced by PLINK `--recode A`. The first six columns are `FID`, `IID`, `PAT`, `MAT`, `SEX` and `PHENOTYPE`. The remaining columns hold one SNP each, coded 0, 1 or 2 copies of the counted allele. Column names have the counted allele appended (for example `_A`, `_C`, `_G` or `_T`). `IID` matches `Sample_Id` in the phenotype files.

### Marker Coding and Relationship Matrices

- **Missing genotypes:** replaced by the marker mean (mean imputation) in the `snpReady` step. The dominance matrix was also checked without imputation, with nearly identical results.
- **Additive coding:** 0/1/2 dosage (`Xa_final`, mean-imputed `M_impute`).
- **Dominance coding:** orthogonal coding of Vitezica et al. (`Xd_final`), where genotypes 0, 1 and 2 are coded `-2p²`, `2pq` and `-2q²` (*p* and *q* are the allele frequencies).
- **Relationship matrices:** additive `Ga` (VanRaden, via `snpReady`) and dominance `Gd` (Vitezica, via `AGHmatrix`), with their inverses `Ginv` and `Gd.inv` in the format ASReml needs.

**CONFIRM:** the number of chromosomes and the SNP physical positions (from the `.bim` files), the reference genome version, and the genotyping rate after QC.

### SNP Panels

| Panel | Description |
| :--- | :--- |
| 4000, 3000, 2000, 1000, 500 | Random SNPs, sampled in proportion to the number of markers per chromosome. Fixed seed (`set.seed(123)`) |
| `top1000_Weight_A`, `top1000_Weight_AD`, `top1000_SurfaceArea_A`, `top1000_SurfaceArea_AD` | The 1,000 SNPs with the largest variance contribution, from GBLUP (additive, or additive plus dominance) on `Adj_Weight` or `Adj_SurfaceArea` |

---

## Files

All files are expected in the working folder of the scripts.

### Raw Input Files

| File | Content |
| :--- | :--- |
| `Pheno_BSF.txt` | Phenotypes of the base population (2,097 individuals). Tab-delimited text |
| `Pheno_BSF_F7.txt` | Phenotypes of the F<sub>7</sub> population (576 individuals). Tab-delimited text |
| `BSF_geno_2097_5562.bed/.bim/.fam` | PLINK binary genotypes of the base population: 2,097 individuals and 5,562 SNPs before QC |
| `BSF_geno_F7_536.bed/.bim/.fam` | PLINK binary genotypes of the F<sub>7</sub> population (536 individuals) |

### Files Created by the Preparation Scripts

| File | Created by | Content |
| :--- | :--- | :--- |
| `BSF_geno_2083_4680.*` | `BasePop_Prep` | Base genotypes after QC: 2,083 individuals and 4,680 SNPs. `.raw` is the additive (0/1/2) coding |
| `BSF_geno_2619_4680.*` | `FwdBwd_Prep` | Base and F<sub>7</sub> genotypes merged: 2,619 individuals and 4,680 SNPs. `.raw` is the additive (0/1/2) coding |
| `Pheno_BSF_2097_576.txt` | `FwdBwd_Prep` | Base and F<sub>7</sub> phenotype files combined, with a `Population` column added |
| `BasePop_prep.RData` | `BasePop_Prep` | Objects for scenarios 1 to 4 (see [R Data Objects](#r-data-objects)) |
| `BasePop_adjusted_phenotypes.csv` | `BasePop_Prep` | The base-population phenotype table including the adjusted phenotypes |
| `GP_panel_input.RData` | `Panels_input` | Marker panels, their relationship matrices and cross-validation folds (scenario 5) |
| `snps_panel_<n>.txt`, `snps_top1000_<trait>_<model>.txt` | `Panels_input` | SNP lists of the reduced panels |
| `BSF_geno_2083_<panel>.raw` | `Panels_input` | Additive genotypes of each reduced panel |
| `Top1000_SNP_<trait>_<model>.csv` | `Panels_input` | The top-1,000 SNPs with allele frequency, estimated effects and ranking scores |
| `FwdBwd_prep.RData` | `FwdBwd_Prep` | Objects for forward and backward prediction (scenario 6) |
| `FwdBwd_adjusted_phenotypes.csv` | `FwdBwd_Prep` | The combined phenotype table including the adjusted phenotypes |

---

## R Data Objects

### `BasePop_prep.RData` (2,083 individuals)

| Object | Description |
| :--- | :--- |
| `pheno_final` | Phenotypes, including all `Adj_` columns. One row per individual, 32 columns |
| `M` | Raw 0/1/2 genotype matrix (individuals × 4,680 SNPs) |
| `M_impute`, `Xa_final` | Mean-imputed additive marker matrix |
| `Xd_final` | Dominance marker matrix |
| `Ga`, `Ginv` | Additive genomic relationship matrix and its inverse |
| `Gd`, `Gd.inv` | Dominance relationship matrix and its inverse |
| `trait_info`, `diet_trait_info` | Tables linking each trait to its adjusted column |
| `diets`, `primary_traits` | Diet codes and trait names |
| `cv_matrix` | Cross-validation fold labels (1 to 5), 2,083 individuals × 20 replications. Folds are stratified by diet, with `set.seed(421)` |

### `GP_panel_input.RData` (2,083 individuals)

`pheno_final`, `panel_mats` (a list with, for each of 9 panels, the dosage matrix `Xa`, `Ginv`, `Gd.inv` and the SNP count), `panel_names`, `cv_matrix` (same folds as the full panel), `trait_info`, `diet_trait_info`, `diets`, `primary_traits`.

### `FwdBwd_prep.RData` (2,619 individuals)

| Object | Description |
| :--- | :--- |
| `pheno_final` | Combined phenotypes, 36 columns. It has `Generation` and `Population`, the diet-specific columns for all four diets (`SYK`, `BSG`, `FVW`, `CSF`), `Adj_Weight`, `Adj_SurfaceArea` and the six within-diet `Adj_` columns. It has no `Adj_Length` or `Adj_Width` |
| `M`, `M_impute`, `Xa_final`, `Xd_final` | Marker matrices for all 2,619 individuals |
| `Ga`, `Ginv`, `Gd`, `Gd.inv` | Relationship matrices built on all 2,619 individuals together, identical for both scenarios |
| `idx_base`, `idx_f7` | Row numbers of the base (2,083) and F<sub>7</sub> (536) individuals |
| `scenarios` | `forward` (train = base, valid = F<sub>7</sub>) and `backward` (train = F<sub>7</sub>, valid = base), as row numbers |
| `trait_info`, `diet_trait_info`, `diets`, `primary_traits`, `adj_traits` | Trait and diet definitions |

In every object, rows of the phenotype table and rows of all matrices are in the same order (the scripts check this).

---

## Linking Phenotypic and Genotypic Data

The `Sample_Id` variable in the phenotypic datasets corresponds to the individual identifier (`IID`) in the PLINK `.fam` and `.raw` files. Phenotypic and genotypic records can be matched using `Sample_Id`. Of the 2,673 phenotyped individuals, 2,619 have genotype data after QC (2,083 base and 536 F<sub>7</sub>).

### Example in R

```R
library(data.table)

pheno <- read.table(
  "Pheno_BSF_2097_576.txt",
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE
)

geno <- fread("BSF_geno_2619_4680.raw")

common_ids  <- intersect(pheno$Sample_Id, geno$IID)
pheno_match <- pheno[pheno$Sample_Id %in% common_ids, ]
geno_match  <- geno[match(pheno_match$Sample_Id, geno$IID), ]

# Marker matrix (columns 1 to 6 are PLINK identifiers)
M <- as.matrix(geno_match[, -(1:6), with = FALSE])
rownames(M) <- pheno_match$Sample_Id
```

---

## Data Use

The dataset was used for genomic prediction of growth traits in six scenarios:

1. Genomic prediction in the base population (`LBW`, `LL`, `LW` and `LSA`)
2. Within-diet genomic prediction (`SYK`, `BSG` and `FVW`)
3. Genotype × diet interactions, with multivariate models
4. Cross-diet prediction
5. SNP marker panel density
6. Forward and backward prediction between the base and F<sub>7</sub> populations

---

## Notes on Use

* **Diets are confounded with population.** The base population has `SYK`, `BSG` and `FVW`; F<sub>7</sub> has `CSF` only. In forward and backward prediction, differences between the populations include both generation and diet.
* **Trays are population-specific.** Tray codes (`B1` to `B3`, `S1` to `S3` and `V1` to `V3` for base, `T1` to `T3` for F<sub>7</sub>) do not overlap, and the adjusted phenotypes are computed separately in each population.
* **`Colour` is numeric.** It should not be converted to a factor.
* **Combined file name.** The combined phenotype file `Pheno_BSF_2097_576.txt` is named for the 2,097 base and 576 F<sub>7</sub> phenotype records. Only 2,619 of them match a genotype after QC.
* **Dominance and relationship matrices** are computed once on the full set of individuals, not within training folds.

---

## Software

* **PLINK v1.9:** Used for SNP genotype processing, quality control and merging.
* **R:** Used for data processing and downstream genetic analyses.
* **ASReml-R:** Used to adjust phenotypes and for GBLUP.
* **snpReady** and **AGHmatrix:** Used to build the additive and dominance genomic relationship matrices.

---

## Data Availability

Raw genotype and phenotype datasets are not included in this repository.

<!-- CONFIRM: add how to request the data, if it can be shared. -->
