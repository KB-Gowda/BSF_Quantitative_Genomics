# Phenotypes and SNP Genotypes of Black Soldier Fly (*Hermetia illucens*) from a Breeding Trial

## Description
This dataset contains phenotypic and SNP genotypic data from a Black Soldier Fly (*Hermetia illucens*) breeding trial. The phenotypic dataset contains 2,097 individuals measured for four growth traits under three dietary treatments.

The dataset supports analyses of genetic parameters, and genotype × diet interactions of growth traits.

---

## Phenotypic Data

- **File:** `Pheno_BSF.txt`
- **Format:** Tab-delimited text file
- **Number of individuals:** 2,097
- **Number of variables:** 22

### Core Variables

| Variable | Description |
| :--- | :--- |
| `Sl_No` | Serial number assigned to each individual |
| `Sample_Id` | Unique identifier for each individual |
| `Tray` | Rearing tray identifier |
| `Diet` | Dietary treatment: `SYK` (soy-okara), `BSG` (brewers' spent grain), or `FVW` (fruit and vegetable waste) |
| `Sampling_Day` | Day on which the individual was phenotyped (13, 14, or 15) |
| `Colour` | Larval colour measured as a grayscale value |
| `Weight` | Larval body weight (LBW), measured in mg |
| `Length` | Larval length (LL), measured in mm |
| `Width` | Larval width (LW), measured in mm |
| `SurfaceArea` | Larval surface area (LSA), measured in mm² |

### Diet-Specific Phenotypes

Phenotypic values are present for the diet to which the individual was assigned and are `NA` for the other diets.

| Variable | Description |
| :--- | :--- |
| `SYK_Weight` | Larval body weight measured under the SYK diet (mg) |
| `BSG_Weight` | Larval body weight measured under the BSG diet (mg) |
| `FVW_Weight` | Larval body weight measured under the FVW diet (mg) |
| `SYK_Length` | Larval length measured under the SYK diet (mm) |
| `BSG_Length` | Larval length measured under the BSG diet (mm) |
| `FVW_Length` | Larval length measured under the FVW diet (mm) |
| `SYK_Width` | Larval width measured under the SYK diet (mm) |
| `BSG_Width` | Larval width measured under the BSG diet (mm) |
| `FVW_Width` | Larval width measured under the FVW diet (mm) |
| `SYK_SurfaceArea` | Larval surface area measured under the SYK diet (mm²) |
| `BSG_SurfaceArea` | Larval surface area measured under the BSG diet (mm²) |
| `FVW_SurfaceArea` | Larval surface area measured under the FVW diet (mm²) |

### Summary of Growth Traits

- **LBW:** Larval body weight (mg)
- **LL:** Larval length (mm)
- **LW:** Larval width (mm)
- **LSA:** Larval surface area (mm²)

---

## SNP Genotypic Data

The corresponding SNP genotype data are provided in PLINK binary format:

- **`BSF_geno_2097_5562.bed`:** Binary SNP genotype data.
- **`BSF_geno_2097_5562.bim`:** SNP marker information, including chromosome, marker ID, position, and allele information.
- **`BSF_geno_2097_5562.fam`:** Individual/sample information corresponding to the genotype dataset.

### Summary Statistics

- **Individuals:** 2,097
- **SNP markers:** 5,562
- **Genotyping rate:** 93.28%
- **Processing tool:** PLINK v1.9

---

## Linking Phenotypic and Genotypic Data

The `Sample_Id` variable in the phenotypic dataset corresponds to the individual identifier in the PLINK `.fam` file. Phenotypic and genotypic records can be matched using `Sample_Id`.

### Example in R

```R
pheno <- read.table(
  "Pheno_BSF.txt",
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE
)

```

---

## Data Use

The dataset was used for analyses of:

1. Genetic parameters of growth traits (`LBW`, `LL`, `LW`, and `LSA`)
2. Genotype × diet interactions across environments (`SYK`, `BSG`, and `FVW`)

---

## Software

* **PLINK v1.9:** Used for SNP genotype processing and quality control.
* **R:** Used for data processing and downstream genetic analyses.
