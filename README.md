# BSF Quantitative Genomics

**Part of the PhD thesis: *Understanding the Genetic Architecture of Black Soldier Fly Production Traits for Enhanced Breeding Applications***

**Author:** Kishor B. Gowda, PhD Candidate, James Cook University (JCU), Townsville, Australia

Code, analysis scripts and results for a quantitative genetics and genomics study of the black soldier fly (*Hermetia illucens*, BSF). The project links phenotypes, and genome-wide SNP data to understand the genetic basis of larval growth traits, and to support genetic improvement in BSF breeding programmes.

> **Status:** work in progress. Chapters 1 and 2 are available. Chapters 3 and 4 are in progress and will be added.

---

## Thesis chapters

| Chapter | Topic | Folder | Status / Publication |
|---|---|---|---|
| 1 | Genetic parameters and genotype-by-diet (G × D) interactions | [Genetic parameters and genotype-by-diet (G × D) interactions](./Genetic%20parameters%20and%20genotype-by-diet%20%28G%20%C3%97%20D%29%20interactions) | *Journal of Insects as Food and Feed* (2026) |
| 2 | Genomic prediction framework | [Genomic prediction](./Genomic%20prediction) | bioRxiv preprint (2026) |
| 3 | Linkage map, QTL and GWAS | n/a | In progress |
| 4 | Long-term breeding evaluations | n/a | In progress |

Chapters 3 and 4 will be added as they are completed. Each chapter folder has its own README with the scripts, methods and run instructions.

---

## Chapter summaries

### Chapter 1: Genetic parameters and G × D interactions
Estimates phenotypic variation, genetic parameters, genetic relationships and population structure for four larval growth traits, and tests for genotype-by-diet interaction across three diets.

- **Data:** base population (2,083 individuals, three diets) genotyped with a custom 6K Allegro targeted SNP panel
- **Traits:** larval body weight (LBW), length (LL), width (LW) and surface area (LSA)
- **Diets:** soy-okara (SYK), brewers' spent grain (BSG), fruit and vegetable waste (FVW)
- **Methods:** descriptive statistics and fixed-effect analysis in SAS. Genotype QC and PCA in PLINK. Additive and additive-dominance genomic mixed models, bivariate models and multivariate G × D models in ASReml-R. GEBV extraction and reaction-norm plots.
- **Outputs:** variance components, heritabilities, genetic correlations between diets, and diet-specific GEBVs

> Gowda, K.B., Septriani, S., Jones, D.B., Jerry, D.R., Tedder, C., & Zenger, K.R. (2026). Genetic parameters and genotype-by-diet interactions for growth traits in Australian black soldier fly larvae: Implications for selective breeding. *Journal of Insects as Food and Feed*. https://doi.org/10.1163/23524588-bja10484

### Chapter 2: Genomic prediction
Develops a genomic prediction framework for BSF growth traits across diets, genetic effect models, marker panel densities and generations.

- **Data:** base population (2,083 individuals, three diets) and an F<sub>7</sub> population (536 individuals), genotyped with a custom 6K Allegro targeted SNP panel
- **Traits:** larval body weight (LBW), length (LL), width (LW) and surface area (LSA)
- **Methods:** GBLUP, BayesA, BayesB, BayesC, Bayesian Lasso, BayesR and random forest, with additive and additive-dominance models
- **Scenarios:** (1) overall base population, (2) within-diet, (3) multivariate G × D, (4) cross-diet, (5) SNP panel density, (6) forward and backward prediction between the base and F<sub>7</sub> populations
- **Validation:** 5-fold cross-validation with 20 replicates. Predictive ability, mean squared error of prediction and prediction scale bias.

> Gowda, K.B., Massault, C., Jerry, D.R. and Zenger, K.R. (2026). A genomic prediction framework for growth traits in black soldier fly (*Hermetia illucens*) breeding programmes. bioRxiv. https://doi.org/10.64898/2026.11.05.764768


---

## Repository structure

```
BSF_Quantitative_Genomics/
├── Genetic parameters and genotype-by-diet (G × D) interactions/   # Chapter 1
├── Genomic prediction/                                             # Chapter 2
└── README.md
```

## Software

- **R** 
- **ASReml-R**
- **SAS** 
- **PLINK** 

Package and software versions are listed in each chapter folder.

## Data availability

Raw genotype and phenotype datasets are not included in this repository.

## Citation

If you use these scripts, please cite the chapter publication listed above, or the thesis: Gowda, K.B. *Understanding the Genetic Architecture of Black Soldier Fly Production Traits for Enhanced Breeding Applications*. PhD thesis, James Cook University (in preparation).

## Contact

**Kishor B. Gowda**, PhD Candidate, James Cook University, Townsville, Australia · [email] · [GitHub](https://github.com/KB-Gowda)
