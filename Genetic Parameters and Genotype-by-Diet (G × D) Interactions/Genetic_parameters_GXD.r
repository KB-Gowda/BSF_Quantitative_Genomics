# ==============================================================================
# Title: Estimation of Genetic Parameters and Genotype-by-Diet (G × D) Interactions
#        for Growth Traits in Australian Black soldier fly (Hermetia illucens)
# Author: Kishor B. Gowda
# Date: 2026-08-21
# ==============================================================================

# ------------------------------------------------------------------------------
# 1. Load Required Packages
# ------------------------------------------------------------------------------
library(data.table)
library(dplyr)
library(readr)
library(ggplot2)
library(readxl)
library(snpReady)
library(MASS)
library(asreml)
library(AGHmatrix)
library(MCMCglmm)

# ------------------------------------------------------------------------------
# 2. Define Paths & Quality Control Filtering (PLINK)
# ------------------------------------------------------------------------------
plink_path <- "C:/Plink/plink.exe"  #(Change to your local path if different)

# Perform SNP and sample quality control using PLINK
system(paste(
  plink_path, "--bfile BSF_geno_2097_5562", "--allow-extra-chr", "--make-bed",
  "--geno 0.1", "--mind 0.10", "--maf 0.005",
  "--recode A", "--out BSF_geno_2083_4680"
))

# ------------------------------------------------------------------------------
# 3. Load Phenotype and Genotype Data
# ------------------------------------------------------------------------------
pheno <- read.table("Pheno_BSF.txt", header = TRUE, sep = "\t")
geno  <- fread("BSF_geno_2083_4680.raw")

# Extract genotype matrix (remove metadata columns)
geno_mat <- geno[, -(1:6), with = FALSE]

# ------------------------------------------------------------------------------
# 4. Match and Align IDs
# ------------------------------------------------------------------------------
common_ids <- intersect(pheno$Sample_Id, geno$IID)

pheno_match <- pheno[pheno$Sample_Id %in% common_ids, ]
geno_match  <- geno[IID %in% common_ids, ]

geno_match <- geno_match[match(pheno_match$Sample_Id, geno_match$IID), ]

M <- as.matrix(geno_match[, -(1:6), with = FALSE])

# ------------------------------------------------------------------------------
# 5. Compute and Invert Genomic Relationship Matrix (GRM)
# ------------------------------------------------------------------------------
geno_ready <- raw.data(
  data = M, frame = "wide", base = FALSE,
  call.rate = 0.9, maf = 0.005,
  imput = TRUE, imput.type = "mean"
)

M_impute <- geno_ready$M.clean

G  <- G.matrix(M = M_impute, method = "VanRaden", format = "wide")
Ga <- G$Ga
colnames(Ga) <- rownames(Ga) <- pheno_match$Sample_Id

Ginv <- ginv(Ga)

attr(Ginv, "rowNames") <- attr(Ginv, "colNames") <- as.character(pheno_match$Sample_Id)
dimnames(Ginv) <- list(as.character(pheno_match$Sample_Id),
                       as.character(pheno_match$Sample_Id))
attr(Ginv, "INVERSE") <- TRUE

# ------------------------------------------------------------------------------
# 6. Construct Dominance Relationship Matrix
# ------------------------------------------------------------------------------
Gd <- Gmatrix(SNPmatrix = as.matrix(M), method = "Vitezica")
colnames(Gd) <- rownames(Gd) <- pheno_match$Sample_Id
iGd <- solve(Gd)
iGd <- as(iGd, "sparseMatrix") 
Gd.inv <- sm2asreml(iGd) 
attr(Gd.inv, "INVERSE") <- TRUE

# ------------------------------------------------------------------------------
# 7. Format Data Types & Define Trait Categories
# ------------------------------------------------------------------------------
pheno_match$Sample_Id    <- as.factor(pheno_match$Sample_Id)
pheno_match$Tray    <- as.factor(pheno_match$Tray)
pheno_match$Sampling_Day <- as.factor(pheno_match$Sampling_Day)
pheno_match$Diet <- as.factor(pheno_match$Diet)


diets          <- c("SYK", "BSG", "FVW")

# Primary (Main) Growth Traits
primary_traits <- c("Weight", "Length", "Width", "SurfaceArea")

# Diet-Specific Growth Traits
diet_traits    <- c(outer(diets, primary_traits, paste, sep = "_"))

all_traits <- c(primary_traits, diet_traits)

# Format all trait variables as numeric
for (tr in all_traits) {
  if (tr %in% colnames(pheno_match)) {
    pheno_match[[tr]] <- as.numeric(pheno_match[[tr]])
  }
}

# ------------------------------------------------------------------------------
# 7A. Principal Component Analysis (PCA) of QC-Filtered Genotypes
# ------------------------------------------------------------------------------
# Run PCA in PLINK
system(paste(plink_path, "--bfile BSF_geno_2083_4680", "--pca 10", "--allow-extra-chr", "--out BSF_pca"))

# Read PCA eigenvalues and calculate variance explained
eigenvals <- read.table("BSF_pca.eigenval", header = FALSE)$V1
var_exp   <- (eigenvals / sum(eigenvals)) * 100
pc1_lab   <- paste0("PC1 (", round(var_exp[1], 2), "%)")
pc2_lab   <- paste0("PC2 (", round(var_exp[2], 2), "%)")

# Read eigenvectors and merge directly with phenotype data
pca_data <- read.table("BSF_pca.eigenvec", header = FALSE, stringsAsFactors = FALSE)
colnames(pca_data) <- c("Fam_ID", "Sample_Id", paste0("PC", 1:(ncol(pca_data) - 2)))

final_diet_pca <- inner_join(pca_data, pheno_match[, c("Sample_Id", "Diet")], by = "Sample_Id") %>%
  filter(!is.na(Diet) & Diet != "")

cat("PCA Individuals:", nrow(pca_data), "| Phenotype Matched:", nrow(final_diet_pca), "\n")

# ------------------------------------------------------------------------------
# PCA Plot by Diet
# ------------------------------------------------------------------------------
plot_diet <- ggplot(final_diet_pca, aes(x = PC1, y = PC2, color = Diet)) +
  geom_point(alpha = 0.6, size = 1.0) +
  scale_color_manual(values = c("BSG" = "#FF7F00", "SYK" = "#D81B60", "FVW" = "#00A88F")) +
  theme_bw() +
  labs(x = pc1_lab, y = pc2_lab, color = "Diet group") +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    legend.position = "right",
    panel.grid.minor = element_blank()
  )

print(plot_diet)
ggsave("BSF_pca_diet_2083.svg", plot = plot_diet, width = 10, height = 5.5, device = "svg")

# ==============================================================================
# PART 1: MAIN GROWTH TRAITS (Weight, Length, Width, Surface Area)
# ==============================================================================

# ------------------------------------------------------------------------------
# 8. Main Growth Traits: Univariate Additive Models (Heritability Estimation)
# ------------------------------------------------------------------------------
models_add_main <- list()

for (tr in primary_traits) {
  cat("\n==========================================\n")
  cat("Univariate Additive Model (Heritability) for Main Trait:", tr, "\n")
  cat("==========================================\n")
  
  fml <- as.formula(paste(tr, "~ Gray_Scale + Sampling_Day + Diet + Tray:Diet"))
  
  mod <- asreml(
    fixed     = fml,
    random    = ~ vm(Sample_Id, Ginv),
    residual  = ~ idv(units),
    data      = pheno_match,
    maxit     = 100,
    workspace = "1GB"
  )
  
  mod <- update.asreml(mod)
  
  summary(mod)
  print(summary(mod)$varcomp)
  print(vpredict(mod, as.formula(paste0("VP_", tr, " ~ V1 + V2"))))
  print(vpredict(mod, as.formula(paste0("h2_", tr, " ~ V1 / (V1 + V2)"))))
  print(vpredict(mod, as.formula(paste0("e2_", tr, " ~ V2 / (V1 + V2)"))))
  wald.asreml(mod)
  plot(mod)
  
  models_add_main[[tr]] <- mod
}

# ------------------------------------------------------------------------------
# 9. Main Growth Traits: Univariate Additive-Dominance Models (Heritability Estimation)
# ------------------------------------------------------------------------------
models_ad_main <- list()

for (tr in primary_traits) {
  cat("\n==========================================\n")
  cat("Univariate Additive-Dominance Model for Main Trait:", tr, "\n")
  cat("==========================================\n")
  
  fml <- as.formula(paste(tr, "~ Gray_Scale + Sampling_Day + Diet + Tray:Diet"))
  
  mod <- asreml(
    fixed     = fml,
    random    = ~ vm(Sample_Id, Ginv) + vm(Sample_Id, Gd.inv),
    residual  = ~ idv(units),
    data      = pheno_match,
    maxit     = 100,
    workspace = "1GB"
  )
  
  mod <- update.asreml(mod)
  
  print(summary(mod)$varcomp)
  print(vpredict(mod, as.formula(paste0("h2_", tr, " ~ V1 / (V1 + V2 + V3)"))))
  print(vpredict(mod, as.formula(paste0("d2_", tr, " ~ V2 / (V1 + V2 + V3)"))))
  print(vpredict(mod, as.formula(paste0("e2_", tr, " ~ V3 / (V1 + V2 + V3)"))))
  print(vpredict(mod, as.formula(paste0("P_", tr, " ~ V1 + V2 + V3"))))
  wald.asreml(mod)
  plot(mod)
  
  models_ad_main[[tr]] <- mod
}

# ------------------------------------------------------------------------------
# 10. Main Growth Traits: Bivariate Models (Genetic & Phenotypic Correlations)
# ------------------------------------------------------------------------------
main_trait_pairs <- combn(primary_traits, 2, simplify = FALSE)
models_bivariate_main <- list()

for (pair in main_trait_pairs) {
  t1 <- pair[1]
  t2 <- pair[2]
  pair_name <- paste(t1, t2, sep = "_vs_")
  
  cat("\n==========================================\n")
  cat("Bivariate Model (Genetic & Phenotypic Correlations) for Main Traits:", t1, "&", t2, "\n")
  cat("==========================================\n")
  
  fml_fixed <- as.formula(
    paste0("cbind(", t1, ", ", t2, ")  ~ trait + trait:Gray_Scale + trait:Diet + trait:Sampling_Day + trait:Tray:Diet")
  )
  
  mod <- asreml(
    fixed     = fml_fixed,
    random    = ~ us(trait, init=c(100,0.1,100)):vm(Sample_Id, Ginv) + us(trait, init=c(100,0.1,100)):vm(Sample_Id, Gd.inv),
    coruh     = ~ id(units):us(trait,init=c(100,0.1,100)),
    data      = pheno_match,
    na.action = na.method(x = "include", y = "include"),
    maxit     = 100,
    ai.sing   = TRUE,
    workspace = "2GB"
  )
  
  mod <- update.asreml(mod)
  
  print(summary(mod)$varcomp)
  print(vpredict(mod, rg ~ V2 / sqrt(V1 * V3)))
  print(vpredict(mod, rD ~ V5 / sqrt(V4 * V6)))
  print(vpredict(mod, re ~ V9 / sqrt(V8 * V10)))
  print(vpredict(mod, rp ~ (V2 + V5 + V9) / sqrt((V1 + V4 + V8) * (V3 + V6 + V10))))
  print(vpredict(mod, h2.1 ~ V1 / (V1 + V4 + V8)))
  print(vpredict(mod, h2.2 ~ V3 / (V3 + V6 + V10)))
  print(vpredict(mod, d2.1 ~ V4 / (V1 + V4 + V8)))
  print(vpredict(mod, d2.2 ~ V6 / (V3 + V6 + V10)))
  wald.asreml(mod)
  plot(mod)
  
  models_bivariate_main[[pair_name]] <- mod
}

# ==============================================================================
# PART 2: DIET-SPECIFIC GROWTH TRAITS (Weight, Length, Width, Surface Area per Diet)
# ==============================================================================

# ------------------------------------------------------------------------------
# 11. Diet-Specific Traits: Univariate Additive Models (Heritability Estimation)
# ------------------------------------------------------------------------------
models_add_diet <- list()

for (tr in diet_traits) {
  cat("\n==========================================\n")
  cat("Univariate Additive Model (Heritability) for Diet-Specific Trait:", tr, "\n")
  cat("==========================================\n")
  
  fml <- as.formula(paste(tr, "~ Gray_Scale + Sampling_Day + Tray:Sampling_Day"))
  
  mod <- asreml(
    fixed     = fml,
    random    = ~ vm(Sample_Id, Ginv),
    residual  = ~ idv(units),
    data      = pheno_match,
    maxit     = 100,
    workspace = "1GB"
  )
  
  mod <- update.asreml(mod)
  
  summary(mod)
  print(summary(mod)$varcomp)
  print(vpredict(mod, as.formula(paste0("VP_", tr, " ~ V1 + V2"))))
  print(vpredict(mod, as.formula(paste0("h2_", tr, " ~ V1 / (V1 + V2)"))))
  print(vpredict(mod, as.formula(paste0("e2_", tr, " ~ V2 / (V1 + V2)"))))
  wald.asreml(mod)
  plot(mod)
  
  models_add_diet[[tr]] <- mod
}

# ------------------------------------------------------------------------------
# 12. Diet-Specific Traits: Univariate Additive-Dominance Models (Heritability Estimation)
# ------------------------------------------------------------------------------
models_ad_diet <- list()

for (tr in diet_traits) {
  cat("\n==========================================\n")
  cat("Univariate Additive-Dominance Model for Diet-Specific Trait:", tr, "\n")
  cat("==========================================\n")
  
  fml <- as.formula(paste(tr, "~ Gray_Scale + Sampling_Day + Tray:Sampling_Day"))
  
  mod <- asreml(
    fixed     = fml,
    random    = ~ vm(Sample_Id, Ginv) + vm(Sample_Id, Gd.inv),
    residual  = ~ idv(units),
    data      = pheno_match,
    maxit     = 100,
    workspace = "1GB"
  )
  
  mod <- update.asreml(mod)
  
  print(summary(mod)$varcomp)
  print(vpredict(mod, as.formula(paste0("h2_", tr, " ~ V1 / (V1 + V2 + V3)"))))
  print(vpredict(mod, as.formula(paste0("d2_", tr, " ~ V2 / (V1 + V2 + V3)"))))
  print(vpredict(mod, as.formula(paste0("e2_", tr, " ~ V3 / (V1 + V2 + V3)"))))
  print(vpredict(mod, as.formula(paste0("P_", tr, " ~ V1 + V2 + V3"))))
  wald.asreml(mod)
  plot(mod)
  
  models_ad_diet[[tr]] <- mod
}

# ==============================================================================
# PART 3: MULTI-VARIATE GxD MODELS (Genotype-by-Diet Interaction)
# ==============================================================================

# ------------------------------------------------------------------------------
# 13. Multivariate GxD Models for Primary Growth Traits across Diets
# ------------------------------------------------------------------------------
models_gxd_main <- list()

for (tr in primary_traits) {
  cat("\n==========================================\n")
  cat("Running GxD Model for:", tr, "\n")
  cat("==========================================\n")
  
  syk_tr <- paste("SYK", tr, sep = "_")
  bsg_tr <- paste("BSG", tr, sep = "_")
  fvw_tr <- paste("FVW", tr, sep = "_")
  
  fml_fixed <- as.formula(
    paste0("cbind(", syk_tr, ", ", bsg_tr, ", ", fvw_tr, ") ~ trait + trait:Gray_Scale + trait:Sampling_Day + trait:Tray:Sampling_Day")
  )
  
  mod <- asreml(
    fixed     = fml_fixed,
    random    = ~ us(trait, init = c(100, 0.1, 100, 0.1, 0.1, 100)):vm(Sample_Id, Ginv),
    coruh     = ~ id(units):us(trait, init = c(100, 0.1, 100, 0.1, 0.1, 100)),
    data      = pheno_match,
    na.action = na.method(x = "include", y = "include"),
    maxit     = 100, 
    ai.sing   = TRUE, 
    workspace = "2GB"
  )
  
  mod <- update.asreml(mod)
  
  summary(mod)
  print(summary(mod)$varcomp)
  
  # Genetic Correlations
  print(vpredict(mod, rg_SYK_BSG ~ V2 / sqrt(V1 * V3)))
  print(vpredict(mod, rg_SYK_FVW ~ V4 / sqrt(V1 * V6)))
  print(vpredict(mod, rg_BSG_FVW ~ V5 / sqrt(V3 * V6)))
  
  # Residual Correlations
  print(vpredict(mod, re_SYK_BSG ~ V9 / sqrt(V8 * V10)))
  print(vpredict(mod, re_SYK_FVW ~ V11 / sqrt(V8 * V13)))
  print(vpredict(mod, re_BSG_FVW ~ V12 / sqrt(V10 * V13)))
  
  # Phenotypic Correlations
  print(vpredict(mod, rp_SYK_BSG ~ (V2 + V9) / sqrt((V1 + V8) * (V3 + V10))))
  print(vpredict(mod, rp_SYK_FVW ~ (V4 + V11) / sqrt((V1 + V8) * (V6 + V13))))
  print(vpredict(mod, rp_BSG_FVW ~ (V5 + V12) / sqrt((V3 + V10) * (V6 + V13))))
  
  # Diet-Specific Heritabilities
  print(vpredict(mod, h2_SYK ~ V1 / (V1 + V8)))
  print(vpredict(mod, h2_BSG ~ V3 / (V3 + V10)))
  print(vpredict(mod, h2_FVW ~ V6 / (V6 + V13)))
  
  wald.asreml(mod)
  plot(mod)
  
  models_gxd_main[[tr]] <- mod

  # ------------------------------------------------------------------------------
  # Extract GEBVs and Export to CSV
  # ------------------------------------------------------------------------------
  coef_mod <- coef(mod)$random
  gebv_rows <- grep("vm\\(Sample_Id, Ginv\\)", rownames(coef_mod))
  gebv_df <- data.frame(
    Term = rownames(coef_mod)[gebv_rows],
    solution = coef_mod[gebv_rows, 1],
    stringsAsFactors = FALSE
  )
  
  gebv_df$Trait <- sub(":vm\\(Sample_Id, Ginv\\).*", "", gebv_df$Term)
  gebv_df$Sample_Id <- sub(".*:vm\\(Sample_Id, Ginv\\)_", "", gebv_df$Term)
  
  write.csv(gebv_df, file = paste0("GEBV_", tr, ".csv"), row.names = FALSE)
}

# ==============================================================================
# PART 4: REACTION NORM PLOTTING (TOP 5% GEBVs PER TRAIT)
# ==============================================================================

# ------------------------------------------------------------------------------
# 14. Reaction norm plotting with GEBVs
# ------------------------------------------------------------------------------

# Read exported CSV files and label traits
GEBV_weight  <- read_csv("GEBV_Weight.csv")      %>% mutate(TraitName = "LBW")
GEBV_length  <- read_csv("GEBV_Length.csv")      %>% mutate(TraitName = "LL")
GEBV_width   <- read_csv("GEBV_Width.csv")       %>% mutate(TraitName = "LW")
GEBV_surface <- read_csv("GEBV_SurfaceArea.csv") %>% mutate(TraitName = "LSA")

# Combine datasets
GEBV_all <- bind_rows(GEBV_weight, GEBV_length, GEBV_width, GEBV_surface)

# Clean and format diet labels
GEBV_clean <- GEBV_all %>%
  mutate(
    Sample_Id = sub("^_", "", Sample_Id),
    diet = case_when(
      grepl("BSG", Trait) ~ "BSG",
      grepl("SYK", Trait) ~ "SYK",
      grepl("FVW|VEG", Trait) ~ "FVW",
      TRUE ~ NA_character_
    ),
    diet = factor(diet, levels = c("BSG", "FVW", "SYK")),
    TraitName = factor(TraitName, levels = c("LBW", "LL", "LW", "LSA"))
  ) %>%
  filter(!is.na(diet))

# Zero-center solutions within each Trait x Diet environment
GEBV_centered <- GEBV_clean %>%
  group_by(TraitName, diet) %>%
  mutate(solution_centered = solution - mean(solution, na.rm = TRUE)) %>%
  ungroup()

# Rank and isolate Top 5% individuals independently PER TRAIT
GEBV_top5 <- GEBV_centered %>%
  group_by(TraitName, Sample_Id) %>%
  mutate(mean_centered_GEBV = mean(solution_centered, na.rm = TRUE)) %>%
  group_by(TraitName) %>%
  filter(dense_rank(-mean_centered_GEBV) <= ceiling(0.05 * n_distinct(Sample_Id))) %>%
  ungroup()

# Plot reaction norms for Top 5% per trait
p_GxD_panel <- ggplot(GEBV_top5, 
                      aes(x = diet, y = solution_centered, group = Sample_Id, colour = Sample_Id)) +
  geom_line(alpha = 0.6) +
  geom_point(size = 1) +
  facet_wrap(~TraitName, scales = "free_y", ncol = 2) +
  theme_bw(base_size = 14) +
  labs(
    x = "Diet", 
    y = "Centered GEBVs (Top 5%)"
  ) +
  theme(
    legend.position = "none",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.text = element_text(face = "bold")
  )

# Display and export figure
print(p_GxD_panel)
ggsave("BSF_GxD_reaction_norms_top5.svg", p_GxD_panel, width = 12, height = 8)
ggsave("BSF_GxD_reaction_norms_top5.jpeg", p_GxD_panel, width = 12, height = 8, dpi = 600)
