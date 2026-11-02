# ==============================================================================
# BAYES R CROSS-VALIDATION FOR DIET-SPECIFIC TRAITS (HIBAYES)
# Models Evaluated: Model A (Additive-Only) vs Model AD (Additive + Dominance)
# ==============================================================================

# ------------------------------------------------------------------------------
# 1. USER SETTINGS, TARGET TRAITS & RUN IDENTIFIER
# ------------------------------------------------------------------------------

run_prefix   <- "BayesR_CV_DietSpecific"
target_traits <- c(
  "Adj_SYK_Weight", "Adj_SYK_SurfaceArea",
  "Adj_BSG_Weight", "Adj_BSG_SurfaceArea",
  "Adj_FVW_Weight", "Adj_FVW_SurfaceArea"
)

# MCMC & Compute Settings
n_iter    <- 12000
burn_in   <- 2000
n_threads <- 4     # Set the number of threads for multi-threading

# ------------------------------------------------------------------------------
# 2. LOAD DATA AND PACKAGES
# ------------------------------------------------------------------------------

library(dplyr)
library(Matrix)
library(hibayes)

load("BasePop_prep.RData")

# Convert to character and move to the first column (required by hibayes)
pheno_final <- pheno_final %>%
  mutate(Sample_Id = as.character(Sample_Id)) %>%
  relocate(Sample_Id, .before = 1)

set.seed(20260825)

# Ensure the cv_matrix row names match the newly formatted Sample_Id
stopifnot(all(rownames(cv_matrix) == pheno_final$Sample_Id))

n_reps  <- ncol(cv_matrix)
n_folds <- max(cv_matrix)

all_results <- list()

# ------------------------------------------------------------------------------
# 3. PREPARE DOMINANCE MATRIX COLUMNS
# ------------------------------------------------------------------------------
# Ensure unique column names before cbinding to avoid hibayes crash
Xd_unique <- Xd_final
colnames(Xd_unique) <- paste0(colnames(Xd_unique), "_dom")
M_AD_combined <- cbind(Xa_final, Xd_unique)

# ------------------------------------------------------------------------------
# 4. CROSS-VALIDATION LOOP
# ------------------------------------------------------------------------------

for (trait in target_traits) {
  cat("\n========================================================================\n")
  cat("STARTING CROSS-VALIDATION FOR TRAIT:", trait, "\n")
  cat("========================================================================\n")
  
  y_obs <- pheno_final[[trait]]
  trait_results <- list()
  
  for (r in seq_len(n_reps)) {
    for (f in seq_len(n_folds)) {
      test_idx <- which(cv_matrix[, r] == f)
      
      # Mask validation set
      y_na <- y_obs
      y_na[test_idx] <- NA
      
      pheno_cv <- pheno_final
      pheno_cv$masked_pheno <- y_na
      
      # ========================================================================
      # MODEL A: BAYES R ADDITIVE-ONLY
      # ========================================================================
      fit_A <- tryCatch({
        set.seed(20260825 + (r * 100) + f)
        ibrm(
          formula = masked_pheno ~ 1,
          data    = pheno_cv,
          M       = Xa_final,
          M.id    = as.character(pheno_cv$Sample_Id),
          method  = "BayesR",
          Pi      = c(0.95, 0.02, 0.02, 0.01),
          fold    = c(0, 0.0001, 0.001, 0.01), 
          niter   = n_iter,
          nburn   = burn_in,
          threads = n_threads
        )
      }, error = function(e) {
        cat("\n[Error in Model A]:", e$message, "\n")
        return(NULL)
      })
      
      pa_A <- NA; bias_A <- NA; msep_A <- NA
      
      if (!is.null(fit_A)) {
        p_A   <- fit_A$g$gebv[test_idx]
        o_val <- y_obs[test_idx]
        
        v_idx_A <- which(!is.na(o_val) & !is.na(p_A))
        if (length(v_idx_A) > 1) {
          p_sub  <- p_A[v_idx_A]
          o_sub  <- o_val[v_idx_A]
          pa_A   <- cor(p_sub, o_sub, use = "complete.obs")
          bias_A <- coef(lm(o_sub ~ p_sub))[2]
          msep_A <- mean((o_sub - p_sub)^2, na.rm = TRUE)
          
          trait_results[[length(trait_results) + 1]] <- data.frame(
            Model = "A", Trait = trait, Rep = r, Fold = f,
            N_validation = length(v_idx_A),
            PA = pa_A, Bias = bias_A, MSEP = msep_A
          )
        }
      }
      rm(fit_A)
      
      # ========================================================================
      # MODEL AD: BAYES R ADDITIVE + DOMINANCE
      # ========================================================================
      fit_AD <- tryCatch({
        set.seed(20260825 + (r * 100) + f)
        ibrm(
          formula = masked_pheno ~ 1,
          data    = pheno_cv,
          M       = M_AD_combined,
          M.id    = as.character(pheno_cv$Sample_Id),
          method  = "BayesR",
          Pi      = c(0.95, 0.02, 0.02, 0.01),
          fold    = c(0, 0.0001, 0.001, 0.01), 
          niter   = n_iter,
          nburn   = burn_in,
          threads = n_threads
        )
      }, error = function(e) {
        cat("\n[Error in Model AD]:", e$message, "\n")
        return(NULL)
      })
      
      pa_AD <- NA; bias_AD <- NA; msep_AD <- NA
      
      if (!is.null(fit_AD)) {
        p_AD  <- fit_AD$g$gebv[test_idx]
        o_val <- y_obs[test_idx]
        
        v_idx_AD <- which(!is.na(o_val) & !is.na(p_AD))
        if (length(v_idx_AD) > 1) {
          p_sub   <- p_AD[v_idx_AD]
          o_sub   <- o_val[v_idx_AD]
          pa_AD   <- cor(p_sub, o_sub, use = "complete.obs")
          bias_AD <- coef(lm(o_sub ~ p_sub))[2]
          msep_AD <- mean((o_sub - p_sub)^2, na.rm = TRUE)
          
          trait_results[[length(trait_results) + 1]] <- data.frame(
            Model = "AD", Trait = trait, Rep = r, Fold = f,
            N_validation = length(v_idx_AD),
            PA = pa_AD, Bias = bias_AD, MSEP = msep_AD
          )
        }
      }
      rm(fit_AD)
      gc(verbose = FALSE)
      
      # Print results after each fold
      cat(sprintf("Trait: %-20s | Rep %2d | Fold %d | [A]  PA: %.4f | Bias: %.4f | MSEP: %.4f\n",
                  trait, r, f, pa_A, bias_A, msep_A))
      cat(sprintf("Trait: %-20s | Rep %2d | Fold %d | [AD] PA: %.4f | Bias: %.4f | MSEP: %.4f\n",
                  trait, r, f, pa_AD, bias_AD, msep_AD))
      flush.console()
    }
  }
  
  if (length(trait_results) > 0) {
    df_trait <- bind_rows(trait_results)
    all_results[[trait]] <- df_trait
    write.csv(df_trait, paste0(run_prefix, "_", trait, "_Raw.csv"), row.names = FALSE)
  }
}

# ------------------------------------------------------------------------------
# 5. COMBINE & SUMMARIZE RESULTS
# ------------------------------------------------------------------------------

if (length(all_results) == 0) {
  stop("No cross-validation results were generated. Check the error prints above.")
}

raw_results <- bind_rows(all_results)

summary_results <- raw_results %>%
  group_by(Model, Trait) %>%
  summarise(
    N_Evaluations = n(),
    Mean_PA   = mean(PA, na.rm = TRUE),
    SD_PA     = sd(PA, na.rm = TRUE),
    SE_PA     = SD_PA / sqrt(N_Evaluations),
    Mean_Bias = mean(Bias, na.rm = TRUE),
    SD_Bias   = sd(Bias, na.rm = TRUE),
    SE_Bias   = SD_Bias / sqrt(N_Evaluations),
    Mean_MSEP = mean(MSEP, na.rm = TRUE),
    SD_MSEP   = sd(MSEP, na.rm = TRUE),
    SE_MSEP   = SD_MSEP / sqrt(N_Evaluations),
    .groups   = "drop"
  )

# ------------------------------------------------------------------------------
# 6. SAVE CSV OUTPUTS & MARKDOWN REPORT
# ------------------------------------------------------------------------------

write.csv(raw_results, paste0(run_prefix, "_Raw_Results.csv"), row.names = FALSE)
write.csv(summary_results, paste0(run_prefix, "_Summary_Results.csv"), row.names = FALSE)

md_file <- paste0(run_prefix, "_Summary_Report.md")
md_content <- c(
  paste0("# ", gsub("_", " ", run_prefix), " Summary Report"),
  "",
  paste0("**Date Generated:** ", Sys.time()),
  paste0("**Replications:** ", n_reps),
  paste0("**Folds:** ", n_folds),
  paste0("**MCMC Iterations:** ", n_iter),
  paste0("**Burn-in:** ", burn_in),
  paste0("**Threads used:** ", n_threads),
  "",
  "## Model Performance Summary Table",
  "",
  "| Model | Trait | N | Mean PA | SD PA | SE PA | Mean Bias | SD Bias | SE Bias | Mean MSEP | SD MSEP | SE MSEP |",
  "|:---|:---|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---:|"
)

for (i in seq_len(nrow(summary_results))) {
  row <- summary_results[i, ]
  md_line <- sprintf(
    "| %s | %s | %d | %.4f | %.4f | %.4f | %.4f | %.4f | %.4f | %.4f | %.4f | %.4f |",
    row$Model, row$Trait, row$N_Evaluations,
    row$Mean_PA, row$SD_PA, row$SE_PA,
    row$Mean_Bias, row$SD_Bias, row$SE_Bias,
    row$Mean_MSEP, row$SD_MSEP, row$SE_MSEP
  )
  md_content <- c(md_content, md_line)
}

writeLines(md_content, con = md_file)
cat("\nMarkdown report saved to:", md_file, "\n")

cat("\n======================================================================\n")
cat("BAYES R CROSS-VALIDATION SUMMARY\n")
cat("======================================================================\n")
print(summary_results)
