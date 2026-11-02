# ==============================================================================
# MULTIVARIATE SPIKE-SLAB CROSS-VALIDATION
# Models Evaluated: Model A (Additive-Only) vs Model AD (Additive + Dominance)
# Target Metrics: PA, Bias, MSEP
# ==============================================================================

# ------------------------------------------------------------------------------
# 1. USER SETTINGS & MULTIVARIATE TRAIT SETS
# ------------------------------------------------------------------------------

run_prefix <- "SpikeSlab_Multivariate_CV"

multivariate_sets <- list(
  Weight      = c("Adj_SYK_Weight", "Adj_BSG_Weight", "Adj_FVW_Weight"),
  SurfaceArea = c("Adj_SYK_SurfaceArea", "Adj_BSG_SurfaceArea", "Adj_FVW_SurfaceArea")
)

# BGLR MCMC Settings
n_iter  <- 12000
burn_in <- 2000
bayes_method <- "SpikeSlab"

# ------------------------------------------------------------------------------
# 2. LOAD DATA AND PACKAGES
# ------------------------------------------------------------------------------

library(dplyr)
library(Matrix)
library(BGLR)

load("BasePop_prep.RData")

set.seed(20260825)

stopifnot(all(rownames(cv_matrix) == as.character(pheno_final$Sample_Id)))

n_reps  <- ncol(cv_matrix)
n_folds <- max(cv_matrix)

all_results <- list()

# ------------------------------------------------------------------------------
# 3. CROSS-VALIDATION LOOP
# ------------------------------------------------------------------------------

for (set_name in names(multivariate_sets)) {
  current_traits <- multivariate_sets[[set_name]]
  
  cat("\n========================================================================\n")
  cat("STARTING MULTIVARIATE CROSS-VALIDATION FOR SET:", set_name, "\n")
  cat("Traits:", paste(current_traits, collapse = ", "), "\n")
  cat("========================================================================\n")
  
  # Extract the full multivariate response matrix
  Y_obs <- as.matrix(pheno_final[, current_traits, drop = FALSE])
  set_results <- list()
  
  for (r in seq_len(n_reps)) {
    for (f in seq_len(n_folds)) {
      test_idx <- which(cv_matrix[, r] == f)
      if (length(test_idx) == 0) next
      
      # Mask validation set with NA for all traits in the multivariate set
      Y_na <- Y_obs
      Y_na[test_idx, ] <- NA
      
      # ========================================================================
      # MODEL A: MULTIVARIATE ADDITIVE-ONLY
      # ========================================================================
      ETA_A <- list(
        Additive = list(X = Xa_final, model = bayes_method)
      )
      
      fit_A <- tryCatch({
        Multitrait(y = Y_na, ETA = ETA_A, nIter = n_iter, burnIn = burn_in, verbose = FALSE)
      }, error = function(e) NULL)
      
      if (!is.null(fit_A)) {
        # Multitrait stores marker effects in $beta (a matrix of markers x traits)
        gebv_all_A <- Xa_final %*% fit_A$ETA[[1]]$beta
        colnames(gebv_all_A) <- current_traits
        
        for (tr in current_traits) {
          p_A   <- gebv_all_A[test_idx, tr]
          o_val <- Y_obs[test_idx, tr]
          
          v_idx_A <- which(!is.na(o_val) & !is.na(p_A))
          if (length(v_idx_A) > 1) {
            p_sub  <- p_A[v_idx_A]
            o_sub  <- o_val[v_idx_A]
            pa_A   <- cor(p_sub, o_sub, use = "complete.obs")
            bias_A <- coef(lm(o_sub ~ p_sub))[2]
            msep_A <- mean((o_sub - p_sub)^2, na.rm = TRUE)
            
            set_results[[length(set_results) + 1]] <- data.frame(
              Model = "A", Set = set_name, Trait = tr, Rep = r, Fold = f,
              N_validation = length(v_idx_A),
              PA = pa_A, Bias = bias_A, MSEP = msep_A
            )
            
            cat(sprintf("Set: %-11s | Trait: %-19s | Rep %2d | Fold %d | [A]  PA: %.4f | Bias: %.4f | MSEP: %.4f\n",
                        set_name, tr, r, f, pa_A, bias_A, msep_A))
          }
        }
      }
      rm(fit_A)
      
      # ========================================================================
      # MODEL AD: MULTIVARIATE ADDITIVE + DOMINANCE
      # ========================================================================
      ETA_AD <- list(
        Additive  = list(X = Xa_final, model = bayes_method),
        Dominance = list(X = Xd_final, model = bayes_method)
      )
      
      fit_AD <- tryCatch({
        Multitrait(y = Y_na, ETA = ETA_AD, nIter = n_iter, burnIn = burn_in, verbose = FALSE)
      }, error = function(e) NULL)
      
      if (!is.null(fit_AD)) {
        gebv_all_AD <- Xa_final %*% fit_AD$ETA[[1]]$beta
        gdbv_all_AD <- Xd_final %*% fit_AD$ETA[[2]]$beta
        gpv_all_AD  <- gebv_all_AD + gdbv_all_AD
        colnames(gpv_all_AD) <- current_traits
        
        for (tr in current_traits) {
          p_AD  <- gpv_all_AD[test_idx, tr]
          o_val <- Y_obs[test_idx, tr]
          
          v_idx_AD <- which(!is.na(o_val) & !is.na(p_AD))
          if (length(v_idx_AD) > 1) {
            p_sub   <- p_AD[v_idx_AD]
            o_sub   <- o_val[v_idx_AD]
            pa_AD   <- cor(p_sub, o_sub, use = "complete.obs")
            bias_AD <- coef(lm(o_sub ~ p_sub))[2]
            msep_AD <- mean((o_sub - p_sub)^2, na.rm = TRUE)
            
            set_results[[length(set_results) + 1]] <- data.frame(
              Model = "AD", Set = set_name, Trait = tr, Rep = r, Fold = f,
              N_validation = length(v_idx_AD),
              PA = pa_AD, Bias = bias_AD, MSEP = msep_AD
            )
            
            cat(sprintf("Set: %-11s | Trait: %-19s | Rep %2d | Fold %d | [AD] PA: %.4f | Bias: %.4f | MSEP: %.4f\n",
                        set_name, tr, r, f, pa_AD, bias_AD, msep_AD))
          }
        }
      }
      rm(fit_AD)
      gc(verbose = FALSE)
      flush.console()
    }
  }
  
  if (length(set_results) > 0) {
    df_set <- bind_rows(set_results)
    all_results[[set_name]] <- df_set
    write.csv(df_set, paste0(run_prefix, "_", set_name, "_Raw.csv"), row.names = FALSE)
  }
}

# ------------------------------------------------------------------------------
# 4. COMBINE & SUMMARIZE RESULTS
# ------------------------------------------------------------------------------

if (length(all_results) == 0) {
  stop("No cross-validation results were generated.")
}

raw_results <- bind_rows(all_results)

summary_results <- raw_results %>%
  group_by(Model, Set, Trait) %>%
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
# 5. SAVE CSV OUTPUTS & MARKDOWN REPORT
# ------------------------------------------------------------------------------

write.csv(raw_results, paste0(run_prefix, "_Raw_Results.csv"), row.names = FALSE)
write.csv(summary_results, paste0(run_prefix, "_Summary_Results.csv"), row.names = FALSE)

# Generate Markdown Summary File (.md)
md_file <- paste0(run_prefix, "_Summary_Report.md")

md_content <- c(
  paste0("# ", gsub("_", " ", run_prefix), " Summary Report"),
  "",
  paste0("**Date Generated:** ", Sys.time()),
  paste0("**Replications:** ", n_reps),
  paste0("**Folds:** ", n_folds),
  paste0("**MCMC Iterations:** ", n_iter),
  paste0("**Burn-in:** ", burn_in),
  "",
  "## Model Performance Summary Table",
  "",
  "| Model | Set | Trait | N | Mean PA | SD PA | SE PA | Mean Bias | SD Bias | SE Bias | Mean MSEP | SD MSEP | SE MSEP |",
  "|:---|:---|:---|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---:|"
)

for (i in seq_len(nrow(summary_results))) {
  row <- summary_results[i, ]
  md_line <- sprintf(
    "| %s | %s | %s | %d | %.4f | %.4f | %.4f | %.4f | %.4f | %.4f | %.4f | %.4f | %.4f |",
    row$Model, row$Set, row$Trait, row$N_Evaluations,
    row$Mean_PA, row$SD_PA, row$SE_PA,
    row$Mean_Bias, row$SD_Bias, row$SE_Bias,
    row$Mean_MSEP, row$SD_MSEP, row$SE_MSEP
  )
  md_content <- c(md_content, md_line)
}

writeLines(md_content, con = md_file)
cat("\nMarkdown report saved to:", md_file, "\n")

cat("\n======================================================================\n")
cat("MULTIVARIATE SPIKE-SLAB CROSS-VALIDATION SUMMARY\n")
cat("======================================================================\n")
print(as.data.frame(summary_results))
