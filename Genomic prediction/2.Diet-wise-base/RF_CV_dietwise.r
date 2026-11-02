# ==============================================================================
# RANDOM FOREST NESTED CROSS-VALIDATION VIA RANGER FOR DIET-SPECIFIC TRAITS
# Models Evaluated: Model A (Additive-Only) vs Model AD (Additive + Dominance)
# Hyperparameter Optimization: Evaluated via Out-Of-Bag (OOB) MSE
# Parallelism: Multi-threading handled directly within ranger (num.threads)
# ==============================================================================

# ------------------------------------------------------------------------------
# 1. USER SETTINGS & COMPUTATIONAL PARALLELISM
# ------------------------------------------------------------------------------

run_prefix   <- "RF_NestedCV_DietSpecific"
target_traits <- c(
  "Adj_SYK_Weight", "Adj_SYK_SurfaceArea",
  "Adj_BSG_Weight", "Adj_BSG_SurfaceArea",
  "Adj_FVW_Weight", "Adj_FVW_SurfaceArea"
)
n_threads     <- 6  # Number of threads assigned directly to ranger

# ------------------------------------------------------------------------------
# 2. LOAD DATA AND PACKAGES
# ------------------------------------------------------------------------------

library(dplyr)
library(Matrix)
library(ranger)

load("BasePop_prep.RData")

# Ensure Sample_Id is character and in column 1
pheno_final <- pheno_final %>%
  mutate(Sample_Id = as.character(Sample_Id)) %>%
  relocate(Sample_Id, .before = 1)

stopifnot(all(rownames(cv_matrix) == pheno_final$Sample_Id))

n_reps  <- ncol(cv_matrix)
n_folds <- max(cv_matrix)

all_results <- list()

# ------------------------------------------------------------------------------
# 3. PREPARE GENOTYPE MATRICES & HYPERPARAMETER GRIDS
# ------------------------------------------------------------------------------

# Additive + Dominance matrix construction
Xd_unique <- Xd_final
colnames(Xd_unique) <- paste0(colnames(Xd_unique), "_dom")
M_AD_combined <- cbind(Xa_final, Xd_unique)

# Function to generate hyperparameter grid based on predictor matrix dimensions
get_grid <- function(p_vars) {
  mtry_vals <- unique(c(
    floor(sqrt(p_vars)),
    floor(p_vars / 10),
    floor(p_vars / 5),
    floor(p_vars / 3)
  ))
  mtry_vals <- mtry_vals[mtry_vals >= 1]
  
  expand.grid(
    num.trees     = c(150, 300),
    mtry          = mtry_vals,
    min.node.size = c(5, 10),
    stringsAsFactors = FALSE
  )
}

grid_A  <- get_grid(ncol(Xa_final))
grid_AD <- get_grid(ncol(M_AD_combined))

# ------------------------------------------------------------------------------
# 4. NESTED CROSS-VALIDATION LOOP
# ------------------------------------------------------------------------------

for (trait in target_traits) {
  cat("\n========================================================================\n")
  cat("STARTING RANDOM FOREST CV FOR TRAIT:", trait, "\n")
  cat("========================================================================\n")
  
  y_obs <- pheno_final[[trait]]
  trait_results <- list()
  
  for (r in seq_len(n_reps)) {
    for (f in seq_len(n_folds)) {
      test_idx  <- which(cv_matrix[, r] == f)
      train_idx <- setdiff(seq_len(nrow(pheno_final)), test_idx)
      
      # Filter non-missing phenotype indices
      train_valid_idx <- train_idx[!is.na(y_obs[train_idx])]
      test_valid_idx  <- test_idx[!is.na(y_obs[test_idx])]
      
      if (length(train_valid_idx) == 0 || length(test_valid_idx) == 0) next
      
      y_train <- y_obs[train_valid_idx]
      y_test  <- y_obs[test_valid_idx]
      
      # ========================================================================
      # MODEL A: ADDITIVE-ONLY
      # ========================================================================
      X_train_A <- Xa_final[train_valid_idx, , drop = FALSE]
      X_test_A  <- Xa_final[test_valid_idx, , drop = FALSE]
      
      fit_A_res <- tryCatch({
        oob_errors <- numeric(nrow(grid_A))
        
        # Grid search using ranger multi-threading
        for (g in seq_len(nrow(grid_A))) {
          set.seed(20260825 + (r * 1000) + (f * 100) + g)
          rf_tune <- ranger(
            x             = X_train_A,
            y             = y_train,
            num.trees     = grid_A$num.trees[g],
            mtry          = grid_A$mtry[g],
            min.node.size = grid_A$min.node.size[g],
            num.threads   = n_threads
          )
          oob_errors[g] <- rf_tune$prediction.error
        }
        
        best_idx_A   <- which.min(oob_errors)
        best_param_A <- grid_A[best_idx_A, ]
        
        # Refit optimal model
        set.seed(20260825 + (r * 100) + f)
        best_fit_A <- ranger(
          x             = X_train_A,
          y             = y_train,
          num.trees     = best_param_A$num.trees,
          mtry          = best_param_A$mtry,
          min.node.size = best_param_A$min.node.size,
          num.threads   = n_threads
        )
        
        preds_A <- predict(best_fit_A, data = X_test_A)$predictions
        
        pa_A   <- cor(preds_A, y_test, use = "complete.obs")
        bias_A <- coef(lm(y_test ~ preds_A))[2]
        msep_A <- mean((y_test - preds_A)^2, na.rm = TRUE)
        
        list(PA = pa_A, Bias = bias_A, MSEP = msep_A, Best_Param = best_param_A)
      }, error = function(e) {
        cat("\n[Error in Model A]:", e$message, "\n")
        return(NULL)
      })
      
      # ========================================================================
      # MODEL AD: ADDITIVE + DOMINANCE
      # ========================================================================
      X_train_AD <- M_AD_combined[train_valid_idx, , drop = FALSE]
      X_test_AD  <- M_AD_combined[test_valid_idx, , drop = FALSE]
      
      fit_AD_res <- tryCatch({
        oob_errors <- numeric(nrow(grid_AD))
        
        # Grid search using ranger multi-threading
        for (g in seq_len(nrow(grid_AD))) {
          set.seed(20260825 + (r * 1000) + (f * 100) + g)
          rf_tune <- ranger(
            x             = X_train_AD,
            y             = y_train,
            num.trees     = grid_AD$num.trees[g],
            mtry          = grid_AD$mtry[g],
            min.node.size = grid_AD$min.node.size[g],
            num.threads   = n_threads
          )
          oob_errors[g] <- rf_tune$prediction.error
        }
        
        best_idx_AD   <- which.min(oob_errors)
        best_param_AD <- grid_AD[best_idx_AD, ]
        
        # Refit optimal model
        set.seed(20260825 + (r * 100) + f)
        best_fit_AD <- ranger(
          x             = X_train_AD,
          y             = y_train,
          num.trees     = best_param_AD$num.trees,
          mtry          = best_param_AD$mtry,
          min.node.size = best_param_AD$min.node.size,
          num.threads   = n_threads
        )
        
        preds_AD <- predict(best_fit_AD, data = X_test_AD)$predictions
        
        pa_AD   <- cor(preds_AD, y_test, use = "complete.obs")
        bias_AD <- coef(lm(y_test ~ preds_AD))[2]
        msep_AD <- mean((y_test - preds_AD)^2, na.rm = TRUE)
        
        list(PA = pa_AD, Bias = bias_AD, MSEP = msep_AD, Best_Param = best_param_AD)
      }, error = function(e) {
        cat("\n[Error in Model AD]:", e$message, "\n")
        return(NULL)
      })
      
      pa_A   <- ifelse(!is.null(fit_A_res), fit_A_res$PA, NA)
      bias_A <- ifelse(!is.null(fit_A_res), fit_A_res$Bias, NA)
      msep_A <- ifelse(!is.null(fit_A_res), fit_A_res$MSEP, NA)
      
      pa_AD   <- ifelse(!is.null(fit_AD_res), fit_AD_res$PA, NA)
      bias_AD <- ifelse(!is.null(fit_AD_res), fit_AD_res$Bias, NA)
      msep_AD <- ifelse(!is.null(fit_AD_res), fit_AD_res$MSEP, NA)
      
      if (!is.null(fit_A_res)) {
        trait_results[[length(trait_results) + 1]] <- data.frame(
          Model = "A", Trait = trait, Rep = r, Fold = f,
          N_validation = length(y_test),
          PA = pa_A, Bias = bias_A, MSEP = msep_A,
          Best_num_trees = fit_A_res$Best_Param$num.trees,
          Best_mtry = fit_A_res$Best_Param$mtry,
          Best_min_node_size = fit_A_res$Best_Param$min.node.size
        )
      }
      
      if (!is.null(fit_AD_res)) {
        trait_results[[length(trait_results) + 1]] <- data.frame(
          Model = "AD", Trait = trait, Rep = r, Fold = f,
          N_validation = length(y_test),
          PA = pa_AD, Bias = bias_AD, MSEP = msep_AD,
          Best_num_trees = fit_AD_res$Best_Param$num.trees,
          Best_mtry = fit_AD_res$Best_Param$mtry,
          Best_min_node_size = fit_AD_res$Best_Param$min.node.size
        )
      }
      
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
  stop("No cross-validation results were generated. Check error prints above.")
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
  paste0("**Threads Used:** ", n_threads),
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
cat("RANDOM FOREST SUMMARY\n")
cat("======================================================================\n")
print(summary_results)
