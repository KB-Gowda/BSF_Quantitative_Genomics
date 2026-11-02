# ==============================================================================
# BACKWARD GENOMIC PREDICTION (F7 -> BASE) - RANGER (RANDOM FOREST)
# Training Set: F7 Generation (n = 536)
# Validation Set: Base Generation (n = 2,083)
# Optimization: Nested Out-of-Bag (OOB) MSE Minimization
# Models Evaluated: Model A (Additive) & Model AD (Additive + Dominance)
# Evaluations:
#   1. Whole Base Population (Adj_Weight / Adj_SurfaceArea)
#   2. Diet-Specific Base Subsets (SYK, BSG, FVW vs Adj_Weight / Adj_SurfaceArea)
# Target Metrics: PA, Bias, MSEP (Using Pre-Adjusted Phenotypes Only)
# ==============================================================================

# ------------------------------------------------------------------------------
# 1. SETUP & LOAD WORKSPACE
# ------------------------------------------------------------------------------
library(ranger)
library(dplyr)
library(Matrix)
library(tidyr)

run_prefix <- "Ranger_Backward_F7_to_Base"
n_threads  <- 6
set.seed(20260825)

load("FwdBwd_prep.RData")

pheno_final <- pheno_final %>% mutate(Sample_Id = trimws(as.character(Sample_Id)))

# Backward scenario indices: Train on F7 (536), Validate on Base (2,083)
idx_train_f7 <- scenarios$backward$train
idx_valid_b  <- scenarios$backward$valid

train_ids <- as.character(pheno_final$Sample_Id[idx_train_f7])
valid_ids <- as.character(pheno_final$Sample_Id[idx_valid_b])

cat(sprintf("Loaded dataset: %d total samples (F7 Train N = %d | Base Target N = %d)\n",
            nrow(pheno_final), length(train_ids), length(valid_ids)))

# Prepare feature matrices directly from Xa_final and Xd_final
df_Xa <- as.data.frame(as.matrix(Xa_final))
colnames(df_Xa) <- paste0("m_A_", seq_len(ncol(df_Xa)))

df_Xd <- as.data.frame(as.matrix(Xd_final))
colnames(df_Xd) <- paste0("m_D_", seq_len(ncol(df_Xd)))

df_Xad <- cbind(df_Xa, df_Xd)

# Hyperparameter tuning grid generator
get_grid <- function(p_vars) {
  mtry_vals <- unique(c(
    floor(sqrt(p_vars)),
    floor(p_vars / 10),
    floor(p_vars / 5),
    floor(p_vars / 3)
  ))
  mtry_vals <- mtry_vals[mtry_vals >= 1]
  
  expand.grid(
    num_trees     = c(150, 300),
    mtry          = mtry_vals,
    min.node.size = c(5, 10),
    stringsAsFactors = FALSE
  )
}

grid_A  <- get_grid(ncol(df_Xa))
grid_AD <- get_grid(ncol(df_Xad))

# ------------------------------------------------------------------------------
# 2. NESTED OOB HYPERPARAMETER TUNING HELPER
# ------------------------------------------------------------------------------
tune_and_fit_ranger <- function(y_vec, X_train, X_test, param_grid, threads) {
  val_rows <- which(!is.na(y_vec))
  y_train  <- y_vec[val_rows]
  X_tr     <- X_train[val_rows, , drop = FALSE]
  
  train_df <- cbind(y = y_train, X_tr)
  
  best_oob_mse <- Inf
  best_param   <- param_grid[1, ]
  
  for (g in seq_len(nrow(param_grid))) {
    fit_tune <- tryCatch({
      ranger(
        y ~ .,
        data          = train_df,
        num.trees     = param_grid$num_trees[g],
        mtry          = param_grid$mtry[g],
        min.node.size = param_grid$min.node.size[g],
        num.threads   = threads,
        seed          = 20260825 + g,
        verbose       = FALSE
      )
    }, error = function(e) NULL)
    
    if (!is.null(fit_tune)) {
      oob_mse <- fit_tune$prediction.error
      if (!is.na(oob_mse) && oob_mse < best_oob_mse) {
        best_oob_mse <- oob_mse
        best_param   <- param_grid[g, ]
      }
    }
  }
  
  final_fit <- ranger(
    y ~ .,
    data          = train_df,
    num.trees     = best_param$num_trees,
    mtry          = best_param$mtry,
    min.node.size = best_param$min.node.size,
    num.threads   = threads,
    seed          = 20260825,
    verbose       = FALSE
  )
  
  preds <- predict(final_fit, data = X_test, num.threads = threads)$predictions
  return(as.numeric(preds))
}

# Metric calculation helper
calc_metrics <- function(pred_df, obs_df, scope_name, diet_subset, predictor_name, trait_name, model_type) {
  d <- inner_join(pred_df, obs_df, by = "ID") %>% filter(!is.na(Obs), !is.na(gEBV))
  if (nrow(d) > 2 && sd(d$gEBV, na.rm = TRUE) > 0) {
    pa   <- cor(d$gEBV, d$Obs, use = "complete.obs")
    bias <- as.numeric(coef(lm(Obs ~ gEBV, data = d))[2])
    p_c  <- d$gEBV - mean(d$gEBV, na.rm = TRUE)
    o_c  <- d$Obs - mean(d$Obs, na.rm = TRUE)
    msep <- mean((o_c - p_c)^2, na.rm = TRUE)
    
    cat(sprintf("  [%-14s | Diet: %-4s | Model %-2s] PA: %.4f | Bias: %.4f | MSEP: %.4f (N = %d)\n",
                scope_name, diet_subset, model_type, pa, bias, msep, nrow(d)))
    
    return(data.frame(
      Trait        = trait_name,
      Scope        = scope_name,
      Diet_Subset  = diet_subset,
      Predictor    = predictor_name,
      Model        = model_type,
      N_Valid      = nrow(d),
      PA           = pa,
      Bias         = bias,
      MSEP         = msep,
      stringsAsFactors = FALSE
    ))
  }
  return(NULL)
}

# ------------------------------------------------------------------------------
# 3. DEFINE TRAIT GROUPS
# ------------------------------------------------------------------------------
trait_groups <- list(
  Weight      = list(train_trait = "Adj_Weight", target_base = "Adj_Weight"),
  SurfaceArea = list(train_trait = "Adj_SurfaceArea", target_base = "Adj_SurfaceArea")
)

all_results <- list()

# ------------------------------------------------------------------------------
# 4. BACKWARD RANGER EXECUTION LOOP
# ------------------------------------------------------------------------------
for (tg in names(trait_groups)) {
  train_tr  <- trait_groups[[tg]]$train_trait
  target_tr <- trait_groups[[tg]]$target_base
  
  cat("\n", paste(rep("=", 85), collapse = ""), "\n", sep = "")
  cat("PROCESSING BACKWARD RANGER FOR TRAIT:", tg, "\n")
  cat("Training on F7:", train_tr, "| Validating on Base:", target_tr, "\n")
  cat(paste(rep("=", 85), collapse = ""), "\n")
  
  # Target observed values in Base (retaining Diet factor for subgroup evaluation)
  base_obs_df <- pheno_final[idx_valid_b, c("Sample_Id", "Diet", target_tr)]
  colnames(base_obs_df) <- c("ID", "Diet", "Obs")
  base_obs_df$ID  <- trimws(as.character(base_obs_df$ID))
  base_obs_df$Obs <- as.numeric(base_obs_df$Obs)
  
  # Training vector from F7
  y_train_f7 <- pheno_final[[train_tr]][idx_train_f7]
  
  # Features partitioned: Train on F7, Predict on Base
  X_A_f7   <- df_Xa[idx_train_f7, , drop = FALSE]
  X_A_base <- df_Xa[idx_valid_b, , drop = FALSE]
  
  X_AD_f7   <- df_Xad[idx_train_f7, , drop = FALSE]
  X_AD_base <- df_Xad[idx_valid_b, , drop = FALSE]
  
  # ============================================================================
  # MODEL A: ADDITIVE RANGER (TRAINED ON F7 -> PREDICT BASE)
  # ============================================================================
  cat("\n--- Running Model A: Additive Ranger (Trained on F7) ---\n")
  pred_A <- tune_and_fit_ranger(y_train_f7, X_A_f7, X_A_base, grid_A, n_threads)
  df_pred_A <- data.frame(ID = valid_ids, gEBV = pred_A, stringsAsFactors = FALSE)
  
  # 4a. Whole Base evaluation
  res <- calc_metrics(df_pred_A, base_obs_df, "Whole_Base", "All", "EBV_F7_Trained", tg, "A")
  if (!is.null(res)) all_results[[length(all_results) + 1]] <- res
  
  # 4b. Diet-specific Base evaluation against overall adjusted phenotype
  for (d in c("SYK", "BSG", "FVW")) {
    sub_obs <- base_obs_df %>% filter(Diet == d)
    res <- calc_metrics(df_pred_A, sub_obs, "Diet_Specific", d, "EBV_F7_Trained", tg, "A")
    if (!is.null(res)) all_results[[length(all_results) + 1]] <- res
  }
  
  # ============================================================================
  # MODEL AD: ADDITIVE + DOMINANCE RANGER (TRAINED ON F7 -> PREDICT BASE)
  # ============================================================================
  cat("\n--- Running Model AD: Additive + Dominance Ranger (Trained on F7) ---\n")
  pred_AD <- tune_and_fit_ranger(y_train_f7, X_AD_f7, X_AD_base, grid_AD, n_threads)
  df_pred_AD <- data.frame(ID = valid_ids, gEBV = pred_AD, stringsAsFactors = FALSE)
  
  # 4a. Whole Base evaluation
  res <- calc_metrics(df_pred_AD, base_obs_df, "Whole_Base", "All", "GPV_F7_Trained", tg, "AD")
  if (!is.null(res)) all_results[[length(all_results) + 1]] <- res
  
  # 4b. Diet-specific Base evaluation against overall adjusted phenotype
  for (d in c("SYK", "BSG", "FVW")) {
    sub_obs <- base_obs_df %>% filter(Diet == d)
    res <- calc_metrics(df_pred_AD, sub_obs, "Diet_Specific", d, "GPV_F7_Trained", tg, "AD")
    if (!is.null(res)) all_results[[length(all_results) + 1]] <- res
  }
  
  flush.console()
}

# ------------------------------------------------------------------------------
# 5. CONSOLIDATE RESULTS & BENCHMARK SUMMARY TABLE
# ------------------------------------------------------------------------------
if (length(all_results) == 0) stop("No results generated.")

results_df <- bind_rows(all_results)

# Benchmark against Whole Base Model A baseline
whole_base_A <- results_df %>%
  filter(Scope == "Whole_Base", Model == "A") %>%
  dplyr::select(Trait, PA_Base_A = PA)

final_table <- results_df %>%
  left_join(whole_base_A, by = "Trait") %>%
  mutate(
    Delta_PA = PA - PA_Base_A,
    Performance = case_when(
      Scope == "Whole_Base" & Model == "A" ~ "Baseline (Whole Base A)",
      Delta_PA > 0.005                     ~ "Higher than Baseline",
      Delta_PA < -0.005                    ~ "Lower than Baseline",
      TRUE                                 ~ "Comparable"
    )
  )

print(as.data.frame(final_table))
write.csv(final_table, paste0(run_prefix, "_Comparison_Results.csv"), row.names = FALSE)

# ------------------------------------------------------------------------------
# 6. GENERATE AND SAVE MARKDOWN REPORT
# ------------------------------------------------------------------------------
md_file <- paste0(run_prefix, "_Comparison_Report.md")

md_lines <- c(
  paste0("# Backward Genomic Prediction: Ranger Random Forest (F7 -> Base)"),
  "",
  paste0("**Generated on:** ", Sys.time()),
  "**Training Generation:** F7 (n = 536) | **Validation Generation:** Base (n = 2,083)",
  "**Optimization:** Nested Out-Of-Bag (OOB) MSE Minimization",
  "- **Pre-Adjusted Phenotypes:** Loaded directly from `FwdBwd_prep.RData`",
  "- **Models:** Model A (Additive markers `Xa`) & Model AD (Additive + Dominance markers `Xa + Xd`)",
  "",
  "### Backward Prediction Accuracy & Validation Metrics",
  "",
  "| Trait | Model | Evaluation Scope | Diet Subset | Predictor | N Valid | PA | Bias | MSEP | Δ PA (vs Whole A) | Performance |",
  "|:---|:---:|:---|:---:|:---|:---:|:---:|:---:|:---:|:---:|:---:|"
)

for (i in seq_len(nrow(final_table))) {
  r <- final_table[i, ]
  delta_str <- if (r$Scope == "Whole_Base" && r$Model == "A") "-" else sprintf("%+.4f", r$Delta_PA)
  md_lines <- c(md_lines, sprintf(
    "| %s | %s | %s | %s | %s | %d | %.4f | %.4f | %.4f | %s | %s |",
    r$Trait, r$Model, r$Scope, r$Diet_Subset, r$Predictor,
    r$N_Valid, r$PA, r$Bias, r$MSEP, delta_str, r$Performance
  ))
}

writeLines(md_lines, con = md_file)
cat("\n✓ Markdown report successfully written to:", md_file, "\n")
cat("✓ CSV results successfully written to:", paste0(run_prefix, "_Comparison_Results.csv"), "\n")
