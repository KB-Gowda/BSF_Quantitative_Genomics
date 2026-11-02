# ==============================================================================
# FORWARD GENOMIC PREDICTION (BASE -> F7) - RANDOM FOREST
# Scenario 1: Univariate Ranger (Fast Pooled Base -> F7 with OOB Tuning)
# Scenario 2: Multivariate randomForestSRC (Multivar(SYK, BSG, FVW) ~ .) ->
#             1. EBV_SYK vs F7
#             2. EBV_BSG vs F7
#             3. EBV_FVW vs F7
#             4. EBV_3Diet_Avg vs F7
# Models: Model A (Additive) & Model AD (Additive + Dominance)
# Target Metrics: PA, Bias, MSEP (Using Pre-Adjusted Phenotypes Only)
# ==============================================================================

# ------------------------------------------------------------------------------
# 1. SETUP & LOAD WORKSPACE
# ------------------------------------------------------------------------------
library(ranger)
library(randomForestSRC)
library(dplyr)
library(Matrix)
library(tidyr)

run_prefix <- "RF_Forward_Ranger_and_SRC"
n_threads  <- 6
options(rf.cores = n_threads)
set.seed(20260825)

load("FwdBwd_prep.RData")

pheno_final <- pheno_final %>% mutate(Sample_Id = trimws(as.character(Sample_Id)))

# Identify Base and F7 row indices
idx_base <- scenarios$forward$train
idx_f7   <- scenarios$forward$valid

target_ids <- as.character(pheno_final$Sample_Id[idx_f7])
train_ids  <- as.character(pheno_final$Sample_Id[idx_base])

cat(sprintf("Loaded dataset: %d total samples (Base Train N = %d | F7 Target N = %d)\n",
            nrow(pheno_final), length(train_ids), length(target_ids)))

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
    num_trees = c(150, 300),
    mtry      = mtry_vals,
    nodesize  = c(5, 10),
    stringsAsFactors = FALSE
  )
}

grid_A  <- get_grid(ncol(df_Xa))
grid_AD <- get_grid(ncol(df_Xad))

# ------------------------------------------------------------------------------
# 2. HELPER FUNCTIONS
# ------------------------------------------------------------------------------

# Ranger tuning and fitting helper (for Scenario 1)
tune_and_fit_ranger <- function(y_vec, X_train, X_test, param_grid, threads) {
  val_rows <- which(!is.na(y_vec))
  train_df <- cbind(y = y_vec[val_rows], X_train[val_rows, , drop = FALSE])
  
  best_oob_mse <- Inf
  best_param   <- param_grid[1, ]
  
  for (g in seq_len(nrow(param_grid))) {
    fit_tune <- tryCatch({
      ranger(
        y ~ ., data = train_df,
        num.trees     = param_grid$num_trees[g],
        mtry          = param_grid$mtry[g],
        min.node.size = param_grid$nodesize[g],
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
    y ~ ., data = train_df,
    num.trees     = best_param$num_trees,
    mtry          = best_param$mtry,
    min.node.size = best_param$nodesize,
    num.threads   = threads,
    seed          = 20260825,
    verbose       = FALSE
  )
  
  preds <- predict(final_fit, data = X_test, num.threads = threads)$predictions
  return(as.numeric(preds))
}

# Metric calculation helper
calc_metrics <- function(pred_df, obs_df, scenario_name, predictor_name, trait_name, model_type) {
  d <- inner_join(pred_df, obs_df, by = "ID") %>% filter(!is.na(Obs), !is.na(gEBV))
  if (nrow(d) > 1 && sd(d$gEBV, na.rm = TRUE) > 0) {
    pa   <- cor(d$gEBV, d$Obs, use = "complete.obs")
    bias <- as.numeric(coef(lm(Obs ~ gEBV, data = d))[2])
    p_c  <- d$gEBV - mean(d$gEBV, na.rm = TRUE)
    o_c  <- d$Obs - mean(d$Obs, na.rm = TRUE)
    msep <- mean((o_c - p_c)^2, na.rm = TRUE)
    
    cat(sprintf("  [%s | %-15s | Model %-2s] PA: %.4f | Bias: %.4f | MSEP: %.4f\n",
                scenario_name, predictor_name, model_type, pa, bias, msep))
    
    return(data.frame(
      Trait       = trait_name,
      Scenario    = scenario_name,
      Predictor   = predictor_name,
      Model       = model_type,
      N_Valid     = nrow(d),
      PA          = pa,
      Bias        = bias,
      MSEP        = msep,
      stringsAsFactors = FALSE
    ))
  }
  return(NULL)
}

# ------------------------------------------------------------------------------
# 3. DEFINE TRAIT GROUPS
# ------------------------------------------------------------------------------
trait_groups <- list(
  Weight = list(
    target   = "Adj_Weight",
    mv_diets = c(SYK = "Adj_SYK_Weight", BSG = "Adj_BSG_Weight", FVW = "Adj_FVW_Weight")
  ),
  SurfaceArea = list(
    target   = "Adj_SurfaceArea",
    mv_diets = c(SYK = "Adj_SYK_SurfaceArea", BSG = "Adj_BSG_SurfaceArea", FVW = "Adj_FVW_SurfaceArea")
  )
)

all_results <- list()

# ------------------------------------------------------------------------------
# 4. FORWARD EXECUTION LOOP
# ------------------------------------------------------------------------------
for (tg in names(trait_groups)) {
  target_tr <- trait_groups[[tg]]$target
  diets_tr  <- trait_groups[[tg]]$mv_diets
  
  cat("\n", paste(rep("=", 85), collapse = ""), "\n", sep = "")
  cat("PROCESSING TRAIT GROUP:", tg, "\n")
  cat(paste(rep("=", 85), collapse = ""), "\n")
  
  # Target observed values in F7
  obs_df <- pheno_final[idx_f7, c("Sample_Id", target_tr)]
  colnames(obs_df) <- c("ID", "Obs")
  obs_df$ID  <- trimws(as.character(obs_df$ID))
  obs_df$Obs <- as.numeric(obs_df$Obs)
  obs_df     <- obs_df[!is.na(obs_df$Obs), ]
  valid_ids  <- intersect(target_ids, obs_df$ID)
  
  # Features partitioned by generation
  X_A_base  <- df_Xa[idx_base, , drop = FALSE]
  X_A_f7    <- df_Xa[idx_f7, , drop = FALSE]
  
  X_AD_base <- df_Xad[idx_base, , drop = FALSE]
  X_AD_f7   <- df_Xad[idx_f7, , drop = FALSE]
  
  # ============================================================================
  # SCENARIO 1: UNIVARIATE RANGER (POOLED BASE -> F7)
  # ============================================================================
  cat("\n--- Running Scenario 1: Univariate Ranger (Pooled Base) ---\n")
  y_base_pooled <- pheno_final[[target_tr]][idx_base]
  
  # --- Model A: Additive Univariate Ranger ---
  pred_s1_A <- tune_and_fit_ranger(y_base_pooled, X_A_base, X_A_f7, grid_A, n_threads)
  df_pred_A <- data.frame(ID = target_ids, gEBV = pred_s1_A, stringsAsFactors = FALSE)
  res_s1_A  <- calc_metrics(df_pred_A, obs_df, "Scenario_1_Ranger", "Pooled_Base", tg, "A")
  if (!is.null(res_s1_A)) all_results[[length(all_results) + 1]] <- res_s1_A
  
  # --- Model AD: Additive + Dominance Univariate Ranger ---
  pred_s1_AD <- tune_and_fit_ranger(y_base_pooled, X_AD_base, X_AD_f7, grid_AD, n_threads)
  df_pred_AD <- data.frame(ID = target_ids, gEBV = pred_s1_AD, stringsAsFactors = FALSE)
  res_s1_AD  <- calc_metrics(df_pred_AD, obs_df, "Scenario_1_Ranger", "Pooled_Base", tg, "AD")
  if (!is.null(res_s1_AD)) all_results[[length(all_results) + 1]] <- res_s1_AD
  
  # ============================================================================
  # SCENARIO 2: MULTIVARIATE randomForestSRC (Multivar(SYK, BSG, FVW) ~ .)
  # ============================================================================
  cat("\n--- Running Scenario 2: Multivariate rfsrc (Multivar(SYK, BSG, FVW) ~ .) ---\n")
  
  train_sub_mv <- which(pheno_final$Population == "Base" & (
    !is.na(pheno_final[[diets_tr[1]]]) | 
    !is.na(pheno_final[[diets_tr[2]]]) | 
    !is.na(pheno_final[[diets_tr[3]]])
  ))
  
  formula_mv <- as.formula(paste("Multivar(", paste(diets_tr, collapse = ", "), ") ~ ."))
  
  # --- Multivariate Model A: Additive rfsrc ---
  train_data_A_mv <- cbind(pheno_final[train_sub_mv, diets_tr, drop = FALSE], df_Xa[train_sub_mv, , drop = FALSE])
  test_data_A     <- df_Xa[idx_f7, , drop = FALSE]
  
  best_msep <- Inf; best_p_A_mv <- grid_A[1, ]
  for (g in seq_len(nrow(grid_A))) {
    set.seed(20260825 + g)
    tune_A_mv <- tryCatch({
      rfsrc(formula = formula_mv, data = train_data_A_mv,
            ntree = grid_A$num_trees[g], mtry = grid_A$mtry[g], nodesize = grid_A$nodesize[g],
            splitrule = "mv", na.action = "na.impute")
    }, error = function(e) NULL)
    
    if (!is.null(tune_A_mv)) {
      msep_sum <- 0
      for (tr in diets_tr) {
        y_obs <- train_data_A_mv[[tr]]
        y_oob <- tune_A_mv$regrOutput[[tr]]$predicted.oob
        val   <- !is.na(y_obs) & !is.na(y_oob)
        if (sum(val) > 1) msep_sum <- msep_sum + mean((y_obs[val] - y_oob[val])^2)
      }
      if (msep_sum < best_msep) {
        best_msep   <- msep_sum
        best_p_A_mv <- grid_A[g, ]
      }
    }
  }
  
  set.seed(20260825)
  fit_A_mv <- tryCatch({
    rfsrc(formula = formula_mv, data = train_data_A_mv,
          ntree = best_p_A_mv$num_trees, mtry = best_p_A_mv$mtry, nodesize = best_p_A_mv$nodesize,
          splitrule = "mv", na.action = "na.impute")
  }, error = function(e) NULL)
  
  if (!is.null(fit_A_mv)) {
    pred_obj_A_mv <- predict(fit_A_mv, newdata = test_data_A, na.action = "na.impute")
    pred_mat_A    <- sapply(diets_tr, function(tr) pred_obj_A_mv$regrOutput[[tr]]$predicted)
    colnames(pred_mat_A) <- names(diets_tr)
    
    # 1. Individual diets: SYK, BSG, FVW
    for (d_name in names(diets_tr)) {
      sub_p <- data.frame(ID = target_ids, gEBV = as.numeric(pred_mat_A[, d_name]), stringsAsFactors = FALSE)
      res   <- calc_metrics(sub_p, obs_df, "Scenario_2_rfsrc", paste0("EBV_", d_name), tg, "A")
      if (!is.null(res)) all_results[[length(all_results) + 1]] <- res
    }
    
    # 2. 4th predictor: Average across all 3 diets per individual
    avg_p_A <- data.frame(ID = target_ids, gEBV = as.numeric(rowMeans(pred_mat_A, na.rm = TRUE)), stringsAsFactors = FALSE)
    res_avg <- calc_metrics(avg_p_A, obs_df, "Scenario_2_rfsrc", "EBV_3Diet_Avg", tg, "A")
    if (!is.null(res_avg)) all_results[[length(all_results) + 1]] <- res_avg
    
    rm(fit_A_mv, pred_obj_A_mv)
  }
  
  # --- Multivariate Model AD: Additive + Dominance rfsrc ---
  train_data_AD_mv <- cbind(pheno_final[train_sub_mv, diets_tr, drop = FALSE], df_Xad[train_sub_mv, , drop = FALSE])
  test_data_AD    <- df_Xad[idx_f7, , drop = FALSE]
  
  best_msep <- Inf; best_p_AD_mv <- grid_AD[1, ]
  for (g in seq_len(nrow(grid_AD))) {
    set.seed(20260825 + g)
    tune_AD_mv <- tryCatch({
      rfsrc(formula = formula_mv, data = train_data_AD_mv,
            ntree = grid_AD$num_trees[g], mtry = grid_AD$mtry[g], nodesize = grid_AD$nodesize[g],
            splitrule = "mv", na.action = "na.impute")
    }, error = function(e) NULL)
    
    if (!is.null(tune_AD_mv)) {
      msep_sum <- 0
      for (tr in diets_tr) {
        y_obs <- train_data_AD_mv[[tr]]
        y_oob <- tune_AD_mv$regrOutput[[tr]]$predicted.oob
        val   <- !is.na(y_obs) & !is.na(y_oob)
        if (sum(val) > 1) msep_sum <- msep_sum + mean((y_obs[val] - y_oob[val])^2)
      }
      if (msep_sum < best_msep) {
        best_msep    <- msep_sum
        best_p_AD_mv <- grid_AD[g, ]
      }
    }
  }
  
  set.seed(20260825)
  fit_AD_mv <- tryCatch({
    rfsrc(formula = formula_mv, data = train_data_AD_mv,
          ntree = best_p_AD_mv$num_trees, mtry = best_p_AD_mv$mtry, nodesize = best_p_AD_mv$nodesize,
          splitrule = "mv", na.action = "na.impute")
  }, error = function(e) NULL)
  
  if (!is.null(fit_AD_mv)) {
    pred_obj_AD_mv <- predict(fit_AD_mv, newdata = test_data_AD, na.action = "na.impute")
    pred_mat_AD    <- sapply(diets_tr, function(tr) pred_obj_AD_mv$regrOutput[[tr]]$predicted)
    colnames(pred_mat_AD) <- names(diets_tr)
    
    # 1. Individual diets: SYK, BSG, FVW
    for (d_name in names(diets_tr)) {
      sub_p <- data.frame(ID = target_ids, gEBV = as.numeric(pred_mat_AD[, d_name]), stringsAsFactors = FALSE)
      res   <- calc_metrics(sub_p, obs_df, "Scenario_2_rfsrc", paste0("EBV_", d_name), tg, "AD")
      if (!is.null(res)) all_results[[length(all_results) + 1]] <- res
    }
    
    # 2. 4th predictor: Average across all 3 diets per individual
    avg_p_AD <- data.frame(ID = target_ids, gEBV = as.numeric(rowMeans(pred_mat_AD, na.rm = TRUE)), stringsAsFactors = FALSE)
    res_avg  <- calc_metrics(avg_p_AD, obs_df, "Scenario_2_rfsrc", "EBV_3Diet_Avg", tg, "AD")
    if (!is.null(res_avg)) all_results[[length(all_results) + 1]] <- res_avg
    
    rm(fit_AD_mv, pred_obj_AD_mv)
  }
  
  gc(verbose = FALSE)
  flush.console()
}

# ------------------------------------------------------------------------------
# 5. CONSOLIDATE RESULTS & BENCHMARK AGAINST SCENARIO 1
# ------------------------------------------------------------------------------
if (length(all_results) == 0) stop("No results generated.")

results_df <- bind_rows(all_results)

# Benchmark against Scenario 1 baseline
s1_base <- results_df %>%
  filter(Scenario == "Scenario_1_Ranger") %>%
  dplyr::select(Trait, Model, PA_Base = PA)

final_table <- results_df %>%
  left_join(s1_base, by = c("Trait", "Model")) %>%
  mutate(
    Delta_PA = PA - PA_Base,
    Verdict  = case_when(
      Scenario == "Scenario_1_Ranger" ~ "Baseline (Ranger Pooled)",
      Delta_PA > 0.005                ~ "Outperforms Pooled",
      Delta_PA < -0.005               ~ "Underperforms Pooled",
      TRUE                            ~ "Comparable"
    )
  )

print(as.data.frame(final_table))
write.csv(final_table, paste0(run_prefix, "_Comparison_Results.csv"), row.names = FALSE)

# ------------------------------------------------------------------------------
# 6. GENERATE AND SAVE MARKDOWN REPORT
# ------------------------------------------------------------------------------
md_file <- paste0(run_prefix, "_Comparison_Report.md")

md_lines <- c(
  paste0("# Forward Genomic Prediction: Hybrid RF (Ranger vs rfsrc)"),
  "",
  paste0("**Generated on:** ", Sys.time()),
  "**Training Generation:** Base (n = 2,083) | **Validation Generation:** F7 (n = 536)",
  "- **Scenario 1:** Univariate `ranger` with OOB MSE tuning (Pooled Base)",
  "- **Scenario 2:** Multivariate `randomForestSRC` (`splitrule = 'mv'`) on SYK, BSG, FVW",
  "",
  "### Comparison against Observed F7 Performance",
  "",
  "| Trait | Model | Scenario | Predictor Used | PA | Bias | MSEP | Δ PA (vs Ranger Pooled) | Verdict |",
  "|:---|:---:|:---|:---|:---:|:---:|:---:|:---:|:---:|"
)

for (i in seq_len(nrow(final_table))) {
  r <- final_table[i, ]
  delta_str <- if (r$Scenario == "Scenario_1_Ranger") "-" else sprintf("%+.4f", r$Delta_PA)
  md_lines <- c(md_lines, sprintf(
    "| %s | %s | %s | %s | %.4f | %.4f | %.4f | %s | %s |",
    r$Trait, r$Model, r$Scenario, r$Predictor,
    r$PA, r$Bias, r$MSEP, delta_str, r$Verdict
  ))
}

writeLines(md_lines, con = md_file)
cat("\n✓ Markdown report successfully written to:", md_file, "\n")
cat("✓ CSV results successfully written to:", paste0(run_prefix, "_Comparison_Results.csv"), "\n")
