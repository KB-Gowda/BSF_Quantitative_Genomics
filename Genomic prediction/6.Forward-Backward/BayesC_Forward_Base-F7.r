# ==============================================================================
# FORWARD GENOMIC PREDICTION (BASE -> F7) - BayesC (BGLR)
# Scenario 1: Univariate BayesC (Pooled Base -> F7)
# Scenario 2: Multitrait SpikeSlab (SYK, BSG, FVW) ->
#             1. EBV_SYK vs F7
#             2. EBV_BSG vs F7
#             3. EBV_FVW vs F7
#             4. EBV_3Diet_Avg vs F7
# Models: Model A (Additive) & Model AD (Additive + Dominance)
# Target Metrics: PA, Bias, MSEP (Pre-Adjusted Phenotypes)
# ==============================================================================

# ------------------------------------------------------------------------------
# 1. SETUP & LOAD WORKSPACE
# ------------------------------------------------------------------------------
library(BGLR)
library(dplyr)
library(Matrix)
library(tidyr)

run_prefix <- "BayesC_Forward_Base_to_F7"
n_iter     <- 12000
burn_in    <- 2000

load("FwdBwd_prep.RData")
set.seed(20260825)

# Verify matrix identity and sample order
sample_ids <- as.character(pheno_final$Sample_Id)
stopifnot(
  all(rownames(Xa_final) == sample_ids),
  all(rownames(Xd_final) == sample_ids)
)

Xa <- as.matrix(Xa_final)
Xd <- as.matrix(Xd_final)

target_ids <- as.character(pheno_final$Sample_Id[scenarios$forward$valid])
train_ids  <- as.character(pheno_final$Sample_Id[scenarios$forward$train])

idx_f7   <- scenarios$forward$valid
idx_base <- scenarios$forward$train

cat(sprintf("Loaded dataset: %d total samples (Base Train N = %d | F7 Target N = %d)\n",
            nrow(pheno_final), length(train_ids), length(target_ids)))

# ------------------------------------------------------------------------------
# 2. METRIC CALCULATION HELPER
# ------------------------------------------------------------------------------
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
# 3. TRAIT CONFIGURATIONS
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
# 4. BayesC FORWARD PREDICTION EXECUTION LOOP
# ------------------------------------------------------------------------------
for (tg in names(trait_groups)) {
  target_tr <- trait_groups[[tg]]$target
  diets_tr  <- trait_groups[[tg]]$mv_diets
  
  cat("\n", paste(rep("=", 85), collapse = ""), "\n", sep = "")
  cat("PROCESSING TRAIT GROUP:", tg, "\n")
  cat(paste(rep("=", 85), collapse = ""), "\n")
  
  # Target observed values in F7
  obs_df <- pheno_final[pheno_final$Population == "F7", c("Sample_Id", target_tr)]
  colnames(obs_df) <- c("ID", "Obs")
  obs_df$ID  <- trimws(as.character(obs_df$ID))
  obs_df$Obs <- as.numeric(obs_df$Obs)
  obs_df     <- obs_df[!is.na(obs_df$Obs), ]
  valid_ids  <- intersect(target_ids, obs_df$ID)
  
  # ============================================================================
  # SCENARIO 1: UNIVARIATE BAYESC (POOLED BASE -> F7)
  # ============================================================================
  cat("\n--- Running Scenario 1: Univariate BayesC (Pooled Base) ---\n")
  y_uni <- as.numeric(pheno_final[[target_tr]])
  y_uni[idx_f7] <- NA
  
  # Model A: Additive BayesC
  fit_bc_A <- tryCatch({
    BGLR(y = y_uni, ETA = list(Additive = list(X = Xa, model = "BayesC")),
         nIter = n_iter, burnIn = burn_in, verbose = FALSE)
  }, error = function(e) NULL)
  
  if (!is.null(fit_bc_A)) {
    gebv_A <- as.vector(Xa %*% fit_bc_A$ETA[[1]]$b)
    pred_A <- data.frame(ID = target_ids, gEBV = gebv_A[idx_f7], stringsAsFactors = FALSE)
    
    res <- calc_metrics(pred_A, obs_df, "Scenario_1_BayesC", "Pooled_Base", tg, "A")
    if (!is.null(res)) all_results[[length(all_results) + 1]] <- res
  }
  
  # Model AD: Additive + Dominance BayesC
  fit_bc_AD <- tryCatch({
    BGLR(y = y_uni, 
         ETA = list(Additive  = list(X = Xa, model = "BayesC"),
                    Dominance = list(X = Xd, model = "BayesC")),
         nIter = n_iter, burnIn = burn_in, verbose = FALSE)
  }, error = function(e) NULL)
  
  if (!is.null(fit_bc_AD)) {
    gebv_AD <- as.vector(Xa %*% fit_bc_AD$ETA[[1]]$b)
    gedv_AD <- as.vector(Xd %*% fit_bc_AD$ETA[[2]]$b)
    pred_AD <- data.frame(ID = target_ids, gEBV = (gebv_AD + gedv_AD)[idx_f7], stringsAsFactors = FALSE)
    
    res <- calc_metrics(pred_AD, obs_df, "Scenario_1_BayesC", "Pooled_Base", tg, "AD")
    if (!is.null(res)) all_results[[length(all_results) + 1]] <- res
  }
  
  # ============================================================================
  # SCENARIO 2: MULTITRAIT SPIKESLAB (SYK, BSG, FVW & 3-DIET AVERAGE vs F7)
  # ============================================================================
  cat("\n--- Running Scenario 2: Multitrait SpikeSlab (Base Diets -> F7) ---\n")
  Y_mv <- as.matrix(pheno_final[, diets_tr])
  Y_mv[idx_f7, ] <- NA
  colnames(Y_mv) <- names(diets_tr)
  
  # Model A: Additive Multitrait SpikeSlab
  fit_ss_A <- tryCatch({
    Multitrait(y = Y_mv, 
               ETA = list(Additive = list(X = Xa, model = "SpikeSlab")),
               nIter = n_iter, burnIn = burn_in, verbose = FALSE)
  }, error = function(e) NULL)
  
  if (!is.null(fit_ss_A)) {
    beta_A     <- fit_ss_A$ETA[[1]]$beta
    pred_mat_A <- Xa %*% beta_A
    colnames(pred_mat_A) <- names(diets_tr)
    pred_f7_A  <- pred_mat_A[idx_f7, , drop = FALSE]
    
    # 1. Individual diets: SYK, BSG, FVW
    for (d_name in names(diets_tr)) {
      sub_p <- data.frame(ID = target_ids, gEBV = pred_f7_A[, d_name], stringsAsFactors = FALSE)
      res   <- calc_metrics(sub_p, obs_df, "Scenario_2_SpikeSlab", paste0("EBV_", d_name), tg, "A")
      if (!is.null(res)) all_results[[length(all_results) + 1]] <- res
    }
    
    # 2. 4th predictor: Average across all 3 diets per individual
    avg_p_A <- data.frame(ID = target_ids, gEBV = rowMeans(pred_f7_A, na.rm = TRUE), stringsAsFactors = FALSE)
    res_avg <- calc_metrics(avg_p_A, obs_df, "Scenario_2_SpikeSlab", "EBV_3Diet_Avg", tg, "A")
    if (!is.null(res_avg)) all_results[[length(all_results) + 1]] <- res_avg
  }
  
  # Model AD: Additive + Dominance Multitrait SpikeSlab
  fit_ss_AD <- tryCatch({
    Multitrait(y = Y_mv,
               ETA = list(Additive  = list(X = Xa, model = "SpikeSlab"),
                          Dominance = list(X = Xd, model = "SpikeSlab")),
               nIter = n_iter, burnIn = burn_in, verbose = FALSE)
  }, error = function(e) NULL)
  
  if (!is.null(fit_ss_AD)) {
    beta_a_AD   <- fit_ss_AD$ETA[[1]]$beta
    beta_d_AD   <- fit_ss_AD$ETA[[2]]$beta
    pred_mat_AD <- (Xa %*% beta_a_AD) + (Xd %*% beta_d_AD)
    colnames(pred_mat_AD) <- names(diets_tr)
    pred_f7_AD  <- pred_mat_AD[idx_f7, , drop = FALSE]
    
    # 1. Individual diets: SYK, BSG, FVW
    for (d_name in names(diets_tr)) {
      sub_p <- data.frame(ID = target_ids, gEBV = pred_f7_AD[, d_name], stringsAsFactors = FALSE)
      res   <- calc_metrics(sub_p, obs_df, "Scenario_2_SpikeSlab", paste0("EBV_", d_name), tg, "AD")
      if (!is.null(res)) all_results[[length(all_results) + 1]] <- res
    }
    
    # 2. 4th predictor: Average across all 3 diets per individual
    avg_p_AD <- data.frame(ID = target_ids, gEBV = rowMeans(pred_f7_AD, na.rm = TRUE), stringsAsFactors = FALSE)
    res_avg  <- calc_metrics(avg_p_AD, obs_df, "Scenario_2_SpikeSlab", "EBV_3Diet_Avg", tg, "AD")
    if (!is.null(res_avg)) all_results[[length(all_results) + 1]] <- res_avg
  }
  flush.console()
}

# ------------------------------------------------------------------------------
# 5. CONSOLIDATE RESULTS & BENCHMARK AGAINST SCENARIO 1
# ------------------------------------------------------------------------------
if (length(all_results) == 0) stop("No results generated.")

results_df <- bind_rows(all_results)

# Benchmark against Scenario 1 (BayesC Pooled) baseline
s1_base <- results_df %>%
  filter(Scenario == "Scenario_1_BayesC") %>%
  dplyr::select(Trait, Model, PA_Base = PA)

final_table <- results_df %>%
  left_join(s1_base, by = c("Trait", "Model")) %>%
  mutate(
    Delta_PA = PA - PA_Base,
    Verdict  = case_when(
      Scenario == "Scenario_1_BayesC" ~ "Baseline (BayesC Pooled)",
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
  paste0("# Forward Genomic Prediction: Univariate BayesC vs Multitrait SpikeSlab"),
  "",
  paste0("**Generated on:** ", Sys.time()),
  "**Training Generation:** Base (n = 2,083) | **Validation Generation:** F7 (n = 536)",
  "",
  "### Comparison against Observed F7 Performance",
  "",
  "| Trait | Model | Scenario | Predictor Used | PA | Bias | MSEP | Δ PA (vs BayesC Pooled) | Verdict |",
  "|:---|:---:|:---|:---|:---:|:---:|:---:|:---:|:---:|"
)

for (i in seq_len(nrow(final_table))) {
  r <- final_table[i, ]
  delta_str <- if (r$Scenario == "Scenario_1_BayesC") "-" else sprintf("%+.4f", r$Delta_PA)
  md_lines <- c(md_lines, sprintf(
    "| %s | %s | %s | %s | %.4f | %.4f | %.4f | %s | %s |",
    r$Trait, r$Model, r$Scenario, r$Predictor,
    r$PA, r$Bias, r$MSEP, delta_str, r$Verdict
  ))
}

writeLines(md_lines, con = md_file)
cat("\n✓ Markdown report successfully written to:", md_file, "\n")
cat("✓ CSV results successfully written to:", paste0(run_prefix, "_Comparison_Results.csv"), "\n")
