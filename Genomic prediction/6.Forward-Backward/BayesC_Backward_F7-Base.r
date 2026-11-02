# ==============================================================================
# BACKWARD GENOMIC PREDICTION (F7 -> BASE) - BAYESC (BGLR)
# Models Evaluated: Model A (Additive) & Model AD (Additive + Dominance)
# Evaluations:
#   1. Whole Base Population (Overall Adjusted Phenotypes)
#   2. Diet-Specific Base Subsets (SYK, BSG, FVW)
# Target Metrics: PA, Bias, MSEP (Using Pre-Adjusted Phenotypes Only)
# ==============================================================================

# ------------------------------------------------------------------------------
# 1. SETUP & LOAD WORKSPACE
# ------------------------------------------------------------------------------
library(BGLR)
library(dplyr)
library(Matrix)
library(tidyr)

run_prefix <- "BayesC_Backward_F7_to_Base"
n_iter     <- 12000
burn_in    <- 2000

load("FwdBwd_prep.RData")
set.seed(20260825)

# Verify identity and matrix alignments
sample_ids <- as.character(pheno_final$Sample_Id)
stopifnot(
  all(rownames(Xa_final) == sample_ids),
  all(rownames(Xd_final) == sample_ids)
)

Xa <- as.matrix(Xa_final)
Xd <- as.matrix(Xd_final)

# Backward scenario indices: Train on F7, Validate on Base
idx_train_f7 <- scenarios$backward$train   # F7 individuals (n = 536)
idx_valid_b  <- scenarios$backward$valid   # Base individuals (n = 2,083)

train_ids <- as.character(pheno_final$Sample_Id[idx_train_f7])
valid_ids <- as.character(pheno_final$Sample_Id[idx_valid_b])

cat(sprintf("Loaded dataset: %d total samples (F7 Train N = %d | Base Target N = %d)\n",
            nrow(pheno_final), length(train_ids), length(valid_ids)))

# ------------------------------------------------------------------------------
# 2. METRIC CALCULATION HELPER
# ------------------------------------------------------------------------------
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
# 3. TRAIT CONFIGURATIONS
# ------------------------------------------------------------------------------
trait_groups <- list(
  Weight = list(
    train_trait = "Adj_Weight",
    target_base = "Adj_Weight",
    diet_traits = c(SYK = "Adj_SYK_Weight", BSG = "Adj_BSG_Weight", FVW = "Adj_FVW_Weight")
  ),
  SurfaceArea = list(
    train_trait = "Adj_SurfaceArea",
    target_base = "Adj_SurfaceArea",
    diet_traits = c(SYK = "Adj_SYK_SurfaceArea", BSG = "Adj_BSG_SurfaceArea", FVW = "Adj_FVW_SurfaceArea")
  )
)

all_results <- list()

# ------------------------------------------------------------------------------
# 4. BACKWARD BAYESC EXECUTION LOOP
# ------------------------------------------------------------------------------
for (tg in names(trait_groups)) {
  train_tr    <- trait_groups[[tg]]$train_trait
  target_tr   <- trait_groups[[tg]]$target_base
  diet_traits <- trait_groups[[tg]]$diet_traits
  
  cat("\n", paste(rep("=", 85), collapse = ""), "\n", sep = "")
  cat("PROCESSING BACKWARD PREDICTION FOR TRAIT:", tg, "\n")
  cat("Training on F7:", train_tr, "| Validating on Base:", target_tr, "\n")
  cat(paste(rep("=", 85), collapse = ""), "\n")
  
  # Base observed phenotypes (with Diet annotation)
  base_obs_df <- pheno_final[idx_valid_b, c("Sample_Id", "Diet", target_tr, unname(diet_traits))]
  colnames(base_obs_df)[colnames(base_obs_df) == "Sample_Id"] <- "ID"
  colnames(base_obs_df)[colnames(base_obs_df) == target_tr]   <- "Obs_Overall"
  base_obs_df$ID          <- trimws(as.character(base_obs_df$ID))
  base_obs_df$Obs_Overall <- as.numeric(base_obs_df$Obs_Overall)
  
  # Phenotype vector: mask Base records to NA so training is strictly on F7
  y_train <- as.numeric(pheno_final[[train_tr]])
  y_train[idx_valid_b] <- NA
  
  # ============================================================================
  # MODEL A: ADDITIVE BAYESC (TRAIN ON F7 -> PREDICT ALL)
  # ============================================================================
  cat("\n--- Fitting Model A (Additive BayesC on F7) ---\n")
  fit_A <- tryCatch({
    BGLR(
      y       = y_train,
      ETA     = list(Additive = list(X = Xa, model = "BayesC")),
      nIter   = n_iter,
      burnIn  = burn_in,
      verbose = FALSE
    )
  }, error = function(e) {
    cat("  [ERROR Model A]:", e$message, "\n")
    NULL
  })
  
  if (!is.null(fit_A)) {
    gebv_A <- as.vector(Xa %*% fit_A$ETA[[1]]$b)
    pred_A_base <- data.frame(ID = valid_ids, gEBV = gebv_A[idx_valid_b], stringsAsFactors = FALSE)
    
    # 4a. Whole Base evaluation (Overall trait)
    obs_overall <- base_obs_df %>% dplyr::select(ID, Obs = Obs_Overall)
    res <- calc_metrics(pred_A_base, obs_overall, "Whole_Base", "All", "EBV_F7_Trained", tg, "A")
    if (!is.null(res)) all_results[[length(all_results) + 1]] <- res
    
    # 4b. Diet-specific Base evaluation against overall adjusted phenotype
    for (d in c("SYK", "BSG", "FVW")) {
      sub_obs <- base_obs_df %>% filter(Diet == d) %>% dplyr::select(ID, Obs = Obs_Overall)
      res <- calc_metrics(pred_A_base, sub_obs, "Diet_Specific_Overall", d, "EBV_F7_Trained", tg, "A")
      if (!is.null(res)) all_results[[length(all_results) + 1]] <- res
    }
    
    # 4c. Diet-specific Base evaluation against within-diet adjusted residuals
    for (d in names(diet_traits)) {
      d_col <- diet_traits[[d]]
      sub_obs_within <- base_obs_df %>% filter(Diet == d) %>% dplyr::select(ID, Obs = all_of(d_col))
      res <- calc_metrics(pred_A_base, sub_obs_within, "Within_Diet_Residual", d, "EBV_F7_Trained", tg, "A")
      if (!is.null(res)) all_results[[length(all_results) + 1]] <- res
    }
  }
  
  # ============================================================================
  # MODEL AD: ADDITIVE + DOMINANCE BAYESC (TRAIN ON F7 -> PREDICT ALL)
  # ============================================================================
  cat("\n--- Fitting Model AD (Additive + Dominance BayesC on F7) ---\n")
  fit_AD <- tryCatch({
    BGLR(
      y       = y_train,
      ETA     = list(
        Additive  = list(X = Xa, model = "BayesC"),
        Dominance = list(X = Xd, model = "BayesC")
      ),
      nIter   = n_iter,
      burnIn  = burn_in,
      verbose = FALSE
    )
  }, error = function(e) {
    cat("  [ERROR Model AD]:", e$message, "\n")
    NULL
  })
  
  if (!is.null(fit_AD)) {
    gebv_AD <- as.vector(Xa %*% fit_AD$ETA[[1]]$b)
    gedv_AD <- as.vector(Xd %*% fit_AD$ETA[[2]]$b)
    gpv_AD  <- gebv_AD + gedv_AD
    pred_AD_base <- data.frame(ID = valid_ids, gEBV = gpv_AD[idx_valid_b], stringsAsFactors = FALSE)
    
    # 4a. Whole Base evaluation (Overall trait)
    obs_overall <- base_obs_df %>% dplyr::select(ID, Obs = Obs_Overall)
    res <- calc_metrics(pred_AD_base, obs_overall, "Whole_Base", "All", "GPV_F7_Trained", tg, "AD")
    if (!is.null(res)) all_results[[length(all_results) + 1]] <- res
    
    # 4b. Diet-specific Base evaluation against overall adjusted phenotype
    for (d in c("SYK", "BSG", "FVW")) {
      sub_obs <- base_obs_df %>% filter(Diet == d) %>% dplyr::select(ID, Obs = Obs_Overall)
      res <- calc_metrics(pred_AD_base, sub_obs, "Diet_Specific_Overall", d, "GPV_F7_Trained", tg, "AD")
      if (!is.null(res)) all_results[[length(all_results) + 1]] <- res
    }
    
    # 4c. Diet-specific Base evaluation against within-diet adjusted residuals
    for (d in names(diet_traits)) {
      d_col <- diet_traits[[d]]
      sub_obs_within <- base_obs_df %>% filter(Diet == d) %>% dplyr::select(ID, Obs = all_of(d_col))
      res <- calc_metrics(pred_AD_base, sub_obs_within, "Within_Diet_Residual", d, "GPV_F7_Trained", tg, "AD")
      if (!is.null(res)) all_results[[length(all_results) + 1]] <- res
    }
  }
  flush.console()
}

# ------------------------------------------------------------------------------
# 5. CONSOLIDATE RESULTS & BENCHMARK SUMMARY TABLE
# ------------------------------------------------------------------------------
if (length(all_results) == 0) stop("No results generated.")

results_df <- bind_rows(all_results)

# Extract Whole Base Model A baseline to compute Delta PA
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
  paste0("# Backward Genomic Prediction: BayesC (F7 -> Base)"),
  "",
  paste0("**Generated on:** ", Sys.time()),
  "**Training Generation:** F7 (n = 536) | **Validation Generation:** Base (n = 2,083)",
  "- **Pre-Adjusted Phenotypes:** Loaded directly from `FwdBwd_prep.RData`",
  "- **Models:** Model A (Additive) & Model AD (Additive + Dominance)",
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
