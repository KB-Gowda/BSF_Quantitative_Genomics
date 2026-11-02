# ==============================================================================
# BACKWARD GENOMIC PREDICTION (F7 -> BASE) - GBLUP (ASREML)
# Training Set: F7 Generation (n = 536)
# Validation Set: Base Generation (n = 2,083)
# Models Evaluated: Model A (Additive) & Model AD (Additive + Dominance)
# Evaluations:
#   1. Whole Base Population (Adj_Weight / Adj_SurfaceArea)
#   2. Diet-Specific Base Subsets (SYK, BSG, FVW vs Adj_Weight / Adj_SurfaceArea)
# Target Metrics: PA, Bias, MSEP (Using Pre-Adjusted Phenotypes Only)
# ==============================================================================

# ------------------------------------------------------------------------------
# 1. SETUP & LOAD WORKSPACE
# ------------------------------------------------------------------------------
library(asreml)
library(dplyr)
library(Matrix)
library(tidyr)

run_prefix     <- "GBLUP_Backward_F7_to_Base"
n_threads      <- 4
workspace_size <- "8GB"

options(sinksToConsole = TRUE)
sink(stdout(), type = "output")
asreml.options(nthreads = n_threads)

load("FwdBwd_prep.RData")
set.seed(20260825)

# Factor alignment to ensure levels match Ginv and Gd.inv dimensions
id_levels <- as.character(pheno_final$Sample_Id)
pheno_final$Sample_Id <- factor(pheno_final$Sample_Id, levels = id_levels)

# Backward scenario indices: Train on F7, Validate on Base
idx_train_f7 <- scenarios$backward$train   # F7 individuals (n = 536)
idx_valid_b  <- scenarios$backward$valid   # Base individuals (n = 2,083)

train_ids <- as.character(pheno_final$Sample_Id[idx_train_f7])
valid_ids <- as.character(pheno_final$Sample_Id[idx_valid_b])

cat(sprintf("Loaded dataset: %d total samples (F7 Train N = %d | Base Target N = %d)\n",
            nrow(pheno_final), length(train_ids), length(valid_ids)))

# Ensure target traits are numeric
eval_traits <- c("Adj_Weight", "Adj_SurfaceArea")
for (tr in eval_traits) {
  if (tr %in% colnames(pheno_final)) pheno_final[[tr]] <- as.numeric(pheno_final[[tr]])
}

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
# 3. DEFINE TRAIT GROUPS
# ------------------------------------------------------------------------------
trait_groups <- list(
  Weight      = list(train_trait = "Adj_Weight", target_base = "Adj_Weight"),
  SurfaceArea = list(train_trait = "Adj_SurfaceArea", target_base = "Adj_SurfaceArea")
)

all_results <- list()

# ------------------------------------------------------------------------------
# 4. BACKWARD GBLUP EXECUTION LOOP
# ------------------------------------------------------------------------------
for (tg in names(trait_groups)) {
  train_tr  <- trait_groups[[tg]]$train_trait
  target_tr <- trait_groups[[tg]]$target_base
  
  cat("\n", paste(rep("=", 85), collapse = ""), "\n", sep = "")
  cat("PROCESSING BACKWARD GBLUP FOR TRAIT:", tg, "\n")
  cat("Training on F7:", train_tr, "| Validating on Base:", target_tr, "\n")
  cat(paste(rep("=", 85), collapse = ""), "\n")
  
  # Base observed phenotypes (retaining Diet factor for subgroup evaluation)
  base_obs_df <- pheno_final[idx_valid_b, c("Sample_Id", "Diet", target_tr)]
  colnames(base_obs_df) <- c("ID", "Diet", "Obs")
  base_obs_df$ID  <- trimws(as.character(base_obs_df$ID))
  base_obs_df$Obs <- as.numeric(base_obs_df$Obs)
  
  # Mask Base population phenotypes so ASReml trains strictly on F7
  df_bwd <- pheno_final
  df_bwd[[train_tr]][df_bwd$Population == "Base"] <- NA
  
  # ============================================================================
  # MODEL A: ADDITIVE GBLUP (TRAIN F7 -> PREDICT ALL)
  # ============================================================================
  cat("\n--- Running Model A: Additive GBLUP (Trained on F7) ---\n")
  mod_A <- tryCatch({
    fit <- asreml(
      fixed     = as.formula(paste(train_tr, "~ 1")),
      random    = ~ vm(Sample_Id, Ginv),
      residual  = ~ idv(units),
      data      = df_bwd,
      na.action = na.method(x = "include", y = "include"),
      maxit     = 100, ai.sing = TRUE, workspace = workspace_size, trace = FALSE
    )
    update.asreml(fit, trace = FALSE)
  }, error = function(e) {
    cat("  [ERROR Model A]:", e$message, "\n")
    NULL
  })
  
  if (!is.null(mod_A) && mod_A$converge) {
    rand_A <- as.data.frame(summary(mod_A, coef = TRUE)$coef.random)
    p_A <- data.frame(Term = rownames(rand_A), gEBV = rand_A[, 1]) %>%
      filter(grepl("vm\\(Sample_Id, Ginv\\)", Term)) %>%
      mutate(ID = trimws(gsub("vm\\(Sample_Id, Ginv\\)_", "", Term))) %>%
      filter(ID %in% valid_ids)
    
    # 4a. Whole Base evaluation (Overall trait)
    res <- calc_metrics(p_A, base_obs_df, "Whole_Base", "All", "EBV_F7_Trained", tg, "A")
    if (!is.null(res)) all_results[[length(all_results) + 1]] <- res
    
    # 4b. Diet-specific Base evaluation against overall adjusted phenotype
    for (d in c("SYK", "BSG", "FVW")) {
      sub_obs <- base_obs_df %>% filter(Diet == d)
      res <- calc_metrics(p_A, sub_obs, "Diet_Specific", d, "EBV_F7_Trained", tg, "A")
      if (!is.null(res)) all_results[[length(all_results) + 1]] <- res
    }
  }
  
  # ============================================================================
  # MODEL AD: ADDITIVE + DOMINANCE GBLUP (TRAIN F7 -> PREDICT ALL)
  # ============================================================================
  cat("\n--- Running Model AD: Additive + Dominance GBLUP (Trained on F7) ---\n")
  mod_AD <- tryCatch({
    fit <- asreml(
      fixed     = as.formula(paste(train_tr, "~ 1")),
      random    = ~ vm(Sample_Id, Ginv) + vm(Sample_Id, Gd.inv),
      residual  = ~ idv(units),
      data      = df_bwd,
      na.action = na.method(x = "include", y = "include"),
      maxit     = 100, ai.sing = TRUE, workspace = workspace_size, trace = FALSE
    )
    update.asreml(fit, trace = FALSE)
  }, error = function(e) {
    cat("  [ERROR Model AD]:", e$message, "\n")
    NULL
  })
  
  if (!is.null(mod_AD) && mod_AD$converge) {
    rand_AD <- as.data.frame(summary(mod_AD, coef = TRUE)$coef.random)
    
    add_df <- data.frame(Term = rownames(rand_AD), Additive = rand_AD[, 1]) %>%
      filter(grepl("vm\\(Sample_Id, Ginv\\)", Term)) %>%
      mutate(ID = trimws(gsub("vm\\(Sample_Id, Ginv\\)_", "", Term))) %>%
      filter(ID %in% valid_ids)
    
    dom_df <- data.frame(Term = rownames(rand_AD), Dominance = rand_AD[, 1]) %>%
      filter(grepl("vm\\(Sample_Id, Gd.inv\\)", Term)) %>%
      mutate(ID = trimws(gsub("vm\\(Sample_Id, Gd.inv\\)_", "", Term))) %>%
      filter(ID %in% valid_ids)
    
    p_AD <- inner_join(add_df, dom_df, by = "ID") %>%
      mutate(gEBV = Additive + Dominance) %>%
      dplyr::select(ID, gEBV)
    
    # 4a. Whole Base evaluation (Overall trait)
    res <- calc_metrics(p_AD, base_obs_df, "Whole_Base", "All", "GPV_F7_Trained", tg, "AD")
    if (!is.null(res)) all_results[[length(all_results) + 1]] <- res
    
    # 4b. Diet-specific Base evaluation against overall adjusted phenotype
    for (d in c("SYK", "BSG", "FVW")) {
      sub_obs <- base_obs_df %>% filter(Diet == d)
      res <- calc_metrics(p_AD, sub_obs, "Diet_Specific", d, "GPV_F7_Trained", tg, "AD")
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
  paste0("# Backward Genomic Prediction: GBLUP (F7 -> Base)"),
  "",
  paste0("**Generated on:** ", Sys.time()),
  "**Training Generation:** F7 (n = 536) | **Validation Generation:** Base (n = 2,083)",
  "- **Pre-Adjusted Phenotypes:** Loaded directly from `FwdBwd_prep.RData`",
  "- **Models:** Model A (Additive: `vm(Sample_Id, Ginv)`) & Model AD (Additive + Dominance: `vm(Sample_Id, Ginv) + vm(Sample_Id, Gd.inv)`)",
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
