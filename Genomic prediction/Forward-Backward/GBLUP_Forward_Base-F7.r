# ==============================================================================
# FORWARD GENOMIC PREDICTION (BASE -> F7) - GBLUP
# Scenario 1: Pooled Base -> F7
# Scenario 2: Multivariate Base (SYK, BSG, FVW) -> 
#             1. EBV_SYK vs F7
#             2. EBV_BSG vs F7
#             3. EBV_FVW vs F7
#             4. EBV_3Diet_Avg (Average of 3 Diets) vs F7
# Models: Model A (Additive) & Model AD (Additive + Dominance)
# ==============================================================================

library(asreml)
library(dplyr)
library(Matrix)
library(tidyr)

run_prefix     <- "GBLUP_Forward_Base_to_F7"
n_threads      <- 4
workspace_size <- "8GB"

options(sinksToConsole = TRUE)
sink(stdout(), type = "output")
asreml.options(nthreads = n_threads)

load("FwdBwd_prep.RData")
set.seed(20260825)

# Factor alignment
id_levels <- as.character(pheno_final$Sample_Id)
pheno_final$Sample_Id <- factor(pheno_final$Sample_Id, levels = id_levels)

target_ids <- as.character(pheno_final$Sample_Id[scenarios$forward$valid])
train_ids  <- as.character(pheno_final$Sample_Id[scenarios$forward$train])

# Ensure numeric traits
eval_traits <- c("Adj_Weight", "Adj_SurfaceArea",
                 "Adj_SYK_Weight", "Adj_BSG_Weight", "Adj_FVW_Weight",
                 "Adj_SYK_SurfaceArea", "Adj_BSG_SurfaceArea", "Adj_FVW_SurfaceArea")
for (tr in eval_traits) {
  if (tr %in% colnames(pheno_final)) pheno_final[[tr]] <- as.numeric(pheno_final[[tr]])
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

trait_groups <- list(
  Weight = list(
    target    = "Adj_Weight",
    mv_diets  = c(SYK = "Adj_SYK_Weight", BSG = "Adj_BSG_Weight", FVW = "Adj_FVW_Weight")
  ),
  SurfaceArea = list(
    target    = "Adj_SurfaceArea",
    mv_diets  = c(SYK = "Adj_SYK_SurfaceArea", BSG = "Adj_BSG_SurfaceArea", FVW = "Adj_FVW_SurfaceArea")
  )
)

all_results <- list()

for (tg in names(trait_groups)) {
  target_tr <- trait_groups[[tg]]$target
  diets_tr  <- trait_groups[[tg]]$mv_diets
  
  cat("\n", paste(rep("=", 85), collapse = ""), "\n", sep = "")
  cat("PROCESSING TRAIT GROUP:", tg, "\n")
  cat(paste(rep("=", 85), collapse = ""), "\n")
  
  # F7 target observed values
  obs_df <- pheno_final[pheno_final$Population == "F7", c("Sample_Id", target_tr)]
  colnames(obs_df) <- c("ID", "Obs")
  obs_df$ID  <- trimws(as.character(obs_df$ID))
  obs_df$Obs <- as.numeric(obs_df$Obs)
  obs_df     <- obs_df[!is.na(obs_df$Obs), ]
  valid_ids  <- intersect(target_ids, obs_df$ID)
  
  # ============================================================================
  # SCENARIO 1: POOLED BASE TRAINED -> PREDICT F7
  # ============================================================================
  cat("\n--- Running Scenario 1 (Pooled Base Training) ---\n")
  df_s1 <- pheno_final
  df_s1[[target_tr]][df_s1$Population == "F7"] <- NA
  
  # Model A
  mod_s1_A <- tryCatch({
    fit <- asreml(fixed = as.formula(paste(target_tr, "~ 1")),
                  random = ~ vm(Sample_Id, Ginv),
                  residual = ~ idv(units),
                  data = df_s1,
                  na.action = na.method(x = "include", y = "include"),
                  maxit = 100, ai.sing = TRUE, workspace = workspace_size, trace = FALSE)
    update.asreml(fit, trace = FALSE)
  }, error = function(e) NULL)
  
  if (!is.null(mod_s1_A) && mod_s1_A$converge) {
    rand_A <- as.data.frame(summary(mod_s1_A, coef = TRUE)$coef.random)
    p_A <- data.frame(Term = rownames(rand_A), gEBV = rand_A[, 1]) %>%
      filter(grepl("vm\\(Sample_Id, Ginv\\)", Term)) %>%
      mutate(ID = trimws(gsub("vm\\(Sample_Id, Ginv\\)_", "", Term))) %>%
      filter(ID %in% valid_ids)
    
    res <- calc_metrics(p_A, obs_df, "Scenario_1", "Pooled_Base", tg, "A")
    if (!is.null(res)) all_results[[length(all_results) + 1]] <- res
  }
  
  # Model AD
  mod_s1_AD <- tryCatch({
    fit <- asreml(fixed = as.formula(paste(target_tr, "~ 1")),
                  random = ~ vm(Sample_Id, Ginv) + vm(Sample_Id, Gd.inv),
                  residual = ~ idv(units),
                  data = df_s1,
                  na.action = na.method(x = "include", y = "include"),
                  maxit = 100, ai.sing = TRUE, workspace = workspace_size, trace = FALSE)
    update.asreml(fit, trace = FALSE)
  }, error = function(e) NULL)
  
  if (!is.null(mod_s1_AD) && mod_s1_AD$converge) {
    rand_AD <- as.data.frame(summary(mod_s1_AD, coef = TRUE)$coef.random)
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
    
    res <- calc_metrics(p_AD, obs_df, "Scenario_1", "Pooled_Base", tg, "AD")
    if (!is.null(res)) all_results[[length(all_results) + 1]] <- res
  }
  
  # ============================================================================
  # SCENARIO 2: MULTIVARIATE BASE -> 3 DIET EBVs + 3-DIET AVERAGE vs F7
  # ============================================================================
  cat("\n--- Running Scenario 2 (Multivariate: SYK, BSG, FVW & 3-Diet Avg) ---\n")
  df_s2 <- pheno_final
  for (tr in diets_tr) {
    df_s2[[tr]][df_s2$Population == "F7"] <- NA
  }
  
  mv_formula <- as.formula(paste("cbind(", paste(diets_tr, collapse = ", "), ") ~ trait"))
  
  # --- Multivariate Model A ---
  mod_s2_A <- tryCatch({
    fit <- asreml(fixed = mv_formula,
                  random = ~ us(trait):vm(Sample_Id, Ginv),
                  residual = ~ id(units):us(trait),
                  data = df_s2,
                  na.action = na.method(x = "include", y = "include"),
                  maxit = 100, ai.sing = TRUE, workspace = workspace_size, trace = FALSE)
    update.asreml(fit, trace = FALSE)
  }, error = function(e) NULL)
  
  if (!is.null(mod_s2_A) && mod_s2_A$converge) {
    rand_A <- as.data.frame(summary(mod_s2_A, coef = TRUE)$coef.random)
    gebv_diet_A <- data.frame(Term = rownames(rand_A), Additive = rand_A[, 1]) %>%
      filter(grepl("vm\\(Sample_Id, Ginv\\)", Term)) %>%
      mutate(Term = gsub("^trait_", "", Term),
             Term = gsub(":vm\\(Sample_Id, Ginv\\)_", "||", Term)) %>%
      separate(Term, into = c("Diet_Trait", "ID"), sep = "\\|\\|") %>%
      mutate(ID = trimws(as.character(ID))) %>%
      filter(ID %in% valid_ids)
    
    # 1. Individual diets: SYK, BSG, FVW
    for (d_name in names(diets_tr)) {
      full_trait <- diets_tr[[d_name]]
      sub_pred <- gebv_diet_A %>%
        filter(Diet_Trait == full_trait) %>%
        rename(gEBV = Additive)
      
      res <- calc_metrics(sub_pred, obs_df, "Scenario_2", paste0("EBV_", d_name), tg, "A")
      if (!is.null(res)) all_results[[length(all_results) + 1]] <- res
    }
    
    # 2. 4th predictor: Average of the 3 diet EBVs per individual
    avg_pred_A <- gebv_diet_A %>%
      group_by(ID) %>%
      summarise(gEBV = mean(Additive, na.rm = TRUE), .groups = "drop")
    
    res_avg_A <- calc_metrics(avg_pred_A, obs_df, "Scenario_2", "EBV_3Diet_Avg", tg, "A")
    if (!is.null(res_avg_A)) all_results[[length(all_results) + 1]] <- res_avg_A
  }
  
  # --- Multivariate Model AD ---
  mod_s2_AD <- tryCatch({
    fit <- asreml(fixed = mv_formula,
                  random = ~ us(trait):vm(Sample_Id, Ginv) + us(trait):vm(Sample_Id, Gd.inv),
                  residual = ~ id(units):us(trait),
                  data = df_s2,
                  na.action = na.method(x = "include", y = "include"),
                  maxit = 100, ai.sing = TRUE, workspace = workspace_size, trace = FALSE)
    update.asreml(fit, trace = FALSE)
  }, error = function(e) NULL)
  
  if (!is.null(mod_s2_AD) && mod_s2_AD$converge) {
    rand_AD <- as.data.frame(summary(mod_s2_AD, coef = TRUE)$coef.random)
    
    add_df <- data.frame(Term = rownames(rand_AD), Additive = rand_AD[, 1]) %>%
      filter(grepl("vm\\(Sample_Id, Ginv\\)", Term)) %>%
      mutate(Term = gsub("^trait_", "", Term),
             Term = gsub(":vm\\(Sample_Id, Ginv\\)_", "||", Term)) %>%
      separate(Term, into = c("Diet_Trait", "ID"), sep = "\\|\\|") %>%
      mutate(ID = trimws(as.character(ID))) %>%
      filter(ID %in% valid_ids)
    
    dom_df <- data.frame(Term = rownames(rand_AD), Dominance = rand_AD[, 1]) %>%
      filter(grepl("vm\\(Sample_Id, Gd.inv\\)", Term)) %>%
      mutate(Term = gsub("^trait_", "", Term),
             Term = gsub(":vm\\(Sample_Id, Gd.inv\\)_", "||", Term)) %>%
      separate(Term, into = c("Diet_Trait", "ID"), sep = "\\|\\|") %>%
      mutate(ID = trimws(as.character(ID))) %>%
      filter(ID %in% valid_ids)
    
    gpv_diet_AD <- inner_join(add_df, dom_df, by = c("Diet_Trait", "ID")) %>%
      mutate(gEBV = Additive + Dominance)
    
    # 1. Individual diets: SYK, BSG, FVW
    for (d_name in names(diets_tr)) {
      full_trait <- diets_tr[[d_name]]
      sub_pred <- gpv_diet_AD %>%
        filter(Diet_Trait == full_trait)
      
      res <- calc_metrics(sub_pred, obs_df, "Scenario_2", paste0("EBV_", d_name), tg, "AD")
      if (!is.null(res)) all_results[[length(all_results) + 1]] <- res
    }
    
    # 2. 4th predictor: Average of the 3 diet genotypic values per individual
    avg_pred_AD <- gpv_diet_AD %>%
      group_by(ID) %>%
      summarise(gEBV = mean(gEBV, na.rm = TRUE), .groups = "drop")
    
    res_avg_AD <- calc_metrics(avg_pred_AD, obs_df, "Scenario_2", "EBV_3Diet_Avg", tg, "AD")
    if (!is.null(res_avg_AD)) all_results[[length(all_results) + 1]] <- res_avg_AD
  }
  flush.console()
}

# ------------------------------------------------------------------------------
# CONSOLIDATE AND PRINT BENCHMARK TABLE
# ------------------------------------------------------------------------------
if (length(all_results) == 0) stop("No results generated.")

results_df <- bind_rows(all_results)

# Benchmark against Scenario 1 baseline
s1_base <- results_df %>%
  filter(Scenario == "Scenario_1") %>%
  dplyr::select(Trait, Model, PA_Base = PA)

final_table <- results_df %>%
  left_join(s1_base, by = c("Trait", "Model")) %>%
  mutate(
    Delta_PA = PA - PA_Base,
    Verdict  = case_when(
      Scenario == "Scenario_1" ~ "Baseline (Pooled)",
      Delta_PA > 0.005         ~ "Outperforms Pooled",
      Delta_PA < -0.005        ~ "Underperforms Pooled",
      TRUE                     ~ "Comparable"
    )
  )

print(as.data.frame(final_table))
write.csv(final_table, paste0(run_prefix, "_Comparison_Results.csv"), row.names = FALSE)

# ------------------------------------------------------------------------------
# GENERATE AND SAVE MARKDOWN REPORT
# ------------------------------------------------------------------------------
md_file <- paste0(run_prefix, "_Comparison_Report.md")

md_lines <- c(
  paste0("# Forward Genomic Prediction: GBLUP (Base -> F7)"),
  "",
  paste0("**Generated on:** ", Sys.time()),
  "**Training Generation:** Base (n = 2,083) | **Validation Generation:** F7 (n = 536)",
  "",
  "### Comparison against Observed F7 Performance",
  "",
  "| Trait | Model | Scenario | Predictor Used | PA | Bias | MSEP | Δ PA (vs Pooled) | Verdict |",
  "|:---|:---:|:---|:---|:---:|:---:|:---:|:---:|:---:|"
)

for (i in seq_len(nrow(final_table))) {
  r <- final_table[i, ]
  delta_str <- if (r$Scenario == "Scenario_1") "-" else sprintf("%+.4f", r$Delta_PA)
  md_lines <- c(md_lines, sprintf(
    "| %s | %s | %s | %s | %.4f | %.4f | %.4f | %s | %s |",
    r$Trait, r$Model, r$Scenario, r$Predictor,
    r$PA, r$Bias, r$MSEP, delta_str, r$Verdict
  ))
}

writeLines(md_lines, con = md_file)
cat("\n✓ Markdown report successfully written to:", md_file, "\n")
cat("✓ CSV results successfully written to:", paste0(run_prefix, "_Comparison_Results.csv"), "\n")
