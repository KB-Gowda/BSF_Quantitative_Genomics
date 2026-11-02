# ==============================================================================
# CROSS-DIET GENOMIC PREDICTION PIPELINE - RANDOM FOREST (BIVARIATE AVERAGING)
# Models Evaluated: Model A (Additive-Only) vs Model AD (Additive + Dominance)
# Target Metrics: PA, Bias, MSEP (Using Pre-Adjusted Phenotypes Only)
# ==============================================================================


library(dplyr)
library(Matrix)
library(randomForestSRC)

run_prefix <- "RF_CrossDiet_Multivariate"
n_threads  <- 6
options(rf.cores = n_threads)
set.seed(20260825)

load("BasePop_prep.RData")

pheno_final <- pheno_final %>% mutate(Sample_Id = trimws(as.character(Sample_Id)))

multivariate_sets <- list(
  Weight      = c("Adj_SYK_Weight", "Adj_BSG_Weight", "Adj_FVW_Weight"),
  SurfaceArea = c("Adj_SYK_SurfaceArea", "Adj_BSG_SurfaceArea", "Adj_FVW_SurfaceArea")
)

# Prepare genotype dataframes directly from Xa_final and Xd_final in BasePop_prep.RData
df_Xa <- as.data.frame(Xa_final)
colnames(df_Xa) <- paste0("m_A_", seq_len(ncol(df_Xa)))

df_Xd <- as.data.frame(Xd_final)
colnames(df_Xd) <- paste0("m_D_", seq_len(ncol(df_Xd)))

df_Xad <- cbind(df_Xa, df_Xd)

get_grid <- function(p_vars) {
  mtry_vals <- unique(c(
    floor(sqrt(p_vars)),
    floor(p_vars / 10),
    floor(p_vars / 5),
    floor(p_vars / 3)
  ))
  mtry_vals <- mtry_vals[mtry_vals >= 1]
  
  expand.grid(
    ntree    = c(150, 300),
    mtry     = mtry_vals,
    nodesize = c(5, 10),
    stringsAsFactors = FALSE
  )
}

grid_A  <- get_grid(ncol(df_Xa))
grid_AD <- get_grid(ncol(df_Xad))

all_results <- list()

# ------------------------------------------------------------------------------
# CROSS-DIET PREDICTION LOOP
# ------------------------------------------------------------------------------
for (set_name in names(multivariate_sets)) {
  current_traits <- multivariate_sets[[set_name]]
  
  cat("\n", paste(rep("=", 80), collapse = ""), "\n", sep = "")
  cat("STARTING CROSS-DIET RANDOM FOREST FOR SET:", set_name, "\n")
  cat("Traits in set:", paste(current_traits, collapse = ", "), "\n")
  cat(paste(rep("=", 80), collapse = ""), "\n")
  
  set_results <- list()
  
  for (target_trait in current_traits) {
    train_traits <- setdiff(current_traits, target_trait)
    
    train_idx <- which(!is.na(pheno_final[[train_traits[1]]]) | !is.na(pheno_final[[train_traits[2]]]))
    test_idx  <- which(!is.na(pheno_final[[target_trait]]))
    
    if (length(test_idx) == 0 || length(train_idx) == 0) next
    
    cat(sprintf("\n--- Target Held-Out Trait: %s | Trained on: %s + %s (Target N = %d) ---\n", 
                target_trait, train_traits[1], train_traits[2], length(test_idx)))
    
    obs_test   <- pheno_final[[target_trait]][test_idx]
    formula_rf <- as.formula(paste("Multivar(", train_traits[1], ",", train_traits[2], ") ~ ."))
    
    # ========================================================================
    # MODEL A: ADDITIVE-ONLY RANDOM FOREST
    # ========================================================================
    train_data_A <- cbind(pheno_final[train_idx, train_traits, drop = FALSE], df_Xa[train_idx, , drop = FALSE])
    test_data_A  <- df_Xa[test_idx, , drop = FALSE]
    
    best_msep <- Inf; best_param_A <- grid_A[1, ]
    
    for (g in seq_len(nrow(grid_A))) {
      set.seed(20260825 + g)
      tune_A <- tryCatch({
        rfsrc(
          formula = formula_rf, data = train_data_A,
          ntree = grid_A$ntree[g], mtry = grid_A$mtry[g], nodesize = grid_A$nodesize[g],
          splitrule = "mv", na.action = "na.impute"
        )
      }, error = function(e) NULL)
      
      if (!is.null(tune_A)) {
        msep_sum <- 0
        for (tr in train_traits) {
          y_obs <- train_data_A[[tr]]
          y_oob <- tune_A$regrOutput[[tr]]$predicted.oob
          valid <- !is.na(y_obs) & !is.na(y_oob)
          if (sum(valid) > 1) msep_sum <- msep_sum + mean((y_obs[valid] - y_oob[valid])^2)
        }
        if (msep_sum < best_msep) {
          best_msep    <- msep_sum
          best_param_A <- grid_A[g, ]
        }
      }
    }
    
    set.seed(20260825)
    fit_A <- tryCatch({
      rfsrc(
        formula = formula_rf, data = train_data_A,
        ntree = best_param_A$ntree, mtry = best_param_A$mtry, nodesize = best_param_A$nodesize,
        splitrule = "mv", na.action = "na.impute"
      )
    }, error = function(e) NULL)
    
    if (!is.null(fit_A)) {
      pred_obj_A <- predict(fit_A, newdata = test_data_A, na.action = "na.impute")
      pred_t1_A  <- pred_obj_A$regrOutput[[train_traits[1]]]$predicted
      pred_t2_A  <- pred_obj_A$regrOutput[[train_traits[2]]]$predicted
      
      pred_target_A <- rowMeans(cbind(pred_t1_A, pred_t2_A), na.rm = TRUE)
      
      valid_A <- which(!is.na(obs_test) & !is.na(pred_target_A))
      if (length(valid_A) > 1 && sd(pred_target_A[valid_A], na.rm = TRUE) > 0) {
        p_val <- pred_target_A[valid_A]; o_val <- obs_test[valid_A]
        pa    <- cor(p_val, o_val, use = "complete.obs")
        bias  <- as.numeric(coef(lm(o_val ~ p_val))[2])
        p_c   <- p_val - mean(p_val, na.rm = TRUE)
        o_c   <- o_val - mean(o_val, na.rm = TRUE)
        msep  <- mean((o_c - p_c)^2, na.rm = TRUE)
        
        set_results[[length(set_results) + 1]] <- data.frame(
          Model = "A", Set = set_name, Target_Trait = target_trait, 
          N_Validation = length(valid_A), PA = pa, Bias = bias, MSEP = msep
        )
        cat(sprintf("  [Model A]  PA: %.4f | Bias: %.4f | MSEP: %.4f\n", pa, bias, msep))
      }
    }
    rm(fit_A)
    
    # ========================================================================
    # MODEL AD: ADDITIVE + DOMINANCE RANDOM FOREST
    # ========================================================================
    train_data_AD <- cbind(pheno_final[train_idx, train_traits, drop = FALSE], df_Xad[train_idx, , drop = FALSE])
    test_data_AD  <- df_Xad[test_idx, , drop = FALSE]
    
    best_msep <- Inf; best_param_AD <- grid_AD[1, ]
    
    for (g in seq_len(nrow(grid_AD))) {
      set.seed(20260825 + g)
      tune_AD <- tryCatch({
        rfsrc(
          formula = formula_rf, data = train_data_AD,
          ntree = grid_AD$ntree[g], mtry = grid_AD$mtry[g], nodesize = grid_AD$nodesize[g],
          splitrule = "mv", na.action = "na.impute"
        )
      }, error = function(e) NULL)
      
      if (!is.null(tune_AD)) {
        msep_sum <- 0
        for (tr in train_traits) {
          y_obs <- train_data_AD[[tr]]
          y_oob <- tune_AD$regrOutput[[tr]]$predicted.oob
          valid <- !is.na(y_obs) & !is.na(y_oob)
          if (sum(valid) > 1) msep_sum <- msep_sum + mean((y_obs[valid] - y_oob[valid])^2)
        }
        if (msep_sum < best_msep) {
          best_msep     <- msep_sum
          best_param_AD <- grid_AD[g, ]
        }
      }
    }
    
    set.seed(20260825)
    fit_AD <- tryCatch({
      rfsrc(
        formula = formula_rf, data = train_data_AD,
        ntree = best_param_AD$ntree, mtry = best_param_AD$mtry, nodesize = best_param_AD$nodesize,
        splitrule = "mv", na.action = "na.impute"
      )
    }, error = function(e) NULL)
    
    if (!is.null(fit_AD)) {
      pred_obj_AD <- predict(fit_AD, newdata = test_data_AD, na.action = "na.impute")
      pred_t1_AD  <- pred_obj_AD$regrOutput[[train_traits[1]]]$predicted
      pred_t2_AD  <- pred_obj_AD$regrOutput[[train_traits[2]]]$predicted
      
      pred_target_AD <- rowMeans(cbind(pred_t1_AD, pred_t2_AD), na.rm = TRUE)
      
      valid_AD <- which(!is.na(obs_test) & !is.na(pred_target_AD))
      if (length(valid_AD) > 1 && sd(pred_target_AD[valid_AD], na.rm = TRUE) > 0) {
        p_val <- pred_target_AD[valid_AD]; o_val <- obs_test[valid_AD]
        pa    <- cor(p_val, o_val, use = "complete.obs")
        bias  <- as.numeric(coef(lm(o_val ~ p_val))[2])
        p_c   <- p_val - mean(p_val, na.rm = TRUE)
        o_c   <- o_val - mean(o_val, na.rm = TRUE)
        msep  <- mean((o_c - p_c)^2, na.rm = TRUE)
        
        set_results[[length(set_results) + 1]] <- data.frame(
          Model = "AD", Set = set_name, Target_Trait = target_trait, 
          N_Validation = length(valid_AD), PA = pa, Bias = bias, MSEP = msep
        )
        cat(sprintf("  [Model AD] PA: %.4f | Bias: %.4f | MSEP: %.4f\n", pa, bias, msep))
      }
    }
    rm(fit_AD)
    gc(verbose = FALSE)
    flush.console()
  }
  
  if (length(set_results) > 0) {
    df_set <- bind_rows(set_results)
    all_results[[set_name]] <- df_set
    write.csv(df_set, paste0(run_prefix, "_", set_name, "_Raw.csv"), row.names = FALSE)
  }
}

# ------------------------------------------------------------------------------
# COMBINE & SUMMARIZE RESULTS
# ------------------------------------------------------------------------------
if (length(all_results) == 0) stop("No results generated.")

raw_results <- bind_rows(all_results)

summary_results <- raw_results %>%
  group_by(Model, Set, Target_Trait) %>%
  summarise(
    N_Evaluations = sum(N_Validation),
    PA   = mean(PA, na.rm = TRUE),
    Bias = mean(Bias, na.rm = TRUE),
    MSEP = mean(MSEP, na.rm = TRUE),
    .groups = "drop"
  )

write.csv(raw_results, paste0(run_prefix, "_Raw_Results.csv"), row.names = FALSE)
write.csv(summary_results, paste0(run_prefix, "_Summary_Results.csv"), row.names = FALSE)

md_file <- paste0(run_prefix, "_Summary_Report.md")
md_content <- c(
  paste0("# ", gsub("_", " ", run_prefix), " Summary Report"),
  "",
  paste0("**Date Generated:** ", Sys.time()),
  "**Cross-Diet Strategy:** Leave-One-Diet-Out (Bivariate Averaging on Pre-Adjusted Traits)",
  "",
  "| Model | Set | Target Trait | N Validation | PA | Bias | MSEP |",
  "|:---|:---|:---|:---:|:---:|:---:|:---:|"
)

for (i in seq_len(nrow(summary_results))) {
  row <- summary_results[i, ]
  md_line <- sprintf("| %s | %s | %s | %d | %.4f | %.4f | %.4f |",
                     row$Model, row$Set, row$Target_Trait, row$N_Evaluations,
                     row$PA, row$Bias, row$MSEP)
  md_content <- c(md_content, md_line)
}

writeLines(md_content, con = md_file)
cat("\nMarkdown report successfully generated:", md_file, "\n")
print(as.data.frame(summary_results))
