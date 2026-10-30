# ==============================================================================
# CROSS-DIET GENOMIC PREDICTION PIPELINE - SPIKE-SLAB
# Models Evaluated: Model A (Additive-Only) vs Model AD (Additive + Dominance)
# Target Metrics: PA, Bias, MSEP 
# ==============================================================================


library(dplyr)
library(Matrix)
library(BGLR)

run_prefix   <- "SpikeSlab_CrossDiet_Multivariate"
n_iter       <- 12000
burn_in      <- 2000
bayes_method <- "SpikeSlab"

load("BasePop_prep.RData")
set.seed(20260825)

pheno_final <- pheno_final %>% mutate(Sample_Id = trimws(as.character(Sample_Id)))

multivariate_sets <- list(
  Weight      = c("Adj_SYK_Weight", "Adj_BSG_Weight", "Adj_FVW_Weight"),
  SurfaceArea = c("Adj_SYK_SurfaceArea", "Adj_BSG_SurfaceArea", "Adj_FVW_SurfaceArea")
)

all_results <- list()

# ------------------------------------------------------------------------------
# CROSS-DIET PREDICTION LOOP
# ------------------------------------------------------------------------------
for (set_name in names(multivariate_sets)) {
  current_traits <- multivariate_sets[[set_name]]
  
  cat("\n", paste(rep("=", 80), collapse = ""), "\n", sep = "")
  cat("STARTING CROSS-DIET SPIKE-SLAB FOR SET:", set_name, "\n")
  cat("Traits in set:", paste(current_traits, collapse = ", "), "\n")
  cat(paste(rep("=", 80), collapse = ""), "\n")
  
  set_results <- list()
  
  for (target_trait in current_traits) {
    train_traits <- setdiff(current_traits, target_trait)
    
    target_idx <- which(!is.na(pheno_final[[target_trait]]))
    target_ids <- pheno_final$Sample_Id[target_idx]
    
    if (length(target_idx) == 0) next
    
    cat(sprintf("\n--- Target Held-Out Trait: %s | Trained on: %s + %s (Target N = %d) ---\n", 
                target_trait, train_traits[1], train_traits[2], length(target_ids)))
    
    # Response matrix using only the two known pre-adjusted traits
    Y_train <- as.matrix(pheno_final[, train_traits, drop = FALSE])
    colnames(Y_train) <- train_traits
    
    obs_df <- pheno_final[target_idx, c("Sample_Id", target_trait)]
    colnames(obs_df) <- c("Sample_Id", "Obs")
    
    # ========================================================================
    # MODEL A: ADDITIVE-ONLY SPIKE-SLAB
    # ========================================================================
    ETA_A <- list(Additive = list(X = Xa_final, model = bayes_method))
    
    fit_A <- tryCatch({
      Multitrait(y = Y_train, ETA = ETA_A, nIter = n_iter, burnIn = burn_in, verbose = FALSE)
    }, error = function(e) NULL)
    
    if (!is.null(fit_A)) {
      gebv_all_A <- Xa_final %*% fit_A$ETA[[1]]$beta
      colnames(gebv_all_A) <- train_traits
      
      pred_A_df <- data.frame(
        Sample_Id = pheno_final$Sample_Id,
        gEBV      = rowMeans(gebv_all_A, na.rm = TRUE)
      ) %>% filter(Sample_Id %in% target_ids)
      
      eval_A <- inner_join(pred_A_df, obs_df, by = "Sample_Id")
      
      if (nrow(eval_A) > 1 && sd(eval_A$gEBV, na.rm = TRUE) > 0) {
        pa   <- cor(eval_A$gEBV, eval_A$Obs, use = "complete.obs")
        bias <- as.numeric(coef(lm(Obs ~ gEBV, data = eval_A))[2])
        p_c  <- eval_A$gEBV - mean(eval_A$gEBV, na.rm = TRUE)
        o_c  <- eval_A$Obs - mean(eval_A$Obs, na.rm = TRUE)
        msep <- mean((o_c - p_c)^2, na.rm = TRUE)
        
        set_results[[length(set_results) + 1]] <- data.frame(
          Model = "A", Set = set_name, Target_Trait = target_trait, 
          N_Validation = nrow(eval_A), PA = pa, Bias = bias, MSEP = msep
        )
        cat(sprintf("  [Model A]  PA: %.4f | Bias: %.4f | MSEP: %.4f\n", pa, bias, msep))
      }
    }
    rm(fit_A)
    
    # ========================================================================
    # MODEL AD: ADDITIVE + DOMINANCE SPIKE-SLAB
    # ========================================================================
    ETA_AD <- list(
      Additive  = list(X = Xa_final, model = bayes_method),
      Dominance = list(X = Xd_final, model = bayes_method)
    )
    
    fit_AD <- tryCatch({
      Multitrait(y = Y_train, ETA = ETA_AD, nIter = n_iter, burnIn = burn_in, verbose = FALSE)
    }, error = function(e) NULL)
    
    if (!is.null(fit_AD)) {
      gebv_all_AD <- Xa_final %*% fit_AD$ETA[[1]]$beta
      gdbv_all_AD <- Xd_final %*% fit_AD$ETA[[2]]$beta
      gpv_all_AD  <- gebv_all_AD + gdbv_all_AD
      colnames(gpv_all_AD) <- train_traits
      
      pred_AD_df <- data.frame(
        Sample_Id = pheno_final$Sample_Id,
        gEBV      = rowMeans(gpv_all_AD, na.rm = TRUE)
      ) %>% filter(Sample_Id %in% target_ids)
      
      eval_AD <- inner_join(pred_AD_df, obs_df, by = "Sample_Id")
      
      if (nrow(eval_AD) > 1 && sd(eval_AD$gEBV, na.rm = TRUE) > 0) {
        pa   <- cor(eval_AD$gEBV, eval_AD$Obs, use = "complete.obs")
        bias <- as.numeric(coef(lm(Obs ~ gEBV, data = eval_AD))[2])
        p_c  <- eval_AD$gEBV - mean(eval_AD$gEBV, na.rm = TRUE)
        o_c  <- eval_AD$Obs - mean(eval_AD$Obs, na.rm = TRUE)
        msep <- mean((o_c - p_c)^2, na.rm = TRUE)
        
        set_results[[length(set_results) + 1]] <- data.frame(
          Model = "AD", Set = set_name, Target_Trait = target_trait, 
          N_Validation = nrow(eval_AD), PA = pa, Bias = bias, MSEP = msep
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
