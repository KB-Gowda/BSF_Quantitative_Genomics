# ==============================================================================
# MULTIVARIATE RANDOM FOREST CROSS-VALIDATION
# Models Evaluated: Model A (Additive-Only) vs Model AD (Additive + Dominance)
# Hyperparameter Optimization: Evaluated via Multivariate Out-Of-Bag (OOB) MSEP
# Traits Evaluated: Weight and Surface Area 
# ==============================================================================

# ------------------------------------------------------------------------------
# 1. USER SETTINGS & SETUP ENVIRONMENT
# ------------------------------------------------------------------------------

run_prefix <- "RF_Multivariate_CV"

library(dplyr)
library(Matrix)
library(randomForestSRC) 

# Set global multithreading for rfsrc
n_threads <- 6
options(rf.cores = n_threads)
set.seed(20260825)

# Load data (Xa_final, Xd_final, pheno_final, cv_matrix are already prepared here)
load("BasePop_prep.RData") 

# # Ensure Sample_Id is character and in column 1
pheno_final <- pheno_final %>%
  mutate(Sample_Id = as.character(Sample_Id)) %>%
  relocate(Sample_Id, .before = 1)

stopifnot(all(rownames(cv_matrix) == pheno_final$Sample_Id))

n_reps  <- ncol(cv_matrix)
n_folds <- max(cv_matrix)

# Define multivariate set groupings
multivariate_sets <- list(
  Weight      = c("Adj_SYK_Weight", "Adj_BSG_Weight", "Adj_FVW_Weight"),
  SurfaceArea = c("Adj_SYK_SurfaceArea", "Adj_BSG_SurfaceArea", "Adj_FVW_SurfaceArea")
)

# ------------------------------------------------------------------------------
# 2. PREPARE GENOTYPE MATRICES & HYPERPARAMETER GRIDS
# ------------------------------------------------------------------------------

# Construct Additive + Dominance matrix directly from workspace matrices
Xd_unique <- Xd_final
colnames(Xd_unique) <- paste0(colnames(Xd_unique), "_dom")
M_AD_combined <- cbind(Xa_final, Xd_unique)

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

grid_A  <- get_grid(ncol(Xa_final))
grid_AD <- get_grid(ncol(M_AD_combined))

all_results <- list()

# ------------------------------------------------------------------------------
# 3. NESTED CROSS-VALIDATION LOOP (MULTIVARIATE)
# ------------------------------------------------------------------------------

for (set_name in names(multivariate_sets)) {
  current_traits <- multivariate_sets[[set_name]]
  
  cat("\n========================================================================\n")
  cat("STARTING MULTIVARIATE RF CV FOR SET:", set_name, "\n")
  cat("Traits:", paste(current_traits, collapse = ", "), "\n")
  cat("========================================================================\n")
  
  formula_mv <- as.formula(paste0("Multivar(", paste(current_traits, collapse = ", "), ") ~ ."))
  set_results <- list()
  
  for (r in seq_len(n_reps)) {
    for (f in seq_len(n_folds)) {
      test_idx  <- which(cv_matrix[, r] == f)
      train_idx <- setdiff(seq_len(nrow(pheno_final)), test_idx)
      
      # Keep rows that have AT LEAST ONE phenotypic observation (drop rows where ALL traits are NA)
      has_pheno <- rowSums(!is.na(pheno_final[, current_traits])) > 0
      train_valid_idx <- intersect(train_idx, which(has_pheno))
      test_valid_idx  <- intersect(test_idx, which(has_pheno))
      
      if (length(train_valid_idx) == 0 || length(test_valid_idx) == 0) next
      
      pheno_train <- pheno_final[train_valid_idx, current_traits, drop = FALSE]
      pheno_test  <- pheno_final[test_valid_idx, current_traits, drop = FALSE]
      
      # ========================================================================
      # MODEL A: MULTIVARIATE ADDITIVE-ONLY
      # ========================================================================
      DataTrain_A <- cbind(pheno_train, as.data.frame(Xa_final[train_valid_idx, , drop = FALSE]))
      DataTest_A  <- cbind(pheno_test, as.data.frame(Xa_final[test_valid_idx, , drop = FALSE]))
      
      fit_A_res <- tryCatch({
        oob_errors <- numeric(nrow(grid_A))
        
        for (g in seq_len(nrow(grid_A))) {
          set.seed(20260825 + (r * 1000) + (f * 100) + g)
          rf_tune <- rfsrc(
            formula = formula_mv,
            data = DataTrain_A,
            ntree = grid_A$ntree[g],
            mtry = grid_A$mtry[g],
            nodesize = grid_A$nodesize[g],
            splitrule = "mv",
            na.action = "na.impute" 
          )
          
          msep_sum <- 0
          for (tr in current_traits) {
            y_obs <- DataTrain_A[[tr]]
            y_oob <- rf_tune$regrOutput[[tr]]$predicted.oob
            
            # Calculate OOB error only on valid (non-NA) observations
            valid <- !is.na(y_obs) & !is.na(y_oob)
            if (sum(valid) > 1) {
              msep_sum <- msep_sum + mean((y_obs[valid] - y_oob[valid])^2)
            }
          }
          oob_errors[g] <- msep_sum
        }
        
        best_param_A <- grid_A[which.min(oob_errors), ]
        
        set.seed(20260825 + (r * 100) + f)
        best_fit_A <- rfsrc(
          formula = formula_mv,
          data = DataTrain_A,
          ntree = best_param_A$ntree,
          mtry = best_param_A$mtry,
          nodesize = best_param_A$nodesize,
          splitrule = "mv",
          na.action = "na.impute"
        )
        
        preds_A <- predict(best_fit_A, newdata = DataTest_A, na.action = "na.impute")$regrOutput
        
        res_list <- lapply(current_traits, function(tr) {
          p_val <- preds_A[[tr]]$predicted
          o_val <- DataTest_A[[tr]]
          
          # Safely calculate metrics, ignoring the NAs for non-target diets
          valid_mask <- !is.na(o_val) & !is.na(p_val)
          if(sum(valid_mask) > 1) {
            pa   <- cor(p_val[valid_mask], o_val[valid_mask])
            bias <- coef(lm(o_val[valid_mask] ~ p_val[valid_mask]))[2]
            msep <- mean((o_val[valid_mask] - p_val[valid_mask])^2)
            n_val <- sum(valid_mask)
          } else {
            pa <- NA; bias <- NA; msep <- NA; n_val <- 0
          }
          
          data.frame(
            Model = "A", Set = set_name, Trait = tr, Rep = r, Fold = f,
            N_validation = n_val, PA = pa, Bias = bias, MSEP = msep,
            Best_ntree = best_param_A$ntree, Best_mtry = best_param_A$mtry, Best_nodesize = best_param_A$nodesize
          )
        })
        bind_rows(res_list)
      }, error = function(e) {
        cat("\n[Error in Model A]:", e$message, "\n")
        return(NULL)
      })
      
      # ========================================================================
      # MODEL AD: MULTIVARIATE ADDITIVE + DOMINANCE
      # ========================================================================
      DataTrain_AD <- cbind(pheno_train, as.data.frame(M_AD_combined[train_valid_idx, , drop = FALSE]))
      DataTest_AD  <- cbind(pheno_test, as.data.frame(M_AD_combined[test_valid_idx, , drop = FALSE]))
      
      fit_AD_res <- tryCatch({
        oob_errors <- numeric(nrow(grid_AD))
        
        for (g in seq_len(nrow(grid_AD))) {
          set.seed(20260825 + (r * 1000) + (f * 100) + g)
          rf_tune <- rfsrc(
            formula = formula_mv,
            data = DataTrain_AD,
            ntree = grid_AD$ntree[g],
            mtry = grid_AD$mtry[g],
            nodesize = grid_AD$nodesize[g],
            splitrule = "mv",
            na.action = "na.impute"
          )
          
          msep_sum <- 0
          for (tr in current_traits) {
            y_obs <- DataTrain_AD[[tr]]
            y_oob <- rf_tune$regrOutput[[tr]]$predicted.oob
            valid <- !is.na(y_obs) & !is.na(y_oob)
            if (sum(valid) > 1) {
              msep_sum <- msep_sum + mean((y_obs[valid] - y_oob[valid])^2)
            }
          }
          oob_errors[g] <- msep_sum
        }
        
        best_param_AD <- grid_AD[which.min(oob_errors), ]
        
        set.seed(20260825 + (r * 100) + f)
        best_fit_AD <- rfsrc(
          formula = formula_mv,
          data = DataTrain_AD,
          ntree = best_param_AD$ntree,
          mtry = best_param_AD$mtry,
          nodesize = best_param_AD$nodesize,
          splitrule = "mv",
          na.action = "na.impute"
        )
        
        preds_AD <- predict(best_fit_AD, newdata = DataTest_AD, na.action = "na.impute")$regrOutput
        
        res_list <- lapply(current_traits, function(tr) {
          p_val <- preds_AD[[tr]]$predicted
          o_val <- DataTest_AD[[tr]]
          
          valid_mask <- !is.na(o_val) & !is.na(p_val)
          if(sum(valid_mask) > 1) {
            pa   <- cor(p_val[valid_mask], o_val[valid_mask])
            bias <- coef(lm(o_val[valid_mask] ~ p_val[valid_mask]))[2]
            msep <- mean((o_val[valid_mask] - p_val[valid_mask])^2)
            n_val <- sum(valid_mask)
          } else {
            pa <- NA; bias <- NA; msep <- NA; n_val <- 0
          }
          
          data.frame(
            Model = "AD", Set = set_name, Trait = tr, Rep = r, Fold = f,
            N_validation = n_val, PA = pa, Bias = bias, MSEP = msep,
            Best_ntree = best_param_AD$ntree, Best_mtry = best_param_AD$mtry, Best_nodesize = best_param_AD$nodesize
          )
        })
        bind_rows(res_list)
      }, error = function(e) {
        cat("\n[Error in Model AD]:", e$message, "\n")
        return(NULL)
      })
      
      # Print results dynamically to console (filtering out rows with 0 validations)
      if (!is.null(fit_A_res)) {
        res_valid <- fit_A_res[fit_A_res$N_validation > 0, ]
        set_results[[length(set_results) + 1]] <- res_valid
        for (i in seq_len(nrow(res_valid))) {
          cat(sprintf("Rep %2d | Fold %d | [A]  Trait: %-19s | PA: %6.4f | Bias: %6.4f | MSEP: %6.4f\n",
                      r, f, res_valid$Trait[i], res_valid$PA[i], res_valid$Bias[i], res_valid$MSEP[i]))
        }
      }
      
      if (!is.null(fit_AD_res)) {
        res_valid <- fit_AD_res[fit_AD_res$N_validation > 0, ]
        set_results[[length(set_results) + 1]] <- res_valid
        for (i in seq_len(nrow(res_valid))) {
          cat(sprintf("Rep %2d | Fold %d | [AD] Trait: %-19s | PA: %6.4f | Bias: %6.4f | MSEP: %6.4f\n",
                      r, f, res_valid$Trait[i], res_valid$PA[i], res_valid$Bias[i], res_valid$MSEP[i]))
        }
      }
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

if (length(all_results) == 0) stop("No results were generated. Check error prints above.")

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
cat("RANDOM FOREST MULTIVARIATE SUMMARY\n")
cat("======================================================================\n")
print(as.data.frame(summary_results))
