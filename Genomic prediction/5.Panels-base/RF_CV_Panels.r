# ==============================================================================
# RANDOM FOREST NESTED CROSS-VALIDATION FOR SNP PANELS (RANGER)
# Models Evaluated: Model A (Additive-Only) vs Model AD (Additive + Dominance)
# Hyperparameter Optimization: Evaluated via Out-Of-Bag (OOB) MSE per Panel
# Parallelism: Multi-threading handled directly within ranger (num.threads)
# ==============================================================================

# ------------------------------------------------------------------------------
# 1. USER SETTINGS & COMPUTATIONAL PARALLELISM
# ------------------------------------------------------------------------------

# Set run prefix
run_prefix    <- "RF_NestedCV_Panel"

# Define target traits to evaluate (pre-adjusted phenotypes in pheno_final)
target_traits <- c("Adj_Weight", "Adj_SurfaceArea")

# Number of threads assigned directly to ranger
n_threads     <- 4

# Master seed
master_seed   <- 20260825

# ------------------------------------------------------------------------------
# 2. LOAD DATA AND PACKAGES
# ------------------------------------------------------------------------------

library(dplyr)
library(Matrix)
library(ranger)

# Enable unbuffered output so logs print live in real-time on HPC/SLURM
options(sinksToConsole = TRUE)
sink(stdout(), type = "output")

# Load prepared panel input containing pheno_final, panel_mats, and cv_matrix
load("GP_panel_input.RData")

# Ensure Sample_Id is character and in column 1
pheno_final <- pheno_final %>%
  mutate(Sample_Id = as.character(Sample_Id)) %>%
  relocate(Sample_Id, .before = 1)

# Verify cv_matrix alignment with pheno_final
stopifnot(all(rownames(cv_matrix) == pheno_final$Sample_Id))

n_reps  <- ncol(cv_matrix)
n_folds <- max(cv_matrix)

# Map trait column names to top-panel filtering codes
panel_trait_map <- c(
  "top1000_Weight_A"       = "Adj_Weight",
  "top1000_Weight_AD"      = "Adj_Weight",
  "top1000_SurfaceArea_A"  = "Adj_SurfaceArea",
  "top1000_SurfaceArea_AD" = "Adj_SurfaceArea"
)

panels <- c("4000", "3000", "2000", "1000", "500",
            "top1000_Weight_A", "top1000_Weight_AD",
            "top1000_SurfaceArea_A", "top1000_SurfaceArea_AD")

all_results <- list()

# Helper function to construct orthogonal dominance genotype matrix (Vitezica coding)
build_dominance_matrix <- function(Xa) {
  p <- colMeans(Xa, na.rm = TRUE) / 2
  q <- 1 - p
  Xd <- Xa
  for (j in seq_len(ncol(Xa))) {
    p_j <- p[j]
    q_j <- q[j]
    geno <- Xa[, j]
    
    Xd[geno == 0, j] <- -2 * (p_j^2)
    Xd[geno == 1, j] <- 2 * p_j * q_j
    Xd[geno == 2, j] <- -2 * (q_j^2)
  }
  return(Xd)
}

# Function to generate hyperparameter grid based on dynamic predictor matrix dimensions
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

# ------------------------------------------------------------------------------
# 3. NESTED CROSS-VALIDATION LOOP ACROSS TRAITS AND PANELS
# ------------------------------------------------------------------------------

for (trait in target_traits) {
  cat("\n", paste(rep("=", 80), collapse = ""), "\n", sep = "")
  cat("STARTING RANDOM FOREST CV FOR TRAIT:", trait, "\n")
  cat(paste(rep("=", 80), collapse = ""), "\n")
  
  y_obs <- pheno_final[[trait]]
  
  for (panel in panels) {
    
    # Skip trait-specific top-1000 panels that do not match current trait
    if (panel %in% names(panel_trait_map) && panel_trait_map[panel] != trait) next
    
    if (!panel %in% names(panel_mats)) {
      cat(sprintf("\n[SKIPPED Panel %s | Trait %s]: Matrix not found in panel_mats\n", panel, trait))
      next
    }
    
    # Extract genotype matrix for current panel from panel_mats
    p_mats <- panel_mats[[panel]]
    
    if (!"Xa" %in% names(p_mats)) {
      cat(sprintf("\n[SKIPPED Panel %s | Trait %s]: Genotype matrix 'Xa' not found in panel_mats[[panel]]\n", panel, trait))
      next
    }
    
    Xa <- p_mats$Xa
    has_AD <- !is.null(p_mats$Xd) || !is.null(p_mats$Gd.inv)
    
    # Generate dynamic hyperparameter grid for Model A
    grid_A <- get_grid(ncol(Xa))
    
    # Retrieve or compute orthogonal dominance matrix if Model AD is supported
    if (has_AD) {
      if ("Xd" %in% names(p_mats)) {
        Xd <- p_mats$Xd
      } else {
        Xd <- build_dominance_matrix(Xa)
      }
      
      Xd_unique <- Xd
      colnames(Xd_unique) <- paste0(colnames(Xd_unique), "_dom")
      M_AD_combined <- cbind(Xa, Xd_unique)
      
      # Generate dynamic hyperparameter grid for Model AD
      grid_AD <- get_grid(ncol(M_AD_combined))
    }
    
    panel_results <- list()
    
    for (r in seq_len(n_reps)) {
      cat(sprintf("\n--- Trait: %s | Panel: %s | Repetition %d/%d ---\n", trait, panel, r, n_reps))
      
      for (f in seq_len(n_folds)) {
        test_idx  <- which(cv_matrix[, r] == f)
        train_idx <- setdiff(seq_len(nrow(pheno_final)), test_idx)
        
        # Filter non-missing phenotype indices
        train_valid_idx <- train_idx[!is.na(y_obs[train_idx])]
        test_valid_idx  <- test_idx[!is.na(y_obs[test_idx])]
        
        if (length(train_valid_idx) == 0 || length(test_valid_idx) == 0) next
        
        y_train <- y_obs[train_valid_idx]
        y_test  <- y_obs[test_valid_idx]
        
        # Calculate a deterministic fold-specific seed
        fold_seed <- master_seed + (r * 1000) + f
        
        # ======================================================================
        # MODEL A: ADDITIVE-ONLY
        # ======================================================================
        X_train_A <- Xa[train_valid_idx, , drop = FALSE]
        X_test_A  <- Xa[test_valid_idx, , drop = FALSE]
        
        fit_A_res <- tryCatch({
          oob_errors <- numeric(nrow(grid_A))
          
          # Grid search using ranger multi-threading
          for (g in seq_len(nrow(grid_A))) {
            set.seed(fold_seed + g)
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
          set.seed(fold_seed)
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
          cat(sprintf("  [ERROR Model A | Panel %s | Rep %d | Fold %d]: %s\n", panel, r, f, e$message))
          return(NULL)
        })
        
        if (!is.null(fit_A_res)) {
          panel_results[[length(panel_results) + 1]] <- data.frame(
            Panel = panel, Model = "A", Trait = trait, Rep = r, Fold = f,
            N_validation = length(y_test),
            PA = fit_A_res$PA, Bias = fit_A_res$Bias, MSEP = fit_A_res$MSEP,
            Best_num_trees = fit_A_res$Best_Param$num.trees,
            Best_mtry = fit_A_res$Best_Param$mtry,
            Best_min_node_size = fit_A_res$Best_Param$min.node.size
          )
          
          cat(sprintf("  Rep %2d | Fold %d | Panel %-22s | Model A  -> PA: %.4f | Bias: %.4f | MSEP: %.4f\n",
                      r, f, panel, fit_A_res$PA, fit_A_res$Bias, fit_A_res$MSEP))
        }
        
        # ======================================================================
        # MODEL AD: ADDITIVE + DOMINANCE
        # ======================================================================
        if (has_AD) {
          X_train_AD <- M_AD_combined[train_valid_idx, , drop = FALSE]
          X_test_AD  <- M_AD_combined[test_valid_idx, , drop = FALSE]
          
          fit_AD_res <- tryCatch({
            oob_errors <- numeric(nrow(grid_AD))
            
            # Grid search using ranger multi-threading
            for (g in seq_len(nrow(grid_AD))) {
              set.seed(fold_seed + 50000 + g)
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
            set.seed(fold_seed + 50000)
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
            cat(sprintf("  [ERROR Model AD | Panel %s | Rep %d | Fold %d]: %s\n", panel, r, f, e$message))
            return(NULL)
          })
          
          if (!is.null(fit_AD_res)) {
            panel_results[[length(panel_results) + 1]] <- data.frame(
              Panel = panel, Model = "AD", Trait = trait, Rep = r, Fold = f,
              N_validation = length(y_test),
              PA = fit_AD_res$PA, Bias = fit_AD_res$Bias, MSEP = fit_AD_res$MSEP,
              Best_num_trees = fit_AD_res$Best_Param$num.trees,
              Best_mtry = fit_AD_res$Best_Param$mtry,
              Best_min_node_size = fit_AD_res$Best_Param$min.node.size
            )
            
            cat(sprintf("  Rep %2d | Fold %d | Panel %-22s | Model AD -> PA: %.4f | Bias: %.4f | MSEP: %.4f\n",
                        r, f, panel, fit_AD_res$PA, fit_AD_res$Bias, fit_AD_res$MSEP))
          }
        }
        
        flush.console()
      }
    }
    
    if (length(panel_results) > 0) {
      all_results[[paste(trait, panel, sep = "_")]] <- bind_rows(panel_results)
    }
  }
}

# ------------------------------------------------------------------------------
# 4. COMBINE & SUMMARIZE RESULTS BY TRAIT
# ------------------------------------------------------------------------------

if (length(all_results) == 0) {
  stop("No cross-validation results were generated. Check for model fit errors printed above.")
}

raw_results <- bind_rows(all_results)

summary_results <- raw_results %>%
  group_by(Trait, Panel, Model) %>%
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
# 5. SAVE SEPARATE TRAIT OUTPUTS & MARKDOWN REPORT
# ------------------------------------------------------------------------------

for (tr in unique(raw_results$Trait)) {
  raw_tr     <- raw_results %>% filter(Trait == tr)
  summary_tr <- summary_results %>% filter(Trait == tr)
  
  write.csv(raw_tr, paste0(run_prefix, "_", tr, "_Raw_Results.csv"), row.names = FALSE)
  write.csv(summary_tr, paste0(run_prefix, "_", tr, "_Summary_Results.csv"), row.names = FALSE)
}

# Generate Markdown Report
md_file <- paste0(run_prefix, "_Summary_Report.md")

md_content <- c(
  paste0("# ", gsub("_", " ", run_prefix), " Summary Report"),
  "",
  paste0("**Date Generated:** ", Sys.time()),
  paste0("**Replications:** ", n_reps),
  paste0("**Folds:** ", n_folds),
  paste0("**Threads Used:** ", n_threads),
  paste0("**Total Samples Evaluated:** ", nrow(pheno_final)),
  ""
)

for (tr in unique(summary_results$Trait)) {
  md_content <- c(
    md_content,
    paste0("## Trait: ", tr),
    "",
    "| Panel | Model | N | Mean PA | SD PA | SE PA | Mean Bias | SD Bias | SE Bias | Mean MSEP | SD MSEP | SE MSEP |",
    "|:---|:---|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---:|"
  )
  
  summary_tr <- summary_results %>% filter(Trait == tr)
  for (i in seq_len(nrow(summary_tr))) {
    row <- summary_tr[i, ]
    md_line <- sprintf(
      "| %s | %s | %d | %.4f | %.4f | %.4f | %.4f | %.4f | %.4f | %.4f | %.4f | %.4f |",
      row$Panel, row$Model, row$N_Evaluations,
      row$Mean_PA, row$SD_PA, row$SE_PA,
      row$Mean_Bias, row$SD_Bias, row$SE_Bias,
      row$Mean_MSEP, row$SD_MSEP, row$SE_MSEP
    )
    md_content <- c(md_content, md_line)
  }
  md_content <- c(md_content, "")
}

md_content <- c(
  md_content,
  "---",
  "### Metric Definitions:",
  "- **PA (Predictive Ability):** Correlation between observed adjusted phenotypes and predicted values.",
  "- **Bias:** Regression slope of observed adjusted phenotypes on predicted values (target = 1.0).",
  "- **MSEP:** Mean Squared Error of Prediction."
)

writeLines(md_content, con = md_file)
cat("\nMarkdown report successfully generated:", md_file, "\n")

# Print overall summary to console
cat("\n======================================================================\n")
cat(toupper(gsub("_", " ", run_prefix)), " SUMMARY\n")
cat("======================================================================\n")
print(summary_results)
