# ==============================================================================
# BAYES A CROSS-VALIDATION FOR SNP PANELS
# Models Evaluated: Model A (Additive-Only) vs Model AD (Additive + Dominance)
# ==============================================================================

# ------------------------------------------------------------------------------
# 1. USER SETTINGS, TARGET TRAITS & RUN IDENTIFIER
# ------------------------------------------------------------------------------

# Set run prefix
run_prefix <- "BayesA_CV_Panel"

# Define target traits to evaluate (pre-adjusted phenotypes in pheno_final)
target_traits <- c("Adj_Weight", "Adj_SurfaceArea")

# Define Bayesian model iterations
n_iter  <- 12000
n_burn  <- 2000

# ------------------------------------------------------------------------------
# 2. LOAD DATA AND PACKAGES
# ------------------------------------------------------------------------------

library(dplyr)
library(Matrix)
library(BGLR)

# Enable unbuffered output so logs print live in real-time on HPC/SLURM
options(sinksToConsole = TRUE)
sink(stdout(), type = "output")

# Load prepared panel input containing pheno_final, panel_mats, and cv_matrix
load("GP_panel_input.RData")

# Verify cv_matrix alignment with pheno_final
stopifnot(all(rownames(cv_matrix) == as.character(pheno_final$Sample_Id)))

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

# Helper function to construct orthogonal dominance genotype matrix
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

# ------------------------------------------------------------------------------
# 3. CROSS-VALIDATION LOOP
# ------------------------------------------------------------------------------

for (trait in target_traits) {
  cat("\n", paste(rep("=", 80), collapse = ""), "\n", sep = "")
  cat("STARTING CROSS-VALIDATION FOR TRAIT:", trait, "\n")
  cat(paste(rep("=", 80), collapse = ""), "\n")
  
  y_obs <- pheno_final[[trait]]
  
  for (panel in panels) {
    
    # Skip trait-specific top-1000 panels that do not match the current trait
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
    
    # Retrieve or compute orthogonal dominance matrix if Model AD is supported
    if (has_AD) {
      if ("Xd" %in% names(p_mats)) {
        Xd <- p_mats$Xd
      } else {
        Xd <- build_dominance_matrix(Xa)
      }
    }
    
    panel_results <- list()
    
    for (r in seq_len(n_reps)) {
      cat(sprintf("\n--- Trait: %s | Panel: %s | Repetition %d/%d ---\n", trait, panel, r, n_reps))
      
      for (f in seq_len(n_folds)) {
        test_idx  <- which(cv_matrix[, r] == f)
        train_idx <- setdiff(seq_len(nrow(pheno_final)), test_idx)
        
        # Mask validation set phenotypes
        y_train <- y_obs
        y_train[test_idx] <- NA
        
        # ======================================================================
        # MODEL A: ADDITIVE-ONLY BAYES A
        # ======================================================================
        ETA_A <- list(
          Additive = list(X = Xa, model = "BayesA")
        )
        
        mod_A <- tryCatch({
          BGLR(
            y = y_train,
            ETA = ETA_A,
            nIter = n_iter,
            burnIn = n_burn,
            verbose = FALSE
          )
        }, error = function(e) {
          cat(sprintf("  [ERROR Model A | Panel %s | Rep %d | Fold %d]: %s\n", panel, r, f, e$message))
          return(NULL)
        })
        
        if (!is.null(mod_A)) {
          p_A   <- (Xa %*% mod_A$ETA[["Additive"]]$b)[test_idx]
          o_val <- y_obs[test_idx]
          
          v_idx_A <- which(!is.na(o_val) & !is.na(p_A))
          if (length(v_idx_A) > 1) {
            p_sub <- p_A[v_idx_A]
            o_sub <- o_val[v_idx_A]
            
            pa_A   <- cor(p_sub, o_sub, use = "complete.obs")
            bias_A <- coef(lm(o_sub ~ p_sub))[2]
            msep_A <- mean((o_sub - p_sub)^2, na.rm = TRUE)
            
            res_A <- data.frame(
              Panel = panel, Model = "A", Trait = trait, Rep = r, Fold = f,
              N_validation = length(v_idx_A),
              PA = pa_A, Bias = bias_A, MSEP = msep_A
            )
            panel_results[[length(panel_results) + 1]] <- res_A
            
            cat(sprintf("  Rep %2d | Fold %d | Panel %-22s | Model A  -> PA: %.4f | Bias: %.4f | MSEP: %.4f\n",
                        r, f, panel, pa_A, bias_A, msep_A))
          }
          rm(mod_A); gc(verbose = FALSE)
        }
        
        # ======================================================================
        # MODEL AD: ADDITIVE + DOMINANCE BAYES A
        # ======================================================================
        if (has_AD) {
          ETA_AD <- list(
            Additive  = list(X = Xa, model = "BayesA"),
            Dominance = list(X = Xd, model = "BayesA")
          )
          
          mod_AD <- tryCatch({
            BGLR(
              y = y_train,
              ETA = ETA_AD,
              nIter = n_iter,
              burnIn = n_burn,
              verbose = FALSE
            )
          }, error = function(e) {
            cat(sprintf("  [ERROR Model AD | Panel %s | Rep %d | Fold %d]: %s\n", panel, r, f, e$message))
            return(NULL)
          })
          
          if (!is.null(mod_AD)) {
            gebv_all <- Xa %*% mod_AD$ETA[["Additive"]]$b
            gdbv_all <- Xd %*% mod_AD$ETA[["Dominance"]]$b
            gpv_all  <- gebv_all + gdbv_all
            
            p_AD  <- gpv_all[test_idx]
            o_val <- y_obs[test_idx]
            
            v_idx_AD <- which(!is.na(o_val) & !is.na(p_AD))
            if (length(v_idx_AD) > 1) {
              p_sub <- p_AD[v_idx_AD]
              o_sub <- o_val[v_idx_AD]
              
              pa_AD   <- cor(p_sub, o_sub, use = "complete.obs")
              bias_AD <- coef(lm(o_sub ~ p_sub))[2]
              msep_AD <- mean((o_sub - p_sub)^2, na.rm = TRUE)
              
              res_AD <- data.frame(
                Panel = panel, Model = "AD", Trait = trait, Rep = r, Fold = f,
                N_validation = length(v_idx_AD),
                PA = pa_AD, Bias = bias_AD, MSEP = msep_AD
              )
              panel_results[[length(panel_results) + 1]] <- res_AD
              
              cat(sprintf("  Rep %2d | Fold %d | Panel %-22s | Model AD -> PA: %.4f | Bias: %.4f | MSEP: %.4f\n",
                          r, f, panel, pa_AD, bias_AD, msep_AD))
            }
            rm(mod_AD); gc(verbose = FALSE)
          }
        }
        
        # Flush log output after each fold
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

# Save combined outputs separately for each trait
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
  paste0("**Bayes A MCMC Iterations:** ", n_iter, " (Burn-in: ", n_burn, ")"),
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
