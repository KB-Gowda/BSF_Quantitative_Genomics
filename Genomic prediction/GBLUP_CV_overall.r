# ==============================================================================
# GBLUP CROSS-VALIDATION FOR BASE POPULATION TRAITS
# Models Evaluated: Model A (Additive-Only) vs Model AD (Additive + Dominance)
# Internal Multithreading Enabled via asreml.options(nthreads = 4)
# ==============================================================================

# ------------------------------------------------------------------------------
# 1. USER SETTINGS, TARGET TRAITS & RUN IDENTIFIER
# ------------------------------------------------------------------------------

# Set run prefix 
run_prefix <- "GBLUP_CV_Overall"

# Define target traits to evaluate
target_traits <- c("Adj_Weight", "Adj_Length", "Adj_Width", "Adj_SurfaceArea")

# Number of CPU threads for ASReml internal calculations
n_threads <- 4

# ------------------------------------------------------------------------------
# 2. LOAD DATA AND PACKAGES
# ------------------------------------------------------------------------------

library(asreml)
library(dplyr)
library(Matrix)

# Enable unbuffered output so logs print live in real-time on HPC/SLURM
options(sinksToConsole = TRUE)
sink(stdout(), type = "output")

# Set ASReml to use 4 CPU threads globally for internal matrix operations
asreml.options(nthreads = n_threads)

# Load prepared input containing pheno_final, Ginv, Gd.inv, and cv_matrix
load("BasePop_prep.RData")

workspace_size <- "8GB"

# Verify cv_matrix alignment with pheno_final
stopifnot(all(rownames(cv_matrix) == as.character(pheno_final$Sample_Id)))

n_reps  <- ncol(cv_matrix)
n_folds <- max(cv_matrix)

all_results <- list()

# ------------------------------------------------------------------------------
# 3. CROSS-VALIDATION LOOP
# ------------------------------------------------------------------------------

for (trait in target_traits) {
  cat("\n", paste(rep("=", 80), collapse = ""), "\n", sep = "")
  cat("STARTING CROSS-VALIDATION FOR TRAIT:", trait, "\n")
  cat(paste(rep("=", 80), collapse = ""), "\n")
  
  y_obs <- pheno_final[[trait]]
  trait_results <- list()
  
  for (r in seq_len(n_reps)) {
    cat(sprintf("\n--- Trait: %s | Repetition %d/%d ---\n", trait, r, n_reps))
    
    for (f in seq_len(n_folds)) {
      test_idx  <- which(cv_matrix[, r] == f)
      train_idx <- setdiff(seq_len(nrow(pheno_final)), test_idx)
      
      # Mask validation set
      temp_data <- pheno_final
      temp_data[[trait]][test_idx] <- NA
      
      # ========================================================================
      # MODEL A: ADDITIVE-ONLY GBLUP
      # ========================================================================
      mod_A <- tryCatch({
        fit <- asreml(
          fixed = as.formula(paste(trait, "~ 1")),
          random = ~ vm(Sample_Id, Ginv),
          residual = ~ idv(units),
          data = temp_data,
          maxit = 50,
          workspace = workspace_size,
          trace = FALSE
        )
        update.asreml(fit, trace = FALSE)
      }, error = function(e) {
        cat(sprintf("  [ERROR Model A | Rep %d | Fold %d]: %s\n", r, f, e$message))
        return(NULL)
      })
      
      if (!is.null(mod_A) && mod_A$converge) {
        blups_A <- as.data.frame(summary(mod_A, coef = TRUE)$coef.random)
        add_rows <- grep("vm\\(Sample_Id, Ginv\\)", rownames(blups_A))
        
        if (length(add_rows) > 0) {
          add_A <- blups_A[add_rows, , drop = FALSE]
          add_A$ID <- gsub(".*vm\\(Sample_Id, Ginv\\)_", "", rownames(add_A))
          
          m_idx_A <- match(as.character(pheno_final$Sample_Id[test_idx]), add_A$ID)
          p_A     <- add_A$solution[m_idx_A]
          o_val   <- y_obs[test_idx]
          
          v_idx_A <- which(!is.na(o_val) & !is.na(p_A))
          if (length(v_idx_A) > 1) {
            p_sub <- p_A[v_idx_A]
            o_sub <- o_val[v_idx_A]
            
            pa_A   <- cor(p_sub, o_sub, use = "complete.obs")
            bias_A <- coef(lm(o_sub ~ p_sub))[2]
            msep_A <- mean((o_sub - p_sub)^2, na.rm = TRUE)
            
            res_A <- data.frame(
              Model = "A", Trait = trait, Rep = r, Fold = f,
              N_validation = length(v_idx_A),
              PA = pa_A, Bias = bias_A, MSEP = msep_A
            )
            trait_results[[length(trait_results) + 1]] <- res_A
            
            cat(sprintf("  Rep %2d | Fold %d | Model A  -> PA: %.4f | Bias: %.4f | MSEP: %.4f\n",
                        r, f, pa_A, bias_A, msep_A))
          }
        }
      }
      
      # ========================================================================
      # MODEL AD: ADDITIVE + DOMINANCE GBLUP
      # ========================================================================
      mod_AD <- tryCatch({
        fit <- asreml(
          fixed = as.formula(paste(trait, "~ 1")),
          random = ~ vm(Sample_Id, Ginv) + vm(Sample_Id, Gd.inv),
          residual = ~ idv(units),
          data = temp_data,
          maxit = 50,
          workspace = workspace_size,
          trace = FALSE
        )
        update.asreml(fit, trace = FALSE)
      }, error = function(e) {
        cat(sprintf("  [ERROR Model AD | Rep %d | Fold %d]: %s\n", r, f, e$message))
        return(NULL)
      })
      
      if (!is.null(mod_AD) && mod_AD$converge) {
        blups_AD <- as.data.frame(summary(mod_AD, coef = TRUE)$coef.random)
        
        add_rows <- grep("vm\\(Sample_Id, Ginv\\)", rownames(blups_AD))
        dom_rows <- grep("vm\\(Sample_Id, Gd.inv\\)", rownames(blups_AD))
        
        if (length(add_rows) > 0 && length(dom_rows) > 0) {
          add_AD <- blups_AD[add_rows, , drop = FALSE]
          dom_AD <- blups_AD[dom_rows, , drop = FALSE]
          
          add_AD$ID <- gsub(".*vm\\(Sample_Id, Ginv\\)_", "", rownames(add_AD))
          dom_AD$ID <- gsub(".*vm\\(Sample_Id, Gd.inv\\)_", "", rownames(dom_AD))
          
          m_idx_add <- match(as.character(pheno_final$Sample_Id[test_idx]), add_AD$ID)
          m_idx_dom <- match(as.character(pheno_final$Sample_Id[test_idx]), dom_AD$ID)
          
          p_AD  <- add_AD$solution[m_idx_add] + dom_AD$solution[m_idx_dom]
          o_val <- y_obs[test_idx]
          
          v_idx_AD <- which(!is.na(o_val) & !is.na(p_AD))
          if (length(v_idx_AD) > 1) {
            p_sub <- p_AD[v_idx_AD]
            o_sub <- o_val[v_idx_AD]
            
            pa_AD   <- cor(p_sub, o_sub, use = "complete.obs")
            bias_AD <- coef(lm(o_sub ~ p_sub))[2]
            msep_AD <- mean((o_sub - p_sub)^2, na.rm = TRUE)
            
            res_AD <- data.frame(
              Model = "AD", Trait = trait, Rep = r, Fold = f,
              N_validation = length(v_idx_AD),
              PA = pa_AD, Bias = bias_AD, MSEP = msep_AD
            )
            trait_results[[length(trait_results) + 1]] <- res_AD
            
            cat(sprintf("  Rep %2d | Fold %d | Model AD -> PA: %.4f | Bias: %.4f | MSEP: %.4f\n",
                        r, f, pa_AD, bias_AD, msep_AD))
          }
        }
      }
      
      # Force flushing output to console/SLURM log immediately after each fold
      flush.console()
    }
  }
  
  if (length(trait_results) > 0) {
    df_trait <- bind_rows(trait_results)
    all_results[[trait]] <- df_trait
    
    # Save individual trait raw output
    write.csv(df_trait, paste0(run_prefix, "_", trait, "_Raw.csv"), row.names = FALSE)
  }
}

# ------------------------------------------------------------------------------
# 4. COMBINE & SUMMARIZE RESULTS
# ------------------------------------------------------------------------------

if (length(all_results) == 0) {
  stop("No cross-validation results were generated. Check for model convergence errors printed above.")
}

raw_results <- bind_rows(all_results)

summary_results <- raw_results %>%
  group_by(Model, Trait) %>%
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

# Generate Markdown Report
md_file <- paste0(run_prefix, "_Summary_Report.md")

md_content <- c(
  paste0("# ", gsub("_", " ", run_prefix), " Summary Report"),
  "",
  paste0("**Date Generated:** ", Sys.time()),
  paste0("**Replications:** ", n_reps),
  paste0("**Folds:** ", n_folds),
  paste0("**ASReml Threads:** ", n_threads),
  paste0("**Total Samples Evaluated:** ", nrow(pheno_final)),
  "",
  "## Overall Model Comparison Table",
  "",
  "| Model | Trait | N | Mean PA | SD PA | SE PA | Mean Bias | SD Bias | SE Bias | Mean MSEP | SD MSEP | SE MSEP |",
  "|:---|:---|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---:|"
)

for (i in seq_len(nrow(summary_results))) {
  row <- summary_results[i, ]
  md_line <- sprintf(
    "| %s | %s | %d | %.4f | %.4f | %.4f | %.4f | %.4f | %.4f | %.4f | %.4f | %.4f |",
    row$Model, row$Trait, row$N_Evaluations,
    row$Mean_PA, row$SD_PA, row$SE_PA,
    row$Mean_Bias, row$SD_Bias, row$SE_Bias,
    row$Mean_MSEP, row$SD_MSEP, row$SE_MSEP
  )
  md_content <- c(md_content, md_line)
}

md_content <- c(
  md_content,
  "",
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
