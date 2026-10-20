# ==============================================================================
# MULTIVARIATE GBLUP CROSS-VALIDATION FOR WEIGHT AND SURFACE AREA (GxD CONTEXT)
# Models Evaluated: Model A (Additive-Only) vs Model AD (Additive + Dominance)
# Target Metrics: PA, Bias, MSEP (Using Pre-Adjusted Phenotypes)
# ==============================================================================

# ------------------------------------------------------------------------------
# 1. USER SETTINGS & TARGET MULTIVARIATE SETS
# ------------------------------------------------------------------------------

# Set run prefix
run_prefix <- "GBLUP_Multivariate_CV"

# Define multivariate set groupings using pre-adjusted phenotypes across environments
multivariate_sets <- list(
  Weight      = c("Adj_SYK_Weight", "Adj_BSG_Weight", "Adj_FVW_Weight"),
  SurfaceArea = c("Adj_SYK_SurfaceArea", "Adj_BSG_SurfaceArea", "Adj_FVW_SurfaceArea")
)

# Number of CPU threads for ASReml internal calculations
n_threads <- 4

# ------------------------------------------------------------------------------
# 2. LOAD DATA AND PACKAGES
# ------------------------------------------------------------------------------

library(asreml)
library(dplyr)
library(Matrix)

# Enable unbuffered console logging for HPC/SLURM environments
options(sinksToConsole = TRUE)
sink(stdout(), type = "output")

# Set global ASReml multithreading
asreml.options(nthreads = n_threads)

# Load prepared input workspace containing pheno_final, Ginv, Gd.inv, and cv_matrix
load("BasePop_prep.RData")

workspace_size <- "8GB"

# Verify cv_matrix alignment with pheno_final
stopifnot(all(rownames(cv_matrix) == as.character(pheno_final$Sample_Id)))

n_reps  <- ncol(cv_matrix)
n_folds <- max(cv_matrix)

all_results <- list()

# ------------------------------------------------------------------------------
# 3. CROSS-VALIDATION LOOP (MULTIVARIATE GxD SETS)
# ------------------------------------------------------------------------------

for (set_name in names(multivariate_sets)) {
  current_traits <- multivariate_sets[[set_name]]

  cat("\n", paste(rep("=", 80), collapse = ""), "\n", sep = "")
  cat("STARTING MULTIVARIATE GxD CROSS-VALIDATION FOR SET:", set_name, "\n")
  cat("Traits in set:", paste(current_traits, collapse = ", "), "\n")
  cat(paste(rep("=", 80), collapse = ""), "\n")

  # Formula for pre-adjusted phenotypes across environments: cbind(Adj_SYK_..., Adj_BSG_..., Adj_FVW_...) ~ trait
  cbind_str    <- paste0("cbind(", paste(current_traits, collapse = ", "), ")")
  multiv_fixed <- as.formula(paste(cbind_str, "~ trait"))

  set_results <- list()

  for (r in seq_len(n_reps)) {
    cat(sprintf("\n--- Set: %s | Repetition %d/%d ---\n", set_name, r, n_reps))

    for (f in seq_len(n_folds)) {
      test_idx <- which(cv_matrix[, r] == f)
      if (length(test_idx) == 0) next

      fold_ids <- as.character(pheno_final$Sample_Id[test_idx])

      # Mask pre-adjusted phenotypes across all 3 environment traits for validation individuals
      temp_data <- pheno_final
      for (tr in current_traits) {
        temp_data[[tr]][test_idx] <- NA
      }

      # ========================================================================
      # MODEL A: MULTIVARIATE ADDITIVE-ONLY GxD GBLUP
      # ========================================================================
      mod_A <- tryCatch({
        fit <- asreml(
          fixed    = multiv_fixed,
          random   = ~ us(trait):vm(Sample_Id, Ginv),
          residual = ~ id(units):us(trait),
          data     = temp_data,
          maxit    = 100,
          ai.sing  = TRUE,
          workspace= workspace_size,
          trace    = FALSE
        )
        update.asreml(fit, trace = FALSE)
      }, error = function(e) {
        cat(sprintf("  [ERROR Model A | Set %s | Rep %d | Fold %d]: %s\n", set_name, r, f, e$message))
        return(NULL)
      })

      if (!is.null(mod_A) && mod_A$converge) {
        blups_A    <- as.data.frame(summary(mod_A, coef = TRUE)$coef.random)
        gebv_all_A <- blups_A[grep("Ginv", rownames(blups_A)), , drop = FALSE]

        # Extract environment-specific GEBVs for each trait in the multivariate set
        for (tr in current_traits) {
          pattern     <- paste0("trait_", tr, ":vm\\(Sample_Id, Ginv\\)_")
          trait_blups <- gebv_all_A[grep(pattern, rownames(gebv_all_A)), , drop = FALSE]

          if (nrow(trait_blups) > 0) {
            trait_blups$ID <- gsub(pattern, "", rownames(trait_blups))
            matched_idx <- match(fold_ids, trait_blups$ID)

            p_sub <- trait_blups$solution[matched_idx] # Environment-specific GEBV
            o_sub <- pheno_final[[tr]][test_idx]       # Pre-adjusted observed phenotype

            valid <- which(!is.na(p_sub) & !is.na(o_sub))
            if (length(valid) > 1) {
              p_val <- p_sub[valid]
              o_val <- o_sub[valid]

              pa   <- cor(p_val, o_val, use = "complete.obs")
              bias <- coef(lm(o_val ~ p_val))[2]

              p_val_c <- p_val - mean(p_val, na.rm = TRUE)
              o_val_c <- o_val - mean(o_val, na.rm = TRUE)
              msep    <- mean((o_val_c - p_val_c)^2, na.rm = TRUE)

              res_A <- data.frame(
                Model = "A", Set = set_name, Trait = tr, Rep = r, Fold = f,
                N_validation = length(valid), PA = pa, Bias = bias, MSEP = msep
              )
              set_results[[length(set_results) + 1]] <- res_A

              cat(sprintf("  Rep %2d | Fold %d | Model A  | Trait: %-19s -> PA: %.4f | Bias: %.4f | MSEP: %.4f\n",
                          r, f, tr, pa, bias, msep))
            }
          }
        }
      }

      # ========================================================================
      # MODEL AD: MULTIVARIATE ADDITIVE + DOMINANCE GxD GBLUP
      # ========================================================================
      mod_AD <- tryCatch({
        fit <- asreml(
          fixed    = multiv_fixed,
          random   = ~ us(trait):vm(Sample_Id, Ginv) + us(trait):vm(Sample_Id, Gd.inv),
          residual = ~ id(units):us(trait),
          data     = temp_data,
          maxit    = 100,
          ai.sing  = TRUE,
          workspace= workspace_size,
          trace    = FALSE
        )
        update.asreml(fit, trace = FALSE)
      }, error = function(e) {
        cat(sprintf("  [ERROR Model AD | Set %s | Rep %d | Fold %d]: %s\n", set_name, r, f, e$message))
        return(NULL)
      })

      if (!is.null(mod_AD) && mod_AD$converge) {
        blups_AD   <- as.data.frame(summary(mod_AD, coef = TRUE)$coef.random)
        add_all_AD <- blups_AD[grep("Ginv", rownames(blups_AD)), , drop = FALSE]
        dom_all_AD <- blups_AD[grep("Gd.inv", rownames(blups_AD)), , drop = FALSE]

        # Extract environment-specific GEBV + GDBV for each trait in the multivariate set
        for (tr in current_traits) {
          pat_add <- paste0("trait_", tr, ":vm\\(Sample_Id, Ginv\\)_")
          pat_dom <- paste0("trait_", tr, ":vm\\(Sample_Id, Gd.inv\\)_")

          trait_add <- add_all_AD[grep(pat_add, rownames(add_all_AD)), , drop = FALSE]
          trait_dom <- dom_all_AD[grep(pat_dom, rownames(dom_all_AD)), , drop = FALSE]

          if (nrow(trait_add) > 0 && nrow(trait_dom) > 0) {
            trait_add$ID <- gsub(pat_add, "", rownames(trait_add))
            trait_dom$ID <- gsub(pat_dom, "", rownames(trait_dom))

            m_idx_add <- match(fold_ids, trait_add$ID)
            m_idx_dom <- match(fold_ids, trait_dom$ID)

            v_add <- trait_add$solution[m_idx_add]
            v_dom <- trait_dom$solution[m_idx_dom]
            p_sub <- v_add + v_dom                     # Total predicted value (GEBV + GDBV)
            o_sub <- pheno_final[[tr]][test_idx]       # Pre-adjusted observed phenotype

            valid <- which(!is.na(p_sub) & !is.na(o_sub))
            if (length(valid) > 1) {
              p_val <- p_sub[valid]
              o_val <- o_sub[valid]

              pa   <- cor(p_val, o_val, use = "complete.obs")
              bias <- coef(lm(o_val ~ p_val))[2]

              p_val_c <- p_val - mean(p_val, na.rm = TRUE)
              o_val_c <- o_val - mean(o_val, na.rm = TRUE)
              msep    <- mean((o_val_c - p_val_c)^2, na.rm = TRUE)

              res_AD <- data.frame(
                Model = "AD", Set = set_name, Trait = tr, Rep = r, Fold = f,
                N_validation = length(valid), PA = pa, Bias = bias, MSEP = msep
              )
              set_results[[length(set_results) + 1]] <- res_AD

              cat(sprintf("  Rep %2d | Fold %d | Model AD | Trait: %-19s -> PA: %.4f | Bias: %.4f | MSEP: %.4f\n",
                          r, f, tr, pa, bias, msep))
            }
          }
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

if (length(all_results) == 0) {
  stop("No cross-validation results were generated. Check for convergence errors above.")
}

raw_results <- bind_rows(all_results)

summary_results <- raw_results %>%
  group_by(Model, Set, Trait) %>%
  summarise(
    N_Evaluations = n(),
    Mean_PA       = mean(PA, na.rm = TRUE),
    SD_PA         = sd(PA, na.rm = TRUE),
    SE_PA         = SD_PA / sqrt(N_Evaluations),
    Mean_Bias     = mean(Bias, na.rm = TRUE),
    SD_Bias       = sd(Bias, na.rm = TRUE),
    SE_Bias       = SD_Bias / sqrt(N_Evaluations),
    Mean_MSEP     = mean(MSEP, na.rm = TRUE),
    SD_MSEP       = sd(MSEP, na.rm = TRUE),
    SE_MSEP       = SD_MSEP / sqrt(N_Evaluations),
    .groups       = "drop"
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

md_content <- c(
  md_content,
  "",
  "---",
  "### Metric Definitions:",
  "- **PA (Predictive Ability):** Pearson correlation between pre-adjusted phenotypes and environment-specific predicted values.",
  "- **Bias:** Regression slope of pre-adjusted phenotypes on predicted values (target = 1.0).",
  "- **MSEP:** Mean Squared Error of Prediction computed on centered values."
)

writeLines(md_content, con = md_file)
cat("\nMarkdown report successfully generated:", md_file, "\n")

# Print overall summary to console
cat("\n======================================================================\n")
cat(toupper(gsub("_", " ", run_prefix)), " SUMMARY\n")
cat("======================================================================\n")
print(as.data.frame(summary_results))
