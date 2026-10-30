# ==============================================================================
# CROSS-DIET GENOMIC PREDICTION PIPELINE - GBLUP (BIVARIATE AVERAGING)
# Models Evaluated: Model A (Additive-Only) vs Model AD (Additive + Dominance)
# Target Metrics: PA, Bias, MSEP (Using Pre-Adjusted Phenotypes Only)
# ==============================================================================

# ------------------------------------------------------------------------------
# 1. SETUP & LOAD WORKSPACE
# ------------------------------------------------------------------------------

library(asreml)
library(dplyr)
library(Matrix)
library(tidyr)

run_prefix     <- "GBLUP_CrossDiet_Multivariate"
n_threads      <- 4
workspace_size <- "8GB"

options(sinksToConsole = TRUE)
sink(stdout(), type = "output")
asreml.options(nthreads = n_threads)

load("BasePop_prep.RData")
set.seed(20260825)

# ------------------------------------------------------------------------------
# 2. DEFINE TRAIT SETS
# ------------------------------------------------------------------------------
multivariate_sets <- list(
  Weight      = c("Adj_SYK_Weight", "Adj_BSG_Weight", "Adj_FVW_Weight"),
  SurfaceArea = c("Adj_SYK_SurfaceArea", "Adj_BSG_SurfaceArea", "Adj_FVW_SurfaceArea")
)

# Ensure all pre-adjusted columns are strictly numeric for ASReml
for (tr in unlist(multivariate_sets)) {
  pheno_final[[tr]] <- as.numeric(pheno_final[[tr]])
}

all_results <- list()

# ------------------------------------------------------------------------------
# 3. CROSS-DIET PREDICTION LOOP
# ------------------------------------------------------------------------------
for (set_name in names(multivariate_sets)) {
  current_traits <- multivariate_sets[[set_name]]
  
  cat("\n", paste(rep("=", 80), collapse = ""), "\n", sep = "")
  cat("STARTING CROSS-DIET GBLUP FOR SET:", set_name, "\n")
  cat("Traits in set:", paste(current_traits, collapse = ", "), "\n")
  cat(paste(rep("=", 80), collapse = ""), "\n")
  
  set_results <- list()
  
  for (target_trait in current_traits) {
    train_traits <- setdiff(current_traits, target_trait)
    
    target_idx <- which(!is.na(pheno_final[[target_trait]]))
    target_ids <- trimws(as.character(pheno_final$Sample_Id[target_idx]))
    
    if (length(target_idx) == 0) next
    
    cat(sprintf("\n--- Target Held-Out Trait: %s | Trained on: %s + %s (Target N = %d) ---\n", 
                target_trait, train_traits[1], train_traits[2], length(target_ids)))
    
    # Bivariate formula using only the two known pre-adjusted traits
    biv_formula <- as.formula(paste("cbind(", train_traits[1], ",", train_traits[2], ") ~ trait"))
    
    obs_df <- pheno_final[target_idx, c("Sample_Id", target_trait)]
    colnames(obs_df) <- c("ID", "Obs")
    obs_df$ID <- trimws(as.character(obs_df$ID))
    
    # ========================================================================
    # MODEL A: ADDITIVE-ONLY GBLUP
    # ========================================================================
    mod_A <- tryCatch({
      fit <- asreml(
        fixed     = biv_formula,
        random    = ~ us(trait):vm(Sample_Id, Ginv),
        residual  = ~ id(units):us(trait),
        data      = pheno_final,
        na.action = na.method(x = "include", y = "include"),
        maxit     = 100, ai.sing = TRUE, workspace = workspace_size, trace = FALSE
      )
      update.asreml(fit, trace = FALSE)
    }, error = function(e) {
      cat(sprintf("  [ERROR Model A | Target %s]: %s\n", target_trait, e$message))
      return(NULL)
    })
    
    if (!is.null(mod_A) && mod_A$converge) {
      rand_A <- as.data.frame(summary(mod_A, coef = TRUE)$coef.random)
      
      gebv_A <- data.frame(Term = rownames(rand_A), Additive = rand_A[, 1]) %>%
        filter(grepl("vm\\(Sample_Id, Ginv\\)", Term)) %>%
        mutate(
          Term = gsub("^trait_", "", Term),
          Term = gsub(":vm\\(Sample_Id, Ginv\\)_", "||", Term)
        ) %>%
        separate(Term, into = c("Trait", "ID"), sep = "\\|\\|") %>%
        mutate(ID = trimws(as.character(ID))) %>%
        filter(ID %in% target_ids) %>%
        group_by(ID) %>%
        summarise(gEBV = mean(Additive, na.rm = TRUE), .groups = "drop")
      
      eval_A <- inner_join(gebv_A, obs_df, by = "ID")
      
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
    
    # ========================================================================
    # MODEL AD: ADDITIVE + DOMINANCE GBLUP
    # ========================================================================
    mod_AD <- tryCatch({
      fit <- asreml(
        fixed     = biv_formula,
        random    = ~ us(trait):vm(Sample_Id, Ginv) + us(trait):vm(Sample_Id, Gd.inv),
        residual  = ~ id(units):us(trait),
        data      = pheno_final,
        na.action = na.method(x = "include", y = "include"),
        maxit     = 100, ai.sing = TRUE, workspace = workspace_size, trace = FALSE
      )
      update.asreml(fit, trace = FALSE)
    }, error = function(e) {
      cat(sprintf("  [ERROR Model AD | Target %s]: %s\n", target_trait, e$message))
      return(NULL)
    })
    
    if (!is.null(mod_AD) && mod_AD$converge) {
      rand_AD <- as.data.frame(summary(mod_AD, coef = TRUE)$coef.random)
      
      add_df <- data.frame(Term = rownames(rand_AD), Additive = rand_AD[, 1]) %>%
        filter(grepl("vm\\(Sample_Id, Ginv\\)", Term)) %>%
        mutate(
          Term = gsub("^trait_", "", Term),
          Term = gsub(":vm\\(Sample_Id, Ginv\\)_", "||", Term)
        ) %>%
        separate(Term, into = c("Trait", "ID"), sep = "\\|\\|") %>%
        mutate(ID = trimws(as.character(ID))) %>%
        filter(ID %in% target_ids)
      
      dom_df <- data.frame(Term = rownames(rand_AD), Dominance = rand_AD[, 1]) %>%
        filter(grepl("vm\\(Sample_Id, Gd.inv\\)", Term)) %>%
        mutate(
          Term = gsub("^trait_", "", Term),
          Term = gsub(":vm\\(Sample_Id, Gd.inv\\)_", "||", Term)
        ) %>%
        separate(Term, into = c("Trait", "ID"), sep = "\\|\\|") %>%
        mutate(ID = trimws(as.character(ID))) %>%
        filter(ID %in% target_ids)
      
      gebv_AD <- inner_join(add_df, dom_df, by = c("Trait", "ID")) %>%
        mutate(Total_gEBV = Additive + Dominance) %>%
        group_by(ID) %>%
        summarise(gEBV = mean(Total_gEBV, na.rm = TRUE), .groups = "drop")
      
      eval_AD <- inner_join(gebv_AD, obs_df, by = "ID")
      
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
    flush.console()
  }
  
  if (length(set_results) > 0) {
    df_set <- bind_rows(set_results)
    all_results[[set_name]] <- df_set
    write.csv(df_set, paste0(run_prefix, "_", set_name, "_Raw.csv"), row.names = FALSE)
  }
}

# ------------------------------------------------------------------------------
# 4. SUMMARIZE & SAVE RESULTS
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
