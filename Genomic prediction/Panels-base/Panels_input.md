It takes the output of `BasePop_prep.RData` and builds the reduced SNP
panels, their genomic relationship matrices and the cross-validation
folds. It saves to `GP_panel_input.RData`.

-   Random panels: 4000, 3000, 2000, 1000, 500 SNPs (proportional
    sampling across chromosomes)
-   Top-1000 panels: ranked from GBLUP (A) and GBLUP (A+D) models fitted
    to `Adj_Weight` and `Adj_SurfaceArea` (no fixed effects, as the
    phenotypes are already adjusted)
-   Every panel gets the same matrices: `Ginv` (additive) and `Gd.inv`
    (dominance), both bent

<!-- -->

    knitr::opts_chunk$set(echo = TRUE, message = FALSE, warning = TRUE)

    # Set root directory relative to the repository location
    knitr::opts_knit$set(root.dir = ".")

    library(data.table); library(dplyr); library(readr); library(snpReady); library(asreml)

    ## 
    ## Attaching package: 'dplyr'

    ## The following objects are masked from 'package:data.table':
    ## 
    ##     between, first, last

    ## The following objects are masked from 'package:stats':
    ## 
    ##     filter, lag

    ## The following objects are masked from 'package:base':
    ## 
    ##     intersect, setdiff, setequal, union

    ## Loading required package: Matrix

    ## Loading required package: matrixcalc

    ## Loading required package: stringr

    ## Loading required package: rgl

    ## Loading required package: impute

    ## Online License checked out Mon Oct 12 11:54:57 2026

    ## Loading ASReml-R version 4.2

    library(AGHmatrix); library(Matrix); library(MCMCglmm); library(MASS); library(caret)

    ## Warning: package 'MCMCglmm' was built under R version 4.4.3

    ## Loading required package: coda

    ## Warning: package 'coda' was built under R version 4.4.3

    ## Loading required package: ape

    ## 
    ## Attaching package: 'ape'

    ## The following object is masked from 'package:dplyr':
    ## 
    ##     where

    ## 
    ## Attaching package: 'MASS'

    ## The following objects are masked from 'package:asreml':
    ## 
    ##     coop, oats

    ## The following object is masked from 'package:dplyr':
    ## 
    ##     select

    ## Warning: package 'caret' was built under R version 4.4.3

    ## Loading required package: ggplot2

    ## Warning: package 'ggplot2' was built under R version 4.4.3

    ## Loading required package: lattice

    options(stringsAsFactors = FALSE)

    workspace_size <- "8GB"
    plink_path     <- "C:/Plink/plink.exe"
    base_file      <- "BasePop_prep.RData"
    out_file       <- "GP_panel_input.RData"

    master_bfile   <- "BSF_geno_2083_4680"      # post-QC PLINK files
    random_sizes   <- c(4000, 3000, 2000, 1000, 500)
    top_n          <- 1000
    sel_traits     <- c("Weight", "SurfaceArea")   # traits used to rank SNPs
    sel_models     <- c("A", "AD")                 # additive only, additive + dominance

## 1. Load base population data

    load(base_file)   # pheno_final, Xa_final, Xd_final, Ga, Ginv, Gd, Gd.inv, cv_matrix, ...

    pheno_final$Sample_Id <- factor(as.character(pheno_final$Sample_Id))
    ids <- as.character(pheno_final$Sample_Id)

    # Marker IDs matching the PLINK .bim file (drops the _A / _C / _G / _T allele suffix)
    marker_ids <- sub("_[ACGT]$", "", colnames(Xa_final))

    stopifnot(all(rownames(Ginv) == ids), all(rownames(Xa_final) == ids))

## 2. Helper functions

    # Proportional random sampling of SNPs per chromosome, topped up / trimmed to exactly n_snps
    pick_random_snps <- function(bim, n_snps) {
      snps <- bim[, .(SNP = sample(SNP, size = min(round(.N / nrow(bim) * n_snps), .N), replace = FALSE)),
                  by = CHR]$SNP
      if (length(snps) > n_snps) snps <- sample(snps, n_snps, FALSE)
      if (length(snps) < n_snps) snps <- c(snps, sample(setdiff(bim$SNP, snps), n_snps - length(snps), FALSE))
      snps
    }

    # Extract a SNP list from the master PLINK file and write a .raw dosage file
    make_panel_raw <- function(snp_file, out_prefix) {
      system(paste(plink_path, "--bfile", master_bfile, "--allow-extra-chr",
                   "--extract", snp_file, "--recode A --out", out_prefix))
    }

    # Make a G matrix positive definite by adding to the diagonal if needed
    bend_G_matrix <- function(G_mat, eps = 1e-4) {
      mat <- as.matrix(G_mat)
      mat <- (mat + t(mat)) / 2
      eig <- eigen(mat, symmetric = TRUE)
      if (any(eig$values <= 1e-6)) diag(mat) <- diag(mat) + abs(min(eig$values)) + eps
      mat
    }

    # Build the additive inverse (Ginv) and dominance inverse (Gd.inv) for one panel
    build_panel_matrices <- function(raw_file) {
      geno_p <- fread(raw_file)
      geno_p <- geno_p[match(ids, geno_p$IID), ]
      stopifnot(!anyNA(geno_p$IID))
      M_p <- as.matrix(geno_p[, -(1:6), with = FALSE])

      # Additive: impute, VanRaden GRM, bend, invert
      ready <- raw.data(data = M_p, frame = "wide", base = FALSE, call.rate = 0.9, maf = 0.005,
                        imput = TRUE, imput.type = "mean")
      Ga_p <- G.matrix(M = ready$M.clean, method = "VanRaden", format = "wide")$Ga
      Ga_p <- bend_G_matrix(Ga_p)
      Ginv_p <- MASS::ginv(Ga_p)
      dimnames(Ginv_p) <- list(ids, ids)
      attr(Ginv_p, "rowNames") <- attr(Ginv_p, "colNames") <- ids
      attr(Ginv_p, "INVERSE") <- TRUE

      # Dominance: Vitezica GRM, bend, invert, ASReml sparse format
      Gd_p <- Gmatrix(SNPmatrix = as.matrix(M_p), method = "Vitezica")
      Gd_p <- bend_G_matrix(Gd_p)
      iGd_p <- as(solve(Gd_p), "sparseMatrix")
      rownames(iGd_p) <- colnames(iGd_p) <- ids
      Gd.inv_p <- MCMCglmm::sm2asreml(iGd_p)
      attr(Gd.inv_p, "INVERSE") <- TRUE

      list(Ginv = Ginv_p, Gd.inv = Gd.inv_p, n_snps = ncol(ready$M.clean))
    }

## 3. Random SNP panels

    set.seed(123)   # fixed seed for reproducible panels

    qc_bim <- fread(paste0(master_bfile, ".bim"), col.names = c("CHR", "SNP", "CM", "BP", "A1", "A2"))

    for (n in random_sizes) {
      snps <- pick_random_snps(qc_bim, n)
      write.table(snps, paste0("snps_panel_", n, ".txt"), quote = FALSE, row.names = FALSE, col.names = FALSE)
      make_panel_raw(paste0("snps_panel_", n, ".txt"), paste0("BSF_geno_2083_", n))
    }

## 4. Top-1000 SNP selection: GBLUP models on adjusted phenotypes

    # Adjusted phenotypes already carry the fixed effects, so the models have no fixed terms
    selection_fits <- list()

    for (tr in sel_traits) {
      for (md in sel_models) {
        random_terms <- if (md == "A") "vm(Sample_Id, Ginv)" else "vm(Sample_Id, Ginv) + vm(Sample_Id, Gd.inv)"

        fit <- asreml(fixed    = as.formula(paste0("Adj_", tr, " ~ 1")),
                      random   = as.formula(paste("~", random_terms)),
                      residual = ~ idv(units),
                      data = pheno_final, maxit = 100, workspace = workspace_size, trace = FALSE)

        selection_fits[[paste(tr, md, sep = "_")]] <- update.asreml(fit)
      }
    }

    ## Notice  : GIV matrix not positive definite: Singular pivots

    ## Notice  : GIV matrix not positive definite: Singular pivots

    ## Notice  : GIV matrix not positive definite: Singular pivots

    ## Notice  : GIV matrix not positive definite: Singular pivots

    ## Notice  : GIV matrix not positive definite: Singular pivots

    ## Notice  : GIV matrix not positive definite: Singular pivots

    ## Notice  : GIV matrix not positive definite: Singular pivots

    ## Notice  : GIV matrix not positive definite: Singular pivots

## 5. Back-solve marker effects and rank SNPs

    # Allele frequencies and marker matrices (full panel)
    p <- colMeans(Xa_final) / 2
    q <- 1 - p

    Za <- sweep(Xa_final, 2, 2 * p, "-")     # centred additive markers
    Zd <- Xd_final                           # Vitezica dominance coding
    scale_G <- 2 * sum(p * q)
    scale_D <- sum((2 * p * q)^2)

    iGd <- solve(Gd)
    Ginv_mat <- as.matrix(Ginv)

    # Pull the solutions for one random term out of an ASReml model
    get_ebv <- function(model, term) {
      sol <- coef(model)$random
      idx <- startsWith(rownames(sol), paste0(term, "_"))
      setNames(sol[idx, 1], substring(rownames(sol)[idx], nchar(term) + 2))
    }

    # Convert breeding values to marker effects: Z' K^-1 u / scale
    back_solve <- function(Z, Kinv, ebv, scale) {
      use <- intersect(rownames(Z), names(ebv))
      out <- as.vector(crossprod(Z[use, , drop = FALSE], Kinv[use, use] %*% ebv[use]) / scale)
      names(out) <- marker_ids
      out
    }

    # Rank SNPs by their variance contribution
    rank_snps <- function(model, dominance = FALSE) {
      alpha <- back_solve(Za, Ginv_mat, get_ebv(model, "vm(Sample_Id, Ginv)"), scale_G)
      res <- data.frame(Marker = marker_ids, Allele_Frequency = as.numeric(p),
                        Additive_Effect = alpha, Additive_Score = 2 * p * q * alpha^2)
      res$Score <- res$Additive_Score

      if (dominance) {
        delta <- back_solve(Zd, iGd, get_ebv(model, "vm(Sample_Id, Gd.inv)"), scale_D)
        res$Dominance_Effect <- delta
        res$Dominance_Score  <- (2 * p * q)^2 * delta^2
        res$Score <- res$Additive_Score + res$Dominance_Score
      }

      res <- res[order(-res$Score), ]
      res$Rank <- seq_len(nrow(res))
      res$Score_Percent <- 100 * res$Score / sum(res$Score)
      res
    }

    top_panels <- character(0)

    for (nm in names(selection_fits)) {
      ranked <- rank_snps(selection_fits[[nm]], dominance = grepl("_AD$", nm))
      top    <- head(ranked, top_n)

      write_csv(top, paste0("Top", top_n, "_SNP_", nm, ".csv"))
      write.table(top$Marker, paste0("snps_top", top_n, "_", nm, ".txt"),
                  quote = FALSE, row.names = FALSE, col.names = FALSE)

      make_panel_raw(paste0("snps_top", top_n, "_", nm, ".txt"),
                     paste0("BSF_geno_2083_top", top_n, "_", nm))
      top_panels <- c(top_panels, paste0("top", top_n, "_", nm))
    }

## 6. Genomic relationship matrices for every panel

    panel_names <- c(as.character(random_sizes), top_panels)
    panel_mats  <- list()

    for (pn in panel_names) {
      panel_mats[[pn]] <- build_panel_matrices(paste0("BSF_geno_2083_", pn, ".raw"))
    }

    print(sapply(panel_mats, function(x) x$n_snps))   # SNPs retained per panel

## 7. Generate cross-validation folds (5-fold, 20 repetitions)

    n_samples <- nrow(pheno_final)
    n_reps    <- 20
    k_folds   <- 5

    cv_matrix <- matrix(NA_integer_, nrow = n_samples, ncol = n_reps,
                        dimnames = list(ids, paste0("Rep_", 1:n_reps)))

    set.seed(421)  # Ensures reproducible fold assignments (same as base population)

    for (r in 1:n_reps) {
      # Stratified sampling by Diet so each fold maintains equal diet representation
      folds <- createFolds(pheno_final$Diet, k = k_folds, list = TRUE, returnTrain = FALSE)
      for (k in 1:k_folds) {
        cv_matrix[folds[[k]], r] <- k
      }
    }

## 8. Save

    save(pheno_final, panel_mats, panel_names, cv_matrix, trait_info, diet_trait_info, diets, primary_traits,
         file = out_file)

    # Quick checks
    print(panel_names)

    ## [1] "4000"                   "3000"                   "2000"                  
    ## [4] "1000"                   "500"                    "top1000_Weight_A"      
    ## [7] "top1000_Weight_AD"      "top1000_SurfaceArea_A"  "top1000_SurfaceArea_AD"

    print(table(cv_matrix[, 1], pheno_final$Diet))

    ##    
    ##     BSG FVW SYK
    ##   1 111 113 193
    ##   2 111 112 193
    ##   3 111 112 193
    ##   4 111 113 193
    ##   5 111 113 193

    cat("\n✓", out_file, "saved with", length(panel_names), "panels and cv_matrix (", n_samples, "samples x", n_reps, "reps)\n")

    ## 
    ## ✓ GP_panel_input.RData saved with 9 panels and cv_matrix ( 2083 samples x 20 reps)

## Session info

    sessionInfo()

    ## R version 4.4.2 (2024-10-31 ucrt)
    ## Platform: x86_64-w64-mingw32/x64
    ## Running under: Windows 11 x64 (build 22631)
    ## 
    ## Matrix products: default
    ## 
    ## 
    ## locale:
    ## [1] LC_COLLATE=English_Australia.utf8  LC_CTYPE=English_Australia.utf8   
    ## [3] LC_MONETARY=English_Australia.utf8 LC_NUMERIC=C                      
    ## [5] LC_TIME=English_Australia.utf8    
    ## 
    ## time zone: Australia/Brisbane
    ## tzcode source: internal
    ## 
    ## attached base packages:
    ## [1] stats     graphics  grDevices utils     datasets  methods   base     
    ## 
    ## other attached packages:
    ##  [1] caret_7.0-1       lattice_0.22-6    ggplot2_4.0.2     MASS_7.3-61      
    ##  [5] MCMCglmm_2.36     ape_5.8-1         coda_0.19-4.1     AGHmatrix_2.1.4  
    ##  [9] asreml_4.2.0.355  snpReady_0.9.6    impute_1.80.0     rgl_1.3.17       
    ## [13] stringr_1.5.1     matrixcalc_1.0-6  Matrix_1.7-1      readr_2.1.5      
    ## [17] dplyr_1.1.4       data.table_1.16.4
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] tidyselect_1.2.1     timeDate_4052.112    farver_2.1.2        
    ##  [4] S7_0.2.0             fastmap_1.2.0        tensorA_0.36.2.1    
    ##  [7] pROC_1.19.1          digest_0.6.37        rpart_4.1.24        
    ## [10] timechange_0.3.0     lifecycle_1.0.4      survival_3.8-3      
    ## [13] magrittr_2.0.3       compiler_4.4.2       rlang_1.1.4         
    ## [16] tools_4.4.2          yaml_2.3.10          knitr_1.49          
    ## [19] htmlwidgets_1.6.4    bit_4.5.0.1          plyr_1.8.9          
    ## [22] RColorBrewer_1.1-3   withr_3.0.2          purrr_1.0.2         
    ## [25] stats4_4.4.2         nnet_7.3-20          grid_4.4.2          
    ## [28] future_1.75.0        globals_0.19.1       scales_1.4.0        
    ## [31] iterators_1.0.14     cli_3.6.3            crayon_1.5.3        
    ## [34] rmarkdown_2.29       generics_0.1.3       rstudioapi_0.17.1   
    ## [37] future.apply_1.20.2  reshape2_1.4.4       tzdb_0.4.0          
    ## [40] splines_4.4.2        parallel_4.4.2       base64enc_0.1-3     
    ## [43] vctrs_0.6.5          hardhat_1.4.3        jsonlite_1.8.9      
    ## [46] hms_1.1.3            bit64_4.6.0-1        listenv_1.0.0       
    ## [49] foreach_1.5.2        gower_1.0.2          recipes_1.4.0       
    ## [52] glue_1.8.0           parallelly_1.48.0    codetools_0.2-20    
    ## [55] lubridate_1.9.4      stringi_1.8.4        cubature_2.1.4      
    ## [58] gtable_0.3.6         tibble_3.2.1         pillar_1.10.1       
    ## [61] htmltools_0.5.8.1    ipred_0.9-16         lava_1.9.3          
    ## [64] R6_2.5.1             vroom_1.6.5          evaluate_1.0.3      
    ## [67] corpcor_1.6.10       class_7.3-23         Rcpp_1.0.14         
    ## [70] nlme_3.1-166         prodlim_2026.03.11   xfun_0.50           
    ## [73] zoo_1.8-12           pkgconfig_2.0.3      ModelMetrics_1.2.2.2
