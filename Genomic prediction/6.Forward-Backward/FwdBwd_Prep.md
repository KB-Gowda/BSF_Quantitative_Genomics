It combines the base (n = 2,083 after QC) and F7 (n = 536) populations,
builds the marker and genomic relationship matrices and the adjusted
phenotypes, defines the forward and backward prediction sets, and saves
to `FwdBwd_prep.RData`.

-   Forward: train on Base, predict F7
-   Backward: train on F7, predict Base
-   Base: overall adjusted phenotypes (Weight, SurfaceArea)
    `Adj_Weight`, `Adj_SurfaceArea`, and within-diet adjusted phenotypes
    (Weight and SurfaceArea in SYK, BSG, FVW): `Adj_SYK_Weight`,
    `Adj_BSG_Weight`, `Adj_FVW_Weight`, `Adj_SYK_SurfaceArea`,
    `Adj_BSG_SurfaceArea`, `Adj_FVW_SurfaceArea`
-   F7: overall adjusted phenotypes only (Weight, SurfaceArea):
    `Adj_Weight`, `Adj_SurfaceArea` (within-diet columns are NA for F7)
-   Objects saved are the same as in `BasePop_prep.RData`:
    `pheno_final`, `M`, `M_impute`, `Xa_final`, `Xd_final`, `Ga`,
    `Ginv`, `Gd`, `Gd.inv`

Phenotypes are adjusted **separately for Base and F7**: Base with
`Colour + Sampling_Day + Diet + Tray:Diet` (overall) and
`Colour + Sampling_Day + Tray + Tray:Sampling_Day` (within diet), as in
`BasePop_prep.Rmd`; F7 (overall only) with `Colour + Tray`.

    knitr::opts_chunk$set(echo = TRUE, message = TRUE, warning = TRUE)

    # Set root directory relative to the repository location
    knitr::opts_knit$set(root.dir = ".")

    library(data.table); library(dplyr); library(snpReady); library(asreml)

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

    ## Online License checked out Fri Oct 30 16:29:35 2026

    ## Loading ASReml-R version 4.2

    library(AGHmatrix); library(Matrix); library(MCMCglmm)   # MCMCglmm provides sm2asreml

    ## Warning: package 'MCMCglmm' was built under R version 4.4.3

    ## Loading required package: coda

    ## Warning: package 'coda' was built under R version 4.4.3

    ## Loading required package: ape

    ## 
    ## Attaching package: 'ape'

    ## The following object is masked from 'package:dplyr':
    ## 
    ##     where

    options(stringsAsFactors = FALSE)

    workspace_size <- "8GB"
    plink_path     <- "C:/Plink/plink.exe"
    out_file       <- "FwdBwd_prep.RData"

    # Fixed effects used to adjust phenotypes, defined separately for each population
    # Fixed-effect models (fitted with asreml; adjusted phenotype = residuals only, no mean added)
    fixed_base_overall <- "Colour + Sampling_Day + Diet + Tray:Diet"           # Base: all diets together
    fixed_base_diet    <- "Colour + Sampling_Day + Tray + Tray:Sampling_Day"   # Base: within ONE diet (Diet is constant)
    fixed_f7_overall   <- "Colour + Tray"                                       # F7: all diets together

    adj_traits <- c("Weight", "SurfaceArea")   # only these two traits are adjusted

## 1. PLINK quality control and merge (skipped if the .raw file already exists)

    # QC of the base population
    if (!file.exists("BSF_geno_2083_4680.bed")) {
      system(paste(plink_path, "--bfile BSF_geno_2097_5562 --allow-extra-chr --make-bed --geno 0.1 --mind 0.1 --maf 0.005 --out BSF_geno_2083_4680"))
    }

    # Merge base and F7, then additive (0/1/2) coding
    if (!file.exists("BSF_geno_2619_4680.raw")) {
      system(paste(plink_path, "--bfile BSF_geno_2083_4680 --allow-extra-chr --bmerge BSF_geno_F7_536 --make-bed --out BSF_geno_2619_4680"))
      system(paste(plink_path, "--bfile BSF_geno_2619_4680 --allow-extra-chr --recode A --out BSF_geno_2619_4680"))
    }

## 2. Load phenotype and genotype data

    # Combine the phenotype files (columns present in only one file are filled with NA)
    pheno_base <- fread("Pheno_BSF.txt")
    pheno_f7   <- fread("Pheno_BSF_F7.txt")
    pheno_base$Population <- "Base"
    pheno_f7$Population   <- "F7"
    pheno <- as.data.frame(rbindlist(list(pheno_base, pheno_f7), fill = TRUE))

    write.table(pheno, "Pheno_BSF_2097_576.txt", sep = "\t", quote = FALSE, row.names = FALSE)

    geno <- fread("BSF_geno_2619_4680.raw")

    common_ids  <- intersect(pheno$Sample_Id, geno$IID)
    pheno_match <- pheno[pheno$Sample_Id %in% common_ids, , drop = FALSE]
    geno_match  <- geno[geno$IID %in% common_ids, ]
    geno_match  <- geno_match[match(pheno_match$Sample_Id, geno_match$IID), ]
    M <- as.matrix(geno_match[, -(1:6), with = FALSE])

## 3. Phenotype preparation

    pheno_match$Sample_Id    <- trimws(as.character(pheno_match$Sample_Id))
    pheno_match$Tray         <- as.factor(pheno_match$Tray)
    pheno_match$Sampling_Day <- as.factor(pheno_match$Sampling_Day)
    pheno_match$Diet         <- as.factor(pheno_match$Diet)
    pheno_match$Population   <- as.factor(pheno_match$Population)
    # Colour is a continuous covariate (numeric), NOT a factor
    if ("Colour" %in% colnames(pheno_match)) pheno_match$Colour <- as.numeric(pheno_match$Colour)

    diets          <- c("SYK", "BSG", "FVW")
    primary_traits <- c("Weight", "Length", "Width", "SurfaceArea")
    diet_traits    <- as.vector(outer(diets, primary_traits, paste, sep = "_"))
    all_traits     <- c(primary_traits, diet_traits)

    for (tr in all_traits) {
      if (tr %in% colnames(pheno_match)) pheno_match[[tr]] <- as.numeric(pheno_match[[tr]])
    }
    stopifnot(all(diets %in% levels(pheno_match$Diet)))
    print(table(pheno_match$Population, pheno_match$Diet))

    ##       
    ##        BSG CSF FVW SYK
    ##   Base 555   0 563 965
    ##   F7     0 536   0   0

## 4. Genomic QC, imputation and genomic relationship matrices

The matrices are built on all 2,619 individuals together, so `Ga`,
`Ginv` and `Gd.inv` are the same for the forward and backward scenarios.

    geno_ready <- raw.data(data = M, frame = "wide", base = FALSE, 
                           imput = TRUE, imput.type = "mean")
    M_impute <- geno_ready$M.clean

    G  <- G.matrix(M = M_impute, method = "VanRaden", format = "wide")
    Ga <- G$Ga
    colnames(Ga) <- rownames(Ga) <- pheno_match$Sample_Id

    Ginv <- MASS::ginv(Ga)
    attr(Ginv, "rowNames") <- attr(Ginv, "colNames") <- as.character(pheno_match$Sample_Id)
    dimnames(Ginv) <- list(as.character(pheno_match$Sample_Id), as.character(pheno_match$Sample_Id))
    attr(Ginv, "INVERSE") <- TRUE

    Gd <- Gmatrix(SNPmatrix = as.matrix(M), method = "Vitezica")
    colnames(Gd) <- rownames(Gd) <- pheno_match$Sample_Id

    # sm2asreml comes from MCMCglmm
    iGd <- solve(Gd)
    iGd <- as(iGd, "sparseMatrix")
    rownames(iGd) <- colnames(iGd) <- as.character(pheno_match$Sample_Id)
    Gd.inv <- MCMCglmm::sm2asreml(iGd)
    attr(Gd.inv, "INVERSE") <- TRUE

## 5. Align phenotype and genotype data

    rownames(M_impute) <- as.character(pheno_match$Sample_Id)
    common_ids <- intersect(pheno_match$Sample_Id, rownames(M_impute))

    pheno_final <- pheno_match[pheno_match$Sample_Id %in% common_ids, , drop = FALSE]
    pheno_final$Sample_Id <- as.character(pheno_final$Sample_Id)
    pheno_final <- pheno_final[match(common_ids, pheno_final$Sample_Id), , drop = FALSE]

    Xa_final <- M_impute[pheno_final$Sample_Id, , drop = FALSE]
    rownames(Xa_final) <- as.character(pheno_final$Sample_Id)

    # Sanity checks: one record per individual, and matrices in the same order as the data
    stopifnot(!anyDuplicated(pheno_final$Sample_Id),
              all(rownames(Ginv)     == pheno_final$Sample_Id),
              all(rownames(Gd)       == pheno_final$Sample_Id),
              all(rownames(Xa_final) == pheno_final$Sample_Id))

## 6. Dominance marker matrix

    p <- colMeans(Xa_final, na.rm = TRUE) / 2
    q <- 1 - p
    Xd_final <- Xa_final

    for (j in seq_len(ncol(Xa_final))) {
      p_j <- p[j]; q_j <- q[j]
      geno_col <- Xa_final[, j]
      Xd_final[geno_col == 0, j] <- -2 * (p_j^2)
      Xd_final[geno_col == 1, j] <-  2 * p_j * q_j
      Xd_final[geno_col == 2, j] <- -2 * (q_j^2)
    }
    colnames(Xd_final) <- paste0(colnames(Xa_final), "_dom")
    rownames(Xd_final) <- as.character(pheno_final$Sample_Id)

## 7. Adjusted phenotypes (separately for Base and F7)

Fitted with `asreml()` on each population’s own data; adjusted phenotype
= residuals only (no mean added). Base uses
`Colour + Sampling_Day + Diet + Tray:Diet` (overall) and
`Colour + Sampling_Day + Tray + Tray:Sampling_Day` (within diet); F7
uses `Colour + Tray` (overall).

    stopifnot(all(c("Base", "F7") %in% levels(pheno_final$Population)))

    # Fit asreml on one data frame; returns residuals for every row (NA where trait or a covariate is missing)
    fit_resid <- function(trait, rhs, dat, label) {
      vars <- unique(c(trait, all.vars(as.formula(paste("~", rhs)))))
      ok   <- complete.cases(dat[, vars, drop = FALSE])
      d    <- droplevels(dat[ok, , drop = FALSE])               # levels present in THIS subset only

      # Design diagnostics (fixed-effect rank and residual df) before fitting
      rk    <- qr(model.matrix(as.formula(paste("~", rhs)), data = d))$rank
      res_df <- nrow(d) - rk
      cat(sprintf("[%s] n = %d, fixed rank = %d, residual df = %d\n", label, nrow(d), rk, res_df))
      if (res_df < 10) {
        cat("\n--- DIAGNOSTICS for", label, "---\n")
        cat("Rows in subset:", nrow(dat), "| rows with trait non-NA:", sum(!is.na(dat[[trait]])),
            "| complete rows used:", nrow(d), "\n")
        for (v in setdiff(vars, trait)) {
          cat("\n", v, ": class ", class(d[[v]])[1], ", ", length(unique(d[[v]])), " distinct values, ",
              sum(is.na(dat[[v]])), " NA in subset\n", sep = "")
          if (!is.numeric(d[[v]])) print(head(sort(table(d[[v]]), decreasing = TRUE), 8))
        }
        stop("[", label, "] only ", res_df, " residual df: the fixed effects use up (almost) all records. ",
             "See DIAGNOSTICS above (look for a term with ~1 record per level, e.g. Tray).", call. = FALSE)
      }

      fit_once <- function(...) asreml(fixed = as.formula(paste(trait, "~", rhs)), data = d,
                                       maxit = 50, workspace = workspace_size, trace = FALSE, ...)
      fit <- tryCatch(fit_once(), error = function(e) {
        message("  [", label, "] asreml failed (", conditionMessage(e), ") - retrying with ai.sing = TRUE")
        fit_once(ai.sing = TRUE)
      })
      fit <- update.asreml(fit)

      out <- rep(NA_real_, nrow(dat))
      out[ok] <- residuals(fit, type = "response")
      out
    }

    idx_b <- which(pheno_final$Population == "Base")
    idx_f <- which(pheno_final$Population == "F7")
    base_dat <- pheno_final[idx_b, ]
    f7_dat   <- pheno_final[idx_f, ]
    cat("Base records:", nrow(base_dat), "| F7 records:", nrow(f7_dat), "\n")

    ## Base records: 2083 | F7 records: 536

    # Non-missing records per population and trait
    print(with(pheno_final, table(Population, Weight_present = !is.na(Weight))))

    ##           Weight_present
    ## Population TRUE
    ##       Base 2083
    ##       F7    536

    print(with(pheno_final, table(Population, SurfaceArea_present = !is.na(SurfaceArea))))

    ##           SurfaceArea_present
    ## Population TRUE
    ##       Base 2083
    ##       F7    536

    # Design check for F7: records per Tray
    cat("F7 Colour range:", range(f7_dat$Colour, na.rm = TRUE), "| F7 Tray levels:", nlevels(droplevels(f7_dat$Tray)), "\n")

    ## F7 Colour range: 135.7566 201.0889 | F7 Tray levels: 3

    print(summary(as.vector(table(droplevels(f7_dat$Tray)))))

    ##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
    ##   171.0   176.0   181.0   178.7   182.5   184.0

    # 7a. Overall (all diets together)
    for (tr in adj_traits) {
      col <- paste0("Adj_", tr)
      pheno_final[[col]] <- NA_real_
      pheno_final[[col]][idx_b] <- fit_resid(tr, fixed_base_overall, base_dat, paste("Base overall", tr))
      pheno_final[[col]][idx_f] <- fit_resid(tr, fixed_f7_overall,   f7_dat,   paste("F7 overall", tr))
    }

    ## [Base overall Weight] n = 2083, fixed rank = 12, residual df = 2071
    ## [F7 overall Weight] n = 536, fixed rank = 4, residual df = 532
    ## [Base overall SurfaceArea] n = 2083, fixed rank = 12, residual df = 2071
    ## [F7 overall SurfaceArea] n = 536, fixed rank = 4, residual df = 532

    # 7b. Within diet (Base data only): Weight and SurfaceArea in SYK, BSG and FVW
    for (d in diets) {
      pos <- which(base_dat$Diet == d)                      # Base animals on this diet
      for (tr in adj_traits) {
        col <- paste0("Adj_", d, "_", tr)
        pheno_final[[col]] <- NA_real_                      # NA for other diets and all F7 animals
        pheno_final[[col]][idx_b[pos]] <- fit_resid(tr, fixed_base_diet, base_dat[pos, ], paste("Base", d, tr))
      }
    }

    ## [Base SYK Weight] n = 965, fixed rank = 10, residual df = 955
    ## [Base SYK SurfaceArea] n = 965, fixed rank = 10, residual df = 955
    ## [Base BSG Weight] n = 555, fixed rank = 10, residual df = 545
    ## [Base BSG SurfaceArea] n = 555, fixed rank = 10, residual df = 545
    ## [Base FVW Weight] n = 563, fixed rank = 10, residual df = 553
    ## [Base FVW SurfaceArea] n = 563, fixed rank = 10, residual df = 553

## 8. Trait tables

    trait_info <- data.frame(Trait = adj_traits, Adjusted = paste0("Adj_", adj_traits))

    diet_trait_info <- expand.grid(Diet = diets, Trait = adj_traits, stringsAsFactors = FALSE)
    diet_trait_info$Adjusted <- paste0("Adj_", diet_trait_info$Diet, "_", diet_trait_info$Trait)

## 9. Forward and backward prediction sets

Row numbers refer to `pheno_final`, `Xa_final`, `Xd_final`, `Ga`,
`Ginv`, `Gd` and `Gd.inv` (all in the same order).

    idx_base <- which(pheno_final$Population == "Base")
    idx_f7   <- which(pheno_final$Population == "F7")

    scenarios <- list(
      forward  = list(train = idx_base, valid = idx_f7),   # train on Base, predict F7
      backward = list(train = idx_f7,   valid = idx_base)  # train on F7, predict Base
    )

    cat("Forward : train =", length(scenarios$forward$train),  " valid =", length(scenarios$forward$valid),  "\n")

    ## Forward : train = 2083  valid = 536

    cat("Backward: train =", length(scenarios$backward$train), " valid =", length(scenarios$backward$valid), "\n")

    ## Backward: train = 536  valid = 2083

## 10. Save

    # Convert Sample_Id to factor before saving so downstream ASReml models recognize it
    pheno_final$Sample_Id <- factor(pheno_final$Sample_Id)

    save(pheno_final, M, M_impute, Xa_final, Xd_final, Ga, Ginv, Gd, Gd.inv,
         trait_info, diet_trait_info, diets, primary_traits, adj_traits,
         idx_base, idx_f7, scenarios, file = out_file)

    write.csv(pheno_final, "FwdBwd_adjusted_phenotypes.csv", row.names = FALSE)

    # Quick checks
    print(table(pheno_final$Population, pheno_final$Diet))

    ##       
    ##        BSG CSF FVW SYK
    ##   Base 555   0 563 965
    ##   F7     0 536   0   0

    print(colSums(!is.na(pheno_final[, grep("^Adj_", colnames(pheno_final))])))

    ##          Adj_Weight     Adj_SurfaceArea      Adj_SYK_Weight Adj_SYK_SurfaceArea 
    ##                2619                2619                 965                 965 
    ##      Adj_BSG_Weight Adj_BSG_SurfaceArea      Adj_FVW_Weight Adj_FVW_SurfaceArea 
    ##                 555                 555                 563                 563

    str(pheno_final)

    ## 'data.frame':    2619 obs. of  36 variables:
    ##  $ Sl_No              : int  1 2 3 4 5 6 7 8 9 10 ...
    ##  $ Sample_Id          : Factor w/ 2619 levels "BSF_F0_B1_19A01",..: 556 557 558 559 568 569 570 571 580 581 ...
    ##  $ Generation         : chr  "Base" "Base" "Base" "Base" ...
    ##  $ Tray               : Factor w/ 12 levels "B1","B2","B3",..: 4 4 4 4 4 4 4 4 4 4 ...
    ##  $ Diet               : Factor w/ 4 levels "BSG","CSF","FVW",..: 4 4 4 4 4 4 4 4 4 4 ...
    ##  $ Sampling_Day       : Factor w/ 4 levels "13","14","15",..: 1 1 1 1 1 1 1 1 1 1 ...
    ##  $ Colour             : num  122 166 145 176 145 ...
    ##  $ Weight             : num  220 197 212 132 220 145 234 149 208 177 ...
    ##  $ Length             : num  18.8 22.2 21.6 17.8 21.3 ...
    ##  $ Width              : num  4.36 5.21 5.36 3.47 5.72 4.21 5.73 4.03 5.78 5.11 ...
    ##  $ SurfaceArea        : num  76 95.7 97.8 47.7 99.6 ...
    ##  $ SYK_Weight         : num  220 197 212 132 220 145 234 149 208 177 ...
    ##  $ BSG_Weight         : num  NA NA NA NA NA NA NA NA NA NA ...
    ##  $ FVW_Weight         : num  NA NA NA NA NA NA NA NA NA NA ...
    ##  $ CSF_Weight         : num  NA NA NA NA NA NA NA NA NA NA ...
    ##  $ SYK_Length         : num  18.8 22.2 21.6 17.8 21.3 ...
    ##  $ BSG_Length         : num  NA NA NA NA NA NA NA NA NA NA ...
    ##  $ FVW_Length         : num  NA NA NA NA NA NA NA NA NA NA ...
    ##  $ CSF_Length         : num  NA NA NA NA NA NA NA NA NA NA ...
    ##  $ SYK_Width          : num  4.36 5.21 5.36 3.47 5.72 4.21 5.73 4.03 5.78 5.11 ...
    ##  $ BSG_Width          : num  NA NA NA NA NA NA NA NA NA NA ...
    ##  $ FVW_Width          : num  NA NA NA NA NA NA NA NA NA NA ...
    ##  $ CSF_Width          : num  NA NA NA NA NA NA NA NA NA NA ...
    ##  $ SYK_SurfaceArea    : num  76 95.7 97.8 47.7 99.6 ...
    ##  $ BSG_SurfaceArea    : num  NA NA NA NA NA NA NA NA NA NA ...
    ##  $ FVW_SurfaceArea    : num  NA NA NA NA NA NA NA NA NA NA ...
    ##  $ CSF_SurfaceArea    : num  NA NA NA NA NA NA NA NA NA NA ...
    ##  $ Population         : Factor w/ 2 levels "Base","F7": 1 1 1 1 1 1 1 1 1 1 ...
    ##  $ Adj_Weight         : num  32.2 27 33.4 -34.1 41.4 ...
    ##  $ Adj_SurfaceArea    : num  -11.8 21.7 17.1 -23.2 18.9 ...
    ##  $ Adj_SYK_Weight     : num  21.7 24.2 26.8 -35.4 34.8 ...
    ##  $ Adj_SYK_SurfaceArea: num  -15 22.2 15.8 -22 17.6 ...
    ##  $ Adj_BSG_Weight     : num  NA NA NA NA NA NA NA NA NA NA ...
    ##  $ Adj_BSG_SurfaceArea: num  NA NA NA NA NA NA NA NA NA NA ...
    ##  $ Adj_FVW_Weight     : num  NA NA NA NA NA NA NA NA NA NA ...
    ##  $ Adj_FVW_SurfaceArea: num  NA NA NA NA NA NA NA NA NA NA ...

    cat("\n✓ FwdBwd_prep.RData saved (", nrow(pheno_final), "samples:",
        length(idx_base), "Base +", length(idx_f7), "F7 )\n")

    ## 
    ## ✓ FwdBwd_prep.RData saved ( 2619 samples: 2083 Base + 536 F7 )

## Session info

    sessionInfo()

    ## R version 4.4.2 (2024-10-31 ucrt)
    ## Platform: x86_64-w64-mingw32/x64
    ## Running under: Windows 11 x64 (build 26100)
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
    ##  [1] MCMCglmm_2.36     ape_5.8-1         coda_0.19-4.1     AGHmatrix_2.1.4  
    ##  [5] asreml_4.2.0.355  snpReady_0.9.6    impute_1.80.0     rgl_1.3.17       
    ##  [9] stringr_1.5.1     matrixcalc_1.0-6  Matrix_1.7-1      dplyr_1.1.4      
    ## [13] data.table_1.16.4
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] tensorA_0.36.2.1   generics_0.1.3     stringi_1.8.4      lattice_0.22-6    
    ##  [5] cubature_2.1.4     digest_0.6.37      magrittr_2.0.3     evaluate_1.0.3    
    ##  [9] grid_4.4.2         RColorBrewer_1.1-3 fastmap_1.2.0      jsonlite_1.8.9    
    ## [13] scales_1.4.0       cli_3.6.3          rlang_1.1.4        base64enc_0.1-3   
    ## [17] yaml_2.3.10        tools_4.4.2        parallel_4.4.2     corpcor_1.6.10    
    ## [21] ggplot2_4.0.2      vctrs_0.6.5        R6_2.5.1           zoo_1.8-12        
    ## [25] lifecycle_1.0.4    htmlwidgets_1.6.4  MASS_7.3-61        pkgconfig_2.0.3   
    ## [29] pillar_1.10.1      gtable_0.3.6       glue_1.8.0         Rcpp_1.0.14       
    ## [33] xfun_0.50          tibble_3.2.1       tidyselect_1.2.1   rstudioapi_0.17.1 
    ## [37] knitr_1.49         farver_2.1.2       htmltools_0.5.8.1  nlme_3.1-166      
    ## [41] rmarkdown_2.29     compiler_4.4.2     S7_0.2.0
