It builds the marker and genomic relationship matrices and the adjusted phenotypes, and saves to `BasePop_prep.RData`.

-   Overall adjusted phenotypes (4 traits): `Adj_Weight`, `Adj_Length`,
    `Adj_Width`, `Adj_SurfaceArea`
-   Within-diet adjusted phenotypes (Weight and SurfaceArea in SYK, BSG,
    FVW): `Adj_SYK_Weight`, `Adj_BSG_Weight`, `Adj_FVW_Weight`,
    `Adj_SYK_SurfaceArea`, `Adj_BSG_SurfaceArea`, `Adj_FVW_SurfaceArea`

<!-- -->

    knitr::opts_chunk$set(echo = TRUE, message = FALSE, warning = TRUE)

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

    ## Online License checked out Mon Oct  5 10:07:21 2026

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
    out_file       <- "BasePop_prep.RData"

    # Fixed effects used to adjust phenotypes
    overall_fixed <- "Colour + Sampling_Day + Diet + Tray:Diet"   # all four traits, all diets
    diet_fixed    <- "Colour + Sampling_Day + Tray + Tray:Sampling_Day"   # within ONE diet (Diet is constant)
    diet_adj_traits <- c("Weight", "SurfaceArea")                 # traits adjusted within each diet

## 1. PLINK quality control (skipped if the .raw file already exists)

    if (!file.exists("BSF_geno_2083_4680.raw")) {
      system(paste(plink_path, "--bfile BSF_geno_2097_5562 --allow-extra-chr --make-bed --geno 0.1 --mind 0.1 --maf 0.005 --recode A --out BSF_geno_2083_4680"))
    }

## 2. Load phenotype and genotype data

    pheno <- read.table("Pheno_BSF.txt", header = TRUE, sep = "\t")
    geno  <- fread("BSF_geno_2083_4680.raw")

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

    diets          <- c("SYK", "BSG", "FVW")
    primary_traits <- c("Weight", "Length", "Width", "SurfaceArea")
    diet_traits    <- as.vector(outer(diets, primary_traits, paste, sep = "_"))
    all_traits     <- c(primary_traits, diet_traits)

    for (tr in all_traits) {
      if (tr %in% colnames(pheno_match)) pheno_match[[tr]] <- as.numeric(pheno_match[[tr]])
    }
    stopifnot(all(diets %in% levels(pheno_match$Diet)))

## 4. Genomic QC, imputation and genomic relationship matrices

    geno_ready <- raw.data(data = M, frame = "wide", base = FALSE, call.rate = 0.9, maf = 0.005,
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

## 7. Adjusted phenotypes

    # 7a. Overall: four traits, all diets together
    for (tr in primary_traits) {
      fit <- asreml(fixed = as.formula(paste(tr, "~", overall_fixed)),
                    data = pheno_final, maxit = 50, workspace = workspace_size, trace = FALSE)
      fit <- update.asreml(fit)
      pheno_final[[paste0("Adj_", tr)]] <- residuals(fit, type = "response")
    }

    # 7b. Within diet: Weight and SurfaceArea in SYK, BSG and FVW
    for (d in diets) {
      for (tr in diet_adj_traits) {
        rows      <- which(pheno_final$Diet == d & !is.na(pheno_final[[tr]]))
        diet_data <- droplevels(pheno_final[rows, ])

        fit <- asreml(fixed = as.formula(paste(tr, "~", diet_fixed)),
                      data = diet_data, maxit = 50, workspace = workspace_size, trace = FALSE)
        fit <- update.asreml(fit)

        col <- paste0("Adj_", d, "_", tr)
        pheno_final[[col]] <- NA_real_                       # NA for animals on the other diets
        pheno_final[[col]][rows] <- residuals(fit, type = "response")
      }
    }

## 8. Trait tables

    trait_info <- data.frame(Trait = primary_traits, Adjusted = paste0("Adj_", primary_traits))

    diet_trait_info <- expand.grid(Diet = diets, Trait = diet_adj_traits, stringsAsFactors = FALSE)
    diet_trait_info$Adjusted <- paste0("Adj_", diet_trait_info$Diet, "_", diet_trait_info$Trait)

## 9. Generate cross-validation folds (5-fold, 20 repetitions)

    library(caret)

    ## Warning: package 'caret' was built under R version 4.4.3

    ## Warning: package 'ggplot2' was built under R version 4.4.3

    n_samples <- nrow(pheno_final)
    n_reps    <- 20
    k_folds   <- 5

    cv_matrix <- matrix(NA_integer_, nrow = n_samples, ncol = n_reps,
                        dimnames = list(pheno_final$Sample_Id, paste0("Rep_", 1:n_reps)))

    set.seed(421)  # Ensures reproducible fold assignments

    for (r in 1:n_reps) {
      # Stratified sampling by Diet so each fold maintains equal diet representation
      folds <- createFolds(pheno_final$Diet, k = k_folds, list = TRUE, returnTrain = FALSE)
      for (k in 1:k_folds) {
        cv_matrix[folds[[k]], r] <- k
      }
    }

## 10. Save

    # Convert Sample_Id to factor before saving so downstream ASReml models recognize it
    pheno_final$Sample_Id <- factor(pheno_final$Sample_Id)

    save(pheno_final, M, M_impute, Xa_final, Xd_final, Ga, Ginv, Gd, Gd.inv,
         trait_info, diet_trait_info, diets, primary_traits, cv_matrix, file = out_file)

    write.csv(pheno_final, "BasePop_adjusted_phenotypes.csv", row.names = FALSE)

    # Quick checks
    print(table(pheno_final$Diet))

    ## 
    ## BSG FVW SYK 
    ## 555 563 965

    print(colSums(!is.na(pheno_final[, grep("^Adj_", colnames(pheno_final))])))

    ##          Adj_Weight          Adj_Length           Adj_Width     Adj_SurfaceArea 
    ##                2083                2083                2083                2083 
    ##      Adj_SYK_Weight Adj_SYK_SurfaceArea      Adj_BSG_Weight Adj_BSG_SurfaceArea 
    ##                 965                 965                 555                 555 
    ##      Adj_FVW_Weight Adj_FVW_SurfaceArea 
    ##                 563                 563

    str(pheno_final)

    ## 'data.frame':    2083 obs. of  32 variables:
    ##  $ Sl_No              : int  1 2 3 4 5 6 7 8 9 10 ...
    ##  $ Sample_Id          : Factor w/ 2083 levels "BSF_F0_B1_19A01",..: 556 557 558 559 568 569 570 571 580 581 ...
    ##  $ Tray               : Factor w/ 9 levels "B1","B2","B3",..: 4 4 4 4 4 4 4 4 4 4 ...
    ##  $ Diet               : Factor w/ 3 levels "BSG","FVW","SYK": 3 3 3 3 3 3 3 3 3 3 ...
    ##  $ Sampling_Day       : Factor w/ 3 levels "13","14","15": 1 1 1 1 1 1 1 1 1 1 ...
    ##  $ Colour             : num  122 166 145 176 145 ...
    ##  $ Weight             : num  220 197 212 132 220 145 234 149 208 177 ...
    ##  $ Length             : num  18.8 22.2 21.6 17.8 21.3 ...
    ##  $ Width              : num  4.36 5.21 5.36 3.47 5.72 4.21 5.73 4.03 5.78 5.11 ...
    ##  $ SurfaceArea        : num  76 95.7 97.8 47.7 99.6 ...
    ##  $ SYK_Weight         : num  220 197 212 132 220 145 234 149 208 177 ...
    ##  $ BSG_Weight         : num  NA NA NA NA NA NA NA NA NA NA ...
    ##  $ FVW_Weight         : num  NA NA NA NA NA NA NA NA NA NA ...
    ##  $ SYK_Length         : num  18.8 22.2 21.6 17.8 21.3 ...
    ##  $ BSG_Length         : num  NA NA NA NA NA NA NA NA NA NA ...
    ##  $ FVW_Length         : num  NA NA NA NA NA NA NA NA NA NA ...
    ##  $ SYK_Width          : num  4.36 5.21 5.36 3.47 5.72 4.21 5.73 4.03 5.78 5.11 ...
    ##  $ BSG_Width          : num  NA NA NA NA NA NA NA NA NA NA ...
    ##  $ FVW_Width          : num  NA NA NA NA NA NA NA NA NA NA ...
    ##  $ SYK_SurfaceArea    : num  76 95.7 97.8 47.7 99.6 ...
    ##  $ BSG_SurfaceArea    : num  NA NA NA NA NA NA NA NA NA NA ...
    ##  $ FVW_SurfaceArea    : num  NA NA NA NA NA NA NA NA NA NA ...
    ##  $ Adj_Weight         : num  32.2 27 33.4 -34.1 41.4 ...
    ##  $ Adj_Length         : num  -0.722 3.735 2.65 -0.484 2.34 ...
    ##  $ Adj_Width          : num  -0.736 0.477 0.45 -1.185 0.81 ...
    ##  $ Adj_SurfaceArea    : num  -11.8 21.7 17.1 -23.2 18.9 ...
    ##  $ Adj_SYK_Weight     : num  21.7 24.2 26.8 -35.4 34.8 ...
    ##  $ Adj_SYK_SurfaceArea: num  -15 22.2 15.8 -22 17.6 ...
    ##  $ Adj_BSG_Weight     : num  NA NA NA NA NA NA NA NA NA NA ...
    ##  $ Adj_BSG_SurfaceArea: num  NA NA NA NA NA NA NA NA NA NA ...
    ##  $ Adj_FVW_Weight     : num  NA NA NA NA NA NA NA NA NA NA ...
    ##  $ Adj_FVW_SurfaceArea: num  NA NA NA NA NA NA NA NA NA NA ...

    cat("\n✓ BasePop_prep.RData saved with cv_matrix (", n_samples, "samples x 20 reps)\n")

    ## 
    ## ✓ BasePop_prep.RData saved with cv_matrix ( 2083 samples x 20 reps)

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
    ##  [1] caret_7.0-1       lattice_0.22-6    ggplot2_4.0.2     MCMCglmm_2.36    
    ##  [5] ape_5.8-1         coda_0.19-4.1     AGHmatrix_2.1.4   asreml_4.2.0.355 
    ##  [9] snpReady_0.9.6    impute_1.80.0     rgl_1.3.17        stringr_1.5.1    
    ## [13] matrixcalc_1.0-6  Matrix_1.7-1      dplyr_1.1.4       data.table_1.16.4
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] tidyselect_1.2.1     timeDate_4052.112    farver_2.1.2        
    ##  [4] S7_0.2.0             fastmap_1.2.0        tensorA_0.36.2.1    
    ##  [7] pROC_1.19.1          digest_0.6.37        rpart_4.1.24        
    ## [10] timechange_0.3.0     lifecycle_1.0.4      survival_3.8-3      
    ## [13] magrittr_2.0.3       compiler_4.4.2       rlang_1.1.4         
    ## [16] tools_4.4.2          yaml_2.3.10          knitr_1.49          
    ## [19] htmlwidgets_1.6.4    plyr_1.8.9           RColorBrewer_1.1-3  
    ## [22] withr_3.0.2          purrr_1.0.2          nnet_7.3-20         
    ## [25] grid_4.4.2           stats4_4.4.2         future_1.75.0       
    ## [28] globals_0.19.1       scales_1.4.0         iterators_1.0.14    
    ## [31] MASS_7.3-61          cli_3.6.3            rmarkdown_2.29      
    ## [34] generics_0.1.3       rstudioapi_0.17.1    future.apply_1.20.2 
    ## [37] reshape2_1.4.4       splines_4.4.2        parallel_4.4.2      
    ## [40] base64enc_0.1-3      vctrs_0.6.5          hardhat_1.4.3       
    ## [43] jsonlite_1.8.9       listenv_1.0.0        foreach_1.5.2       
    ## [46] gower_1.0.2          recipes_1.4.0        glue_1.8.0          
    ## [49] parallelly_1.48.0    codetools_0.2-20     lubridate_1.9.4     
    ## [52] stringi_1.8.4        cubature_2.1.4       gtable_0.3.6        
    ## [55] tibble_3.2.1         pillar_1.10.1        htmltools_0.5.8.1   
    ## [58] ipred_0.9-16         lava_1.9.3           R6_2.5.1            
    ## [61] evaluate_1.0.3       corpcor_1.6.10       class_7.3-23        
    ## [64] Rcpp_1.0.14          nlme_3.1-166         prodlim_2026.03.11  
    ## [67] xfun_0.50            zoo_1.8-12           pkgconfig_2.0.3     
    ## [70] ModelMetrics_1.2.2.2
