/* ==============================================================================
   Title: Descriptive Statistics and GLM Analysis of Black Soldier Fly Growth Traits
   ============================================================================== */


/* Step 1: Import Phenotype Data */

proc import datafile="/home/u51116922/RKG/BSF/Pheno_BSF.txt"
    out=work.pheno
    dbms=dlm
    replace;
    delimiter='09'x;
    getnames=yes;
    guessingrows=max;
run;


/* Step 2: Overall Descriptive Statistics */

proc means data=pheno
    n min max mean stddev stderr cv maxdec=4;

    var Weight Length Width SurfaceArea;

run;


/* Step 3: Subgroup Descriptive Statistics */

proc means data=pheno
    n min max mean stddev stderr cv maxdec=4;

    var Weight Length Width SurfaceArea;

    class Diet Sampling_Day Tray;

run;


/* Step 4: Distribution and Normality Assessment */

proc univariate data=pheno normal;

    var Weight Length Width SurfaceArea;

    histogram / normal;

    qqplot / normal(mu=est sigma=est);

run;


/* Step 5: Fixed-Effects General Linear Model (GLM) */

proc glm data=pheno;

    class Diet Tray Sampling_Day;

    model Weight Length Width SurfaceArea =
        Gray_Scale
        Diet
        Tray(Diet)
        Sampling_Day;

    lsmeans Diet Sampling_Day Tray(Diet)
        / pdiff=all
          adjust=tukey
          stderr;

run;

quit;
