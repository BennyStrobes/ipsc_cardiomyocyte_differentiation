# Dynamic eQTL Calling
Run linear and non-linear dynamic eQTLs on ipsc-cardiomyocyte differentiation data



## Running the code

Dynamic eqtl calling can be run within 'driver_key.sh'. 
'driver_key.sh' splits pipeline into 5 steps that can be run in series (all the user has to do is remove the if false; then code arround the desired part and submit jobs corresponding to the current section). 
Comments within 'driver_key'.sh' explain what each section does. But briefly:

PART 1: Preprocess data for dynamic eQTL calling

PART 2: Run GLM/GLMM dynamic qtl modeling over a spectrum of parameters

PART 3:  Run Downstream analysis on eQTL results

PART 4: Perform simulated dynamic qtl power analysis

PART 5: Visualize results from dynamic qtl analysis



## Computer cluster

This pipeline was written to run on midway2 rcc

## Authors

* **Ben Strober** -- [BennyStrobes](https://github.com/BennyStrobes)
