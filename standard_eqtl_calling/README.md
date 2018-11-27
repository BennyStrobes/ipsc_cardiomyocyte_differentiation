# Standard eQTL Calling
Run the WASP Combined Haplotype Test (CHT; https://github.com/bmvdgeijn/WASP) in each time step, independently.



## Running the code

Standard eqtl calling can be run within 'driver_key.sh'. 
'driver_key.sh' splits pipeline into 6 steps that can be run in series (all the user has to do is remove the if false; then code arround the desired part and submit jobs corresponding to the current section). 
Comments within 'driver_key'.sh' explain what each section does. But briefly:

PART 1: Extract list of variant-gene pairs to be tested by WASP-CHT (script based on orignial WASP code 'get_target_regions.py') in each time step seperately.

PART 2: Merge list of variant-gene pairs (from PART 1) across all 16 time steps. And created merged list of variant gene pairs that pass filters in all 16 time steps.

PART 3: Make CHT input file for each RNA-seq sample

PART 4: Update total read depth (directly from https://github.com/bmvdgeijn/WASP)

PART 5: Fit hyper-parameters for WASP, offline. Also, compute PCs to control for

PART 6: Run WASP combined haplotype test (very computationally expensive b/c we have to run it for all time steps, over a range of PCs)

PART 7: Organize and perform multiple testing correction on WASP CHT results

PART 8: Run downstream analysis/visualization on WASP CHT results.



## Computer cluster

This pipeline was written to run on midway2 rcc

## Authors

* **Ben Strober** -- [BennyStrobes](https://github.com/BennyStrobes)
