# ipsc preprocess pipeline

This pipeline processes/quantifies the fastq files to organized allelic count matrices & quantile normalized expression matrices. 

## Running the code

Preprocessing can be run within 'preprocess_driver.key.sh'. 
'preprocess_driver.key.sh' splits pipeline into 6 steps that can be run in series (all the user has to do is remove the if false; then code arround the desired part and submit jobs corresponding to the current section). 
Comments within 'preprocess_driver_key'.sh' explain what each section does. But briefly:

PART 1: Concatenate all fastq files for each sample into one "merged" fastq file per sample. Samples have more than one fastq file initially because there are multiple sequencing rounds (to increase read depth).

PART 2: Runs fastqc on each of the fastq files, as well as downloads the reference genome

PART 3: Aligns fastq files and generates bam files using Rsubread.

PART 4: Quantify and normalize bam files, prepares covariate files, and visualizes processed total expression data

PART 5: Part 1 of WASP Mapping pipeline (for allele specific expression; https://github.com/bmvdgeijn/WASP)

Part 6: Part 2 of WASP Mapping pipeline


## Computer cluster

This pipeline was written to run on midway rcc

## Authors

* **Ben Strober** -- [BennyStrobes](https://github.com/BennyStrobes)
