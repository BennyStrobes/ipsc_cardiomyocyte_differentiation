#!/bin/bash
#SBATCH --time=24:00:00 --mem=12G --partition=broadwl
fastq_input_dir="$1"
fastqc_dir="$2"
genome_dir="$3"

# Script provided by John Blischak (https://github.com/jdblischak/midway-subread-pipeline)
Rscript download-genome.R $genome_dir
date

# Script provided by John Blischak (https://github.com/jdblischak/midway-subread-pipeline)
Rscript download-exons.R $genome_dir


module load fastqc/0.11.5

# Run fastqc
fastqc $fastq_input_dir*fastq.gz -o $fastqc_dir
date

# Multiqc combines all of our fastqc files into one organized one!
multiqc --force --outdir $fastqc_dir $fastqc_dir $fastq_input_dir
