#!/usr/bin/env Rscript

# Run Subread to map reads to the genome.
#
# To submit on RCC Midway:
#
#   sbatch --mem=12g --partition=broadwl run-subread.R <file.fastq.gz>

# Input ------------------------------------------------------------------------

# Obtain name of fastq file passed at the command line
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 3)
fastq <- args[1]
bam_dir <- args[2]
genome_dir <- args[3]
stopifnot(file.exists(fastq))

# Prefix of indexed genome files

genome_prefix <- paste0(genome_dir,"GRCh37")
# Output directory
outdir <- bam_dir

# Setup ------------------------------------------------------------------------

library("Rsubread")

dir.create(outdir, showWarnings = FALSE)

# Create name of output BAM file from input fastq file
outfile <- paste0(outdir,
                  sub(".fastq.gz$", ".bam", basename(fastq)))

# Map reads with Subread -------------------------------------------------------

align(index = genome_prefix, readfile1 = fastq, output_file = outfile)
