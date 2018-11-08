#!/usr/bin/env Rscript

# Download and index human genome for mapping with Subread.
#
# To submit on RCC Midway:
#
#   sbatch --mem=12G --partition=broadwl download-genome.R
#
# Notes:
#
#  * The indexing step requires at least 8 GB of memory.
#
#  * The available releases on the Ensembl FTP site can be viewed at
#    ftp://ftp.ensembl.org/pub/

# Input ------------------------------------------------------------------------


# Script provided by John Blischak (https://github.com/jdblischak/midway-subread-pipeline)


args = commandArgs(trailingOnly=TRUE)

# Directory to save genome fasta files
outdir=args[1]


# Chromsomes for mapping
chroms <- c(1:22, "X", "Y", "MT")


# Ensembl release to use for annotation
ensembl_rel <- "75"

# The genome build used by the Ensembl release
ensembl_genome <- "GRCh37"

# Setup ------------------------------------------------------------------------

library("Rsubread")

dir.create(outdir, showWarnings = FALSE)

# Construct the URL to the Ensembl FTP site
ensembl_ftp <- paste0("ftp://ftp.ensembl.org/pub/release-",
                      ensembl_rel, "/fasta/homo_sapiens/dna/")

# Download chromosome fasta files ----------------------------------------------

for (chr in chroms) {
  cat("Downloading chromesome", chr, "...\n")
  chr_url <- paste0(ensembl_ftp, "Homo_sapiens.",
                    ensembl_genome, ".",ensembl_rel,".dna_sm.chromosome.",
                    chr, ".fa.gz")
  chr_file <- paste0(outdir, "Homo_sapiens.",
                     ensembl_genome, ".dna_sm.chromosome.",
                     chr, ".fa.gz")
  download.file(url = chr_url, destfile = chr_file)
}

# Decompress and combine fasta files -------------------------------------------

# Name of temporary fasta file
fa_tmp <- paste0("/tmp/tmp.fa")

cmd <- sprintf("zcat %s*fa.gz > %s", outdir, fa_tmp)
system(cmd)

# Index genome for mappig with Subread -----------------------------------------

buildindex(basename = paste0(outdir, ensembl_genome),
           reference = fa_tmp)

print("Done downloading genome")
