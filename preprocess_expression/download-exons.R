#!/usr/bin/env Rscript

# Download exons for counting reads per gene with Subread featureCounts.
#
# To submit on RCC Midway:
#
#   sbatch --mem=2G --partition=broadwl download-exons.R
#
# Notes:
#
#  * View available Ensembl archives at http://www.ensembl.org/info/website/archives/index.html
#  * Output is in Simplified Annotation Format (SAF)
#     * Columns: GeneID, Chr, Start, End, Strand
#     * Coordinates are 1-based, inclusive on both ends
#  * Contains duplicate and overlapping exons (featureCounts handles this)
#  * See ?Rsubread::featureCounts for more details

# Input ------------------------------------------------------------------------


# Script provided by John Blischak (https://github.com/jdblischak/midway-subread-pipeline)


args = commandArgs(trailingOnly=TRUE)

# Directory to save genome fasta files
outdir=args[1]

# The Ensembl archive
archive <- "feb2014.archive.ensembl.org"


# The file to save the exon coordinates in SAF format
outfile <- paste0(outdir, "exons.saf")

# Setup ------------------------------------------------------------------------

suppressMessages(library("biomaRt"))

ensembl <- useMart(host = archive,
                   biomart = "ENSEMBL_MART_ENSEMBL",
                   dataset = "hsapiens_gene_ensembl")

# Download exons ---------------------------------------------------------------

exons_all <- getBM(attributes = c("ensembl_gene_id", "ensembl_exon_id",
                                  "chromosome_name", "exon_chrom_start",
                                  "exon_chrom_end", "strand", "gene_biotype"),
                   mart = ensembl)

exons_final <- exons_all[exons_all$chromosome_name %in% c(1:22, "X", "Y", "MT"),
                         c("ensembl_gene_id", "chromosome_name",
                           "exon_chrom_start", "exon_chrom_end",
                           "strand", "gene_biotype")]

colnames(exons_final) <- c("GeneID", "Chr", "Start", "End", "Strand", "Biotype")

# Sort by Ensembl gene ID, then start and end positins
exons_final <- exons_final[order(exons_final$GeneID,
                                 exons_final$Start,
                                 exons_final$End), ]
# Fix strand
exons_final$Strand <- ifelse(exons_final$Strand == 1, "+", "-")

# Save as tab-separated file in Simplified Annotation Format (SAF)
write.table(exons_final, outfile, quote = FALSE, sep = "\t",
            row.names = FALSE)

print("Done downloading exons")
