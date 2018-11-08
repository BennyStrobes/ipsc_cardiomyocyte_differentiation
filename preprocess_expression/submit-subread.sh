#!/bin/bash
fastq_dir="$1"
bam_dir="$2"
genome_dir="$3"


for fq in $fastq_dir*fastq.gz
do
  echo "Submitting" $fq
  sbatch --mem=30G --partition=broadwl run-subread.R $fq $bam_dir $genome_dir
done

