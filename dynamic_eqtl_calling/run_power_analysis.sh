#!/bin/bash
#SBATCH --time=10:00:00 --partition=broadwl --mem=5GB

power_analysis_dir="$1"

module unload R
module load R/3.5.1

Rscript power_analysis.R $power_analysis_dir
