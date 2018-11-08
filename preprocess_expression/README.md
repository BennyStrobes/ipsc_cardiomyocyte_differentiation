# ipsc preprocess pipeline

This pipeline processes/quantifies the fastq files to organized allelic count matrices & quantile normalized expression matrices. It can be run through 'preprocess_driver.key.sh'
'preprocess_driver.key.sh' splits pipeline into 7 steps that can be run in series (all the user has to do is remove the if false; then code arround the desired part). The code is well-commented so it should be clear from 'preprocess_driver_key'.sh' what each section does.

## Deliverables

As for important output files from this pipeline (using dir_names defined in 'preprocess_driver_key.sh':

	1. $preprocess_total_expression_dir"quantile_normalized.txt" contains quantile normalized expression matrices for all samples across all protein-coding autosomal genes that have at least 10 samples such that RPKM >= .1 and counts >= 6.

	2. $preprocess_total_expression_dir"rpkm.txt" contains raw rpkm expression matrices for all samples across alll protein-coding autosomal genes that have at least 10 samples such that RPKM >= .1 and counts >= 6

	3. $covariate_dir"principal_components_10.txt" contains matrix for loadings of first 10 PCs for all samples

	4. $covariate_dir"processed_covariates_categorical.txt" contains an organized matrix of covariates

	5. $covariate_dir"sva_loadings.txt" contains results of running sva (while inputting time to allow sva to make factors independent of time)

	6. $visualize_total_expression_dir: a lot of plots describing total expression quantification

	7. $genotype_dir"YRI_genotype.vcf": imputed-dosage based genotypes for our 10 cell lines

	8. $genotype_dir"YRI_het_prob_genotype.vcf": heterozygous probabilities for our 10 cell lines (based on impute2)

	9. $processed_allelic_counts_dir"allelic_counts_gene_mapped_het_prob_*.txt" where * is the heterozygous probability threshold used to call heterozygous sites. I computed this matrices for a bunch of different thresholds. This matrix contains refAlleleCounts_totalCounts in each cell. Currently no filtering on sites (this will be changed)

	10. $visualize_allelic_counts_dir contains plots describing allelic count quantification

	11. $fastq_input_dir"fastq_mapping.txt" contains list of mapping from sample_ids (cellLine_timeStep) to fastq files composing those samples

	12. $preprocess_total_expression_dir"quant_expr_sva_corrected.txt contains residual quantile normalized expression matrix after regressing out sva latent factors.


Let me know if something isn't clear or you think you found a mistake!

## Computer cluster

This pipeline was written to run on midway rcc

## Authors

* **Ben Strober** -- [BennyStrobes](https://github.com/BennyStrobes)
