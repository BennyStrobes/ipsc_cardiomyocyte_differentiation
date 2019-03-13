args = commandArgs(trailingOnly=TRUE)
library(lme4)
library(lmtest)




simulate_haplotype <- function(num_individuals, num_minor_alleles) {
	haplotype <- numeric(num_individuals)
	haplotype[sample(1:num_individuals, num_minor_alleles)] = numeric(num_minor_alleles)+1
	if (sum(haplotype) != num_minor_alleles) {
		print("Haplotype simulation error")
		quit(status=1)
	}

	return(haplotype)
}


simulate_genotypes <- function(maf, num_individuals) {
	num_minor_alleles = maf*num_individuals
	h1 = simulate_haplotype(num_individuals, num_minor_alleles)
	h2 = simulate_haplotype(num_individuals, num_minor_alleles)
	genotype = h1 + h2
	return(genotype)
}



simulate_linear_model_expression <- function(genotype, genotype_beta, time_beta, interaction_beta, sdev, num_time_steps, num_individuals) {
	# Get number of samples
	num_samples = num_time_steps*num_individuals

	# Create genotype and time step vectors ACROSS samples
	sample_genotype <- rep(genotype, num_time_steps)
	sample_time <- rep(1:num_time_steps, each=num_individuals)

	# Create vector of predicted means based on linear model
	predicted_means = sample_genotype*genotype_beta + sample_time*time_beta + sample_genotype*sample_time*interaction_beta
	# Simulate vector
	simulated_expression = rnorm(predicted_means, mean=predicted_means, sd=sdev)

	return(simulated_expression)
}


run_linear_model <- function(expr, genotype, num_time_steps, num_individuals) {
	# Get number of samples
	num_samples = num_time_steps*num_individuals
	# Create genotype and time step vectors ACROSS samples
	sample_genotype <- rep(genotype, num_time_steps)
	sample_time <- rep(1:num_time_steps, each=num_individuals)

	# Fit the linear model using R's LM
	fit <- lm(expr ~ sample_genotype + sample_time + sample_genotype:sample_time)
	pvalue <- summary(fit)$coefficients[4,4]
	return(pvalue)
}

run_simulation <- function(maf, num_simulations, num_individuals, num_time_steps, interaction_effect_size, sdev, linear_effects_sdev, fraction_of_positives, output_file) {
	num_positive_simulations <- num_simulations*fraction_of_positives
	num_negative_simulations <- num_simulations - num_simulations*fraction_of_positives


	sink(output_file)
	header_line <- paste0("alternate_model\tpvalue\n")
	cat(header_line)
	# Run alternate model simulation
	for (sim_num in 1:num_positive_simulations) {
		genotype = simulate_genotypes(maf, num_individuals)
		expr <- simulate_linear_model_expression(genotype, rnorm(1, 0, linear_effects_sdev), rnorm(1, 0, linear_effects_sdev), interaction_effect_size, sdev, num_time_steps, num_individuals)
		pvalue <- run_linear_model(expr, genotype, num_time_steps, num_individuals)
		cat(paste0("1\t",pvalue,"\n"))
	}
	# Run null model simulation
	for (sim_num in 1:num_negative_simulations) {
		genotype = simulate_genotypes(maf, num_individuals)
		expr <- simulate_linear_model_expression(genotype, rnorm(1, 0, linear_effects_sdev), rnorm(1, 0, linear_effects_sdev), 0, sdev, num_time_steps, num_individuals)
		pvalue <- run_linear_model(expr, genotype, num_time_steps, num_individuals)
		cat(paste0("0\t",pvalue,"\n"))
	}
	sink()
}





#####################
# Command Line args
#####################
power_analysis_dir = args[1]  # Output dir



######################
# Simulation Settings
######################
maf = .2
num_simulations = 10000
num_time_steps = 16


num_individuals_arr = c(10, 20, 30, 40, 50, 60, 70, 80, 90, 100)

interaction_effect_size_arr = c(.002, .003, .005, .0075, .01)

sdev_arr = c(.1)
linear_effects_sdev_arr = c(.1)
fraction_of_positives_arr = c(.1, .2, .3)

counter <- 0
for (num_individual_iter in 1:length(num_individuals_arr)) {
	for (interaction_effect_size_iter in 1:length(interaction_effect_size_arr)) {
		for (sdev_iter in 1:length(sdev_arr)) {
			for (linear_effects_sdev_iter in 1:length(linear_effects_sdev_arr)) {
				for (fraction_of_positives_iter in 1:length(fraction_of_positives_arr)) {
					num_individuals <- num_individuals_arr[num_individual_iter]
					interaction_effect_size <- interaction_effect_size_arr[interaction_effect_size_iter]
					sdev <- sdev_arr[sdev_iter]
					linear_effects_sdev <- linear_effects_sdev_arr[linear_effects_sdev_iter]
					fraction_of_positives <- fraction_of_positives_arr[fraction_of_positives_iter]
					counter <- counter + 1
					output_file <- paste0(power_analysis_dir, "simulation_results_", fraction_of_positives, "_fraction_positives_", num_individuals, "_individuals_", interaction_effect_size, "_interaction_beta_", sdev, "_sdev_", linear_effects_sdev, "_linear_effects_sd_", num_simulations,"_simulations_", num_time_steps, "_time_steps_", maf, "_maf.txt")
					run_simulation(maf, num_simulations, num_individuals, num_time_steps, interaction_effect_size, sdev, linear_effects_sdev, fraction_of_positives, output_file)
				}
			}
		}
	}
}


