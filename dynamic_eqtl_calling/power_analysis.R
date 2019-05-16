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

	sample_time <- rep(1:num_time_steps, each=num_individuals) -1

	# Create vector of predicted means based on linear model
	predicted_means = sample_genotype*genotype_beta + sample_time*time_beta + sample_genotype*sample_time*interaction_beta
	# Simulate vector
	simulated_expression = rnorm(predicted_means, mean=predicted_means, sd=sdev)

	total_sd = sd(simulated_expression)
	interaction_sd = sd(sample_genotype*sample_time*interaction_beta)
	genotype_sd = sd(sample_genotype*genotype_beta)
	time_sd = sd(sample_time*time_beta)
	pve = (interaction_sd*interaction_sd)/var(simulated_expression)

	return(simulated_expression)
}

simulate_step_linear_model_expression <- function(genotype, genotype_beta, time_beta, interaction_beta, sdev, num_time_steps, num_individuals) {
	# Get number of samples
	num_samples = num_time_steps*num_individuals

	# Create genotype and time step vectors ACROSS samples
	sample_genotype <- rep(genotype, num_time_steps)

	sample_time <- rep(1:num_time_steps, each=num_individuals) - 1

	sample_time_step = 1.0*(sample_time > 7)

	sample_time_step_normalized = (sample_time_step/sd(sample_time_step))*sd(sample_time)


	# Create vector of predicted means based on linear model
	predicted_means = sample_genotype*genotype_beta + sample_time_step_normalized*time_beta + sample_genotype*sample_time_step_normalized*interaction_beta
	# Simulate vector
	simulated_expression = rnorm(predicted_means, mean=predicted_means, sd=sdev)

	return(simulated_expression)
}


simulate_cubic_linear_model_expression <- function(genotype, genotype_beta, time_beta, interaction_beta, sdev, num_time_steps, num_individuals) {
	# Get number of samples
	num_samples = num_time_steps*num_individuals

	# Create genotype and time step vectors ACROSS samples
	sample_genotype <- rep(genotype, num_time_steps)

	sample_time <- rep(1:num_time_steps, each=num_individuals) - 1

	sample_time_cubic = sample_time*(sample_time-7)*(sample_time-15)

	sample_time_cubic_normalized = (sample_time_cubic/sd(sample_time_cubic))*sd(sample_time)


	# Create vector of predicted means based on linear model
	predicted_means = sample_genotype*genotype_beta + sample_time_cubic_normalized*time_beta + sample_genotype*sample_time_cubic_normalized*interaction_beta
	# Simulate vector
	simulated_expression = rnorm(predicted_means, mean=predicted_means, sd=sdev)

	return(simulated_expression)
}


simulate_sinusoidal_linear_model_expression <- function(genotype, genotype_beta, time_beta, interaction_beta, sdev, num_time_steps, num_individuals) {
	# Get number of samples
	num_samples = num_time_steps*num_individuals

	# Create genotype and time step vectors ACROSS samples
	sample_genotype <- rep(genotype, num_time_steps)
	sample_time <- rep(1:num_time_steps, each=num_individuals) - 1

	sinusoidal_time <- sin(.2*pi*sample_time)

	sinusoidal_time_normalized = (sinusoidal_time/sd(sinusoidal_time))*sd(sample_time)


	# Create vector of predicted means based on linear model
	predicted_means = sample_genotype*genotype_beta + sinusoidal_time_normalized*time_beta + sample_genotype*sinusoidal_time_normalized*interaction_beta
	# Simulate vector
	simulated_expression = rnorm(predicted_means, mean=predicted_means, sd=sdev)

	return(simulated_expression)
}

simulate_quadratic_linear_model_expression <- function(genotype, genotype_beta, time_beta, interaction_beta, sdev, num_time_steps, num_individuals, intercept1, intercept2) {
	# Get number of samples
	num_samples = num_time_steps*num_individuals

	# Create genotype and time step vectors ACROSS samples
	sample_genotype <- rep(genotype, num_time_steps)

	sample_time <- rep(1:num_time_steps, each=num_individuals) - 1

	sample_time_quadratic = (sample_time-intercept1)*(sample_time-intercept2)

	sample_time_quadratic_normalized = (sample_time_quadratic/sd(sample_time_quadratic))*sd(sample_time)

	# Create vector of predicted means based on linear model
	predicted_means = sample_genotype*genotype_beta + sample_time_quadratic_normalized*time_beta + sample_genotype*sample_time_quadratic_normalized*interaction_beta
	# Simulate vector
	simulated_expression = rnorm(predicted_means, mean=predicted_means, sd=sdev)

	return(simulated_expression)
}



compute_variance_explained <- function(observed, predicted, mean) {
	ve <- (sum((predicted - mean)^2))/(sum((observed-mean)^2))
	return(ve)
}

run_linear_model <- function(expr, genotype, num_time_steps, num_individuals) {
	# Get number of samples
	num_samples = num_time_steps*num_individuals
	# Create genotype and time step vectors ACROSS samples
	sample_genotype <- rep(genotype, num_time_steps)
	sample_time <- rep(1:num_time_steps, each=num_individuals)

	# Fit the linear model using R's LM
	fit <- lm(expr ~ sample_genotype + sample_time + sample_genotype:sample_time)

	intercept_term <- summary(fit)$coefficients[1,1]
	genotype_term <- summary(fit)$coefficients[2,1]
	time_term <- summary(fit)$coefficients[3,1]
	interaction_term <- summary(fit)$coefficients[4,1]

	#mean <- mean(expr)
	#predicted <- sample_genotype*genotype_term + sample_time*time_term + sample_genotype*sample_time*interaction_term + intercept_term


	# total_ve <- compute_variance_explained(expr, predicted, mean)
	pvalue <- summary(fit)$coefficients[4,4]
	return(pvalue)
}

run_linear_quadratic_model <- function(expr, genotype, num_time_steps, num_individuals) {
	# Get number of samples
	num_samples = num_time_steps*num_individuals
	# Create genotype and time step vectors ACROSS samples
	sample_genotype <- rep(genotype, num_time_steps)
	sample_time <- rep(1:num_time_steps, each=num_individuals)

	sample_time_squared <- sample_time*sample_time
	fit_full <- lm(expr ~  sample_genotype + sample_time + sample_time_squared + sample_genotype:sample_time + sample_genotype:sample_time_squared)
	fit_null <- lm(expr ~  sample_genotype + sample_time + sample_time_squared)
	
	obj <- lrtest(fit_null, fit_full)
	pvalue <- obj[[5]][2]

	return(pvalue)
}

run_anova_model <- function(expr, genotype, num_time_steps, num_individuals) {
	# Get number of samples
	num_samples = num_time_steps*num_individuals
	# Create genotype and time step vectors ACROSS samples
	sample_genotype_factor <- factor(rep(genotype, num_time_steps))
	sample_time_factor <- factor(rep(1:num_time_steps, each=num_individuals))

	fit <- aov(expr ~ sample_genotype_factor + sample_time_factor + sample_genotype_factor:sample_time_factor)

	#fit <- aov(factor(expr) ~ factor(sample_genotype) + factor(sample_time) + factor(sample_genotype):factor(sample_time))

	pvalue <- summary(fit)[[1]][["Pr(>F)"]][3]
	return(pvalue)

}

run_simulation <- function(maf, num_simulations, num_individuals, num_time_steps, interaction_effect_size, sdev, linear_effects_sdev, fraction_of_positives, output_file) {
	num_positive_simulations <- num_simulations*fraction_of_positives
	num_negative_simulations <- num_simulations - num_simulations*fraction_of_positives


	sink(output_file)
	header_line <- paste0("alternate_model\tlm_pvalue\tlm_quad_pvalue\tanova_pvalue\n")
	cat(header_line)
	# Run alternate model simulation
	for (sim_num in 1:num_positive_simulations) {
		genotype = simulate_genotypes(maf, num_individuals)
		expr <- simulate_linear_model_expression(genotype, rnorm(1, 0, linear_effects_sdev), rnorm(1, 0, linear_effects_sdev), interaction_effect_size, sdev, num_time_steps, num_individuals)
		valz <- c(valz, expr)
		lm_pvalue <- run_linear_model(expr, genotype, num_time_steps, num_individuals)
		lm_quad_pvalue <- run_linear_quadratic_model(expr, genotype, num_time_steps, num_individuals)
		anova_pvalue <- run_anova_model(expr, genotype, num_time_steps, num_individuals)
		cat(paste0("1\t",lm_pvalue,"\t", lm_quad_pvalue, "\t",anova_pvalue,"\n"))
	}

	# Run null model simulation
	for (sim_num in 1:num_negative_simulations) {
		genotype = simulate_genotypes(maf, num_individuals)
		expr <- simulate_linear_model_expression(genotype, rnorm(1, 0, linear_effects_sdev), rnorm(1, 0, linear_effects_sdev), 0, sdev, num_time_steps, num_individuals)
		lm_pvalue <- run_linear_model(expr, genotype, num_time_steps, num_individuals)
		lm_quad_pvalue <- run_linear_quadratic_model(expr, genotype, num_time_steps, num_individuals)
		anova_pvalue <- run_anova_model(expr, genotype, num_time_steps, num_individuals)
		cat(paste0("0\t",lm_pvalue,"\t", lm_quad_pvalue, "\t", anova_pvalue,"\n"))
	}
}

run_non_linear_simulation <- function(maf, num_simulations, num_individuals, num_time_steps, interaction_effect_size, sdev, linear_effects_sdev, fraction_of_positives, non_linear_version, output_file) {
	num_positive_simulations <- num_simulations*fraction_of_positives
	num_negative_simulations <- num_simulations - num_simulations*fraction_of_positives


	sink(output_file)
	header_line <- paste0("alternate_model\tlm_pvalue\tlm_quad_pvalue\tanova_pvalue\n")
	cat(header_line)
	# Run alternate model simulation
	for (sim_num in 1:num_positive_simulations) {
		genotype = simulate_genotypes(maf, num_individuals)
		if (non_linear_version == "quadratic_ends_unequal") {
			expr <- simulate_quadratic_linear_model_expression(genotype, rnorm(1, 0, linear_effects_sdev), rnorm(1, 0, linear_effects_sdev), interaction_effect_size, sdev, num_time_steps, num_individuals, 0, 10)
		} else if (non_linear_version == "quadratic_ends_equal") {
			expr <- simulate_quadratic_linear_model_expression(genotype, rnorm(1, 0, linear_effects_sdev), rnorm(1, 0, linear_effects_sdev), interaction_effect_size, sdev, num_time_steps, num_individuals, 0, 15)
		} else if (non_linear_version == "sinusoidal") {
			expr <- simulate_sinusoidal_linear_model_expression(genotype, rnorm(1, 0, linear_effects_sdev), rnorm(1, 0, linear_effects_sdev), interaction_effect_size, sdev, num_time_steps, num_individuals)
		} else if (non_linear_version == "step_function") {
			expr <- simulate_step_linear_model_expression(genotype, rnorm(1, 0, linear_effects_sdev), rnorm(1, 0, linear_effects_sdev), interaction_effect_size, sdev, num_time_steps, num_individuals)
		} else if (non_linear_version == "cubic") {
			expr <- simulate_cubic_linear_model_expression(genotype, rnorm(1, 0, linear_effects_sdev), rnorm(1, 0, linear_effects_sdev), interaction_effect_size, sdev, num_time_steps, num_individuals)
		} else if (non_linear_version == "linear") {
			expr <- simulate_linear_model_expression(genotype, rnorm(1, 0, linear_effects_sdev), rnorm(1, 0, linear_effects_sdev), interaction_effect_size, sdev, num_time_steps, num_individuals)
		}
		lm_pvalue <- run_linear_model(expr, genotype, num_time_steps, num_individuals)
		lm_quad_pvalue <- run_linear_quadratic_model(expr, genotype, num_time_steps, num_individuals)
		anova_pvalue <- run_anova_model(expr, genotype, num_time_steps, num_individuals)
		cat(paste0("1\t",lm_pvalue, "\t", lm_quad_pvalue,"\t",anova_pvalue,"\n"))
	}
	# Run null model simulation
	for (sim_num in 1:num_negative_simulations) {
		genotype = simulate_genotypes(maf, num_individuals)
		if (non_linear_version == "quadratic_ends_unequal") {
			expr <- simulate_quadratic_linear_model_expression(genotype, rnorm(1, 0, linear_effects_sdev), rnorm(1, 0, linear_effects_sdev), 0, sdev, num_time_steps, num_individuals, 0, 10)
		} else if (non_linear_version == "quadratic_ends_equal") {
			expr <- simulate_quadratic_linear_model_expression(genotype, rnorm(1, 0, linear_effects_sdev), rnorm(1, 0, linear_effects_sdev), 0, sdev, num_time_steps, num_individuals, 0, 15)
		} else if (non_linear_version == "sinusoidal") {
			expr <- simulate_sinusoidal_linear_model_expression(genotype, rnorm(1, 0, linear_effects_sdev), rnorm(1, 0, linear_effects_sdev), 0, sdev, num_time_steps, num_individuals)
		} else if (non_linear_version == "step_function") {
			expr <- simulate_step_linear_model_expression(genotype, rnorm(1, 0, linear_effects_sdev), rnorm(1, 0, linear_effects_sdev), 0, sdev, num_time_steps, num_individuals)
		} else if (non_linear_version == "cubic") {
			expr <- simulate_cubic_linear_model_expression(genotype, rnorm(1, 0, linear_effects_sdev), rnorm(1, 0, linear_effects_sdev), 0, sdev, num_time_steps, num_individuals)
		} else if (non_linear_version == "linear") {
			expr <- simulate_linear_model_expression(genotype, rnorm(1, 0, linear_effects_sdev), rnorm(1, 0, linear_effects_sdev), 0, sdev, num_time_steps, num_individuals)
		}
		lm_pvalue <- run_linear_model(expr, genotype, num_time_steps, num_individuals)
		lm_quad_pvalue <- run_linear_quadratic_model(expr, genotype, num_time_steps, num_individuals)
		anova_pvalue <- run_anova_model(expr, genotype, num_time_steps, num_individuals)
		cat(paste0("0\t",lm_pvalue, "\t", lm_quad_pvalue,"\t", anova_pvalue,"\n"))
	}
	sink()
}





#####################
# Command Line args
#####################
power_analysis_dir = args[1]  # Output dir



######################
# Simulation linear effect of genotype, time, and genotypeXtime
######################

num_simulations = 10000
num_time_steps = 16


num_individuals_arr = c(10, 20, 30, 40, 50, 60, 70, 80, 90, 100)

interaction_effect_size_arr = c(.002, .003, .005, .0075, .01)

sdev_arr = c(.1)
linear_effects_sdev_arr = c(.1)
fraction_of_positives_arr = c(.1, .2, .3)

mafs = c(.1, .4)

counter <- 0
for (num_individual_iter in 1:length(num_individuals_arr)) {
	for (interaction_effect_size_iter in 1:length(interaction_effect_size_arr)) {
		for (sdev_iter in 1:length(sdev_arr)) {
			for (linear_effects_sdev_iter in 1:length(linear_effects_sdev_arr)) {
				for (fraction_of_positives_iter in 1:length(fraction_of_positives_arr)) {
					for (maf_iter in 1:length(mafs)) {
						print(counter)
						maf <- mafs[maf_iter]
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
}


######################
# Simulation non-linear effect of genotype, time, and genotypeXtime
######################


num_simulations = 10000
num_time_steps = 16


num_individuals_arr = c(10, 20, 30, 40, 50, 60, 70, 80, 90, 100)

interaction_effect_size_arr = c(.002, .003, .005, .0075, .01)

non_linear_versions_arr = c( "linear", "quadratic_ends_equal", "step_function", "sinusoidal", "cubic", "quadratic_ends_unequal")


sdev_arr = c(.1)
linear_effects_sdev_arr = c(.1)
fraction_of_positives_arr = c(.2)
mafs = c(.4)

counter <- 0
for (num_individual_iter in 1:length(num_individuals_arr)) {
	for (interaction_effect_size_iter in 1:length(interaction_effect_size_arr)) {
		for (sdev_iter in 1:length(sdev_arr)) {
			for (linear_effects_sdev_iter in 1:length(linear_effects_sdev_arr)) {
				for (fraction_of_positives_iter in 1:length(fraction_of_positives_arr)) {
					for (maf_iter in 1:length(mafs)) {
						for (version_iter in 1:length(non_linear_versions_arr)) {
							print(counter)
							maf <- mafs[maf_iter]
							num_individuals <- num_individuals_arr[num_individual_iter]
							interaction_effect_size <- interaction_effect_size_arr[interaction_effect_size_iter]
							sdev <- sdev_arr[sdev_iter]
							linear_effects_sdev <- linear_effects_sdev_arr[linear_effects_sdev_iter]
							fraction_of_positives <- fraction_of_positives_arr[fraction_of_positives_iter]
							non_linear_version <- non_linear_versions_arr[version_iter]
							counter <- counter + 1
							output_file <- paste0(power_analysis_dir, "non_linear_", non_linear_version, "_simulation_results_", fraction_of_positives, "_fraction_positives_", num_individuals, "_individuals_", interaction_effect_size, "_interaction_beta_", sdev, "_sdev_", linear_effects_sdev, "_linear_effects_sd_", num_simulations,"_simulations_", num_time_steps, "_time_steps_", maf, "_maf.txt")
							run_non_linear_simulation(maf, num_simulations, num_individuals, num_time_steps, interaction_effect_size, sdev, linear_effects_sdev, fraction_of_positives,non_linear_version, output_file)
						}
					}
				}
			}
		}
	}
}


