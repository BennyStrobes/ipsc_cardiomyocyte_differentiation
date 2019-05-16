args = commandArgs(trailingOnly=TRUE)
library(ggplot2)
library(ggthemes)
library(cowplot)
library(reshape)
library(lme4)
library(lmtest)

run_linear_model_five_cov <- function(model_version, expr, genotype, genotype_interaction, time_steps, time_steps_interaction, cell_lines, cov1, cov2, cov3, cov4, cov5) {
    # Run the models
    if (model_version == "glm") {  # Regular linear model
        fit <- lm(expr ~ cov1 + cov1:time_steps + cov2 + cov2:time_steps + cov3 + cov3:time_steps + cov4 + cov4:time_steps + cov5 + cov5:time_steps + genotype + time_steps + genotype_interaction:time_steps_interaction)
        pvalue <- summary(fit)$coefficients[14,4]
        af <- anova(fit)
        afss <- af$"Sum Sq"
        percent_explained <- afss/sum(afss)*100
        genotype_pve <- percent_explained[6]
        time_pve <- percent_explained[7]
        interaction_pve <- percent_explained[13]
        pve_vec <- c(genotype_pve, time_pve, interaction_pve)
    } else if (model_version == "glmm") { # Linear mixed model
        fit_full <- lmer(expr ~  genotype + time_steps + cov1 + cov1:time_steps + cov2 + cov2:time_steps + cov3 + cov3:time_steps + cov4 + cov4:time_steps + cov5 + cov5:time_steps + genotype_interaction:time_steps_interaction + (1|cell_lines), REML=FALSE)
        fit_null <- lmer(expr ~  genotype + time_steps + cov1 + cov1:time_steps + cov2 + cov2:time_steps + cov3 + cov3:time_steps + cov4 + cov4:time_steps + cov5 + cov5:time_steps + (1|cell_lines), REML=FALSE)
        coefs <- data.frame(coef(summary(fit_full)))
        coef <- coefs[14,1]
        #tvalue <- coefs[4,3]
        #pvalue <- 2*pt(-abs(tvalue),df=num_samp-1)
        lrt <- anova(fit_null,fit_full)
        pvalue <- lrt[[8]][2]
        corrected_expr <- expr - (coefs[4,1])*cov1 - (coefs[5,1])*cov2 - (coefs[6,1])*cov3 - (coefs[7,1])*cov4 - (coefs[8,1])*cov5 - (coefs[9,1])*cov1*time_steps - (coefs[10,1])*cov2*time_steps - (coefs[11,1])*cov3*time_steps - (coefs[12,1])*cov4*time_steps - (coefs[13,1])*cov5*time_steps

    } else if (model_version == "glm_quadratic") {
        squared_time_steps <- time_steps*time_steps
        squared_time_steps_interaction <- time_steps_interaction*time_steps_interaction
        fit_full <- lm(expr ~  cov1 + cov1:time_steps + cov1:squared_time_steps + cov2 + cov2:time_steps + cov2:squared_time_steps + cov3 + cov3:time_steps + cov3:squared_time_steps + cov4 + cov4:time_steps + cov4:squared_time_steps + cov5 + cov5:time_steps + cov5:squared_time_steps + genotype + time_steps + squared_time_steps + genotype_interaction:time_steps_interaction + genotype_interaction:squared_time_steps_interaction)
        fit_null <- lm(expr ~  genotype + time_steps + squared_time_steps + cov1 + cov1:time_steps + cov1:squared_time_steps + cov2 + cov2:time_steps + cov2:squared_time_steps + cov3 + cov3:time_steps + cov3:squared_time_steps + cov4 + cov4:time_steps + cov4:squared_time_steps + cov5 + cov5:time_steps + cov5:squared_time_steps)
        
        af <- anova(fit_full)
        afss <- af$"Sum Sq"
        percent_explained <- afss/sum(afss)*100
        genotype_pve <- percent_explained[6]
        time_pve <- percent_explained[7]
        time_squared_pve <- percent_explained[8]
        interaction_pve <- percent_explained[19]
        interaction_squared_pve <- percent_explained[20]
  		pve_vec <- c(genotype_pve, time_pve, time_squared_pve, interaction_pve, interaction_squared_pve)
    }

    return(pve_vec)
}



input_file <- args[1]
output_file <- args[2]
model_version <- args[3]
covariate_method <- args[4]

f = file(input_file, "r")
stop = FALSE
# Read current line
next_line = readLines(f, n = 1)
sink(output_file)

if (model_version == "glm") {
	# print header to output file
	header <- paste0("rs_id\tensamble_id\tgenotype_pve\ttime_pve\tinteraction_pve\n")
    cat(header)
} else if (model_version == "glm_quadratic") {
	header <- paste0("rs_id\tensamble_id\tgenotype_pve\ttime_pve\ttime_squared_pve\tinteraction_pve\tinteraction_squared_pve\n")
    cat(header)
}

while(!stop) {

	# Read current line

    # Parse the line
        data = strsplit(next_line,'\t')[[1]]    
    	rs_id = data[1]
    	ensamble_id = data[2]


    	time_steps = as.numeric(strsplit(data[3],';')[[1]])
    	genotype = as.numeric(strsplit(data[4],';')[[1]])
    	expr = as.numeric(strsplit(data[5],';')[[1]])
    	cell_lines = as.factor(strsplit(data[6],';')[[1]])
    	num_samp = length(genotype)
    	time_steps_interaction = as.numeric(strsplit(data[3],';')[[1]])
    	genotype_interaction = as.numeric(strsplit(data[4],';')[[1]])
    	pc1 = as.numeric(strsplit(data[7],';')[[1]])
    	pc2 = as.numeric(strsplit(data[8],';')[[1]])
    	pc3 = as.numeric(strsplit(data[9],';')[[1]])
    	pc4 = as.numeric(strsplit(data[10],';')[[1]])
    	pc5 = as.numeric(strsplit(data[11],';')[[1]])
    	pc6 = as.numeric(strsplit(data[12],';')[[1]])
    	pc7 = as.numeric(strsplit(data[13],';')[[1]])

    	if (covariate_method == "none") {
      		lm_results <- run_linear_model(model_version, expr, genotype, genotype_interaction, time_steps, time_steps_interaction, cell_lines)
    	} else if (covariate_method == "pc1") {
      		lm_results <- run_linear_model_one_cov(model_version, expr, genotype, genotype_interaction, time_steps, time_steps_interaction, cell_lines, pc1)
    	} else if (covariate_method == "pc1_2") {
      		lm_results <- run_linear_model_two_cov(model_version, expr, genotype, genotype_interaction, time_steps, time_steps_interaction, cell_lines, pc1, pc2)
    	} else if (covariate_method == "pc1_3") {
      		lm_results <- run_linear_model_three_cov(model_version, expr, genotype, genotype_interaction, time_steps, time_steps_interaction, cell_lines, pc1, pc2, pc3)
    	} else if (covariate_method == "pc1_4") {
      		lm_results <- run_linear_model_four_cov(model_version, expr, genotype, genotype_interaction, time_steps, time_steps_interaction, cell_lines, pc1, pc2, pc3, pc4)
    	} else if (covariate_method == "pc1_5") {
      		lm_results <- run_linear_model_five_cov(model_version, expr, genotype, genotype_interaction, time_steps, time_steps_interaction, cell_lines, pc1, pc2, pc3, pc4, pc5)
    	} else if (covariate_method == "pc1_6") {
      		lm_results <- run_linear_model_six_cov(model_version, expr, genotype, genotype_interaction, time_steps, time_steps_interaction, cell_lines, pc1, pc2, pc3, pc4, pc5, pc6)
    	} 

    	if (model_version == "glm") {
    		liner <- paste0(rs_id, "\t", ensamble_id, "\t", lm_results[1], "\t", lm_results[2], "\t", lm_results[3], "\n")
    		cat(liner)
    	} else if (model_version == "glm_quadratic") {
    		liner <- paste0(rs_id, "\t", ensamble_id, "\t", lm_results[1], "\t", lm_results[2], "\t", lm_results[3], "\t", lm_results[4], "\t", lm_results[5], "\n")
    		cat(liner)
    	}

    	next_line = readLines(f, n = 1)
    	if(length(next_line) == 0) {
        	stop = TRUE
        	close(f)
    	}
}
sink()