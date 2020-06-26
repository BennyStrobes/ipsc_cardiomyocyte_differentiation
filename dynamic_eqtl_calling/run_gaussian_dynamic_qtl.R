args = commandArgs(trailingOnly=TRUE)
library(lme4)
library(lmtest)



run_linear_model <- function(model_version, expr, genotype, genotype_interaction, time_steps, time_steps_interaction, cell_lines) {
    # Run the models
    if (model_version == "glm") {  # Regular linear model
        fit <- lm(expr ~ genotype + time_steps + genotype_interaction:time_steps_interaction)
        pvalue <- summary(fit)$coefficients[4,4]
        coef <- paste(summary(fit)$coefficients[,1],collapse=',')
    } else if (model_version == "glmm") { # Linear mixed model
        fit_full <- lmer(expr ~  genotype + time_steps + genotype_interaction:time_steps_interaction + (1|cell_lines), REML=FALSE)
        fit_null <- lmer(expr ~  genotype + time_steps + (1|cell_lines), REML=FALSE)
        coefs <- data.frame(coef(summary(fit_full)))
        coef <- paste(coefs[,1],collapse=',')
        #tvalue <- coefs[4,3]
        #pvalue <- 2*pt(-abs(tvalue),df=num_samp-1)
        lrt <- anova(fit_null,fit_full)
        pvalue <- lrt[[8]][2]
    } else if (model_version == "glm_quadratic") {
        squared_time_steps <- time_steps*time_steps
        squared_time_steps_interaction <- time_steps_interaction*time_steps_interaction
        fit_full <- lm(expr ~  genotype + time_steps + squared_time_steps + genotype_interaction:time_steps_interaction + genotype_interaction:squared_time_steps_interaction)
        fit_null <- lm(expr ~  genotype + time_steps + squared_time_steps)
        
        # lrt <- anova(fit_null, fit_full)
        # pvalue2 <- lrt[[6]][2]
        coef <- paste(summary(fit_full)$coefficients[,1],collapse=',')

        obj <- lrtest(fit_null, fit_full)
        pvalue <- obj[[5]][2]
    }

    return(list(coef=coef, pvalue=pvalue))
}


run_linear_model_one_cov <- function(model_version, expr, genotype, genotype_interaction, time_steps, time_steps_interaction, cell_lines, cov1) {
    # Run the models
    if (model_version == "glm") {  # Regular linear model
        fit <- lm(expr ~ genotype + time_steps + cov1 + cov1:time_steps + genotype_interaction:time_steps_interaction)
        pvalue <- summary(fit)$coefficients[6,4]
        coef <- paste(summary(fit)$coefficients[,1],collapse=',')
    } else if (model_version == "glmm") { # Linear mixed model
        fit_full <- lmer(expr ~  genotype + time_steps + cov1 + cov1:time_steps + genotype_interaction:time_steps_interaction + (1|cell_lines), REML=FALSE)
        fit_null <- lmer(expr ~  genotype + time_steps + cov1 + cov1:time_steps + (1|cell_lines), REML=FALSE)
        coefs <- data.frame(coef(summary(fit_full)))
        coef <- paste(coefs[,1],collapse=',')
        #tvalue <- coefs[4,3]
        #pvalue <- 2*pt(-abs(tvalue),df=num_samp-1)
        lrt <- anova(fit_null,fit_full)
        pvalue <- lrt[[8]][2]
    } else if (model_version == "glmm_time") { # Linear mixed model
        fit_full <- lmer(expr ~  genotype + time_steps + cov1 + cov1:time_steps + genotype_interaction:time_steps_interaction + (1|factor(time_steps)), REML=FALSE)
        fit_null <- lmer(expr ~  genotype + time_steps + cov1 + cov1:time_steps + (1|factor(time_steps)), REML=FALSE)
        coefs <- data.frame(coef(summary(fit_full)))
        coef <- coefs[6,1]
        #tvalue <- coefs[4,3]
        #pvalue <- 2*pt(-abs(tvalue),df=num_samp-1)
        lrt <- anova(fit_null,fit_full)
        pvalue <- lrt[[8]][2]
    } 
    return(list(coef=coef, pvalue=pvalue))
}

run_linear_model_two_cov <- function(model_version, expr, genotype, genotype_interaction, time_steps, time_steps_interaction, cell_lines, cov1, cov2) {
    # Run the models
    if (model_version == "glm") {  # Regular linear model
        fit <- lm(expr ~ genotype + time_steps + cov1 + cov1:time_steps + cov2 + cov2:time_steps + genotype_interaction:time_steps_interaction)
        pvalue <- summary(fit)$coefficients[8,4]
        coef <- paste(summary(fit)$coefficients[,1],collapse=',')
    } else if (model_version == "glmm") { # Linear mixed model
        fit_full <- lmer(expr ~  genotype + time_steps + cov1 + cov1:time_steps + cov2 + cov2:time_steps + genotype_interaction:time_steps_interaction + (1|cell_lines), REML=FALSE)
        fit_null <- lmer(expr ~  genotype + time_steps + cov1 + cov1:time_steps + cov2 + cov2:time_steps + (1|cell_lines), REML=FALSE)
        coefs <- data.frame(coef(summary(fit_full)))
        coef <- paste(coefs[,1],collapse=',')
        #tvalue <- coefs[4,3]
        #pvalue <- 2*pt(-abs(tvalue),df=num_samp-1)
        lrt <- anova(fit_null,fit_full)
        pvalue <- lrt[[8]][2]
    } else if (model_version == "glmm_time") { # Linear mixed model
        fit_full <- lmer(expr ~  genotype + time_steps + cov1 + cov1:time_steps + cov2 + cov2:time_steps + genotype_interaction:time_steps_interaction + (1|factor(time_steps)), REML=FALSE)
        fit_null <- lmer(expr ~  genotype + time_steps + cov1 + cov1:time_steps + cov2 + cov2:time_steps + (1|factor(time_steps)), REML=FALSE)
        coefs <- data.frame(coef(summary(fit_full)))
        coef <- coefs[8,1]
        #tvalue <- coefs[4,3]
        #pvalue <- 2*pt(-abs(tvalue),df=num_samp-1)
        lrt <- anova(fit_null,fit_full)
        pvalue <- lrt[[8]][2]
    } 
    return(list(coef=coef, pvalue=pvalue))
}

run_linear_model_three_cov <- function(model_version, expr, genotype, genotype_interaction, time_steps, time_steps_interaction, cell_lines, cov1, cov2, cov3) {
    # Run the models
    if (model_version == "glm") {  # Regular linear model
        fit <- lm(expr ~ genotype + time_steps + cov1 + cov1:time_steps + cov2 + cov2:time_steps + cov3 + cov3:time_steps + genotype_interaction:time_steps_interaction)
        pvalue <- summary(fit)$coefficients[10,4]
        coef <- paste(summary(fit)$coefficients[,1],collapse=',')
    } else if (model_version == "glmm") { # Linear mixed model
        fit_full <- lmer(expr ~  genotype + time_steps + cov1 + cov1:time_steps + cov2 + cov2:time_steps + cov3 + cov3:time_steps + genotype_interaction:time_steps_interaction + (1|cell_lines), REML=FALSE)
        fit_null <- lmer(expr ~  genotype + time_steps + cov1 + cov1:time_steps + cov2 + cov2:time_steps + cov3 + cov3:time_steps + (1|cell_lines), REML=FALSE)
        coefs <- data.frame(coef(summary(fit_full)))
        coef <- paste(coefs[,1],collapse=',')
        #tvalue <- coefs[4,3]
        #pvalue <- 2*pt(-abs(tvalue),df=num_samp-1)
        lrt <- anova(fit_null,fit_full)
        pvalue <- lrt[[8]][2]
    } else if (model_version == "glmm_time") { # Linear mixed model
        fit_full <- lmer(expr ~  genotype + time_steps + cov1 + cov1:time_steps + cov2 + cov2:time_steps + cov3 + cov3:time_steps + genotype_interaction:time_steps_interaction + (1|factor(time_steps)), REML=FALSE)
        fit_null <- lmer(expr ~  genotype + time_steps + cov1 + cov1:time_steps + cov2 + cov2:time_steps + cov3 + cov3:time_steps + (1|factor(time_steps)), REML=FALSE)
        coefs <- data.frame(coef(summary(fit_full)))
        coef <- coefs[10,1]
        #tvalue <- coefs[4,3]
        #pvalue <- 2*pt(-abs(tvalue),df=num_samp-1)
        lrt <- anova(fit_null,fit_full)
        pvalue <- lrt[[8]][2]
    } 

    return(list(coef=coef, pvalue=pvalue))
}

run_linear_model_four_cov <- function(model_version, expr, genotype, genotype_interaction, time_steps, time_steps_interaction, cell_lines, cov1, cov2, cov3, cov4) {
    # Run the models
    if (model_version == "glm") {  # Regular linear model
        fit <- lm(expr ~ genotype + time_steps + cov1 + cov1:time_steps + cov2 + cov2:time_steps + cov3 + cov3:time_steps + cov4 + cov4:time_steps + genotype_interaction:time_steps_interaction)
        pvalue <- summary(fit)$coefficients[12,4]
        coef <- paste(summary(fit)$coefficients[,1],collapse=',')
    } else if (model_version == "glmm") { # Linear mixed model
        fit_full <- lmer(expr ~  genotype + time_steps + cov1 + cov1:time_steps + cov2 + cov2:time_steps + cov3 + cov3:time_steps + cov4 + cov4:time_steps + genotype_interaction:time_steps_interaction + (1|cell_lines), REML=FALSE)
        fit_null <- lmer(expr ~  genotype + time_steps + cov1 + cov1:time_steps + cov2 + cov2:time_steps + cov3 + cov3:time_steps + cov4 + cov4:time_steps + (1|cell_lines), REML=FALSE)
        coefs <- data.frame(coef(summary(fit_full)))
        coef <- paste(coefs[,1],collapse=',')
        #tvalue <- coefs[4,3]
        #pvalue <- 2*pt(-abs(tvalue),df=num_samp-1)
        lrt <- anova(fit_null,fit_full)
        pvalue <- lrt[[8]][2]
    } else if (model_version == "glmm_time") { # Linear mixed model
        fit_full <- lmer(expr ~  genotype + time_steps + cov1 + cov1:time_steps + cov2 + cov2:time_steps + cov3 + cov3:time_steps + cov4 + cov4:time_steps + genotype_interaction:time_steps_interaction + (1|factor(time_steps)), REML=FALSE)
        fit_null <- lmer(expr ~  genotype + time_steps + cov1 + cov1:time_steps + cov2 + cov2:time_steps + cov3 + cov3:time_steps + cov4 + cov4:time_steps + (1|factor(time_steps)), REML=FALSE)
        coefs <- data.frame(coef(summary(fit_full)))
        coef <- coefs[12,1]
        #tvalue <- coefs[4,3]
        #pvalue <- 2*pt(-abs(tvalue),df=num_samp-1)
        lrt <- anova(fit_null,fit_full)
        pvalue <- lrt[[8]][2]
    } 

    return(list(coef=coef, pvalue=pvalue))
}

run_linear_model_five_cov <- function(model_version, expr, genotype, genotype_interaction, time_steps, time_steps_interaction, cell_lines, cov1, cov2, cov3, cov4, cov5) {
    # Run the models
    if (model_version == "glm") {  # Regular linear model
        fit <- lm(expr ~ genotype + time_steps + cov1 + cov1:time_steps + cov2 + cov2:time_steps + cov3 + cov3:time_steps + cov4 + cov4:time_steps + cov5 + cov5:time_steps + genotype_interaction:time_steps_interaction)
        pvalue <- summary(fit)$coefficients[14,4]
        coef <- paste(summary(fit)$coefficients[,1],collapse=',')
    } else if (model_version == "glmm") { # Linear mixed model
        fit_full <- lmer(expr ~  genotype + time_steps + cov1 + cov1:time_steps + cov2 + cov2:time_steps + cov3 + cov3:time_steps + cov4 + cov4:time_steps + cov5 + cov5:time_steps + genotype_interaction:time_steps_interaction + (1|cell_lines), REML=FALSE)
        fit_null <- lmer(expr ~  genotype + time_steps + cov1 + cov1:time_steps + cov2 + cov2:time_steps + cov3 + cov3:time_steps + cov4 + cov4:time_steps + cov5 + cov5:time_steps + (1|cell_lines), REML=FALSE)
        coefs <- data.frame(coef(summary(fit_full)))
        coef <- paste(coefs[,1],collapse=',')
        #tvalue <- coefs[4,3]
        #pvalue <- 2*pt(-abs(tvalue),df=num_samp-1)
        lrt <- anova(fit_null,fit_full)
        pvalue <- lrt[[8]][2]
    } else if (model_version == "glmm_time") { # Linear mixed model
        fit_full <- lmer(expr ~  genotype + time_steps + cov1 + cov1:time_steps + cov2 + cov2:time_steps + cov3 + cov3:time_steps + cov4 + cov4:time_steps + cov5 + cov5:time_steps + genotype_interaction:time_steps_interaction + (1|factor(time_steps)), REML=FALSE)
        fit_null <- lmer(expr ~  genotype + time_steps + cov1 + cov1:time_steps + cov2 + cov2:time_steps + cov3 + cov3:time_steps + cov4 + cov4:time_steps + cov5 + cov5:time_steps + (1|factor(time_steps)), REML=FALSE)
        coefs <- data.frame(coef(summary(fit_full)))
        coef <- coefs[14,1]
        #tvalue <- coefs[4,3]
        #pvalue <- 2*pt(-abs(tvalue),df=num_samp-1)
        lrt <- anova(fit_null,fit_full)
        pvalue <- lrt[[8]][2]
    }  else if (model_version == "glm_quadratic") {
        squared_time_steps <- time_steps*time_steps
        squared_time_steps_interaction <- time_steps_interaction*time_steps_interaction
        fit_full <- lm(expr ~  genotype + time_steps + squared_time_steps + cov1 + cov1:time_steps + cov1:squared_time_steps + cov2 + cov2:time_steps + cov2:squared_time_steps + cov3 + cov3:time_steps + cov3:squared_time_steps + cov4 + cov4:time_steps + cov4:squared_time_steps + cov5 + cov5:time_steps + cov5:squared_time_steps + genotype_interaction:time_steps_interaction + genotype_interaction:squared_time_steps_interaction)
        fit_null <- lm(expr ~  genotype + time_steps + squared_time_steps + cov1 + cov1:time_steps + cov1:squared_time_steps + cov2 + cov2:time_steps + cov2:squared_time_steps + cov3 + cov3:time_steps + cov3:squared_time_steps + cov4 + cov4:time_steps + cov4:squared_time_steps + cov5 + cov5:time_steps + cov5:squared_time_steps)
        
        # lrt <- anova(fit_null, fit_full)
        # pvalue2 <- lrt[[6]][2]
        coef <- paste(summary(fit_full)$coefficients[,1],collapse=',')

        obj <- lrtest(fit_null, fit_full)
        pvalue <- obj[[5]][2]

    } else if (model_version == "anova") {
        genotype <- factor(genotype)
        time_steps <- factor(time_steps)

        genotype <- factor(genotype)
        genotype_interaction <- factor(genotype_interaction)
        time_steps_factor <- factor(time_steps)
        time_steps_interaction <- factor(time_steps_interaction)

        fit <- aov(expr ~ genotype + time_steps_factor + genotype_interaction:time_steps_interaction)
        coef <- summary(fit)[[1]][["F value"]][3]
        pvalue <- summary(fit)[[1]][["Pr(>F)"]][3]

    }

    return(list(coef=coef, pvalue=pvalue))
}

run_linear_model_six_cov <- function(model_version, expr, genotype, genotype_interaction, time_steps, time_steps_interaction, cell_lines, cov1, cov2, cov3, cov4, cov5, cov6) {
    # Run the models
    if (model_version == "glm") {  # Regular linear model
        fit <- lm(expr ~ genotype + time_steps + cov1 + cov1:time_steps + cov2 + cov2:time_steps + cov3 + cov3:time_steps + cov4 + cov4:time_steps + cov5 + cov5:time_steps + cov6 + cov6:time_steps + genotype_interaction:time_steps_interaction)
        pvalue <- summary(fit)$coefficients[16,4]
        coef <- paste(summary(fit)$coefficients[,1],collapse=',')
    } else if (model_version == "glmm") { # Linear mixed model
        fit_full <- lmer(expr ~  genotype + time_steps + cov1 + cov1:time_steps + cov2 + cov2:time_steps + cov3 + cov3:time_steps + cov4 + cov4:time_steps + cov5 + cov5:time_steps + cov6 + cov6:time_steps + genotype_interaction:time_steps_interaction + (1|cell_lines), REML=FALSE)
        fit_null <- lmer(expr ~  genotype + time_steps + cov1 + cov1:time_steps + cov2 + cov2:time_steps + cov3 + cov3:time_steps + cov4 + cov4:time_steps + cov5 + cov5:time_steps + cov6 + cov6:time_steps + (1|cell_lines), REML=FALSE)
        coefs <- data.frame(coef(summary(fit_full)))
        coef <- paste(coefs[,1],collapse=',')
        #tvalue <- coefs[4,3]
        #pvalue <- 2*pt(-abs(tvalue),df=num_samp-1)
        lrt <- anova(fit_null,fit_full)
        pvalue <- lrt[[8]][2]
    } else if (model_version == "glmm_time") { # Linear mixed model
        fit_full <- lmer(expr ~  genotype + time_steps + cov1 + cov1:time_steps + cov2 + cov2:time_steps + cov3 + cov3:time_steps + cov4 + cov4:time_steps + cov5 + cov5:time_steps + cov6 + cov6:time_steps + genotype_interaction:time_steps_interaction + (1|factor(time_steps)), REML=FALSE)
        fit_null <- lmer(expr ~  genotype + time_steps + cov1 + cov1:time_steps + cov2 + cov2:time_steps + cov3 + cov3:time_steps + cov4 + cov4:time_steps + cov5 + cov5:time_steps + cov6 + cov6:time_steps + (1|factor(time_steps)), REML=FALSE)
        coefs <- data.frame(coef(summary(fit_full)))
        coef <- coefs[16,1]
        #tvalue <- coefs[4,3]
        #pvalue <- 2*pt(-abs(tvalue),df=num_samp-1)
        lrt <- anova(fit_null,fit_full)
        pvalue <- lrt[[8]][2]
    } 

    return(list(coef=coef, pvalue=pvalue))
}

run_linear_model_seven_cov <- function(model_version, expr, genotype, genotype_interaction, time_steps, time_steps_interaction, cell_lines, cov1, cov2, cov3, cov4, cov5, cov6, cov7) {
    # Run the models
    if (model_version == "glm") {  # Regular linear model
        fit <- lm(expr ~ genotype + time_steps + cov1 + cov1:time_steps + cov2 + cov2:time_steps + cov3 + cov3:time_steps + cov4 + cov4:time_steps + cov5 + cov5:time_steps + cov6 + cov6:time_steps + cov7 + cov7:time_steps + genotype_interaction:time_steps_interaction)
        pvalue <- summary(fit)$coefficients[18,4]
        coef <- paste(summary(fit)$coefficients[,1],collapse=',')
    } else if (model_version == "glmm") { # Linear mixed model
        fit_full <- lmer(expr ~  genotype + time_steps + cov1 + cov1:time_steps + cov2 + cov2:time_steps + cov3 + cov3:time_steps + cov4 + cov4:time_steps + cov5 + cov5:time_steps + cov6 + cov6:time_steps + cov7 + cov7:time_steps + genotype_interaction:time_steps_interaction + (1|cell_lines), REML=FALSE)
        fit_null <- lmer(expr ~  genotype + time_steps + cov1 + cov1:time_steps + cov2 + cov2:time_steps + cov3 + cov3:time_steps + cov4 + cov4:time_steps + cov5 + cov5:time_steps + cov6 + cov6:time_steps + cov7 + cov7:time_steps + (1|cell_lines), REML=FALSE)
        coefs <- data.frame(coef(summary(fit_full)))
        coef <- paste(coefs[,1],collapse=',')
        #tvalue <- coefs[4,3]
        #pvalue <- 2*pt(-abs(tvalue),df=num_samp-1)
        lrt <- anova(fit_null,fit_full)
        pvalue <- lrt[[8]][2]
    } else if (model_version == "glmm_time") { # Linear mixed model
        fit_full <- lm(expr ~  genotype + time_steps + cov1 + cov1:time_steps + cov2 + cov2:time_steps + cov3 + cov3:time_steps + cov4 + cov4:time_steps + cov5 + cov5:time_steps + cov6 + cov6:time_steps + cov7 + cov7:time_steps + genotype_interaction:time_steps_interaction + (1|factor(time_steps)), REML=FALSE)
        fit_null <- lm(expr ~  genotype + time_steps + cov1 + cov1:time_steps + cov2 + cov2:time_steps + cov3 + cov3:time_steps + cov4 + cov4:time_steps + cov5 + cov5:time_steps + cov6 + cov6:time_steps + cov7 + cov7:time_steps + (1|factor(time_steps)), REML=FALSE)
        coefs <- data.frame(coef(summary(fit_full)))
        coef <- coefs[18,1]
        #tvalue <- coefs[4,3]
        #pvalue <- 2*pt(-abs(tvalue),df=num_samp-1)
        lrt <- anova(fit_null,fit_full)
        pvalue <- lrt[[8]][2]
    } 

    return(list(coef=coef, pvalue=pvalue))
}


run_linear_model_eight_cov <- function(model_version, expr, genotype, genotype_interaction, time_steps, time_steps_interaction, cell_lines, cov1, cov2, cov3, cov4, cov5, cov6, cov7,cov8) {
    # Run the models
    if (model_version == "glm") {  # Regular linear model
        fit <- lm(expr ~ genotype + time_steps + cov1 + cov1:time_steps + cov2 + cov2:time_steps + cov3 + cov3:time_steps + cov4 + cov4:time_steps + cov5 + cov5:time_steps + cov6 + cov6:time_steps + cov7 + cov7:time_steps + cov8 + cov8:time_steps + genotype_interaction:time_steps_interaction)
        pvalue <- summary(fit)$coefficients[20,4]
        coef <- paste(summary(fit)$coefficients[,1],collapse=',')
    } else if (model_version == "glmm") { # Linear mixed model
        fit_full <- lmer(expr ~  genotype + time_steps + cov1 + cov1:time_steps + cov2 + cov2:time_steps + cov3 + cov3:time_steps + cov4 + cov4:time_steps + cov5 + cov5:time_steps + cov6 + cov6:time_steps + cov7 + cov7:time_steps + genotype_interaction:time_steps_interaction + (1|cell_lines), REML=FALSE)
        fit_null <- lmer(expr ~  genotype + time_steps + cov1 + cov1:time_steps + cov2 + cov2:time_steps + cov3 + cov3:time_steps + cov4 + cov4:time_steps + cov5 + cov5:time_steps + cov6 + cov6:time_steps + cov7 + cov7:time_steps + (1|cell_lines), REML=FALSE)
        coefs <- data.frame(coef(summary(fit_full)))
        coef <- paste(coefs[,1],collapse=',')
        #tvalue <- coefs[4,3]
        #pvalue <- 2*pt(-abs(tvalue),df=num_samp-1)
        lrt <- anova(fit_null,fit_full)
        pvalue <- lrt[[8]][2]
    } else if (model_version == "glmm_time") { # Linear mixed model
        fit_full <- lm(expr ~  genotype + time_steps + cov1 + cov1:time_steps + cov2 + cov2:time_steps + cov3 + cov3:time_steps + cov4 + cov4:time_steps + cov5 + cov5:time_steps + cov6 + cov6:time_steps + cov7 + cov7:time_steps + genotype_interaction:time_steps_interaction + (1|factor(time_steps)), REML=FALSE)
        fit_null <- lm(expr ~  genotype + time_steps + cov1 + cov1:time_steps + cov2 + cov2:time_steps + cov3 + cov3:time_steps + cov4 + cov4:time_steps + cov5 + cov5:time_steps + cov6 + cov6:time_steps + cov7 + cov7:time_steps + (1|factor(time_steps)), REML=FALSE)
        coefs <- data.frame(coef(summary(fit_full)))
        coef <- coefs[18,1]
        #tvalue <- coefs[4,3]
        #pvalue <- 2*pt(-abs(tvalue),df=num_samp-1)
        lrt <- anova(fit_null,fit_full)
        pvalue <- lrt[[8]][2]
    } 

    return(list(coef=coef, pvalue=pvalue))
}

run_linear_model_nine_cov <- function(model_version, expr, genotype, genotype_interaction, time_steps, time_steps_interaction, cell_lines, cov1, cov2, cov3, cov4, cov5, cov6, cov7,cov8,cov9) {
    # Run the models
    if (model_version == "glm") {  # Regular linear model
        fit <- lm(expr ~ genotype + time_steps + cov1 + cov1:time_steps + cov2 + cov2:time_steps + cov3 + cov3:time_steps + cov4 + cov4:time_steps + cov5 + cov5:time_steps + cov6 + cov6:time_steps + cov7 + cov7:time_steps + cov8 + cov8:time_steps + cov9 + cov9:time_steps + genotype_interaction:time_steps_interaction)
        pvalue <- summary(fit)$coefficients[22,4]
        coef <- paste(summary(fit)$coefficients[,1],collapse=',')
    } else if (model_version == "glmm") { # Linear mixed model
        fit_full <- lmer(expr ~  genotype + time_steps + cov1 + cov1:time_steps + cov2 + cov2:time_steps + cov3 + cov3:time_steps + cov4 + cov4:time_steps + cov5 + cov5:time_steps + cov6 + cov6:time_steps + cov7 + cov7:time_steps + genotype_interaction:time_steps_interaction + (1|cell_lines), REML=FALSE)
        fit_null <- lmer(expr ~  genotype + time_steps + cov1 + cov1:time_steps + cov2 + cov2:time_steps + cov3 + cov3:time_steps + cov4 + cov4:time_steps + cov5 + cov5:time_steps + cov6 + cov6:time_steps + cov7 + cov7:time_steps + (1|cell_lines), REML=FALSE)
        coefs <- data.frame(coef(summary(fit_full)))
        coef <- paste(coefs[,1],collapse=',')
        #tvalue <- coefs[4,3]
        #pvalue <- 2*pt(-abs(tvalue),df=num_samp-1)
        lrt <- anova(fit_null,fit_full)
        pvalue <- lrt[[8]][2]
    } else if (model_version == "glmm_time") { # Linear mixed model
        fit_full <- lm(expr ~  genotype + time_steps + cov1 + cov1:time_steps + cov2 + cov2:time_steps + cov3 + cov3:time_steps + cov4 + cov4:time_steps + cov5 + cov5:time_steps + cov6 + cov6:time_steps + cov7 + cov7:time_steps + genotype_interaction:time_steps_interaction + (1|factor(time_steps)), REML=FALSE)
        fit_null <- lm(expr ~  genotype + time_steps + cov1 + cov1:time_steps + cov2 + cov2:time_steps + cov3 + cov3:time_steps + cov4 + cov4:time_steps + cov5 + cov5:time_steps + cov6 + cov6:time_steps + cov7 + cov7:time_steps + (1|factor(time_steps)), REML=FALSE)
        coefs <- data.frame(coef(summary(fit_full)))
        coef <- coefs[18,1]
        #tvalue <- coefs[4,3]
        #pvalue <- 2*pt(-abs(tvalue),df=num_samp-1)
        lrt <- anova(fit_null,fit_full)
        pvalue <- lrt[[8]][2]
    } 

    return(list(coef=coef, pvalue=pvalue))
}


run_linear_model_ten_cov <- function(model_version, expr, genotype, genotype_interaction, time_steps, time_steps_interaction, cell_lines, cov1, cov2, cov3, cov4, cov5, cov6, cov7,cov8,cov9,cov10) {
    # Run the models
    if (model_version == "glm") {  # Regular linear model
        fit <- lm(expr ~ genotype + time_steps + cov1 + cov1:time_steps + cov2 + cov2:time_steps + cov3 + cov3:time_steps + cov4 + cov4:time_steps + cov5 + cov5:time_steps + cov6 + cov6:time_steps + cov7 + cov7:time_steps + cov8 + cov8:time_steps + cov9 + cov9:time_steps + cov10 + cov10:time_steps + genotype_interaction:time_steps_interaction)
        pvalue <- summary(fit)$coefficients[24,4]
        coef <- paste(summary(fit)$coefficients[,1],collapse=',')
    } else if (model_version == "glmm") { # Linear mixed model
        fit_full <- lmer(expr ~  genotype + time_steps + cov1 + cov1:time_steps + cov2 + cov2:time_steps + cov3 + cov3:time_steps + cov4 + cov4:time_steps + cov5 + cov5:time_steps + cov6 + cov6:time_steps + cov7 + cov7:time_steps + genotype_interaction:time_steps_interaction + (1|cell_lines), REML=FALSE)
        fit_null <- lmer(expr ~  genotype + time_steps + cov1 + cov1:time_steps + cov2 + cov2:time_steps + cov3 + cov3:time_steps + cov4 + cov4:time_steps + cov5 + cov5:time_steps + cov6 + cov6:time_steps + cov7 + cov7:time_steps + (1|cell_lines), REML=FALSE)
        coefs <- data.frame(coef(summary(fit_full)))
        coef <- paste(coefs[,1],collapse=',')
        #tvalue <- coefs[4,3]
        #pvalue <- 2*pt(-abs(tvalue),df=num_samp-1)
        lrt <- anova(fit_null,fit_full)
        pvalue <- lrt[[8]][2]
    } else if (model_version == "glmm_time") { # Linear mixed model
        fit_full <- lm(expr ~  genotype + time_steps + cov1 + cov1:time_steps + cov2 + cov2:time_steps + cov3 + cov3:time_steps + cov4 + cov4:time_steps + cov5 + cov5:time_steps + cov6 + cov6:time_steps + cov7 + cov7:time_steps + genotype_interaction:time_steps_interaction + (1|factor(time_steps)), REML=FALSE)
        fit_null <- lm(expr ~  genotype + time_steps + cov1 + cov1:time_steps + cov2 + cov2:time_steps + cov3 + cov3:time_steps + cov4 + cov4:time_steps + cov5 + cov5:time_steps + cov6 + cov6:time_steps + cov7 + cov7:time_steps + (1|factor(time_steps)), REML=FALSE)
        coefs <- data.frame(coef(summary(fit_full)))
        coef <- coefs[18,1]
        #tvalue <- coefs[4,3]
        #pvalue <- 2*pt(-abs(tvalue),df=num_samp-1)
        lrt <- anova(fit_null,fit_full)
        pvalue <- lrt[[8]][2]
    } 

    return(list(coef=coef, pvalue=pvalue))
}



null_response <- function(){
    return(list(coef=0,pvalue=1))
}

#####################
# Command line args
#####################

input_data_file = args[1]
output_file = args[2]
model_version = args[3]
permute = args[4]
covariate_method = args[5]
job_number = as.numeric(args[6])
num_jobs = as.numeric(args[7])
num_lines_string = args[8]




# Extract total number of lines in file
total_lines = as.numeric(strsplit(num_lines_string," ")[[1]][1])

# Determine number of lines each parrallelized job will complete
lines_per_job = ceiling(total_lines/num_jobs)
start_num = job_number*lines_per_job
end_num = (job_number + 1)*lines_per_job


# Stream input file
stop = FALSE
count = 0
f = file(input_data_file, "r")

next_line = readLines(f, n = 1)

sink(output_file)

while(!stop) {
    # Only consider lines between start_num and end_num (for parallelization purposes)
    if (count >= start_num & count < end_num) {
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
    pc8 = as.numeric(strsplit(data[14],';')[[1]])
    pc9 = as.numeric(strsplit(data[15],';')[[1]])
    pc10 = as.numeric(strsplit(data[16],';')[[1]])

    # Permute the data
    if (permute == "True") {
        time_steps_interaction <- sample(time_steps_interaction)
    } 
    # Run the model in try-catch statement
    tryCatch(
    {
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
        } else if (covariate_method == "pc1_7") {
            lm_results <- run_linear_model_seven_cov(model_version, expr, genotype, genotype_interaction, time_steps, time_steps_interaction, cell_lines, pc1, pc2, pc3, pc4, pc5, pc6, pc7)
        } else if (covariate_method == "pc1_8") {
            lm_results <- run_linear_model_eight_cov(model_version, expr, genotype, genotype_interaction, time_steps, time_steps_interaction, cell_lines, pc1, pc2, pc3, pc4, pc5, pc6, pc7, pc8)
        } else if (covariate_method == "pc1_9") {
            lm_results <- run_linear_model_nine_cov(model_version, expr, genotype, genotype_interaction, time_steps, time_steps_interaction, cell_lines, pc1, pc2, pc3, pc4, pc5, pc6, pc7, pc8, pc9)
        } else if (covariate_method == "pc1_10") {
            lm_results <- run_linear_model_ten_cov(model_version, expr, genotype, genotype_interaction, time_steps, time_steps_interaction, cell_lines, pc1, pc2, pc3, pc4, pc5, pc6, pc7, pc8, pc9, pc10)
        }

        # print result to output file!!
        new_line <- paste0(rs_id, "\t", ensamble_id ,"\t",lm_results$coef,"\t", lm_results$pvalue,"\n")
        cat(new_line)
    },
    error = function(e){
        new_line <- paste0(rs_id, "\t", ensamble_id,"\t",0.0,"\t", 1.0,"\n")
        cat(new_line)
    }
    )
    }
    count = count + 1

    next_line = readLines(f, n = 1)

    if(length(next_line) == 0) {
        stop = TRUE
        close(f)
    }

}
# close output file handle
sink()