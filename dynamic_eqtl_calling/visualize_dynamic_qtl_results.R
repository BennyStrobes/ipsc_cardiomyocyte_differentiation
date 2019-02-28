args = commandArgs(trailingOnly=TRUE)
library(ggplot2)
library(ggthemes)
library(cowplot)
library(reshape)
library(lme4)
library(lmtest)


make_te_plot <- function(ensamble_id, rs_id, te_df, lower, upper, ref_allele,alt_allele) {
    homo_ref <- paste0(ref_allele,ref_allele)
    het <- paste0(ref_allele,alt_allele)
    homo_alt <- paste0(alt_allele,alt_allele)
    for (index in 1:length(te_df$genotype)) {
        if (te_df$genotype[index] == 0) {
            te_df$genotype[index] <- homo_ref
        }
        if (te_df$genotype[index] == 1) {
            te_df$genotype[index] <- het
        }
        if (te_df$genotype[index] == 2) {
            te_df$genotype[index] <- homo_alt
        }
    }
    # Change box plot line colors by groups
    p <- ggplot(te_df, aes(x=time_step, y=gene_counts, fill=factor(genotype, levels=c(homo_ref,het,homo_alt)))) + geom_boxplot(width=.6, outlier.shape=NA) + scale_fill_manual(values=c("#56B4E9","#E69F00","plum2")) +
       theme(text = element_text(size=8),axis.text=element_text(size=7), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.text = element_text(size=7), legend.title = element_text(size=8))  +
       labs(x = "Day", y = paste0(ensamble_id), fill=rs_id, title=paste0(rs_id, "-", ensamble_id)) +
        theme(legend.position="bottom") + ylim(lower,upper)
    return(p)
}

make_te_plot_alternate <- function(ensamble_id, rs_id, te_df) {
    # Change box plot line colors by groups

   # line_plot <- ggplot(te_df, aes(x=time, y=expression, group=cell_line)) + geom_line(aes(color=cell_line)) +
     #           geom_point(aes(color=cell_line)) 
    p <-ggplot(te_df,aes(x=time_step,y=gene_counts,group=cell_line, color=genotype)) + geom_line()  + scale_color_manual(values=c("plum2", "#56B4E9","#E69F00")) +
       theme(text = element_text(size=12),plot.title = element_text(size=12), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black")) +
       labs(x = "Time Step", y = "Expression", title=paste0(rs_id, "-", ensamble_id)) +
        theme(legend.position="bottom") + theme(legend.text = element_text(size=10))+ theme(legend.title = element_text(size=10))
    return(p)
}

run_linear_model_five_cov <- function(model_version, expr, genotype, genotype_interaction, time_steps, time_steps_interaction, cell_lines, cov1, cov2, cov3, cov4, cov5) {
    # Run the models
    if (model_version == "glm") {  # Regular linear model
        fit <- lm(expr ~ genotype + time_steps + cov1 + cov1:time_steps + cov2 + cov2:time_steps + cov3 + cov3:time_steps + cov4 + cov4:time_steps + cov5 + cov5:time_steps + genotype_interaction:time_steps_interaction)
        pvalue <- summary(fit)$coefficients[14,4]
        coef <- summary(fit)$coefficients[14,1]
        coefs <- summary(fit)$coefficients
        corrected_expr <- expr - (coefs[4,1])*cov1 - (coefs[5,1])*cov2 - (coefs[6,1])*cov3 - (coefs[7,1])*cov4 - (coefs[8,1])*cov5 - (coefs[9,1])*cov1*time_steps - (coefs[10,1])*cov2*time_steps - (coefs[11,1])*cov3*time_steps - (coefs[12,1])*cov4*time_steps - (coefs[13,1])*cov5*time_steps

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

    } else if (model_version == "glmm_time") { # Linear mixed model
        fit_full <- lmer(expr ~  genotype + time_steps + cov1 + cov1:time_steps + cov2 + cov2:time_steps + cov3 + cov3:time_steps + cov4 + cov4:time_steps + cov5 + cov5:time_steps + genotype_interaction:time_steps_interaction + (1|factor(time_steps)), REML=FALSE)
        fit_null <- lmer(expr ~  genotype + time_steps + cov1 + cov1:time_steps + cov2 + cov2:time_steps + cov3 + cov3:time_steps + cov4 + cov4:time_steps + cov5 + cov5:time_steps + (1|factor(time_steps)), REML=FALSE)
        coefs <- data.frame(coef(summary(fit_full)))
        coef <- coefs[14,1]
        #tvalue <- coefs[4,3]
        #pvalue <- 2*pt(-abs(tvalue),df=num_samp-1)
        lrt <- anova(fit_null,fit_full)
        pvalue <- lrt[[8]][2]
    } else if (model_version == "glm_quadratic") {
        squared_time_steps <- time_steps*time_steps
        squared_time_steps_interaction <- time_steps_interaction*time_steps_interaction
        fit_full <- lm(expr ~  genotype + time_steps + squared_time_steps + cov1 + cov1:time_steps + cov1:squared_time_steps + cov2 + cov2:time_steps + cov2:squared_time_steps + cov3 + cov3:time_steps + cov3:squared_time_steps + cov4 + cov4:time_steps + cov4:squared_time_steps + cov5 + cov5:time_steps + cov5:squared_time_steps + genotype_interaction:time_steps_interaction + genotype_interaction:squared_time_steps_interaction)
        fit_null <- lm(expr ~  genotype + time_steps + squared_time_steps + cov1 + cov1:time_steps + cov1:squared_time_steps + cov2 + cov2:time_steps + cov2:squared_time_steps + cov3 + cov3:time_steps + cov3:squared_time_steps + cov4 + cov4:time_steps + cov4:squared_time_steps + cov5 + cov5:time_steps + cov5:squared_time_steps)
        
        # lrt <- anova(fit_null, fit_full)
        # pvalue2 <- lrt[[6]][2]
        coef <- paste(summary(fit_full)$coefficients[,1],collapse=',')

        obj <- lrtest(fit_null, fit_full)
        pvalue <- obj[[5]][2]

        coefs <- summary(fit_full)$coefficients

        corrected_expr <- expr - (coefs[5,1])*cov1 - (coefs[6,1])*cov2 - (coefs[7,1])*cov3 - (coefs[8,1])*cov4 - (coefs[9,1])*cov5 - (coefs[10,1])*cov1*time_steps - (coefs[11,1])*cov1*squared_time_steps - (coefs[12,1])*cov2*time_steps - (coefs[13,1])*cov2*squared_time_steps - (coefs[14,1])*cov3*time_steps - (coefs[15,1])*cov3*squared_time_steps - (coefs[16,1])*cov4*time_steps - (coefs[17,1])*cov4*squared_time_steps - (coefs[18,1])*cov5*time_steps - (coefs[19,1])*cov5*squared_time_steps 


    }

    return(list(coef=coef, pvalue=pvalue, expr=corrected_expr))
}

make_dynamic_qtl_plot <- function(input_file, desired_rsid, desired_ensamble_id, gene_name, model_version, covariate_method, lower, upper, ref_allele,alt_allele) {
	f = file(input_file, "r")
	stop = FALSE
	# Read current line
    next_line = readLines(f, n = 1)

	while(!stop) {

    	# Read current line

    	# Parse the line
    	if (startsWith(next_line, desired_rsid) == TRUE) {
    	data = strsplit(next_line,'\t')[[1]]    
    	rs_id = data[1]
    	ensamble_id = data[2]


    	if (rs_id == desired_rsid & ensamble_id == desired_ensamble_id) {
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


    		corrected_expr <- lm_results$expr
    		pvalue <- lm_results$pvalue
    		beta <- lm_results$coef

    		te_df <- data.frame(gene_counts=corrected_expr, time_step=factor(time_steps), genotype=round(genotype))

    		te_plot <- make_te_plot(gene_name, rs_id, te_df, lower, upper, ref_allele, alt_allele)
    	}
    	}
    	next_line = readLines(f, n = 1)
    	if(length(next_line) == 0) {
        	stop = TRUE
        	close(f)
    	}
	}
    return(te_plot)
}

load_in_odds_ratios <- function(file_name, adding_constant) {
    aa <- read.table(file_name,header=TRUE)

    real_overlaps <- as.numeric(aa$real_overlaps) + adding_constant
    real_misses <- as.numeric(aa$real_misses) + adding_constant
    perm_overlaps <- as.numeric(aa$perm_overlaps) + adding_constant
    perm_misses <- as.numeric(aa$perm_misses) + adding_constant
    odds_ratios <- (real_overlaps/real_misses)/(perm_overlaps/perm_misses)
    return(odds_ratios)
}
cre_enrichments_boxplot <- function(cre, cell_lines, hits_versions, adding_constant, num_permutations, threshold, input_root, title) {
    odds_ratios <- c()
    roadmap_cell_types <- c()
    dynamic_qtl_versions <- c()
    for (cell_line_counter in 1:length(cell_lines)) {
        for (hits_version_counter in 1:length(hits_versions)) {
            cell_line <- cell_lines[cell_line_counter]
            hits_version <- hits_versions[hits_version_counter]
            input_file <- paste0(input_root, cre, "_", cell_line, "_cell_lines_", hits_version,"_hits_", num_permutations, "_",threshold, ".txt" )
            or <- load_in_odds_ratios(input_file, adding_constant)
            odds_ratios <- c(odds_ratios, or)
            if (cell_line == "ipsc_only") {
            	cell_line <- "iPSC"
            }
            if (cell_line == "heart_only") {
            	cell_line <- "heart"
            }
            if (hits_version == "early_time_step") {
            	hits_version <- "early"
            }
            if (hits_version == "late_time_step") {
            	hits_version <- "late"
            }
            roadmap_cell_types <- c(roadmap_cell_types, rep(cell_line, length(or)))
            dynamic_qtl_versions <- c(dynamic_qtl_versions, rep(hits_version, length(or)))
        }
    }
    df <- data.frame(odds_ratio=log10(odds_ratios),roadmap_cell_type=factor(paste0(roadmap_cell_types, " enhancer"),levels=c("iPSC enhancer", "heart enhancer")), qtl_version=factor(dynamic_qtl_versions))
    # PLOT
    boxplot <- ggplot(df, aes(x=roadmap_cell_type,y=odds_ratio,fill=qtl_version)) + geom_boxplot(outlier.shape=NA) + labs(x = "Roadmap Cell Type", y = expression(log[10]("odds ratio")), fill="") + scale_fill_manual(values=c("darkgrey", "firebrick"))
    boxplot <- boxplot + theme(text = element_text(size=8),axis.text=element_text(size=7), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.text = element_text(size=7), legend.title = element_text(size=8)) 
    boxplot <- boxplot + geom_hline(yintercept = 0.0) 
    if (title != "none") {
        boxplot <- boxplot + ggtitle(title) + theme(plot.title =element_text(size=8, face='plain'))
    }
	return(boxplot)
}

make_manhatten_plot <- function(qtl_df, miny, maxy, coloring,variant_pos, reverse_boolean, variant_chrom) {
    print(variant_pos)
    qtl_df$pvalue <- -log10(qtl_df$pvalue)
    p <- ggplot(qtl_df, aes(x=variant_position,y=pvalue,color=phenotype))+ geom_vline(xintercept = variant_pos,size=.15) + geom_point(size=.25) + scale_color_manual(values=c(coloring)) +
        theme(text = element_text(size=8),axis.text=element_text(size=7), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.text = element_text(size=7), legend.title = element_text(size=8)) +
        labs(x = "", y = expression(-log[10]("p-value")),color="") + xlim(miny,maxy) + theme(legend.title=element_blank()) +
        scale_x_continuous(breaks=c(75460000,75500000,75540000), labels=c("", "", ""))
    if (reverse_boolean == TRUE) {
        p <- p + scale_y_reverse() + scale_x_continuous(position="top") + labs(x="")
    }
    return(p)
}

make_joint_manhatten_plot <- function(qtl_df,miny,maxy, variant_pos, variant_chrom) {
    qtl_df$pvalue <- -log10(qtl_df$pvalue)
    p <- ggplot(qtl_df, aes(x=variant_position,y=pvalue,color=factor(phenotype))) + geom_vline(xintercept = variant_pos,size=.15)+ geom_point(size=.25) +scale_color_manual(values=c("dodgerblue3")) +
        theme(text = element_text(size=8),axis.text=element_text(size=7), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.text = element_text(size=7), legend.title = element_text(size=8)) +
        labs(x = "", y = expression(-log[10]("p-value")),color="") + xlim(miny,maxy)
    p <- p + scale_y_reverse()  + labs(x="") + theme(legend.title=element_blank())
    p <- p + scale_x_continuous(position="top",breaks=c(75460000,75500000,75540000), labels=c("75.46", "75.50", "75.54"))
    return(p)

}


make_joint_miami_plot_one_phen <- function(dynamic_qtl_file_name, gwas_file_name1, gwas_file_name2, gwas_file_name3, phenotype1, phenotype2, phenotype3, rs_id, ensamble_id, gene_symbol,variant_pos,variant_chrom) {
    dynamic_qtls <- read.table(dynamic_qtl_file_name,header=TRUE)
    gwas_stats1 <- read.table(gwas_file_name1, header=TRUE)



    legible_pheno_name1 = paste0("GWAS: ",paste(strsplit(phenotype1,"_")[[1]], collapse=" "),"  ")
    legible_pheno_name1 <-"GWAS: BMI"
    gwas_stats1$phenotype <- rep(legible_pheno_name1,dim(gwas_stats1)[1])


    gwas_stats <- rbind(gwas_stats1)


    dynamic_qtls$phenotype <- rep(paste0("Dynamic eQTL for ", gene_symbol), dim(dynamic_qtls)[1])

    miny <- min(min(dynamic_qtls$variant_position), min(gwas_stats1$variant_position))
    maxy <- max(max(dynamic_qtls$variant_position), max(gwas_stats1$variant_position))

    dynamic_qtl_title <- paste0("Dynamic eqtl for ", gene_symbol)
    dynamic_qtl_plot <- make_manhatten_plot(dynamic_qtls, miny, maxy, "chartreuse4", variant_pos,FALSE, variant_chrom) + theme(text = element_text(size=8),axis.text=element_text(size=7), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.text = element_text(size=7), legend.title = element_text(size=8)) 
    #legible_pheno_name = paste(strsplit(phenotype_name,"_")[[1]], collapse=" ")
    #gwas_title <- paste0("GWAS: ", legible_pheno_name)
    gwas_qtl_plot <- make_joint_manhatten_plot(gwas_stats,miny,maxy, variant_pos, variant_chrom) + theme(text = element_text(size=8),axis.text=element_text(size=7), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.text = element_text(size=7), legend.title = element_text(size=8)) 


    combined <- ggdraw() + 
        draw_plot(dynamic_qtl_plot+ theme(legend.position="none") ,.05,.43,.95,.787) + 
        draw_plot(gwas_qtl_plot + theme(legend.position="none") ,.05,0,.95,.77) 
    return(combined)

}

make_miami_plot_one_phen <- function(gwas_overlap_dir) {
    rs_id = 'rs28818910'
    ensamble_id = 'ENSG00000167173'
    gene_symbol <- "C15orf39"
    variant_chrom = '15'
    variant_pos = 75440669
    phenotype1 = 'UKB_21001_Body_mass_index_BMI'
    phenotype2 = 'Astle_et_al_2016_Red_blood_cell_count'
    phenotype3 = 'UKB_23099_Body_fat_percentage'

    dynamic_qtl_file_name <- paste0(gwas_overlap_dir,"gaussian_dynamic_qtl_input_file_environmental_variable_time_steps_genotype_version_dosage_model_type_glm_quadratic_covariate_method_pc1_5_",rs_id,"_",ensamble_id,"_",phenotype1,"_nearby_dynamic_qtl_pvalues.txt")
    gwas_file_name1 <- paste0(gwas_overlap_dir,"gaussian_dynamic_qtl_input_file_environmental_variable_time_steps_genotype_version_dosage_model_type_glm_quadratic_covariate_method_pc1_5_",rs_id,"_",ensamble_id,"_",phenotype1,"_nearby_gwas_pvalues.txt")
    gwas_file_name2 <- paste0(gwas_overlap_dir,"gaussian_dynamic_qtl_input_file_environmental_variable_time_steps_genotype_version_dosage_model_type_glm_quadratic_covariate_method_pc1_5_",rs_id,"_",ensamble_id,"_",phenotype2,"_nearby_gwas_pvalues.txt")
    gwas_file_name3 <- paste0(gwas_overlap_dir,"gaussian_dynamic_qtl_input_file_environmental_variable_time_steps_genotype_version_dosage_model_type_glm_quadratic_covariate_method_pc1_5_",rs_id,"_",ensamble_id,"_",phenotype3,"_nearby_gwas_pvalues.txt")

    miami_plot <- make_joint_miami_plot_one_phen(dynamic_qtl_file_name, gwas_file_name1, gwas_file_name2, gwas_file_name3, phenotype1, phenotype2, phenotype3, rs_id, ensamble_id, gene_symbol,variant_pos,variant_chrom)
    return(miami_plot)
}

make_joint_manhatten_plot_v2 <- function(qtl_df,miny,maxy, variant_pos, variant_chrom) {
    qtl_df$pvalue <- -log10(qtl_df$pvalue)
    p <- ggplot(qtl_df, aes(x=variant_position,y=pvalue,color=factor(phenotype))) + geom_vline(xintercept = variant_pos,size=.3)+ geom_point(size=.95) + scale_color_manual(values=c("salmon2","dodgerblue3","purple3")) +
        theme(text = element_text(size=8),axis.text=element_text(size=7), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.text = element_text(size=7), legend.title = element_text(size=8)) +
        labs(x = paste0("Position on chromosome ", variant_chrom), y = expression(-log[10]("p-value")),color="") + xlim(miny,maxy)
    p <- p + scale_y_reverse() + scale_x_continuous(position="top") + labs(x="") + theme(legend.title=element_blank())
    return(p)

}

make_manhatten_plot_v2 <- function(qtl_df, miny, maxy, coloring,variant_pos, reverse_boolean, variant_chrom) {
    print(variant_pos)
    qtl_df$pvalue <- -log10(qtl_df$pvalue)
    p <- ggplot(qtl_df, aes(x=variant_position,y=pvalue,color=phenotype))+ geom_vline(xintercept = variant_pos,size=.3) + geom_point(size=.95) + scale_color_manual(values=c(coloring)) +
        theme(text = element_text(size=8),axis.text=element_text(size=7), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.text = element_text(size=7), legend.title = element_text(size=8)) +
        labs(x = paste0("Position on chromosome ", variant_chrom), y = expression(-log[10]("p-value")),color="") + xlim(miny,maxy) + theme(legend.title=element_blank())
    if (reverse_boolean == TRUE) {
        p <- p + scale_y_reverse() + scale_x_continuous(position="top") + labs(x="")
    }
    return(p)
}

make_joint_miami_plot <- function(dynamic_qtl_file_name, gwas_file_name1, gwas_file_name2, gwas_file_name3, phenotype1, phenotype2, phenotype3, rs_id, ensamble_id, gene_symbol,variant_pos,variant_chrom) {
    dynamic_qtls <- read.table(dynamic_qtl_file_name,header=TRUE)
    gwas_stats1 <- read.table(gwas_file_name1, header=TRUE)
    gwas_stats2 <- read.table(gwas_file_name2, header=TRUE)
    gwas_stats3 <- read.table(gwas_file_name3, header=TRUE)


    legible_pheno_name1 = paste0("GWAS: ",paste(strsplit(phenotype1,"_")[[1]], collapse=" "),"  ")
    gwas_stats1$phenotype <- rep(legible_pheno_name1,dim(gwas_stats1)[1])
    legible_pheno_name2 = paste0("GWAS: ",paste(strsplit(phenotype2,"_")[[1]], collapse=" "),"  ")
    gwas_stats2$phenotype <- rep(legible_pheno_name2,dim(gwas_stats2)[1])
    legible_pheno_name3 = paste0("GWAS: ",paste(strsplit(phenotype3,"_")[[1]], collapse=" "),"  ")
    gwas_stats3$phenotype <- rep(legible_pheno_name3,dim(gwas_stats3)[1])

    gwas_stats <- rbind(gwas_stats1,gwas_stats2, gwas_stats3)


    dynamic_qtls$phenotype <- rep(paste0("Dynamic eQTL for ", gene_symbol), dim(dynamic_qtls)[1])

    miny <- min(min(dynamic_qtls$variant_position), min(gwas_stats1$variant_position), min(gwas_stats2$variant_position), min(gwas_stats3$variant_position))
    maxy <- max(max(dynamic_qtls$variant_position), max(gwas_stats1$variant_position), max(gwas_stats2$variant_position), max(gwas_stats3$variant_position))

    dynamic_qtl_title <- paste0("Dynamic eqtl for ", gene_symbol)
    dynamic_qtl_plot <- make_manhatten_plot_v2(dynamic_qtls, miny, maxy, "chartreuse4", variant_pos,FALSE, variant_chrom)

    #legible_pheno_name = paste(strsplit(phenotype_name,"_")[[1]], collapse=" ")
    #gwas_title <- paste0("GWAS: ", legible_pheno_name)
    gwas_qtl_plot <- make_joint_manhatten_plot_v2(gwas_stats,miny,maxy, variant_pos, variant_chrom)

    gwas_legend <- get_legend(gwas_qtl_plot+ theme(legend.background = element_rect(size=0.5, linetype="solid", colour ="black")))

    dynamic_qtl_legend <- get_legend(dynamic_qtl_plot+ theme(legend.background = element_rect(size=0.5, linetype="solid", colour ="black")))

    combined <- ggdraw() + 
        draw_plot(dynamic_qtl_plot+ theme(legend.position="none"),0,.47,1,.54) + 
        draw_plot(gwas_qtl_plot + theme(legend.position="none"),0,0,1,.53) +
        draw_plot(gwas_legend,.33,-.42,1,1) +
        draw_plot(dynamic_qtl_legend,.405,.45,1,1) + 
        draw_plot_label(c("A","B"),c(.02,.02),c(1,.48),size=12)
        #draw_plot_label(c(dynamic_qtl_title, gwas_title), c(.24,.2),c(.98,.1),size=10)
       # draw_plot_label(c("a","b","c", "d"),c(.01,.01,.66,.01),c(1,.68,.68,.36),size=12)
    ggsave(combined, file=output_file, width=7.2, height=5.5,units="in")

}

make_miami_plot <- function(gwas_overlap_dir, output_file) {
    rs_id = 'rs28818910'
    ensamble_id = 'ENSG00000167173'
    gene_symbol <- "C15orf39"
    variant_chrom = '15'
    variant_pos = 75440669
    phenotype1 = 'UKB_21001_Body_mass_index_BMI'
    phenotype2 = 'Astle_et_al_2016_Red_blood_cell_count'
    phenotype3 = 'UKB_23099_Body_fat_percentage'

    dynamic_qtl_file_name <- paste0(gwas_overlap_dir,"gaussian_dynamic_qtl_input_file_environmental_variable_time_steps_genotype_version_dosage_model_type_glm_quadratic_covariate_method_pc1_5_",rs_id,"_",ensamble_id,"_",phenotype1,"_nearby_dynamic_qtl_pvalues.txt")
    gwas_file_name1 <- paste0(gwas_overlap_dir,"gaussian_dynamic_qtl_input_file_environmental_variable_time_steps_genotype_version_dosage_model_type_glm_quadratic_covariate_method_pc1_5_",rs_id,"_",ensamble_id,"_",phenotype1,"_nearby_gwas_pvalues.txt")
    gwas_file_name2 <- paste0(gwas_overlap_dir,"gaussian_dynamic_qtl_input_file_environmental_variable_time_steps_genotype_version_dosage_model_type_glm_quadratic_covariate_method_pc1_5_",rs_id,"_",ensamble_id,"_",phenotype2,"_nearby_gwas_pvalues.txt")
    gwas_file_name3 <- paste0(gwas_overlap_dir,"gaussian_dynamic_qtl_input_file_environmental_variable_time_steps_genotype_version_dosage_model_type_glm_quadratic_covariate_method_pc1_5_",rs_id,"_",ensamble_id,"_",phenotype3,"_nearby_gwas_pvalues.txt")

    miami_plot <- make_joint_miami_plot(dynamic_qtl_file_name, gwas_file_name1, gwas_file_name2, gwas_file_name3, phenotype1, phenotype2, phenotype3, rs_id, ensamble_id, gene_symbol,variant_pos,variant_chrom)
    ggsave(miami_plot, file=output_file, width=7.2, height=5.5,units="in")

}

cre_enrichment_over_range_of_pcs_boxplot <- function(tissue_specific_chrom_hmm_enrichment_dir, output_file) {
    cell_lines <- c("ipsc_only", "heart_only")
    hits_versions <- c("early_time_step", "late_time_step")
    cre <- "enhancer"
    adding_constant <- 1
    num_permutations="1000"
    threshold = "1.0"

    #Range of PCs
    input_root <- paste0(tissue_specific_chrom_hmm_enrichment_dir,"gaussian_dynamic_qtl_input_file_environmental_variable_time_steps_genotype_version_dosage_model_type_glm_covariate_method_none_")
    plot_0 <- cre_enrichments_boxplot(cre, cell_lines, hits_versions, adding_constant, num_permutations, threshold, input_root, "0 PC")
    legend <- get_legend(plot_0 + theme(legend.position="bottom"))
    plot_0 <- plot_0 + theme(legend.position="none")

    input_root <- paste0(tissue_specific_chrom_hmm_enrichment_dir,"gaussian_dynamic_qtl_input_file_environmental_variable_time_steps_genotype_version_dosage_model_type_glm_covariate_method_pc1_")
    plot_1 <- cre_enrichments_boxplot(cre, cell_lines, hits_versions, adding_constant, num_permutations, threshold, input_root, "1 PC") + theme(legend.position="none")

    input_root <- paste0(tissue_specific_chrom_hmm_enrichment_dir,"gaussian_dynamic_qtl_input_file_environmental_variable_time_steps_genotype_version_dosage_model_type_glm_covariate_method_pc1_2_")
    plot_2 <- cre_enrichments_boxplot(cre, cell_lines, hits_versions, adding_constant, num_permutations, threshold, input_root, "2 PC") + theme(legend.position="none")

    input_root <- paste0(tissue_specific_chrom_hmm_enrichment_dir,"gaussian_dynamic_qtl_input_file_environmental_variable_time_steps_genotype_version_dosage_model_type_glm_covariate_method_pc1_3_")
    plot_3 <- cre_enrichments_boxplot(cre, cell_lines, hits_versions, adding_constant, num_permutations, threshold, input_root, "3 PC") + theme(legend.position="none")

    input_root <- paste0(tissue_specific_chrom_hmm_enrichment_dir,"gaussian_dynamic_qtl_input_file_environmental_variable_time_steps_genotype_version_dosage_model_type_glm_covariate_method_pc1_4_")
    plot_4 <- cre_enrichments_boxplot(cre, cell_lines, hits_versions, adding_constant, num_permutations, threshold, input_root, "4 PC") + theme(legend.position="none")

    input_root <- paste0(tissue_specific_chrom_hmm_enrichment_dir,"gaussian_dynamic_qtl_input_file_environmental_variable_time_steps_genotype_version_dosage_model_type_glm_covariate_method_pc1_5_")
    plot_5 <- cre_enrichments_boxplot(cre, cell_lines, hits_versions, adding_constant, num_permutations, threshold, input_root, "5 PC") + theme(legend.position="none")

    input_root <- paste0(tissue_specific_chrom_hmm_enrichment_dir,"gaussian_dynamic_qtl_input_file_environmental_variable_time_steps_genotype_version_dosage_model_type_glm_covariate_method_pc1_6_")
    plot_6 <- cre_enrichments_boxplot(cre, cell_lines, hits_versions, adding_constant, num_permutations, threshold, input_root, "6 PC") + theme(legend.position="none")

    input_root <- paste0(tissue_specific_chrom_hmm_enrichment_dir,"gaussian_dynamic_qtl_input_file_environmental_variable_time_steps_genotype_version_dosage_model_type_glm_covariate_method_pc1_7_")
    plot_7 <- cre_enrichments_boxplot(cre, cell_lines, hits_versions, adding_constant, num_permutations, threshold, input_root, "7 PC") + theme(legend.position="none")

    input_root <- paste0(tissue_specific_chrom_hmm_enrichment_dir,"gaussian_dynamic_qtl_input_file_environmental_variable_time_steps_genotype_version_dosage_model_type_glm_covariate_method_pc1_8_")
    plot_8 <- cre_enrichments_boxplot(cre, cell_lines, hits_versions, adding_constant, num_permutations, threshold, input_root, "8 PC") + theme(legend.position="none")

    input_root <- paste0(tissue_specific_chrom_hmm_enrichment_dir,"gaussian_dynamic_qtl_input_file_environmental_variable_time_steps_genotype_version_dosage_model_type_glm_covariate_method_pc1_9_")
    plot_9 <- cre_enrichments_boxplot(cre, cell_lines, hits_versions, adding_constant, num_permutations, threshold, input_root, "9 PC") + theme(legend.position="none")

    input_root <- paste0(tissue_specific_chrom_hmm_enrichment_dir,"gaussian_dynamic_qtl_input_file_environmental_variable_time_steps_genotype_version_dosage_model_type_glm_covariate_method_pc1_10_")
    plot_10 <- cre_enrichments_boxplot(cre, cell_lines, hits_versions, adding_constant, num_permutations, threshold, input_root, "10 PC") + theme(legend.position="none")



    combined <- plot_grid(plot_0, plot_1, plot_2, plot_3, plot_4, plot_5, plot_6, plot_7, plot_8, plot_9, plot_10, legend, labels = c('A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', ''), ncol = 3)
    ggsave(combined, file=output_file, width=7.2,height=7.0,units="in")
}



produce_figure_3 <- function(qtl_results_dir, time_step_comparison_dir, tissue_specific_chrom_hmm_enrichment_dir, gwas_overlap_dir, visualization_input_dir, output_file) {
    #############################
    # Figure 3a
    ###############################
    dynamic_qtl_input_file <- paste0(visualization_input_dir, "gaussian_dynamic_qtl_input_file_environmental_variable_time_steps_genotype_version_dosage_model_type_glm_covariate_method_pc1_5_dynamic_qtl_efdr_05_visualization_input.txt")
    fig3a <- make_dynamic_qtl_plot(dynamic_qtl_input_file, "rs11124033", "ENSG00000115641", "FHL2", "glm", "pc1_5", -1.5,2.5, "G","A")


    #############################
    # Figure 3b
    ###############################
    input_root <- paste0(tissue_specific_chrom_hmm_enrichment_dir,"gaussian_dynamic_qtl_input_file_environmental_variable_time_steps_genotype_version_dosage_model_type_glm_covariate_method_pc1_5_")
    cell_lines <- c("ipsc_only", "heart_only")
    hits_versions <- c("early_time_step", "late_time_step")
    cre <- "enhancer"
    adding_constant <- 1
    num_permutations="1000"
    threshold = "1.0"
    fig3b <- cre_enrichments_boxplot(cre, cell_lines, hits_versions, adding_constant, num_permutations, threshold, input_root, "none")

    #############################
    # Figure 3c
    ###############################
    dynamic_qtl_input_file <- paste0(visualization_input_dir, "gaussian_dynamic_qtl_input_file_environmental_variable_time_steps_genotype_version_dosage_model_type_glm_quadratic_covariate_method_pc1_5_dynamic_qtl_efdr_05_visualization_input.txt")
    fig3c <- make_dynamic_qtl_plot(dynamic_qtl_input_file, "rs28818910", "ENSG00000167173", "C15orf39", "glm_quadratic", "pc1_5", -2,3.5, "C","T")

    #############################
    # Figure 3c
    ###############################

    fig3d <- make_miami_plot_one_phen(gwas_overlap_dir)

    legend_3a <- get_legend(fig3a)
    legend_3b <- get_legend(fig3b + theme(legend.position="top"))
    legend_3c <- get_legend(fig3c)

    combined <- ggdraw() + 
        draw_plot(fig3a + theme(legend.position='none') + labs(title=""),-.01,.5,.605,.5) + 
        draw_plot(legend_3a,.08,.45,1,1) +
        draw_plot(fig3b + theme(legend.position='none'),.59,.5,.4,.5) + 
        draw_plot(legend_3b,.755,.45,1,1) +
        draw_plot(fig3c + theme(legend.position='none')+ labs(title=""),-.01,0,.605,.5) + 
        draw_plot(legend_3c,.08,-.05,1,1) +
        draw_plot(fig3d + theme(legend.position='none'),.565,0,.434,.4) +
        draw_plot_label(c("A","B","C", "D"),c(.01,.58,.01,.58),c(1,1.0,.5,.5),size=12) + 
        draw_plot_label("Dynamic eQTL for C15orf39", c(.63),c(.47),size=8,fontface="plain") +
        draw_plot_label("GWAS: body mass index", c(.65),c(.05),size=8,fontface="plain") +
        draw_plot_label("Mb",c(.65),c(.258),size=7,fontface="plain")
    ggsave(combined, file=output_file, width=7.2, height=4.5,units="in")


}


two_dynamic_qtls_that_are_known_gwas_variants <- function(output_file, visualization_input_dir) {

    dynamic_qtl_input_file <- paste0(visualization_input_dir, "gaussian_dynamic_qtl_input_file_environmental_variable_time_steps_genotype_version_dosage_model_type_glm_covariate_method_pc1_5_dynamic_qtl_efdr_05_visualization_input.txt")


    dynamic_qtl_plot1 <- make_dynamic_qtl_plot(dynamic_qtl_input_file, "rs7633988", "ENSG00000183873", "SCN5A", "glm", "pc1_5", -2,3.5, "T","A")

    dynamic_qtl_plot2 <- make_dynamic_qtl_plot(dynamic_qtl_input_file, "rs6599234", "ENSG00000183873", "SCN5A", "glm", "pc1_5", -2,3.5, "A","T")


    combined <- plot_grid(dynamic_qtl_plot1+ labs(title="")+ theme(legend.position='right'), dynamic_qtl_plot2+ labs(title="")+ theme(legend.position='right'), labels = c("A", "B"), ncol=1)


    ggsave(combined, file=output_file, width=7.2, height=5.5,units="in")
}



get_vec_from_matrix <- function(mat1) {
    num_lines <- dim(mat1)[1]
    diff <- c()
    for (ii in 1:(num_lines-1)) {
        for (jj in (ii+1):num_lines) {
            value <- mat1[ii,jj]
            diff <- c(diff, value)
        }
    }
    return(diff)
}

load_in_overlap_matrix <- function(real_file_name) {
    aa <- read.table(real_file_name,header=TRUE)
    data <- aa[,2:(1+dim(aa)[1])]
    return(data)
}

# Plot absolute difference between real and observed data for each of the covariate methods
violin_plot_top_n_genes <- function(cell_line_overlap_analysis_dir, output_file, nn, covariate_methods, model_options, covariate_method_names) {
    # Initialize arrays
    frequency <- c()
    methods <- c()
    data_type <- c()
    used_names <- c()
    for (itera in 1:length(covariate_methods)) {
        for (itera_model in 1:length(model_options)) {
        covariate_method <- covariate_methods[itera]
        curr_model <- model_options[itera_model]

        real_file_name <- paste0(cell_line_overlap_analysis_dir, "gaussian_dynamic_qtl_input_file_environmental_variable_time_steps_genotype_version_dosage_model_type_", curr_model, "_covariate_method_", covariate_method, "_", nn, "_genes_real_overlap_matrix.txt")
        perm_file_name <- paste0(cell_line_overlap_analysis_dir, "gaussian_dynamic_qtl_input_file_environmental_variable_time_steps_genotype_version_dosage_model_type_", curr_model, "_covariate_method_", covariate_method, "_", nn, "_genes_perm_overlap_matrix.txt")

        real_overlap_mat <- load_in_overlap_matrix(real_file_name)
        perm_overlap_mat <- load_in_overlap_matrix(perm_file_name)

        real_vec <- get_vec_from_matrix(real_overlap_mat)
        perm_vec <- get_vec_from_matrix(perm_overlap_mat)
        frequency <- c(frequency, real_vec)
        frequency <- c(frequency, perm_vec)

        namer <- covariate_method_names[itera]
        methods <- c(methods, rep(namer, length(real_vec)))
        methods <- c(methods, rep(namer, length(perm_vec)))

        data_type <- c(data_type, rep("eQTL", length(real_vec)))
        data_type <- c(data_type, rep("Background", length(perm_vec)))


        used_names <- c(used_names, namer)
        }
    }


    real_file_name <- paste0(cell_line_overlap_analysis_dir, "gaussian_dynamic_qtl_input_file_environmental_variable_time_steps_genotype_version_dosage_model_type_", curr_model, "_covariate_method_", covariate_method, "_", nn, "_genes_t0_real_overlap_matrix.txt")
    perm_file_name <- paste0(cell_line_overlap_analysis_dir, "gaussian_dynamic_qtl_input_file_environmental_variable_time_steps_genotype_version_dosage_model_type_", curr_model, "_covariate_method_", covariate_method, "_", nn, "_genes_t0_perm_overlap_matrix.txt")

    real_overlap_mat <- load_in_overlap_matrix(real_file_name)
    perm_overlap_mat <- load_in_overlap_matrix(perm_file_name)

    real_vec <- get_vec_from_matrix(real_overlap_mat)
    perm_vec <- get_vec_from_matrix(perm_overlap_mat)
    frequency <- c(frequency, real_vec)
    frequency <- c(frequency, perm_vec)

    methods <- c(methods, rep("Day 0 eQTL", length(real_vec)))
    methods <- c(methods, rep("Day 0 eQTL", length(perm_vec)))

    data_type <- c(data_type, rep("eQTL", length(real_vec)))
    data_type <- c(data_type, rep("Background", length(perm_vec)))



    df <- data.frame(frequency=frequency, covariate_method=factor(methods,levels= c("Day 0 eQTL",used_names)), data_type=factor(data_type,levels=c("eQTL","Background")))

    #t0_indices <- df$covariate_method=="time_step_0"
    #none_indices <- df$covariate_method=="none_glm"
    #pc1_glm <- df$covariate_method=="pc1_glm"
    #pc1_2_glm <- df$covariate_method=="pc1_2_glm"
    #pc1_3_glm <- df$covariate_method=="pc1_3_glm"
    #pc1_4_glm <- df$covariate_method=="pc1_4_glm"
    #pc1_5_glm <- df$covariate_method=="pc1_5_glm"

    #print(wilcox.test(df$absolute_difference[t0_indices], df$absolute_difference[none_indices]))
    #print(wilcox.test(df$absolute_difference[t0_indices], df$absolute_difference[pc1_glm]))
    #print(wilcox.test(df$absolute_difference[t0_indices], df$absolute_difference[pc1_2_glm]))
    #print(wilcox.test(df$absolute_difference[t0_indices], df$absolute_difference[pc1_3_glm]))
    #print(wilcox.test(df$absolute_difference[t0_indices], df$absolute_difference[pc1_4_glm]))
    #print(wilcox.test(df$absolute_difference[t0_indices], df$absolute_difference[pc1_5_glm]))
    p <- ggplot(df, aes(x=covariate_method, y=frequency, fill=data_type))
    p <- p + theme(legend.text = element_text(size=10))+ theme(legend.title = element_text(size=10)) + scale_fill_manual(values=c("goldenrod3","steelblue3"))
    p <- p + labs(x = "", y = "Cell line overalap frequency", fill="")
    p <- p + theme(text = element_text(size=8),axis.text=element_text(size=7), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.text = element_text(size=7), legend.title = element_text(size=8),axis.text.x = element_text(angle = 90,hjust=1, vjust=.5)) 
    p <- p + theme(legend.position="bottom")
    p <- p + geom_violin()

    ggsave(p, file=output_file, width=7.2, height=5.5,units="in")


}


cmp_glm_glmm <- function(input_file, output_file) {
    df <- data.frame(-log10(read.table(input_file,header=TRUE) + .000001))
    cory <- cor(df$glm_pvalue, df$glmm_pvalue)
    maxy <- max(max(df$glm_pvalue), max(df$glmm_pvalue))
    print(cor(df$glm_pvalue, df$glmm_pvalue))
    p <- ggplot(df, aes(glm_pvalue, glmm_pvalue)) + geom_point(size=.1) + scale_color_manual(values=c("chartreuse4"))
    p <- p + theme(text = element_text(size=8),axis.text=element_text(size=7), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.text = element_text(size=7), legend.title = element_text(size=8)) 
    p <- p + labs(x = expression(-log[10]("LM p-value")), y = expression(-log[10]("LMM p-value")))
    p <- p + geom_abline()
    p <- p + scale_x_continuous(limits = c(0, maxy)) +  scale_y_continuous(limits = c(0, maxy))


    ggsave(p, file=output_file, width=7.2, height=5.5,units="in")
}

histogram_showing_nominal_time_steps <- function(dynamic_egenes, pvalue_threshold) {
    time_steps <- c()
    num_genes <- c()
    for (time_step in 0:15) {
        time_steps <- c(time_steps, time_step)
        num_hits <- sum(dynamic_egenes[,(time_step+2)] <= pvalue_threshold)
        num_genes <- c(num_genes, num_hits)
    }
    df <- data.frame(time_steps=time_steps, num_dynamic_qtls=num_genes)
    barplot <- ggplot(df, aes(time_steps, num_dynamic_qtls)) + geom_bar(stat = "identity",aes(fill=time_steps))
    barplot <- barplot + theme(text = element_text(size=8),axis.text=element_text(size=7), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.text = element_text(size=7), legend.title = element_text(size=8)) 
    barplot <- barplot + labs(x = "Day", y = paste0("# dynamic eQTLs"))
    barplot <- barplot + scale_fill_gradient(low="darkgrey",high="firebrick")
             
    return(barplot + theme(legend.position="none"))

}

histogram_showing_number_nominal_time_steps <- function(dynamic_egenes, pvalue_threshold) {
    time_step_pvalz <- dynamic_egenes[,2:17]
    num_time_steps <- rowSums(time_step_pvalz <= pvalue_threshold)
    df <- data.frame(num_time_steps=num_time_steps)
    histo <- ggplot(data=df, aes(df$num_time_steps)) +
            geom_histogram(breaks=seq(-.5, 16.5, by = 1),col="grey", fill="dodgerblue3") +
            labs(x=paste0("# days dynamic eQTL is significant (p <=", pvalue_threshold,") in"), y="Count") +
            theme(text = element_text(size=8),axis.text=element_text(size=7), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.text = element_text(size=7), legend.title = element_text(size=8)) 
    return(histo)
}


joint_plot_summarizing_time_step_independent_comparison <- function(pvalue_threshold, time_step_independent_comparison_file, output_file) {
    dynamic_egenes <- read.table(time_step_independent_comparison_file,header=TRUE)
    time_step_position_bar <- histogram_showing_nominal_time_steps(dynamic_egenes, pvalue_threshold)
    time_step_number_histo <- histogram_showing_number_nominal_time_steps(dynamic_egenes, pvalue_threshold)
    
    combined <- plot_grid(time_step_number_histo, time_step_position_bar, labels = c("A", "B"), ncol=1)
    ggsave(combined, file=output_file, width=7.2, height=4.5,units="in")
}

boxplot_comparing_time_steps_grouped_by_dynamic_qtl_classes <- function(data, output_file) {
    time_step_pvalz <- c()
    time_step_arr <- c()
    dynamic_qtl_grouping <- c()
    for (time_step in 0:15) {
        time_step_pvalz <- c(time_step_pvalz, -log10(data[,(time_step+2)]+ .000000000001))
        dynamic_qtl_grouping <- c(dynamic_qtl_grouping, as.character(data[,19]))

        time_step_arr <- c(time_step_arr, rep(time_step, length(data[,18])))
    }
    dynamic_qtl_grouping2 <- c()
    for (index in 1:length(dynamic_qtl_grouping)) {
        if (dynamic_qtl_grouping[index] == "change") {
            dynamic_qtl_grouping2 <- c(dynamic_qtl_grouping2, "switch")
        } else {
            dynamic_qtl_grouping2 <- c(dynamic_qtl_grouping2, dynamic_qtl_grouping[index])
        }
    }

    if (length(unique(dynamic_qtl_grouping2)) == 3) {
        df <- data.frame(pvalue=time_step_pvalz,time_step=factor(time_step_arr),int_time_step=time_step_arr, dynamic_qtl_class=factor(dynamic_qtl_grouping2, levels = c("early", "switch", "late")))
    } else {
        df <- data.frame(pvalue=time_step_pvalz,time_step=factor(time_step_arr),int_time_step=time_step_arr, dynamic_qtl_class=factor(dynamic_qtl_grouping2, levels = c("early", "switch", "middle", "late")))
    }

    # PLOT
    boxplot <- ggplot(df, aes(x=time_step,y=pvalue,fill=dynamic_qtl_class)) + geom_boxplot(outlier.shape=NA) + labs(x = "Day", y = expression(-log[10]("p-value")),fill="Dynamic QTL Class") + scale_x_discrete(name="Day") + scale_fill_manual(values=c("darkgrey","steelblue3","firebrick"))
    boxplot <- boxplot + theme(text = element_text(size=8),axis.text=element_text(size=7), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.text = element_text(size=7), legend.title = element_text(size=8))
    boxplot <- boxplot + ylim(0,5.0)
    boxplot <- boxplot + theme(legend.position="bottom")

    ggsave(boxplot, file=output_file, width=7.2, height=4.5,units="in")

}

boxplot_comparing_time_steps_grouped_by_dynamic_qtl_classes <- function(input_file, output_file) {
    data <- read.table(input_file,header=TRUE)
    time_step_pvalz <- c()
    time_step_arr <- c()
    dynamic_qtl_grouping <- c()
    for (time_step in 0:15) {
        time_step_pvalz <- c(time_step_pvalz, -log10(data[,(time_step+2)]+ .000000000001))
        dynamic_qtl_grouping <- c(dynamic_qtl_grouping, as.character(data[,19]))

        time_step_arr <- c(time_step_arr, rep(time_step, length(data[,18])))
    }
    dynamic_qtl_grouping2 <- c()
    for (index in 1:length(dynamic_qtl_grouping)) {
        if (dynamic_qtl_grouping[index] == "change") {
            dynamic_qtl_grouping2 <- c(dynamic_qtl_grouping2, "switch")
        } else {
            dynamic_qtl_grouping2 <- c(dynamic_qtl_grouping2, dynamic_qtl_grouping[index])
        }
    }

    if (length(unique(dynamic_qtl_grouping2)) == 3) {
        df <- data.frame(pvalue=time_step_pvalz,time_step=factor(time_step_arr),int_time_step=time_step_arr, dynamic_qtl_class=factor(dynamic_qtl_grouping2, levels = c("early", "switch", "late")))
    } else {
        df <- data.frame(pvalue=time_step_pvalz,time_step=factor(time_step_arr),int_time_step=time_step_arr, dynamic_qtl_class=factor(dynamic_qtl_grouping2, levels = c("early", "switch", "middle", "late")))
    }

    # PLOT
    boxplot <- ggplot(df, aes(x=time_step,y=pvalue,fill=dynamic_qtl_class)) + geom_boxplot(outlier.shape=NA) + labs(x = "Day", y = expression(-log[10]("p-value")),fill="Dynamic QTL Class") + scale_x_discrete(name="Day") + scale_fill_manual(values=c("darkgrey","steelblue3","firebrick"))
    boxplot <- boxplot + theme(text = element_text(size=8),axis.text=element_text(size=7), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.text = element_text(size=7), legend.title = element_text(size=8))
    boxplot <- boxplot + ylim(0,5.0)
    boxplot <- boxplot + theme(legend.position="bottom")

    ggsave(boxplot, file=output_file, width=7.2, height=4.5,units="in")

}
boxplot_comparing_time_steps_grouped_by_two_dynamic_qtl_classes <- function(input_file, output_file) {
    data <- read.table(input_file,header=TRUE)
    time_step_pvalz <- c()
    time_step_arr <- c()
    dynamic_qtl_grouping <- c()
    num_rows <- dim(data)[1]
    for (time_step in 0:15) {
        for (row_num in 1:num_rows) {
            grouper <- as.character(data[row_num, 19])
            if (grouper == "change") {
                word <- "Do Nothing"
            } else {
                time_step_pvalz <- c(time_step_pvalz, -log10(data[row_num,(time_step+2)]+ .000000000001))
                dynamic_qtl_grouping <- c(dynamic_qtl_grouping, as.character(data[row_num,19]))

                time_step_arr <- c(time_step_arr, time_step)
            }
        }
    }
    if (length(unique(dynamic_qtl_grouping)) == 2) {
    df <- data.frame(pvalue=time_step_pvalz,time_step=factor(time_step_arr),int_time_step=time_step_arr, dynamic_qtl_class=factor(dynamic_qtl_grouping))
    } else {
        df <- data.frame(pvalue=time_step_pvalz,time_step=factor(time_step_arr),int_time_step=time_step_arr, dynamic_qtl_class=factor(dynamic_qtl_grouping, levels = c("early", "middle", "late")))
    }

    # PLOT
    boxplot <- ggplot(df, aes(x=time_step,y=pvalue,fill=dynamic_qtl_class)) + geom_boxplot(outlier.shape=NA) + labs(x = "Day", y = expression(-log[10]("p-value")),fill= "") + scale_x_discrete(name="Day") + scale_fill_manual(values=c("darkgrey","salmon","firebrick"))
    boxplot <- boxplot + theme(text = element_text(size=8),axis.text=element_text(size=7), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.text = element_text(size=7), legend.title = element_text(size=8))
    boxplot <- boxplot + ylim(0,5.0)
    boxplot <- boxplot + theme(legend.position="bottom")
    ggsave(boxplot, file=output_file, width=7.2, height=4.5,units="in")

}

# Helper method to qq_plot_vs_uniform
# Makes a qqplot for a specific time step
make_qq_plot_vs_uniform_one_time_step <- function(real_eqtl_file, null_eqtl_file,title_name) {
    # Read in real data
    all_eqtl_nominal_pvalues <- read.table(real_eqtl_file, header=FALSE)

    # Extract pvalues
    pvalues <- sort(all_eqtl_nominal_pvalues$V4)
    # Simulate uniform distribution
    uniform_1 <- sort(runif(length(pvalues)))
    # Sample points (because impossible to plot ALL hits)
    pvalues <- sample_non_significant_hits(pvalues)
    
    # Read in null data
    null_data <- read.table(null_eqtl_file,header=FALSE)
    # Extract pvalues
    null_pvalues <- sort(null_data$V4)
    # Simulate uniform distribution
    uniform_2 <- sort(runif(length(null_pvalues)))
    # Sample points (because impossible to plot ALL hits)
    null_pvalues <- sample_non_significant_hits(null_pvalues)


    # Organize into data frame
    all_pvalues <- c(pvalues, null_pvalues)
    type <- c(rep("real",length(pvalues)),rep("permuted",length(null_pvalues)))

    uniform <- c(sample_non_significant_hits(uniform_1), sample_non_significant_hits(uniform_2))

    df <- data.frame(pvalues=-log10(all_pvalues + 1e-17), expected_pvalues=-log10(uniform + 1e-17), type=factor(type, levels=c("real", "permuted")))

    # PLOT!
    max_val <-max(max(-log10(uniform + 1e-17)), max(-log10(all_pvalues + 1e-17)))
    #PLOT!
    scatter <- ggplot(df, aes(x = expected_pvalues, y = pvalues, colour = type)) + geom_point() 
    scatter <- scatter + theme(text = element_text(size=14), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black")) 
    scatter <- scatter + labs(title=title_name, colour="",x = expression(log[10]("expected p-value")), y = expression(log[10]("observed p-value")))
    scatter <- scatter + geom_abline() 
    scatter <- scatter + theme(legend.position="bottom")
    scatter <- scatter + scale_colour_manual(values=c("dodgerblue3","chartreuse4"))
    scatter <- scatter + scale_x_continuous(limits = c(-.1, max_val + .1), breaks = round(seq(0, max_val, by = 5),1)) + scale_y_continuous(limits = c(-.1,max_val+.1), breaks = round(seq(0, max_val, by = 5),1))
    scatter <- scatter + theme(plot.title=element_text(size=8, face="plain"), text = element_text(size=8),axis.text=element_text(size=7), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.text = element_text(size=7), legend.title = element_text(size=8)) 
    return(scatter)
}

sample_non_significant_hits <- function(pvalues, fraction_kept=.01,fraction_sampled=.001) {
    index <- floor(length(pvalues)*fraction_kept)
    to_keep <- pvalues[1:index]
    to_filter <- pvalues[(index+1):length(pvalues)]
    filtered <- sort(sample(to_filter,floor(length(to_filter)*fraction_sampled)))
    return(c(to_keep,filtered))
}

###############################################################################
# QQPlot showing both:
## 1. Real data pvalues compared to uniform
## 2. Permuted data pvalues compared to pvalues
# One plot for linear dynamic eqtls and 1 plot for nonlinear dynamic qtls
#################################################################################
qq_plot_for_linear_and_non_linear_dynamic_qtls <- function(qtl_results_dir, output_file) {
    linear_real_file <- paste0(qtl_results_dir, "gaussian_dynamic_qtl_input_file_environmental_variable_time_steps_genotype_version_dosage_model_type_glm_covariate_method_pc1_5_permute_False_results.txt")
    linear_perm_file <- paste0(qtl_results_dir, "gaussian_dynamic_qtl_input_file_environmental_variable_time_steps_genotype_version_dosage_model_type_glm_covariate_method_pc1_5_permute_True_results.txt")
    nonlinear_real_file <- paste0(qtl_results_dir, "gaussian_dynamic_qtl_input_file_environmental_variable_time_steps_genotype_version_dosage_model_type_glm_quadratic_covariate_method_pc1_5_permute_False_results.txt")
    nonlinear_perm_file <- paste0(qtl_results_dir, "gaussian_dynamic_qtl_input_file_environmental_variable_time_steps_genotype_version_dosage_model_type_glm_quadratic_covariate_method_pc1_5_permute_True_results.txt")

    linear_qq <- make_qq_plot_vs_uniform_one_time_step(linear_real_file,linear_perm_file, "Linear dynamic eQTL")
    nonlinear_qq <- make_qq_plot_vs_uniform_one_time_step(nonlinear_real_file, nonlinear_perm_file, "Nonlinear dynamic eQTL")

    legend <- get_legend(linear_qq)

    gg <- plot_grid(linear_qq + theme(legend.position="none"),nonlinear_qq + theme(legend.position="none"), legend, ncol=1, label_size=8, labels = c('A', 'B'), rel_heights=c(1,1,.06))

    ggsave(gg, file=output_file, width=7.2,height=5, units="in")

}

compare_eqtl_results_to_banovich_eqtls <- function(eqtl_data_set_comparison_dir, output_file) {
    ipsc_file <- paste0(eqtl_data_set_comparison_dir, "gaussian_dynamic_qtl_input_file_environmental_variable_time_steps_genotype_version_dosage_model_type_glm_covariate_method_pc1_5_1.0_compare_dynamic_qtls_to_ipsc.txt")
    ipsc_cm_file <- paste0(eqtl_data_set_comparison_dir, "gaussian_dynamic_qtl_input_file_environmental_variable_time_steps_genotype_version_dosage_model_type_glm_covariate_method_pc1_5_1.0_compare_dynamic_qtls_to_ipsc_cm.txt")

    ipsc_data <- read.table(ipsc_file, header=TRUE)
    ipsc_cm_data <- read.table(ipsc_cm_file, header=TRUE)
    ipsc_data <- ipsc_data[as.character(ipsc_data$dynamic_eqtl_class) != "change",]
    ipsc_cm_data <- ipsc_cm_data[as.character(ipsc_cm_data$dynamic_eqtl_class) != "change",]

    pvalue_arr <- c()
    dynamic_eqtl_class_arr <- c()
    cell_type_arr <- c()

    pvalue_arr <- c(pvalue_arr, ipsc_data$eqtl_pvalue + 1e-9)
    dynamic_eqtl_class_arr <- c(dynamic_eqtl_class_arr, as.character(ipsc_data$dynamic_eqtl_class))
    cell_type_arr <- c(cell_type_arr, rep("ipsc", length(ipsc_data$dynamic_eqtl_class)))

    pvalue_arr <- c(pvalue_arr, ipsc_cm_data$eqtl_pvalue + 1e-9)
    dynamic_eqtl_class_arr <- c(dynamic_eqtl_class_arr, as.character(ipsc_cm_data$dynamic_eqtl_class))
    cell_type_arr <- c(cell_type_arr, rep("ipsc_cm", length(ipsc_cm_data$dynamic_eqtl_class)))

    df <- data.frame(pvalue=-log10(pvalue_arr), cell_type = as.factor(cell_type_arr), dynamic_eqtl_class = as.factor(dynamic_eqtl_class_arr))
    print(summary(df))
    # PLOT
    boxplot <- ggplot(df, aes(x=cell_type,y=pvalue,fill=dynamic_eqtl_class)) + geom_boxplot() + labs(x = "Cell type", y = expression(-log[10]("p-value")),fill= "") + scale_fill_manual(values=c("darkgrey", "firebrick"))
    boxplot <- boxplot + theme(text = element_text(size=8),axis.text=element_text(size=7), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.text = element_text(size=7), legend.title = element_text(size=8))
    boxplot <- boxplot + theme(legend.position="bottom")
    ggsave(boxplot, file=output_file, width=7.2, height=4.5,units="in")
}




#############################
## Command Line Args
#############################
qtl_results_dir = args[1]
cell_line_overlap_analysis_dir = args[2]
tissue_specific_chrom_hmm_enrichment_dir = args[3]
time_step_independent_comparison_dir = args[4]
gwas_overlap_dir = args[5]
eqtl_data_set_comparison_dir = args[6]
visualization_input_dir = args[7]
visualization_dir = args[8]



###############################################################################
# Make Plot comparing dynamic eQTLs with Banovich eqtl results
#################################################################################
output_file <- paste0(visualization_dir, "dynamic_eqtl_comparison_to_banovich_eqtls.png")
compare_eqtl_results_to_banovich_eqtls(eqtl_data_set_comparison_dir, output_file)

print("DONE")
###############################################################################
# Make Manuscript Figure 3
#################################################################################
output_file <- paste0(visualization_dir, "figure3.png")
# produce_figure_3(qtl_results_dir, time_step_independent_comparison_dir, tissue_specific_chrom_hmm_enrichment_dir, gwas_overlap_dir, visualization_input_dir, output_file)


###############################################################################
# Make CRE enrichment boxplot over a range of number of PCs
#################################################################################
output_file <- paste0(visualization_dir, "cre_enrichment_boxplot_over_a_range_of_pcs.png")
#cre_enrichment_over_range_of_pcs_boxplot(tissue_specific_chrom_hmm_enrichment_dir, output_file)


###############################################################################
# QQPlot showing both:
## 1. Real data pvalues compared to uniform
## 2. Permuted data pvalues compared to pvalues
# One plot for linear dynamic eqtls and 1 plot for nonlinear dynamic qtls
#################################################################################
output_file <- paste0(visualization_dir, "linear_and_nonlinear_pc1_5_qq_plots.png")
#qq_plot_for_linear_and_non_linear_dynamic_qtls(qtl_results_dir, output_file)

###############################################################################
# Make plot showing two dynamic QTLs for SCN5A that are known GWAS variants
#################################################################################
output_file <- paste0(visualization_dir, "two_gwas_dynamic_qtls.png")
#two_dynamic_qtls_that_are_known_gwas_variants(output_file, visualization_input_dir)


###############################################################################
# Make plot showing middle dynamic QTL for rs8107849-ZNF606
#################################################################################
output_file <- paste0(visualization_dir, "rs8107849_ENSG00000166704_nonlinear_viz.png")
dynamic_qtl_input_file <- paste0(visualization_input_dir, "gaussian_dynamic_qtl_input_file_environmental_variable_time_steps_genotype_version_dosage_model_type_glm_quadratic_covariate_method_pc1_5_dynamic_qtl_efdr_05_visualization_input.txt")
#non_linear_dynamic_qtl_plot <- make_dynamic_qtl_plot(dynamic_qtl_input_file, "rs8107849", "ENSG00000166704", "ZNF606", "glm_quadratic", "pc1_5", -2.5,4, 'T','C')
output_file <- paste0(visualization_dir, "rs8107849_ENSG00000166704_nonlinear_viz.png")
#ggsave(non_linear_dynamic_qtl_plot + labs(title=""), file=output_file, width=7.2, height=4.5,units="in")


###############################################################################
# Plot frequency distributions for real and background for each of the covariate methods using top n genes
###############################################################################
num_genes <- "200"
covariate_methods <- c("none", "pc1", "pc1_2", "pc1_3", "pc1_4", "pc1_5", "pc1_6", "pc1_7", "pc1_8", "pc1_9", "pc1_10")
covariate_method_names <- c("Dynamic eQTL (0 PC)", "Dynamic eQTL (1 PC)", "Dynamic eQTL (2 PC)", "Dynamic eQTL (3 PC)", "Dynamic eQTL (4 PC)", "Dynamic eQTL (5 PC)", "Dynamic eQTL (6 PC)", "Dynamic eQTL (7 PC)", "Dynamic eQTL (8 PC)", "Dynamic eQTL (9 PC)", "Dynamic eQTL (10 PC)")
model_options <- c("glm")
output_file <- paste0(visualization_dir, "real_and_observed_in_cell_line_overlap_glm_", num_genes,"_violin_plot.png")
# violin_plot_top_n_genes(cell_line_overlap_analysis_dir, output_file, num_genes, covariate_methods, model_options, covariate_method_names)

###############################################################################
# Scatter plot comparing pvalues from glm dynamic QTLs and glmm dynamic qtls
###############################################################################
output_file <- paste0(visualization_dir, "compare_glm_glmm_pc1_5_scatter.png")
input_file <- paste0(visualization_input_dir, "gaussian_dynamic_qtl_input_file_environmental_variable_time_steps_genotype_version_dosage_covariate_method_pc1_5_merge_glm_glmm.txt")
#cmp_glm_glmm(input_file, output_file)


############################################################################
# Make cowplot combined of:
# Histogram showing number of time-steps have pvalue < $pvalue_threshold for a given dynamic qtl
# Histogram showing which time-steps have pvalue < $pvalue_threshold for our dynamic qtls 
############################################################################
pvalue_threshold <- .05
time_step_independent_comparison_file <- paste0(time_step_independent_comparison_dir, "gaussian_dynamic_qtl_input_file_environmental_variable_time_steps_genotype_version_dosage_model_type_glm_covariate_method_pc1_5_1.0_dynamic_standard_egenes_comparison.txt")
output_file <- paste0(visualization_dir, "joint_plot_summarizing_per_time_step_eqtl_comparison_",pvalue_threshold,".png")
#joint_plot_summarizing_time_step_independent_comparison(pvalue_threshold, time_step_independent_comparison_file, output_file)


############################################################################
# Boxplot showing Standard eQTL p- values (y-axis) in all 16 time steps (x-axis) of linear dynamic eQTLs (most significant variant per dynamic eQTL gene) stratified by linear dynamic eQTL classifications (early, switch, and late) 
############################################################################
output_file <- paste0(visualization_dir, "dynamic_egenes_glm_pc1_5_boxplot_comparing_per_time_step_qtls.png")
time_step_independent_comparison_file <- paste0(time_step_independent_comparison_dir, "gaussian_dynamic_qtl_input_file_environmental_variable_time_steps_genotype_version_dosage_model_type_glm_covariate_method_pc1_5_1.0_dynamic_standard_egenes_comparison.txt")
# boxplot_comparing_time_steps_grouped_by_dynamic_qtl_classes(time_step_independent_comparison_file, output_file)

############################################################################
# Boxplot showing Standard eQTL p- values (y-axis) in all 16 time steps (x-axis) of non-linear dynamic eQTLs (most significant variant per dynamic eQTL gene) stratified by linear dynamic eQTL classifications (early, switch, and late) 
############################################################################
output_file <- paste0(visualization_dir, "nonlinear_dynamic_egenes_glm_quadratic_pc1_5_boxplot_comparing_per_time_step_qtls.png")
time_step_independent_comparison_file <- paste0(time_step_independent_comparison_dir, "gaussian_dynamic_qtl_input_file_environmental_variable_time_steps_genotype_version_dosage_model_type_glm_quadratic_covariate_method_pc1_5_1.0_dynamic_standard_egenes_comparison.txt")
# boxplot_comparing_time_steps_grouped_by_two_dynamic_qtl_classes(time_step_independent_comparison_file, output_file)

############################################################################
# Miami Plot with all three significant phenotypes
############################################################################
output_file <- paste0(visualization_dir, "joint_miami_plot_rs28818910_ENSG00000167173.png")
# make_miami_plot(gwas_overlap_dir, output_file)