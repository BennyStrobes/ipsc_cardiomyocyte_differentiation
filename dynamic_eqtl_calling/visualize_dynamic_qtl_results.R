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
cre_enrichments_boxplot <- function(cre, cell_lines, hits_versions, adding_constant, num_permutations, threshold, input_root) {
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


make_joint_miami_plot <- function(dynamic_qtl_file_name, gwas_file_name1, gwas_file_name2, gwas_file_name3, phenotype1, phenotype2, phenotype3, rs_id, ensamble_id, gene_symbol,variant_pos,variant_chrom) {
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

make_miami_plot <- function(gwas_overlap_dir) {
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
    return(miami_plot)
}

produce_figure_3 <- function(qtl_results_dir, time_step_comparison_dir, tissue_specific_chrom_hmm_enrichment_dir, gwas_overlap_dir, dynamic_qtl_input_file, output_file) {
    #############################
    # Figure 3a
    ###############################
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
    fig3b <- cre_enrichments_boxplot(cre, cell_lines, hits_versions, adding_constant, num_permutations, threshold, input_root)

    #############################
    # Figure 3c
    ###############################
    fig3c <- make_dynamic_qtl_plot(dynamic_qtl_input_file, "rs28818910", "ENSG00000167173", "C15orf39", "glm_quadratic", "pc1_5", -2,3.5, "C","T")

    #############################
    # Figure 3c
    ###############################

    fig3d <- make_miami_plot(gwas_overlap_dir)

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












#############################
## Command Line Args
#############################
qtl_results_dir = args[1]
cell_line_overlap_analysis_dir = args[2]
tissue_specific_chrom_hmm_enrichment_dir = args[3]
time_step_independent_comparison_dir = args[4]
gwas_overlap_dir = args[5]
visualization_dir = args[6]
dynamic_qtl_input_file = args[7]


###############################################################################
# Make Manuscript Figure 3
#################################################################################
output_file <- paste0(visualization_dir, "figure3.png")
produce_figure_3(qtl_results_dir, time_step_independent_comparison_dir, tissue_specific_chrom_hmm_enrichment_dir, gwas_overlap_dir, dynamic_qtl_input_file, output_file)

