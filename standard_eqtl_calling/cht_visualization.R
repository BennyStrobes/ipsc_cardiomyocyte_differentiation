args = commandArgs(trailingOnly=TRUE)
library(ggplot2)
library(ggthemes)
library(reshape)
library(cowplot)






###############################################
#  Boxplot of number of egenes as a function of number of pcs
#  Each point in boxplot is a time step
egenes_as_function_of_pcs_boxplot <- function(cht_output_dir, parameter_string, output_file) {
    # Extract data
    egenes <- c()
    time_step <- c()
    num_pcs <- c()

    # loop through time steps and pcs
    for (temp_time_step in 0:15) {
        for (temp_num_pc in 0:5) {
            # egene file for this time step and this pc
            egene_file <- paste0(cht_output_dir,"cht_results_",parameter_string,"_num_pc_", temp_num_pc,"_time_",temp_time_step,"_efdr_thresh_.05_significant_egenes.txt")
            egene_data <- read.table(egene_file, header=TRUE)
            # Get number of egenes at this time step and pc num
            num_egenes <- dim(egene_data)[1]
            # Now store data 
            egenes <- c(egenes,num_egenes)
            time_step <- c(time_step, temp_time_step)
            num_pcs <- c(num_pcs, temp_num_pc)
        }
    }
    # Put data in organized data frame
    df <- data.frame(egenes=as.numeric(egenes),time_step = time_step, num_pcs = factor(num_pcs))
    # PLOT!!
    box_plot <- ggplot(df, aes(x=num_pcs, y=egenes)) + geom_boxplot(width=.54)
    box_plot <- box_plot + theme(text = element_text(size=18), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black")) 
    box_plot <- box_plot + labs(x = "Number of PCs", y = "Number of eGenes")
    ggsave(box_plot, file=output_file,width = 20,height=10.5,units="cm")
}

boxplot_showing_top_n_variant_gene_pairs_per_time_step_with_banovich_results <- function(input_file, output_file, num_genes, pc_num) {
    data <- read.table(input_file, header=TRUE)

    df_ipsc <- data.frame(time_step=factor(data[data$data_type=="iPSC",]$time_step), pvalue=-log10(data[data$data_type=="iPSC",]$pvalue + .00000000001))
    box_plot_ipsc <- ggplot(df_ipsc, aes(x=time_step, y=pvalue)) + geom_boxplot()
    box_plot_ipsc <- box_plot_ipsc + theme(text = element_text(size=18), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black")) 
    box_plot_ipsc <- box_plot_ipsc+ labs(x = "Time Step", y = "-log10(pvalue)")

    df_ipsc_cm <- data.frame(time_step=factor(data[data$data_type=="iPSC_CM",]$time_step), pvalue=-log10(data[data$data_type=="iPSC_CM",]$pvalue + .00000000001))
    box_plot_ipsc_cm <- ggplot(df_ipsc_cm, aes(x=time_step, y=pvalue)) + geom_boxplot()
    box_plot_ipsc_cm <- box_plot_ipsc_cm + theme(text = element_text(size=18), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black")) 
    box_plot_ipsc_cm <- box_plot_ipsc_cm + labs(x = "Time Step", y = "-log10(pvalue)")


    pdf(output_file)
    gg <- plot_grid(box_plot_ipsc, box_plot_ipsc_cm,nrow=2,ncol=1,label_size=8)
    combined_gg <- ggdraw() + draw_plot(gg,0,0,1,1)
    print(combined_gg)
    dev.off()

}


###############################################
#  Boxplot of number of egenes as a function of number of pcs
#  Each point in boxplot is a time step
egenes_as_function_of_pcs_violinplot <- function(cht_output_dir, parameter_string, output_file) {
    # Extract data
    egenes <- c()
    time_step <- c()
    num_pcs <- c()

    # loop through time steps and pcs
    for (temp_time_step in 0:15) {
        for (temp_num_pc in 0:5) {
            # egene file for this time step and this pc
            egene_file <- paste0(cht_output_dir,"cht_results_",parameter_string,"_num_pc_", temp_num_pc,"_time_",temp_time_step,"_efdr_thresh_.05_significant_egenes.txt")
            egene_data <- read.table(egene_file, header=TRUE)
            # Get number of egenes at this time step and pc num
            num_egenes <- dim(egene_data)[1]
            # Now store data 
            egenes <- c(egenes,num_egenes)
            time_step <- c(time_step, temp_time_step)
            num_pcs <- c(num_pcs, temp_num_pc)
        }
    }
    # Put data in organized data frame
    df <- data.frame(egenes=as.numeric(egenes),time_step = time_step, num_pcs = factor(num_pcs))
    # PLOT!!
    box_plot <- ggplot(df, aes(x=num_pcs, y=egenes,colour=time_step)) + geom_violin()
    box_plot <- box_plot  + geom_jitter(aes(colour=time_step),shape=16, position=position_jitter(0.06))
    box_plot <- box_plot + scale_color_gradient(low="darkgrey",high="firebrick")
    box_plot <- box_plot + labs(x = "Number of PCs", y = "# significant genes", colour= "Day") 
    box_plot <- box_plot +  theme(text = element_text(size=8),axis.text=element_text(size=7), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.text = element_text(size=7), legend.title = element_text(size=8)) 
    ggsave(box_plot, file=output_file, width=7.2, height=4.0, units="in")
    return(box_plot)
}



read_in_summary_statistic_data <- function(file_name) {
    xx <- read.table(file_name)
    xxx <- xx[,2:dim(xx)[2]]
    return(xxx)
}

# Plot correlation histogram for summary stat
correlation_heatmap <- function(correlation_matrix, output_file, version,num_pc) {
    colnames(correlation_matrix) <- paste0(0:15)
    rownames(correlation_matrix) <- paste0(0:15)
    nn <- dim(correlation_matrix)[1]
    vec <- c()
    time_step_diff <- c()
    for (i in 1:nn){
        for (j in 1:nn) {
            if (i != j) {
                if (i > j) {
                    vec <- c(vec,correlation_matrix[i,j])
                    time_step_diff <- c(time_step_diff, abs(i-j))
                }
            }
        }
    }
    lm_result <- lm(vec ~ time_step_diff)
    print(summary(lm_result))
    maxy <- max(vec)
    for (i in 1:nn) {
        correlation_matrix[i,i] <- maxy
    }
    melted_corr <- melt(correlation_matrix)

    # Axis labels are factors
    melted_corr$X1 <- factor(melted_corr$X1)
    melted_corr$X2 <- factor(melted_corr$X2)

    #  PLOT!
    heatmap <- ggplot(data=melted_corr, aes(x=X1, y=X2)) + geom_tile(aes(fill=value)) #+ scale_fill_gradient(low="grey",high="plum2")
    #heatmap <- heatmap + scale_fill_distiller()
    #heatmap <- heatmap + scale_fill_brewer(values = brewer.pal(3,"RdPu"))
    heatmap <- heatmap + scale_fill_distiller(palette = "Blues", direction=1)
    heatmap <- heatmap + theme(text = element_text(size=8),axis.text=element_text(size=7), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.text = element_text(size=7), legend.title = element_text(size=8), axis.text.x = element_text(angle = 0, vjust=.5)) 
    heatmap <- heatmap + labs(x = "Day", y = "Day", fill= expression(paste("  ",rho)))
    heatmap <- heatmap + scale_x_discrete(breaks=c("0","1","2","3","4","5","6","7","8","9","10","11","12","13","14","15"),labels=c("0","","","","","5","","","","","10","","","","","15"))
    heatmap <- heatmap + scale_y_discrete(breaks=c("0","1","2","3","4","5","6","7","8","9","10","11","12","13","14","15"),labels=c("0","","","","","5","","","","","10","","","","","15"))

    #hetamap <- heatmap + scale_fill_gradient2(low = "darkgreen", mid = "white", high = "darkred")
    print("good bye")

    return(heatmap)


}




###############################################
# Heatmap of correlation of summary statistics between time steps
# Ie. heatmap is of dimension number of time steps by number of time steps
# Do independently for each number of PCs
summary_statistic_correlation_heatmap <- function(parameter_string, pc_num, cht_output_dir, cht_visualization_dir) {
    # All summary statistic files have the following file stem
    file_stem <- paste0(cht_output_dir,"best_variant_per_egene_",parameter_string, "_num_pc_",pc_num)

    # Specific summmary statistic files
    alpha_file <- paste0(file_stem,"_fdr_.05_alpha.txt")
    beta_file <- paste0(file_stem, "_fdr_.05_beta.txt")
    pvalue_file <- paste0(file_stem, "_fdr_.05_pvalues.txt")

    # Load in summary stat data
    # Of dimension (number of eGenes) X (number of time steps)
    alpha <- read_in_summary_statistic_data(alpha_file)
    beta <- read_in_summary_statistic_data(beta_file)
    pvalue <- read_in_summary_statistic_data(pvalue_file)
    # Compute allelic fraction (or p as it is referred to in the WASP paper)
    allelic_fraction <- alpha/(alpha + beta)

    # All output files have the following stem
    output_file_stem <- paste0(cht_visualization_dir, parameter_string,"_num_pc_",pc_num,"_correlation_heatmap_")

    # Plot correlation histogram for pvalue summary stat
    heatmap <- correlation_heatmap(cor(pvalue,method="spearman"), "Pvalue", pc_num)

    # Plot correlation histogram for pvalue summary stat
    #correlation_heatmap(cor(alpha,method="spearman"), paste0(output_file_stem, "alpha.png"), "Alpha", pc_num)

    # Plot correlation histogram for pvalue summary stat
    #correlation_heatmap(cor(beta,method="spearman"), paste0(output_file_stem, "beta.png"), "Beta", pc_num)

    # Plot correlation histogram for pvalue summary stat
    #correlation_heatmap(cor(allelic_fraction,method="spearman"), paste0(output_file_stem, "allelic_fraction.png"), "Allelic Fraction",pc_num)
    return(heatmap)
}


factor_matrix_heatmap <- function(factor_matrix_file) {
    factor_matrix <- t(as.matrix(read.table(factor_matrix_file)))

    row.names(factor_matrix) = as.character(0:15)

    melted_corr <- melt(factor_matrix)

    # Axis labels are factors
    melted_corr$X1 <- factor(melted_corr$X1)
    melted_corr$X2 <- factor(melted_corr$X2)

    #  PLOT!
    heatmap <- ggplot(data=melted_corr, aes(x=X1, y=X2)) + geom_tile(aes(fill=value)) 
    #heatmap <- heatmap + scale_fill_distiller()
    #heatmap <- heatmap + scale_fill_brewer(values = brewer.pal(3,"RdPu"))
    heatmap <- heatmap + scale_fill_distiller(palette = "Greens", direction=1)
    heatmap <- heatmap + theme(text = element_text(size=8),axis.text=element_text(size=7), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.text = element_text(size=7), legend.title = element_text(size=8), axis.text.x = element_text(angle = 0, vjust=.5)) 
    heatmap <- heatmap + labs(x = "Day", y = "Latent Factor", fill= " ")
    heatmap <- heatmap + scale_x_discrete(breaks=c("0","1","2","3","4","5","6","7","8","9","10","11","12","13","14","15"),labels=c("0","","","","","5","","","","","10","","","","","15"))

    return(heatmap)


}
make_grid_of_factor_matrices_old <- function(factor_matrix_root, output_file) {
    ###########
    alpha <- "0"
    num_factor <- "3"
    factor_matrix_file <- paste0(factor_matrix_root, alpha,"_",num_factor, "_loading_matrix.txt")
    heatmap_0_3 <- factor_matrix_heatmap(factor_matrix_file) + theme(legend.position="bottom") + theme(legend.text = element_text(size=7)) + theme(legend.title = element_text(size=8))+ theme(legend.key.size =  unit(0.25, "in"))
    heatmap_0_3_legend <- get_legend(heatmap_0_3)
    ###########
    alpha <- "0.5"
    num_factor <- "3"
    factor_matrix_file <- paste0(factor_matrix_root, alpha,"_",num_factor, "_loading_matrix.txt")
    heatmap_.5_3 <- factor_matrix_heatmap(factor_matrix_file)+ theme(legend.position="bottom")+ theme(legend.text = element_text(size=7)) + theme(legend.title = element_text(size=8)) + theme(legend.key.size =  unit(0.25, "in"))
    heatmap_.5_3_legend <- get_legend(heatmap_.5_3)
    ###########
    alpha <- "1"
    num_factor <- "3"
    factor_matrix_file <- paste0(factor_matrix_root, alpha,"_",num_factor, "_loading_matrix.txt")
    heatmap_1_3 <- factor_matrix_heatmap(factor_matrix_file)+ theme(legend.position="bottom")+ theme(legend.text = element_text(size=7)) + theme(legend.title = element_text(size=8))+ theme(legend.key.size =  unit(0.25, "in"))
    heatmap_1_3_legend <- get_legend(heatmap_1_3)
    ###########
    alpha <- "0"
    num_factor <- "4"
    factor_matrix_file <- paste0(factor_matrix_root, alpha,"_",num_factor, "_loading_matrix.txt")
    heatmap_0_4 <- factor_matrix_heatmap(factor_matrix_file) + theme(legend.position="bottom") + theme(legend.text = element_text(size=7)) + theme(legend.title = element_text(size=8))+ theme(legend.key.size =  unit(0.25, "in"))
    heatmap_0_4_legend <- get_legend(heatmap_0_4)
    ###########
    alpha <- "0.5"
    num_factor <- "4"
    factor_matrix_file <- paste0(factor_matrix_root, alpha,"_",num_factor, "_loading_matrix.txt")
    heatmap_.5_4 <- factor_matrix_heatmap(factor_matrix_file)+ theme(legend.position="bottom")+ theme(legend.text = element_text(size=7)) + theme(legend.title = element_text(size=8)) + theme(legend.key.size =  unit(0.25, "in"))
    heatmap_.5_4_legend <- get_legend(heatmap_.5_4)
    ###########
    alpha <- "1"
    num_factor <- "4"
    factor_matrix_file <- paste0(factor_matrix_root, alpha,"_",num_factor, "_loading_matrix.txt")
    heatmap_1_4 <- factor_matrix_heatmap(factor_matrix_file)+ theme(legend.position="bottom")+ theme(legend.text = element_text(size=7)) + theme(legend.title = element_text(size=8))+ theme(legend.key.size =  unit(0.25, "in"))
    heatmap_1_4_legend <- get_legend(heatmap_1_4)
    ###########
    alpha <- "0"
    num_factor <- "5"
    factor_matrix_file <- paste0(factor_matrix_root, alpha,"_",num_factor, "_loading_matrix.txt")
    heatmap_0_5 <- factor_matrix_heatmap(factor_matrix_file)+ theme(legend.position="bottom")+ theme(legend.text = element_text(size=7)) + theme(legend.title = element_text(size=8))+ theme(legend.key.size =  unit(0.25, "in"))
    heatmap_0_5_legend <- get_legend(heatmap_0_5)
    ###########
    alpha <- "0.5"
    num_factor <- "5"
    factor_matrix_file <- paste0(factor_matrix_root, alpha,"_",num_factor, "_loading_matrix.txt")
    heatmap_.5_5 <- factor_matrix_heatmap(factor_matrix_file)+ theme(legend.position="bottom")+ theme(legend.text = element_text(size=7)) + theme(legend.title = element_text(size=8))+ theme(legend.key.size =  unit(0.25, "in"))
    heatmap_.5_5_legend <- get_legend(heatmap_.5_5)
    ###########
    alpha <- "1"
    num_factor <- "5"
    factor_matrix_file <- paste0(factor_matrix_root, alpha,"_",num_factor, "_loading_matrix.txt")
    heatmap_1_5 <- factor_matrix_heatmap(factor_matrix_file)+ theme(legend.position="bottom")+ theme(legend.text = element_text(size=7)) + theme(legend.title = element_text(size=8))+ theme(legend.key.size =  unit(0.25, "in"))
    heatmap_1_5_legend <- get_legend(heatmap_1_5)

    combined_heatmap <- plot_grid(heatmap_0_3, heatmap_.5_3, heatmap_1_3, heatmap_0_4, heatmap_.5_4,heatmap_1_4, heatmap_0_5, heatmap_.5_5, heatmap_1_5, labels = c("A","B","C","D","E","F", "G","H",""), ncol=3)
    
    combined_heatmap2 <- ggdraw() + 
                draw_plot(combined_heatmap, 0,0,1,.97) +
                draw_plot_label(c("alpha=0","alpha=.5","alpha=1"),c(.13,.46,.8),c(.98,.98,.98),size=8,fontface="plain")


    ggsave(combined_heatmap2,file=output_file, width=7.2, height=5,units="in")

}


make_grid_of_factor_matrices <- function(factor_matrix_root, output_file) {
    ###########
    alpha <- "0"
    num_factor <- "3"
    factor_matrix_file <- paste0(factor_matrix_root, alpha,"_",num_factor, "_loading_matrix.txt")
    heatmap_0_3 <- factor_matrix_heatmap(factor_matrix_file) + theme(legend.position="top") + theme(legend.text = element_text(size=7)) + theme(legend.title = element_text(size=8))+ theme(legend.key.size =  unit(0.25, "in"))
    heatmap_0_3_legend <- get_legend(heatmap_0_3)
    ###########
    alpha <- "0.5"
    num_factor <- "3"
    factor_matrix_file <- paste0(factor_matrix_root, alpha,"_",num_factor, "_loading_matrix.txt")
    heatmap_.5_3 <- factor_matrix_heatmap(factor_matrix_file)+ theme(legend.position="top")+ theme(legend.text = element_text(size=7)) + theme(legend.title = element_text(size=8)) + theme(legend.key.size =  unit(0.25, "in"))
    heatmap_.5_3_legend <- get_legend(heatmap_.5_3)
    ###########
    alpha <- "1"
    num_factor <- "3"
    factor_matrix_file <- paste0(factor_matrix_root, alpha,"_",num_factor, "_loading_matrix.txt")
    heatmap_1_3 <- factor_matrix_heatmap(factor_matrix_file)+ theme(legend.position="top")+ theme(legend.text = element_text(size=7)) + theme(legend.title = element_text(size=8))+ theme(legend.key.size =  unit(0.25, "in"))
    heatmap_1_3_legend <- get_legend(heatmap_1_3)
    ###########
    alpha <- "0"
    num_factor <- "4"
    factor_matrix_file <- paste0(factor_matrix_root, alpha,"_",num_factor, "_loading_matrix.txt")
    heatmap_0_4 <- factor_matrix_heatmap(factor_matrix_file) + theme(legend.position="top") + theme(legend.text = element_text(size=7)) + theme(legend.title = element_text(size=8))+ theme(legend.key.size =  unit(0.25, "in"))
    heatmap_0_4_legend <- get_legend(heatmap_0_4)
    ###########
    alpha <- "0.5"
    num_factor <- "4"
    factor_matrix_file <- paste0(factor_matrix_root, alpha,"_",num_factor, "_loading_matrix.txt")
    heatmap_.5_4 <- factor_matrix_heatmap(factor_matrix_file)+ theme(legend.position="top")+ theme(legend.text = element_text(size=7)) + theme(legend.title = element_text(size=8)) + theme(legend.key.size =  unit(0.25, "in"))
    heatmap_.5_4_legend <- get_legend(heatmap_.5_4)
    ###########
    alpha <- "1"
    num_factor <- "4"
    factor_matrix_file <- paste0(factor_matrix_root, alpha,"_",num_factor, "_loading_matrix.txt")
    heatmap_1_4 <- factor_matrix_heatmap(factor_matrix_file)+ theme(legend.position="top")+ theme(legend.text = element_text(size=7)) + theme(legend.title = element_text(size=8))+ theme(legend.key.size =  unit(0.25, "in"))
    heatmap_1_4_legend <- get_legend(heatmap_1_4)
    ###########
    alpha <- "0"
    num_factor <- "5"
    factor_matrix_file <- paste0(factor_matrix_root, alpha,"_",num_factor, "_loading_matrix.txt")
    heatmap_0_5 <- factor_matrix_heatmap(factor_matrix_file)+ theme(legend.position="top")+ theme(legend.text = element_text(size=7)) + theme(legend.title = element_text(size=8))+ theme(legend.key.size =  unit(0.25, "in"))
    heatmap_0_5_legend <- get_legend(heatmap_0_5)
    ###########
    alpha <- "0.5"
    num_factor <- "5"
    factor_matrix_file <- paste0(factor_matrix_root, alpha,"_",num_factor, "_loading_matrix.txt")
    heatmap_.5_5 <- factor_matrix_heatmap(factor_matrix_file)+ theme(legend.position="top")+ theme(legend.text = element_text(size=7)) + theme(legend.title = element_text(size=8))+ theme(legend.key.size =  unit(0.25, "in"))
    heatmap_.5_5_legend <- get_legend(heatmap_.5_5)
    ###########
    alpha <- "1"
    num_factor <- "5"
    factor_matrix_file <- paste0(factor_matrix_root, alpha,"_",num_factor, "_loading_matrix.txt")
    heatmap_1_5 <- factor_matrix_heatmap(factor_matrix_file)+ theme(legend.position="top")+ theme(legend.text = element_text(size=7)) + theme(legend.title = element_text(size=8))+ theme(legend.key.size =  unit(0.25, "in"))
    heatmap_1_5_legend <- get_legend(heatmap_1_5)

    combined_heatmap <- plot_grid(heatmap_0_3, heatmap_.5_3, heatmap_1_3, heatmap_0_4, heatmap_.5_4,heatmap_1_4, heatmap_0_5, heatmap_.5_5, heatmap_1_5, labels = c("A","B","C","D","E","F", "G","H","I"), ncol=3)
    
    combined_heatmap2 <- ggdraw() + 
                draw_plot(combined_heatmap, 0,0,1,.97) +
                draw_plot_label(c("alpha=0","alpha=.5","alpha=1"),c(.13,.46,.8),c(.98,.98,.98),size=8,fontface="plain")


    ggsave(combined_heatmap2,file=output_file, width=7.2, height=7,units="in")

}


# Plot correlation histogram for summary stat
symmetric_correlation_heatmap_general <- function(correlation_matrix, output_file, version,num_pc) {
    nn <- dim(correlation_matrix)[1]
    vec <- c()
    for (i in 1:nn){
        for (j in 1:nn) {
            if (i != j) {
                vec <- c(vec,correlation_matrix[i,j])
            }
        }
    }
    maxy <- max(vec)
    for (i in 1:nn) {
        correlation_matrix[i,i] <- maxy
    }
    melted_corr <- melt(correlation_matrix)

    # Axis labels are factors
    melted_corr$X1 <- factor(melted_corr$X1, levels = rownames(correlation_matrix))
    melted_corr$X2 <- factor(melted_corr$X2, levels = colnames(correlation_matrix))

    #  PLOT!
    heatmap <- ggplot(data=melted_corr, aes(x=X1, y=X2)) + geom_tile(aes(fill=value)) 
    #heatmap <- heatmap + scale_fill_distiller()
    #heatmap <- heatmap + scale_fill_brewer(values = brewer.pal(3,"RdPu"))
    heatmap <- heatmap + scale_fill_distiller(palette = "RdPu", direction=1)
    heatmap <- heatmap + theme(text = element_text(size=18), panel.background = element_blank(), axis.text.x = element_text(angle = 90))
    heatmap <- heatmap + labs(x = "Cell type", y = "Cell type", title=paste0("Num PCs = ",num_pc," / version = ", version), fill= "Spearman Rho")

    ggsave(heatmap, file=output_file,width = 25,height=18,units="cm")


}


# Plot correlation histogram for summary stat
assymetric_correlation_heatmap_general <- function(correlation_matrix, output_file, version,num_pc) {
    nn1 <- dim(correlation_matrix)[1]
    nn2 <- dim(correlation_matrix)[2]
    vec <- c()
    for (i in 1:nn1){
        for (j in 1:nn2) {
            vec <- c(vec,correlation_matrix[i,j])
        }
    }
    maxy <- max(vec)
    melted_corr <- melt(correlation_matrix)

    # Axis labels are factors
    melted_corr$X1 <- factor(melted_corr$X1, levels = rownames(correlation_matrix))
    melted_corr$X2 <- factor(melted_corr$X2, levels = colnames(correlation_matrix))

    #  PLOT!
    heatmap <- ggplot(data=melted_corr, aes(x=X1, y=X2)) + geom_tile(aes(fill=value)) 
    #heatmap <- heatmap + scale_fill_distiller()
    #heatmap <- heatmap + scale_fill_brewer(values = brewer.pal(3,"RdPu"))
    heatmap <- heatmap + scale_fill_distiller(palette = "RdPu", direction=1)
    heatmap <- heatmap + theme(text = element_text(size=18), panel.background = element_blank(), axis.text.x = element_text(angle = 90))
    heatmap <- heatmap + labs(x = "Cell type", y = "Cell type", title=paste0("Num PCs = ",num_pc," / version = ", version), fill= "Spearman Rho")

    ggsave(heatmap, file=output_file,width = 25,height=18,units="cm")


}

line_plot_of_spearman_correlation_with_banovich_results_across_time <- function(parameter_string, pc_num, cht_output_dir, cht_visualization_dir) {
    input_file <- paste0(cht_output_dir, parameter_string, "_num_pc_", pc_num, "_eqtls_across_time_steps_and_banovich_studies_geometric_mean_05.txt")
    
    raw_data <- read.table(input_file, header=TRUE)
    pvalues <- raw_data[,2:(dim(raw_data)[2])]
    corr_mat <- cor(pvalues,method="spearman")
    ipsc_corr_vec <- corr_mat[(dim(raw_data)[2] -2),]
    ipsc_cm_corr_vec <- corr_mat[(dim(raw_data)[2]-1),]
    # Remove non-time steps
    ipsc_corr_vec <- ipsc_corr_vec[1:(length(ipsc_corr_vec)-2)]
    ipsc_cm_corr_vec <- ipsc_cm_corr_vec[1:(length(ipsc_cm_corr_vec)-2)]


    correlation_vec <- c()
    time_steps <- c()
    banovich_data <- c()

    correlation_vec <- c(correlation_vec, ipsc_corr_vec)
    time_steps <- c(time_steps,0:15)
    banovich_data <- c(banovich_data, rep("iPSC", length(0:15)))

    correlation_vec <- c(correlation_vec, ipsc_cm_corr_vec)
    time_steps <- c(time_steps,0:15)
    banovich_data <- c(banovich_data, rep("iPSC-derived cardiomyocytes", length(0:15)))


    df <- data.frame(correlation=correlation_vec, time_step=time_steps, data_type=factor(banovich_data))

    #PLOT!
    line_plot <- ggplot(df, aes(x = time_steps, y = correlation, colour = data_type)) + geom_line(size=2.2) 
    line_plot <- line_plot + labs(colour="",x = "Day", y = expression(paste("Spearman's ", rho))) + scale_colour_manual(values=c("darkgrey", "firebrick"))
    line_plot <- line_plot + theme(text = element_text(size=8),axis.text=element_text(size=7), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.text = element_text(size=7), legend.title = element_text(size=8)) 

    return(line_plot)
}


summary_statistic_correlation_heatmap_include_gtex <- function(parameter_string, pc_num, cht_output_dir, cht_visualization_dir) {
    
    # Run when variant gene pairs are selected by those significant egenes (efdr <= .05). Best variant per gene selected by geometric mean
    input_file <- paste0(cht_output_dir, parameter_string, "_num_pc_", pc_num, "_eqtls_across_time_steps_and_gtex_v7_geometric_mean_05.txt")
    output_file <- paste0(cht_visualization_dir, parameter_string, "_num_pc_", pc_num, "_eqtls_across_time_steps_and_gtex_v7_geometric_mean_05_asymmetric_heatmap.png")
    raw_data <- read.table(input_file,header=TRUE)
    pvalues <- raw_data[,2:(dim(raw_data)[2])]
    corr_mat <- cor(pvalues,method="spearman")
    corr_mat_asymetric <- corr_mat[1:16,]
    corr_mat_asymetric <- corr_mat_asymetric[,17:(dim(corr_mat_asymetric)[2])]
    assymetric_correlation_heatmap_general(corr_mat_asymetric, output_file, "geometric_mean_05", pc_num)
    
    # Run when variant gene pairs are selected by those significant egenes (efdr <= .1). Best variant per gene selected by geometric mean
    input_file <- paste0(cht_output_dir, parameter_string, "_num_pc_", pc_num, "_eqtls_across_time_steps_and_gtex_v7_geometric_mean_1.txt")
    output_file <- paste0(cht_visualization_dir, parameter_string, "_num_pc_", pc_num, "_eqtls_across_time_steps_and_gtex_v7_geometric_mean_1_asymmetric_heatmap.png")
    raw_data <- read.table(input_file,header=TRUE)
    pvalues <- raw_data[,2:(dim(raw_data)[2])]
    corr_mat <- cor(pvalues,method="spearman")
    corr_mat_asymetric <- corr_mat[1:16,]
    corr_mat_asymetric <- corr_mat_asymetric[,17:(dim(corr_mat_asymetric)[2])]
    assymetric_correlation_heatmap_general(corr_mat_asymetric, output_file, "geometric_mean_1", pc_num)

    # Run when variant gene pairs are selected by those significant egenes (efdr <= .1). Best variant per gene selected by geometric mean
    input_file <- paste0(cht_output_dir, parameter_string, "_num_pc_", pc_num, "_eqtls_across_time_steps_and_gtex_v7_all_egenes.txt")
    output_file <- paste0(cht_visualization_dir, parameter_string, "_num_pc_", pc_num, "_eqtls_across_time_steps_and_gtex_v7_all_egenes_asymmetric_heatmap.png")
    raw_data <- read.table(input_file,header=TRUE)
    pvalues <- raw_data[,2:(dim(raw_data)[2])]
    corr_mat <- cor(pvalues,method="spearman")
    corr_mat_asymetric <- corr_mat[1:16,]
    corr_mat_asymetric <- corr_mat_asymetric[,17:(dim(corr_mat_asymetric)[2])]
    assymetric_correlation_heatmap_general(corr_mat_asymetric, output_file, "all_egenes", pc_num)



    # Run when variant gene pairs are selected by those significant egenes (efdr <= .05). Best variant per gene selected by geometric mean
    input_file <- paste0(cht_output_dir, parameter_string, "_num_pc_", pc_num, "_eqtls_across_time_steps_and_gtex_v7_geometric_mean_05.txt")
    output_file <- paste0(cht_visualization_dir, parameter_string, "_num_pc_", pc_num, "_eqtls_across_time_steps_and_gtex_v7_geometric_mean_05_symmetric_heatmap.png")
    raw_data <- read.table(input_file,header=TRUE)
    pvalues <- raw_data[,2:(dim(raw_data)[2])]
    corr_mat <- cor(pvalues,method="spearman")
    symmetric_correlation_heatmap_general(corr_mat, output_file, "geometric_mean_05", pc_num)


    # Run when variant gene pairs are selected by those significant egenes (efdr <= .1). Best variant per gene selected by geometric mean
    input_file <- paste0(cht_output_dir, parameter_string, "_num_pc_", pc_num, "_eqtls_across_time_steps_and_gtex_v7_geometric_mean_1.txt")
    output_file <- paste0(cht_visualization_dir, parameter_string, "_num_pc_", pc_num, "_eqtls_across_time_steps_and_gtex_v7_geometric_mean_1_symmetric_heatmap.png")
    raw_data <- read.table(input_file,header=TRUE)
    pvalues <- raw_data[,2:(dim(raw_data)[2])]
    corr_mat <- cor(pvalues,method="spearman")
    symmetric_correlation_heatmap_general(corr_mat, output_file, "geometric_mean_1", pc_num)

    # Run when variant gene pairs are selected by those significant egenes (efdr <= .1). Best variant per gene selected by geometric mean
    input_file <- paste0(cht_output_dir, parameter_string, "_num_pc_", pc_num, "_eqtls_across_time_steps_and_gtex_v7_all_egenes.txt")
    output_file <- paste0(cht_visualization_dir, parameter_string, "_num_pc_", pc_num, "_eqtls_across_time_steps_and_gtex_v7_all_egenes_symmetric_heatmap.png")
    raw_data <- read.table(input_file,header=TRUE)
    pvalues <- raw_data[,2:(dim(raw_data)[2])]
    corr_mat <- cor(pvalues,method="spearman")
    symmetric_correlation_heatmap_general(corr_mat, output_file, "all_egenes", pc_num)

}


sample_non_significant_hits <- function(pvalues, fraction_kept=.01,fraction_sampled=.001) {
    index <- floor(length(pvalues)*fraction_kept)
    to_keep <- pvalues[1:index]
    to_filter <- pvalues[(index+1):length(pvalues)]
    filtered <- sort(sample(to_filter,floor(length(to_filter)*fraction_sampled)))
    return(c(to_keep,filtered))
}

###############################################
# QQPlot showing both:
## 1. Real-data pvalues compared to uniform
## 2. Permuted-data pvalues compared to uniform
## Done independently for each PC 
## Each output image has 16 subplots (1 for each time step)

qq_plot_vs_uniform <- function(input_stem, null_stem, output_file) {
    # Make qq-plot for each of 16 time steps

    time_step <- 0
    p0 <- make_qq_plot_vs_uniform_one_time_step(paste0(input_stem,time_step,"_eqtl_results.txt"),paste0(null_stem,time_step,"_eqtl_results.txt"), time_step)
    
    time_step <- 1
    p1 <- make_qq_plot_vs_uniform_one_time_step(paste0(input_stem,time_step,"_eqtl_results.txt"),paste0(null_stem,time_step,"_eqtl_results.txt"), time_step)
    
    time_step <- 2
    p2 <- make_qq_plot_vs_uniform_one_time_step(paste0(input_stem,time_step,"_eqtl_results.txt"),paste0(null_stem,time_step,"_eqtl_results.txt"), time_step)
        
    time_step <- 4
    p3 <- make_qq_plot_vs_uniform_one_time_step(paste0(input_stem,time_step,"_eqtl_results.txt"),paste0(null_stem,time_step,"_eqtl_results.txt"), time_step)
    
    time_step <- 4
    p4 <- make_qq_plot_vs_uniform_one_time_step(paste0(input_stem,time_step,"_eqtl_results.txt"),paste0(null_stem,time_step,"_eqtl_results.txt"), time_step)
        
    time_step <- 5
    p5 <- make_qq_plot_vs_uniform_one_time_step(paste0(input_stem,time_step,"_eqtl_results.txt"),paste0(null_stem,time_step,"_eqtl_results.txt"), time_step)
    
    time_step <- 6
    p6 <- make_qq_plot_vs_uniform_one_time_step(paste0(input_stem,time_step,"_eqtl_results.txt"),paste0(null_stem,time_step,"_eqtl_results.txt"), time_step)
        
    time_step <- 7
    p7 <- make_qq_plot_vs_uniform_one_time_step(paste0(input_stem,time_step,"_eqtl_results.txt"),paste0(null_stem,time_step,"_eqtl_results.txt"), time_step)
    
    time_step <- 8
    p8 <- make_qq_plot_vs_uniform_one_time_step(paste0(input_stem,time_step,"_eqtl_results.txt"),paste0(null_stem,time_step,"_eqtl_results.txt"), time_step)
    
    time_step <- 9
    p9 <- make_qq_plot_vs_uniform_one_time_step(paste0(input_stem,time_step,"_eqtl_results.txt"),paste0(null_stem,time_step,"_eqtl_results.txt"), time_step)
    
    time_step <- 10
    p10 <- make_qq_plot_vs_uniform_one_time_step(paste0(input_stem,time_step,"_eqtl_results.txt"),paste0(null_stem,time_step,"_eqtl_results.txt"), time_step)
    
    time_step <- 11
    p11 <- make_qq_plot_vs_uniform_one_time_step(paste0(input_stem,time_step,"_eqtl_results.txt"),paste0(null_stem,time_step,"_eqtl_results.txt"), time_step)
    
    time_step <- 12
    p12 <- make_qq_plot_vs_uniform_one_time_step(paste0(input_stem,time_step,"_eqtl_results.txt"),paste0(null_stem,time_step,"_eqtl_results.txt"), time_step)
    
    time_step <- 13
    p13 <- make_qq_plot_vs_uniform_one_time_step(paste0(input_stem,time_step,"_eqtl_results.txt"),paste0(null_stem,time_step,"_eqtl_results.txt"), time_step)
    
    time_step <- 14
    p14 <- make_qq_plot_vs_uniform_one_time_step(paste0(input_stem,time_step,"_eqtl_results.txt"),paste0(null_stem,time_step,"_eqtl_results.txt"), time_step)

    time_step <- 15
    p15 <- make_qq_plot_vs_uniform_one_time_step(paste0(input_stem,time_step,"_eqtl_results.txt"),paste0(null_stem,time_step,"_eqtl_results.txt"), time_step)
    
    # Merge those 16 plots with cowplot!
    legend = get_legend(p0)
    gg <- plot_grid(p0 + theme(legend.position="none"),p1 + theme(legend.position="none"),p2 + theme(legend.position="none"),p3 + theme(legend.position="none"),p4 + theme(legend.position="none"),p5 + theme(legend.position="none"),p6 + theme(legend.position="none"),p7 + theme(legend.position="none"),p8 + theme(legend.position="none"),p9 + theme(legend.position="none"),p10 + theme(legend.position="none"),p11 + theme(legend.position="none"),p12 + theme(legend.position="none"),p13 + theme(legend.position="none"),p14 + theme(legend.position="none"),p15 + theme(legend.position="none"),nrow=4,ncol=4,label_size=8)

    gg_combined <- plot_grid(gg, legend,ncol=1,rel_heights = c(1, .04))

    ggsave(gg_combined,file=output_file, width=7.2, height=8,units="in")

}

# Helper method to qq_plot_vs_uniform
# Makes a qqplot for a specific time step
make_qq_plot_vs_uniform_one_time_step <- function(real_eqtl_file, null_eqtl_file, time_step) {
    # Read in real data
    all_eqtl_nominal_pvalues <- read.table(real_eqtl_file, header=TRUE)
    # Extract pvalues
    pvalues <- sort(all_eqtl_nominal_pvalues$pvalue)
    # Simulate uniform distribution
    uniform_1 <- sort(runif(length(pvalues)))
    # Sample points (because impossible to plot ALL hits)
    pvalues <- sample_non_significant_hits(pvalues)
    
    # Read in null data
    null_data <- read.table(null_eqtl_file,header=TRUE)
    # Extract pvalues
    null_pvalues <- sort(null_data$pvalue)
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
    scatter <- scatter + labs(colour="",x = expression(log[10]("expected p-value")), y = expression(log[10]("observed p-value")), title = paste0("Day ", time_step))
    scatter <- scatter + geom_abline() 
    scatter <- scatter + theme(legend.position="bottom")
    scatter <- scatter + scale_colour_manual(values=c("dodgerblue3","chartreuse4"))
    scatter <- scatter + scale_x_continuous(limits = c(-.1, max_val + .1), breaks = round(seq(0, max_val, by = 5),1)) + scale_y_continuous(limits = c(-.1,max_val+.1), breaks = round(seq(0, max_val, by = 5),1))
    scatter <- scatter + theme(plot.title=element_text(size=8, face="plain"), text = element_text(size=8),axis.text=element_text(size=7), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.text = element_text(size=7), legend.title = element_text(size=8)) 
    return(scatter)
}

###############################################
# QQPlot showing both:
## 1. Real-data pvalues compared to permuted
## Done independently for each PC 
## Each output image has 16 subplots (1 for each time step)
qq_plot_vs_permuted <- function(input_stem, null_stem, output_file) {
    # Make qq-plot for each of 16 time steps
    time_step <- 0
    p0 <- make_qq_plot_vs_permuted_one_time_step(paste0(input_stem,time_step,"_eqtl_results.txt"),paste0(null_stem,time_step,"_eqtl_results.txt"), time_step)

    time_step <- 1
    p1 <- make_qq_plot_vs_permuted_one_time_step(paste0(input_stem,time_step,"_eqtl_results.txt"),paste0(null_stem,time_step,"_eqtl_results.txt"), time_step)

    time_step <- 2
    p2 <- make_qq_plot_vs_permuted_one_time_step(paste0(input_stem,time_step,"_eqtl_results.txt"),paste0(null_stem,time_step,"_eqtl_results.txt"), time_step)

    time_step <- 3
    p3 <- make_qq_plot_vs_permuted_one_time_step(paste0(input_stem,time_step,"_eqtl_results.txt"),paste0(null_stem,time_step,"_eqtl_results.txt"), time_step)

    time_step <- 4
    p4 <- make_qq_plot_vs_permuted_one_time_step(paste0(input_stem,time_step,"_eqtl_results.txt"),paste0(null_stem,time_step,"_eqtl_results.txt"), time_step)

    time_step <- 5
    p5 <- make_qq_plot_vs_permuted_one_time_step(paste0(input_stem,time_step,"_eqtl_results.txt"),paste0(null_stem,time_step,"_eqtl_results.txt"), time_step)

    time_step <- 6
    p6 <- make_qq_plot_vs_permuted_one_time_step(paste0(input_stem,time_step,"_eqtl_results.txt"),paste0(null_stem,time_step,"_eqtl_results.txt"), time_step)

    time_step <- 7
    p7 <- make_qq_plot_vs_permuted_one_time_step(paste0(input_stem,time_step,"_eqtl_results.txt"),paste0(null_stem,time_step,"_eqtl_results.txt"), time_step)

    time_step <- 8
    p8 <- make_qq_plot_vs_permuted_one_time_step(paste0(input_stem,time_step,"_eqtl_results.txt"),paste0(null_stem,time_step,"_eqtl_results.txt"), time_step)

    time_step <- 9
    p9 <- make_qq_plot_vs_permuted_one_time_step(paste0(input_stem,time_step,"_eqtl_results.txt"),paste0(null_stem,time_step,"_eqtl_results.txt"), time_step)

    time_step <- 10
    p10 <- make_qq_plot_vs_permuted_one_time_step(paste0(input_stem,time_step,"_eqtl_results.txt"),paste0(null_stem,time_step,"_eqtl_results.txt"), time_step)

    time_step <- 11
    p11 <- make_qq_plot_vs_permuted_one_time_step(paste0(input_stem,time_step,"_eqtl_results.txt"),paste0(null_stem,time_step,"_eqtl_results.txt"), time_step)

    time_step <- 12
    p12 <- make_qq_plot_vs_permuted_one_time_step(paste0(input_stem,time_step,"_eqtl_results.txt"),paste0(null_stem,time_step,"_eqtl_results.txt"), time_step)

    time_step <- 13
    p13 <- make_qq_plot_vs_permuted_one_time_step(paste0(input_stem,time_step,"_eqtl_results.txt"),paste0(null_stem,time_step,"_eqtl_results.txt"), time_step)

    time_step <- 14
    p14 <- make_qq_plot_vs_permuted_one_time_step(paste0(input_stem,time_step,"_eqtl_results.txt"),paste0(null_stem,time_step,"_eqtl_results.txt"), time_step)

    time_step <- 15
    p15 <- make_qq_plot_vs_permuted_one_time_step(paste0(input_stem,time_step,"_eqtl_results.txt"),paste0(null_stem,time_step,"_eqtl_results.txt"), time_step)

    # Merge those 16 plots with cowplot!
    pdf(output_file)
    gg <- plot_grid(p0,p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12,p13,p14,p15,nrow=4,ncol=4,label_size=8)
    combined_gg <- ggdraw() + draw_plot(gg,0,0,1,1)
    print(combined_gg)
    dev.off()

}

# Helper method to qq_plot_vs_permuted
# Makes a qqplot for a specific time step
make_qq_plot_vs_permuted_one_time_step <- function(real_eqtl_file, null_eqtl_file, time_step) {
    # Extract pvalues from real data
    all_eqtl_nominal_pvalues <- read.table(real_eqtl_file, header=TRUE)
    # Subselect hits (can't plot all hits)
    pvalues <- sample_non_significant_hits(sort(all_eqtl_nominal_pvalues$pvalue))

    # Extract pvalues from permuted data
    null_data <- read.table(null_eqtl_file,header=TRUE)
    # Subselect hits (can't plot all hits)
    null_pvalues <- sample_non_significant_hits(sort(null_data$pvalue))

    # put into convenient data frame
    df <- data.frame(real_pvalues=-log10(pvalues + .000000000001), expected_pvalues=-log10(null_pvalues + .000000000001))


    max_val <-max(max(-log10(pvalues + .000000000001)), max(-log10(null_pvalues + .000000000001)))
    #PLOT!
    scatter <- ggplot(df, aes(x = expected_pvalues, y = real_pvalues)) + geom_point() 
    scatter <- scatter + theme(text = element_text(size=14), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black")) 
    scatter <- scatter + labs(x = "Permuted", y = "Real", title = paste0("Time step ", time_step))
    scatter <- scatter + geom_abline() +  theme(legend.position="none")
    scatter <- scatter + scale_x_continuous(limits = c(-.1, max_val + .1), breaks = round(seq(0, max_val, by = 5),1)) + scale_y_continuous(limits = c(-.1,max_val+.1), breaks = round(seq(0, max_val, by = 5),1))
    return(scatter)
}




eqtl_comparison_to_reference_bar_plot <- function(cht_enrichment_dir,parameter_string, data_set_name, pc_num, plot_file) {
    # First extract data. And get into nice data format
    pvalues <- c()
    version <- c()
    time_step <- c()
    # Loop through time steps
    for (temp_time_step in 0:15) {
        # Get and parse enrichment file for this time step
        ipsc_file <- paste0(cht_enrichment_dir,"enrichment_results_",parameter_string,"_num_pc_",pc_num,"_time_",temp_time_step,"_", data_set_name,"_real_v_matched_controls.txt")
        data <- read.table(ipsc_file,header=TRUE)
        pvalues <- c(pvalues,data$real_pvalue)
        time_step <- c(time_step, rep(temp_time_step, length(data$real_pvalue)))
        version <- c(version, as.character(rep("real", length(data$real_pvalue))))
        pvalues <- c(pvalues,data$matched_pvalue)
        time_step <- c(time_step, rep(temp_time_step, length(data$matched_pvalue)))
        version <- c(version, as.character(rep("matched", length(data$matched_pvalue))))
        #print(temp_time_step)
        #print(wilcox.test(data$real_pvalue,data$matched_pvalue))
    }
    df <- data.frame(pvalues = as.numeric(pvalues), version = factor(version,c("real","matched")), time_step = factor(time_step))

    box_plot <- ggplot(df, aes(x=time_step, y=pvalues, fill=version)) + geom_boxplot(width=.54)
    box_plot <- box_plot + theme(text = element_text(size=18), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black")) 
    box_plot <- box_plot + labs(fill= "Test version",x = "time step", y = "pvalue")
    ggsave(box_plot, file=plot_file,width = 28,height=10.5,units="cm")
}

###############################################
# Boxplot showing pvalues found in our data for only the eqtls in a specific data set
# 1 plot per data set
# 1 plot for pc
eqtl_comparison_to_background_shell <- function(pc_num, cht_visualization_dir, parameter_string, cht_enrichment_dir, eqtl_data_set_file) {
    # Load in data set data
    data_sets <- read.table(eqtl_data_set_file)
    num_data_sets <- dim(data_sets)[1]
    # Loop through each of data sets
    for (data_set_num in 1:num_data_sets) {
        # Name of data set of this line
        data_set_name <- paste0(data_sets[data_set_num, 1])
        # Output file
        plot_file <- paste0(cht_visualization_dir,parameter_string,"_num_pc_",pc_num,"_data_set_comparison_",data_set_name,"_boxplot.png")
        # Make boxplot
        eqtl_comparison_to_reference_bar_plot(cht_enrichment_dir, parameter_string, data_set_name, pc_num, plot_file)
    }
}


eqtl_comparison_to_reference_bar_plot_elegant <- function(cht_enrichment_dir,parameter_string, data_set_name, pc_num, plot_file) {
    # First extract data. And get into nice data format
    pvalues <- c()
    version <- c()
    time_step <- c()
    median_matched_pvalues <- c()
    # Loop through time steps
    for (temp_time_step in 0:15) {
        # Get and parse enrichment file for this time step
        ipsc_file <- paste0(cht_enrichment_dir,"enrichment_results_",parameter_string,"_num_pc_",pc_num,"_time_",temp_time_step,"_", data_set_name,"_real_v_matched_controls.txt")
        data <- read.table(ipsc_file,header=TRUE)
        pvalues <- c(pvalues,data$real_pvalue)
        time_step <- c(time_step, rep(temp_time_step, length(data$real_pvalue)))
        version <- c(version, as.character(rep("real", length(data$real_pvalue))))
        #pvalues <- c(pvalues,data$matched_pvalue)
        median_matched_pvalues <- c(median_matched_pvalues, data$matched_pvalue)
        #time_step <- c(time_step, rep(temp_time_step, length(data$matched_pvalue)))
        #version <- c(version, as.character(rep("matched", length(data$matched_pvalue))))
        #print(temp_time_step)
        #print(wilcox.test(data$real_pvalue,data$matched_pvalue))
    }

    matched_pvalz <- -log10(median_matched_pvalues +.000000001)


    df <- data.frame(pvalues = as.numeric(-log10(pvalues+.000000001)), time_step_int=as.numeric(time_step), time_step = factor(time_step))

    #df_line <- data.frame(time=0:15, pvalues=median_matched_pvalues)

    box_plot <- ggplot() + geom_boxplot(data=df,aes(x=time_step, y=pvalues), fill="deepskyblue2",width=.54, outlier.colour = NULL, outlier.fill=NULL, outlier.shape=NA)
    #box_plot <- boxplot + geom_line(data = df_line, aes(x =time, y = pvalues, group=1))
    box_plot <- box_plot + theme(text = element_text(size=18), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black")) 
    box_plot <- box_plot + labs(fill= "Test version",x = "time step", y = "-log10(pvalue)") + theme(legend.position="none") #+ scale_fill_gradient(low="pink",high="blue")
    box_plot <- box_plot + geom_hline(yintercept = median(matched_pvalz),colour="red",size=.8)
    box_plot <- box_plot + coord_cartesian(ylim=c(0,5.0))
    ggsave(box_plot, file=plot_file,width = 28,height=10.5,units="cm")
}


eqtl_comparison_reference_shell <- function(pc_num, cht_visualization_dir, parameter_string, cht_enrichment_dir, eqtl_data_set_file) {
    # Load in data set data
    data_sets <- read.table(eqtl_data_set_file)
    num_data_sets <- dim(data_sets)[1]
    # Loop through each of data sets
    for (data_set_num in 1:num_data_sets) {
        print(data_set_num)
        # Name of data set of this line
        data_set_name <- paste0(data_sets[data_set_num, 1])
        # Output file
        plot_file <- paste0(cht_visualization_dir,parameter_string,"_num_pc_",pc_num,"_data_set_comparison_",data_set_name,"_boxplot_elegant.png")
        # Make boxplot
        eqtl_comparison_to_reference_bar_plot_elegant(cht_enrichment_dir, parameter_string, data_set_name, pc_num, plot_file)
    }  
}

###############################################
# Boxplot showing pvalues found in our data based on eqtls found in GTEx v7:
## 1. HLV
## 2. skin_not_sun_exposed
## 3. stomach
## 4. Breast
# 1 plot for pc
gtex_tissue_comparison_bar_plot <- function(tissue_comparison_plot_file, pc_num, parameter_string, cht_enrichment_dir) {
    # First extract data. And get into nice data format
    pvalues <- c()
    version <- c()
    time_step <- c()
    for (temp_time_step in 0:15) {
        # HLV
        hlv_file <- paste0(cht_enrichment_dir,"enrichment_results_",parameter_string,"_num_pc_",pc_num, "_time_", temp_time_step, "_gtex_v7_hlv_beta_filter_real_v_matched_controls.txt")
        data <- read.table(hlv_file,header=TRUE)
        pvalues <- c(pvalues,data$real_pvalue)
        time_step <- c(time_step, rep(temp_time_step, length(data$real_pvalue)))
        version <- c(version, as.character(rep("heart_left", length(data$real_pvalue))))
        # stomach
        stomach_file <- paste0(cht_enrichment_dir,"enrichment_results_",parameter_string,"_num_pc_",pc_num, "_time_", temp_time_step, "_gtex_v7_stomach_beta_filter_real_v_matched_controls.txt")
        data <- read.table(stomach_file, header=TRUE)
        pvalues <- c(pvalues,data$real_pvalue)
        time_step <- c(time_step, rep(temp_time_step, length(data$real_pvalue)))
        version <- c(version, as.character(rep("stomach", length(data$real_pvalue))))
        # breast
        breast_file <- paste0(cht_enrichment_dir,"enrichment_results_",parameter_string,"_num_pc_",pc_num, "_time_", temp_time_step, "_gtex_v7_breast_mammary_beta_filter_real_v_matched_controls.txt")
        data <- read.table(breast_file,header=TRUE)
        pvalues <- c(pvalues,data$real_pvalue)
        time_step <- c(time_step, rep(temp_time_step, length(data$real_pvalue)))
        version <- c(version, as.character(rep("breast", length(data$real_pvalue))))
        # Skin not sun exposed
        skin_no_sun_file <- paste0(cht_enrichment_dir,"enrichment_results_",parameter_string,"_num_pc_",pc_num, "_time_", temp_time_step, "_gtex_v7_skin_not_sun_exposed_beta_filter_real_v_matched_controls.txt")
        data <- read.table(skin_no_sun_file,header=TRUE)
        pvalues <- c(pvalues,data$real_pvalue)
        time_step <- c(time_step, rep(temp_time_step, length(data$real_pvalue)))
        version <- c(version, as.character(rep("skin_not_sun", length(data$real_pvalue))))  
        
    }
    # Put everything into data frame
    df <- data.frame(pvalues = -log10(as.numeric(pvalues + .000001)), version = factor(version,c("heart_left","breast","skin_not_sun","stomach")), time_step = factor(time_step))

    # PLOT!
    box_plot <- ggplot(df, aes(x=time_step, y=pvalues, fill=version)) + geom_boxplot(width=.7,outlier.size = 0.1)
    box_plot <- box_plot + theme(text = element_text(size=18), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black")) 
    box_plot <- box_plot + labs(fill= "Test version",x = "time step", y = "-log10(pvalue)")
    box_plot <- box_plot + theme(legend.position="bottom") 
    ggsave(box_plot, file=tissue_comparison_plot_file,width = 20,height=10.5,units="cm")


}


visualize_trajectories_for_one_cluster <- function(filtered_statistic_matrix, cluster_center, statistic_type, k) {
    filtered_statistic_matrix <- t(scale(t(filtered_statistic_matrix)))
    #filtered_statistic_matrix <- t(apply(filtered_statistic_matrix, 1, function(x)(x-min(x))/(max(x)-min(x))))

    group <- c()
    time_step <- c()
    statistics <- c()
    nrow <- dim(filtered_statistic_matrix)[1]
    for (row_num in 1:nrow) {
        time_step <- c(time_step,1:16)
        statistics <- c(statistics, as.numeric(filtered_statistic_matrix[row_num,]))
        group <- c(group, rep(row_num, 16))
    }
    df <- data.frame(statistic = statistics, time_step = time_step, group = factor(group))

    line_plot <- ggplot(df,aes(time_step,statistic,group=group))
    line_plot <- line_plot + geom_line(alpha=0.15)
    line_plot <- line_plot + theme(text = element_text(size=18), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black")) 
    line_plot <- line_plot + labs(x = "time step", y = statistic_type, title = paste0("Cluster ", k))


    return(line_plot)
}

kmeans_cluster_of_summary_statistics <- function(statistic_matrix, output_file, pc_num, k, statistic_type) {
    # Run Kmeans clustering
    kmeans_obj <- kmeans(statistic_matrix, k,nstart=30,iter.max=40)
    # Get centers of kmeans clusters
    # Matrix is of dimension k X (time_steps)
    centers <- kmeans_obj$centers
    # Get assignments of each sample to 1 of k clusters
    cluster_assignments <- kmeans_obj$cluster




    # Loop through clusters
    cluster_num <- 1
    # filter matrix to only samples belonging to current cluster
    filtered_statistic_matrix <- statistic_matrix[cluster_assignments == cluster_num,] 
    # get center of this cluster
    cluster_center <- centers[cluster_num,]
    # Make Line plot
    p1 <- visualize_trajectories_for_one_cluster(filtered_statistic_matrix, cluster_center, statistic_type, cluster_num)

    # Loop through clusters
    cluster_num <- 2
    # filter matrix to only samples belonging to current cluster
    filtered_statistic_matrix <- statistic_matrix[cluster_assignments == cluster_num,] 
    # get center of this cluster
    cluster_center <- centers[cluster_num,]
    # Make Line plot
    p2 <- visualize_trajectories_for_one_cluster(filtered_statistic_matrix, cluster_center, statistic_type, cluster_num)

    # Loop through clusters
    cluster_num <- 3
    # filter matrix to only samples belonging to current cluster
    filtered_statistic_matrix <- statistic_matrix[cluster_assignments == cluster_num,] 
    # get center of this cluster
    cluster_center <- centers[cluster_num,]
    # Make Line plot
    p3 <- visualize_trajectories_for_one_cluster(filtered_statistic_matrix, cluster_center, statistic_type, cluster_num)

    # Loop through clusters
    cluster_num <- 4
    # filter matrix to only samples belonging to current cluster
    filtered_statistic_matrix <- statistic_matrix[cluster_assignments == cluster_num,] 
    # get center of this cluster
    cluster_center <- centers[cluster_num,]
    # Make Line plot
    p4 <- visualize_trajectories_for_one_cluster(filtered_statistic_matrix, cluster_center, statistic_type, cluster_num)



    pdf(output_file)
    gg <- plot_grid(p1,p2,p3,p4,nrow=2,ncol=2,label_size=8)
    combined_gg <- ggdraw() + draw_plot(gg,0,0,1,1)
    print(combined_gg)
    dev.off()

}

visualize_number_of_genome_wide_significant_egenes <- function(input_stem, output_file) {
    num_genes <- c()
    for (temp_time_step in 0:15) {
        sig_egene_file <- paste0(input_stem, temp_time_step, "_efdr_thresh_.05_significant_egenes.txt")
        data <- read.table(sig_egene_file,header=TRUE)
        time_step_num_egenes <- dim(data)[1]
        num_genes <- c(num_genes, time_step_num_egenes)
    }
    print("MEAN NUM GENES")
    print(mean(num_genes))
    df <- data.frame(time_step=0:15, num_egenes=num_genes)
    p <- ggplot(df, aes(time_step, num_genes)) 
    p <- p + geom_bar(stat = "identity",aes(fill=time_step)) 
    p <- p + theme(text = element_text(size=8),axis.text=element_text(size=7), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.text = element_text(size=7), legend.title = element_text(size=8)) 
    p <- p + labs(x = "Day", y = "# significant genes")
    p <- p + scale_fill_gradient(low="darkgrey",high="firebrick")
    p <- p + theme(legend.position="none")
    ggsave(p, file=output_file, width=7.2, height=4.0, units="in")
    return(p)

}

eqtl_sharing_plot <- function(input_file, output_file, pc_num) {
    data = read.table(input_file)
    num_hits <- data$V3

    df <- data.frame(num_time_steps=num_hits)


    histo <- ggplot(data=df, aes(df$num_time_steps)) +
    geom_histogram(breaks=seq(.5, 16.5, by = 1),col="grey", fill="dodgerblue3") +
            labs(x="Number of days eQTL is significant (eFDR <= .05) in", y="Count") +
             theme(text = element_text(size=8),axis.text=element_text(size=7), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.text = element_text(size=7), legend.title = element_text(size=8)) 
    ggsave(histo, file=output_file, width=7.2, height=4.0, units="in")

}

make_figure_2 <- function(figure_2a, figure_2b, figure_2c, output_file) {
    figure_2a_legend <- get_legend(figure_2a)
    figure_2b_legend <- get_legend(figure_2b)
    figure_2c_legend <- get_legend(figure_2c)

    combined <- ggdraw() + 
                draw_plot(figure_2a + theme(legend.position='none'), -.0098,.58,1,.4) +
                draw_plot(figure_2a_legend,.72,.46,1,1) +
                draw_plot(figure_2b + theme(legend.position='none'),-.0028,0,.43,.58) + 
                draw_plot(figure_2b_legend,.41,-.14,1,1) +
                draw_plot(figure_2c + theme(legend.position='none'),.5,0,.43,.58) +
                draw_plot(figure_2c_legend,.91,-.14,1,1) +
                draw_plot_label(c("A","B","C"),c(.01,.01,.52),c(1,.6,.6),size=12)
    ggsave(combined, file=output_file, width=7.2, height=5.0,units="in")


}

number_of_significant_genes_merge_plots <- function(num_qtl_violin, num_qtl_bar_plot, pc_pve_line_plot, output_file) {

    combined <- plot_grid(pc_pve_line_plot, num_qtl_violin, num_qtl_bar_plot, labels=c("A","B","C"), ncol=1)

    ggsave(combined,file=output_file, width=7.2, height=6.3,units="in")

}

###############################################
# Line Plot showing PVE of PCs in each time step
per_time_step_pcs_pve_line_plot <- function(input_stem, output_file, n) {
    pve_arr <- c()
    time_step_arr <- c()
    pc_num_arr <- c()
    # Loop through time steps
    for (time_step in 0:15) {
        pve_file <- paste0(input_stem, time_step, "_variance_explained.txt")
        data <- read.table(pve_file,header=FALSE)$V1
        for (pc_num in 1:n) {
            pve_arr <- c(pve_arr, data[pc_num])
            time_step_arr <- c(time_step_arr, time_step)
            pc_num_arr <- c(pc_num_arr, pc_num)
        }
    }
    # Put data into nice clean data frame
    df <- data.frame(pve=pve_arr, time=time_step_arr,pc=pc_num_arr)

    # PLOT AWAY
    line_plot <- ggplot(data=df, aes(x=pc, y=pve, group=factor(time), color=time)) +
                geom_line() +
                geom_point() +
                scale_color_gradient(low="darkgrey",high="firebrick") + 
                ylim(0,max(pve_arr) + .01) + 
                scale_x_continuous(breaks=1:n) +
                labs(x = "PC number", y = "Variance Explained", color="Day") + 
                 theme(plot.title=element_text(size=8,face="plain"),text = element_text(size=8),axis.text=element_text(size=7), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.text = element_text(size=7), legend.title = element_text(size=8), axis.text.x = element_text(vjust=.5)) 

    ggsave(line_plot, file=output_file, width=7.2, height=5.0,units="in")
    return(line_plot)

}



parameter_string = args[1]  # string used to keep track of files used with specified parameter settting
cht_output_dir = args[2]  # input directory with cht test results
cht_visualization_dir = args[3]  # output directory to save images
matrix_factorization_dir = args[4]  # Directory containing results from matrix factorization analysis
cht_input_file_dir = args[5]  # Directory containing PCs from time step independent analysis


pc_num <- 3


###############################################
# Violinplot of number of egenes as a function of number of pcs
#  Each point in violin is a time step
output_file <- paste0(cht_visualization_dir, parameter_string, "egenes_as_function_of_pcs_violinplot.png")
#num_qtl_violin <- egenes_as_function_of_pcs_violinplot(cht_output_dir, parameter_string, output_file)

###############################################
# Bar plot showing number of genome wide significant egenes at each of the 16 time steps
# Do this for each of the pc_nums in dependently
output_file <- paste0(cht_visualization_dir, parameter_string,"_num_pc_",pc_num,"_number_of_significant_egenes_per_time_step_bar_plot.png")
input_stem <- paste0(cht_output_dir, "cht_results_",parameter_string,"_num_pc_",pc_num,"_time_")
#num_qtl_bar_plot <- visualize_number_of_genome_wide_significant_egenes(input_stem, output_file)


###############################################
# Line Plot showing PVE of PCs in each time step
n <- 10
output_file <- paste0(cht_visualization_dir, parameter_string,"_per_time_step_pcs_pve_line_plot.png")
input_stem <- paste0(cht_input_file_dir, "pcs_cis_distance_50000_maf_cutoff_0.1_min_reads_100_min_as_reads_25_min_het_counts_5_time_")
#pc_pve_line_plot <- per_time_step_pcs_pve_line_plot(input_stem, output_file, n)


###############################################
# Merged plot showing num_qtl_violin and num_qtl_bar_plot side by side
output_file <- paste0(cht_visualization_dir, parameter_string,"_num_pc_",pc_num,"_number_of_significant_genes_merged_plot.pdf")
#number_of_significant_genes_merge_plots(num_qtl_violin, num_qtl_bar_plot, pc_pve_line_plot, output_file)


###############################################
# Histogram showing the number of time points each per-time step eqtl is significant in 
input_file <- paste0(cht_output_dir,parameter_string,"_num_pc_",pc_num,"_fdr_.05_eqtl_sharing.txt")
output_file <- paste0(cht_visualization_dir, parameter_string,"_num_pc_",pc_num,"_eqtl_sharing_histogram.pdf")
#eqtl_sharing_plot(input_file, output_file, pc_num)


###############################################
# QQPlot showing both:
## 1. Real-data pvalues compared to uniform
## 2. Permuted-data pvalues compared to uniform
## Done independently for each PC 
## Each output image has 16 subplots (1 for each time step)

# Output file
output_file <- paste0(cht_visualization_dir, parameter_string,"_num_pc_",pc_num,"_qq_plot_vs_uniform.pdf")
# Input file stems
input_stem <- paste0(cht_output_dir, "cht_results_",parameter_string,"_num_pc_",pc_num,"_time_")
null_stem <- paste0(cht_output_dir,"cht_perm1_results_",parameter_string,"_num_pc_",pc_num,"_time_")
#qq_plot_vs_uniform(input_stem,null_stem,output_file)


###############################################
# QQPlot showing both:
## 1. Real-data pvalues compared to permuted
## Done independently for each PC 
## Each output image has 16 subplots (1 for each time step)
    
# Output file
output_file <- paste0(cht_visualization_dir, parameter_string,"_num_pc_",pc_num,"_qq_plot_vs_permuted.pdf")
# Input file stems
input_stem <- paste0(cht_output_dir, "cht_results_",parameter_string,"_num_pc_",pc_num,"_time_")
null_stem <- paste0(cht_output_dir,"cht_perm1_results_",parameter_string,"_num_pc_",pc_num,"_time_")
#qq_plot_vs_permuted(input_stem,null_stem,output_file)




###############################################
# Line plot showing spearman correlation between banovich ipsc/cm results (two colors) and each of the 16 time points (x-axis)

#figure_2a <- line_plot_of_spearman_correlation_with_banovich_results_across_time(parameter_string, pc_num, cht_output_dir)
output_file <- paste0(cht_visualization_dir, parameter_string, "_num_pc_", pc_num,"_spearman_correlation_with_banovich_results_line_plot.png")
#ggsave(figure_2a, file=output_file, width=7.2, height=3.0, units="in")




###############################################
# Heatmap of correlation of summary statistics between time steps
# Ie. heatmap is of dimension number of time steps by number of time steps
# Do independently for each number of PCs

#figure_2b <- summary_statistic_correlation_heatmap(parameter_string, pc_num, cht_output_dir, cht_visualization_dir)
output_file <- paste0(cht_visualization_dir, parameter_string, "_num_pc_", pc_num,"_pvalue_correlation_heatmap.png")
#ggsave(figure_2b, file=output_file, width=7.2, height=5.3, units="in")




###############################################
# Heatmap showing factor matrix from sparse matrix factor analysis

sparsity_parameter = "0.5"
num_factors = "3"
factor_matrix_file <- paste0(matrix_factorization_dir,parameter_string,"_num_pc_", pc_num, "_fdr_.05_log_pvalue_factorization_alpha_", sparsity_parameter, "_", num_factors,"_loading_matrix.txt")
#figure_2c <- factor_matrix_heatmap(factor_matrix_file)
output_file <- paste0(cht_visualization_dir, parameter_string, "_num_pc_", pc_num,"_pvalue_factor_matrix_",sparsity_parameter,"_",num_factors,".png")
#ggsave(figure_2c, file=output_file, width=7.2, height=5.3, units="in")





################################################
# Make combined plot for figure 2
output_file <- paste0(cht_visualization_dir, "figure2.pdf")
#make_figure_2(figure_2a, figure_2b, figure_2c, output_file)



###############################################
# Grid of heatmap showing factor matrix from sparse matrix factor analysis
# Heatmaps will span number of latent factors and sparse prior choice
    
factor_matrix_root <- paste0(matrix_factorization_dir, parameter_string, "_num_pc_", pc_num,"_fdr_.05_log_pvalue_factorization_alpha_")
output_file <- paste0(cht_visualization_dir, parameter_string, "_num_pc_", pc_num,"_log_pvalue_factor_matrices.pdf")
make_grid_of_factor_matrices(factor_matrix_root, output_file)



###############################################
# Grid of heatmap showing factor matrix from sparse matrix factor analysis
# Heatmaps will span number of latent factors and sparse prior choice
    
factor_matrix_root <- paste0(matrix_factorization_dir, parameter_string, "_num_pc_", pc_num,"_fdr_.05_allelic_fraction_factorization_alpha_")
output_file <- paste0(cht_visualization_dir, parameter_string, "_num_pc_", pc_num,"_allelic_fraction_factor_matrices.pdf")
make_grid_of_factor_matrices(factor_matrix_root, output_file)

























###########################################
# Old (retired scripts)
###########################################

###############################################
# Heatmap of correlation of summary statistics between time steps and between gtex tissues
# Ie. heatmap is of dimension number of time steps+num_gtex_tissues by number of time steps+num_gtex_tissues
# Do independently for each number of PCs
for (pc_num in 0:5) {
    #summary_statistic_correlation_heatmap_include_gtex(parameter_string, pc_num, cht_output_dir, cht_visualization_dir)
}



###############################################
# Boxplot taking significant top n variant gene pairs in each time step, and plotting -log10 pvalue of those variant gene pairs for banovich ipscs and ipsc-cms
num_genes <- 100
for (pc_num in 3:3) {
    input_file <- paste0(cht_output_dir, parameter_string, "_num_pc_", pc_num, "_top_", num_genes, "_eqtls_in_time_steps_for_banovich_results.txt")
    output_file <- paste0(cht_visualization_dir, parameter_string, "_top_", num_genes, "_eqtls_in_time_steps_for_banovich_results_boxplot.pdf")
    #boxplot_showing_top_n_variant_gene_pairs_per_time_step_with_banovich_results(input_file, output_file, num_genes, pc_num)
}



###############################################
# Boxplot of number of egenes as a function of number of pcs
#  Each point in boxplot is a time step
# output_file <- paste0(cht_visualization_dir, parameter_string, "egenes_as_function_of_pcs_boxplot.png")
# egenes_as_function_of_pcs_boxplot(cht_output_dir, parameter_string, output_file)


###############################################
# Boxplot showing pvalues found in our data for only the eqtls in a specific data set
# 1 plot per data set
# 1 plot for pc
for (pc_num in 0:3) {
    #eqtl_comparison_to_background_shell(pc_num, cht_visualization_dir, parameter_string, cht_enrichment_dir, eqtl_data_set_file)
}

###############################################
# Boxplot showing pvalues found in our data for only the eqtls in a specific data set
# 1 plot per data set
# 1 plot for pc
for (pc_num in 0:3) {
    #eqtl_comparison_reference_shell(pc_num, cht_visualization_dir, parameter_string, cht_enrichment_dir, eqtl_data_set_file)
}



###############################################
# Boxplot showing pvalues found in our data based on eqtls found in GTEx v7:
## 1. HLV
## 2. skin_not_sun_exposed
## 3. stomach
## 4. Breast
# 1 plot for pc
for (pc_num in 0:3) {
    tissue_comparison_plot_file <- paste0(cht_visualization_dir,parameter_string,"_num_pc_",pc_num, "_gtex_tissue_beta_filter_comparison.png")
    #gtex_tissue_comparison_bar_plot(tissue_comparison_plot_file, pc_num, parameter_string, cht_enrichment_dir)
}




