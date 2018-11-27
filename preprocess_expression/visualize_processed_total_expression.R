args = commandArgs(trailingOnly=TRUE)
library(ggplot2)
library(ggthemes)
library(glmnet)
library(reshape)
library(reshape2)
library(cowplot)
library(RColorBrewer)
library(ica)
library(gplots)


#  BETABINOMIAL_GLM=stan_model(file="sparse_betabinomial_glm.stan", save_dso=T, auto_write=T)

#  Plot first two PC's. Color points by time step
plot_pca_time_step <- function(sample_info, quant_expr, output_file) {

    #  Compute singular value decomposition
    svd1 <- svd(as.matrix(quant_expr))

    #  Scores of first 2 pc's across all samples
    pc1 <- svd1$v[,1]
    pc2 <- svd1$v[,2]

    # Put all information into data structure
    df <- data.frame(pc1 = pc1, pc2 = pc2, time_step = sample_info$time)

    #PLOT!
    pca_scatter <-  ggplot(df,aes(pc1,pc2)) + geom_point(aes(colour=time_step)) + theme(text = element_text(size=18), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black")) 
    pca_scatter <- pca_scatter + scale_color_gradient(low="pink",high="blue")


    ggsave(pca_scatter, file=output_file,width = 15,height=10.5,units="cm",dpi=600)

}

#  Plot first two PC's. Color points by time step
plot_pca_time_step_modular <- function(sample_info, quant_expr) {

    #  Compute singular value decomposition
    svd1 <- svd(as.matrix(quant_expr))

    #  Scores of first 2 pc's across all samples
    pc1 <- svd1$v[,1]
    pc2 <- svd1$v[,2]

    # Put all information into data structure
    df <- data.frame(pc1 = pc1, pc2 = pc2, time_step = sample_info$time)

    #PLOT!
    pca_scatter <-  ggplot(df,aes(pc1,pc2)) + geom_point(aes(colour=time_step),size=.7) 
    pca_scatter <- pca_scatter + scale_color_gradient(low="darkgrey",high="firebrick")#scale_color_gradient(low="darkgrey",high="firebrick")
    pca_scatter <- pca_scatter +  labs(colour="Day",x = "PC1", y = "PC2") + theme(legend.position="left") + theme(legend.key.size =  unit(0.1, "in"))
    pca_scatter <- pca_scatter + theme(text = element_text(size=8),axis.text=element_text(size=7), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.text = element_text(size=7), legend.title = element_text(size=8)) 
    return(pca_scatter)
}

plot_pca_old_vs_new <- function(sample_info, quant_expr, output_file) {
    #  Compute singular value decomposition
    svd1 <- svd(as.matrix(quant_expr))

    #  Scores of first 2 pc's across all samples
    pc1 <- svd1$v[,1]
    pc2 <- svd1$v[,2]

    # Put all information into data structure
    df <- data.frame(pc1 = pc1, pc2 = pc2, old_vs_new=factor(sample_info$old_vs_new))

    #PLOT!
    pca_scatter <-  ggplot(df,aes(pc1,pc2)) + geom_point(aes(colour=old_vs_new)) + theme(text = element_text(size=18), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black")) 



    ggsave(pca_scatter, file=output_file,width = 15,height=10.5,units="cm",dpi=600)

}




#  Plot first two PC's. Color points by time step
plot_ica_time_step <- function(sample_info, quant_expr, output_file) {

    #  Compute singular value decomposition
    imod <- icafast(as.matrix(quant_expr),nc=5)

    #  Scores of first 2 pc's across all samples
    ic1 <- imod$M[,1]
    ic2 <- imod$M[,2]

    # Put all information into data structure
    df <- data.frame(ic1 = ic1, ic2 = ic2, time_step = sample_info$time)

    #PLOT!
    ica_scatter <-  ggplot(df,aes(ic1,ic2)) + geom_point(aes(colour=time_step)) + theme(text = element_text(size=18), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black")) 
    ica_scatter <- ica_scatter + scale_color_gradient(low="pink",high="blue")


    ggsave(ica_scatter, file=output_file,width = 15,height=10.5,units="cm")

}

make_one_eigenvector_plot <- function(quant_expr, actual_pc_num, sample_info) {
    #  Compute singular value decomposition
    svd1 <- svd(as.matrix(quant_expr))

    #  Scores of first 2 pc's across all samples
    pc <- svd1$v[,(actual_pc_num + 1)] 


    df <- data.frame(pc = pc, time_step = sample_info$time)

    pca_scatter <-  ggplot(df,aes(time_step,pc)) + geom_point(aes(colour=time_step)) + theme(text = element_text(size=18), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black")) 
    pca_scatter <- pca_scatter + scale_color_gradient(low="pink",high="blue")
    pca_scatter <- pca_scatter + labs(colour="time step",x = "time step", y = paste0("PC",(actual_pc_num+1)), title = paste0("PC",(actual_pc_num+1)))

    return(pca_scatter)
}

make_one_eigenvector_plot2 <- function(quant_expr, actual_pc_num, sample_info) {
    #  Compute singular value decomposition
    svd1 <- svd(as.matrix(quant_expr))

    #  Scores of first 2 pc's across all samples
    pc <- svd1$v[,(actual_pc_num + 1)] 


    df <- data.frame(pc = pc, time_step = sample_info$time, cell_line = factor(sample_info$cell_line))

    colourCount = length(unique(df$cell_line))
    getPalette = colorRampPalette(brewer.pal(9, "Set1"))


    pca_scatter <-  ggplot(df,aes(time_step,pc)) + geom_point(aes(colour=cell_line)) + theme(text = element_text(size=18), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black")) 
    pca_scatter <- pca_scatter + scale_colour_manual(values = getPalette(colourCount))
    pca_scatter <- pca_scatter + labs(colour="cell line",x = "time step", y = paste0("PC",(actual_pc_num+1)), title = paste0("PC",(actual_pc_num+1)))

    return(pca_scatter)
}

plot_pca_eigenvectors <- function(sample_info, quant_expr, eigenvectors_output_file) {
    t0 <- make_one_eigenvector_plot(quant_expr, 0, sample_info)
    t1 <- make_one_eigenvector_plot(quant_expr, 1, sample_info)+ theme(legend.position="none")
    t2 <- make_one_eigenvector_plot(quant_expr, 2, sample_info)+ theme(legend.position="none")
    t3 <- make_one_eigenvector_plot(quant_expr, 3, sample_info)+ theme(legend.position="none")


    legend <- get_legend(t0)
    pdf(eigenvectors_output_file)
    gg <- plot_grid(t0+ theme(legend.position="none"),t1,t2,t3,nrow=2,ncol=2,label_size=10)
    combined_gg <- ggdraw() + draw_plot(gg,0,0,1,1) 
    print(combined_gg)
    dev.off()


}

plot_pca_eigenvectors_by_line <- function(sample_info, quant_expr, eigenvectors_output_file) {
    t0 <- make_one_eigenvector_plot2(quant_expr, 0, sample_info)
    t1 <- make_one_eigenvector_plot2(quant_expr, 1, sample_info)+ theme(legend.position="none")
    t2 <- make_one_eigenvector_plot2(quant_expr, 2, sample_info)+ theme(legend.position="none")
    t3 <- make_one_eigenvector_plot2(quant_expr, 3, sample_info)+ theme(legend.position="none")


    legend <- get_legend(t0)
    pdf(eigenvectors_output_file)
    gg <- plot_grid(t0+ theme(legend.position="none"),t1,t2,t3,nrow=2,ncol=2,label_size=10)
    combined_gg <- ggdraw() + draw_plot(gg,0,0,1,1) 
    print(combined_gg)
    dev.off()


}

make_eigenvector_plot <- function(i, cell_lines, pc, sample_info) {
    # Get name of ith cell line
    i_cell_line <- cell_lines[i]

    # Get indices  of all samples that are from ith cell line
    i_indices <- sample_info$cell_line == i_cell_line

    #  Extract 1st two pcs of all samples that belong to the ith cell lein
    i_pc <- pc[i_indices]

    # Get time steps of these samples
    i_time <- sample_info$time[i_indices]

    # Put into compact data frame for plotting
    df <- data.frame(pc= i_pc, time_step = i_time)

    #PLOT!
    pca_scatter <-  ggplot(df,aes(time_step,pc)) + geom_point(aes(colour=time_step)) + theme(text = element_text(size=12), axis.text = element_text(size=8), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black")) 
    pca_scatter <- pca_scatter + scale_color_gradient(low="pink",high="blue") + ggtitle(i_cell_line) + ylim(min(pc)-.01,max(pc) + .01)

    return(pca_scatter)
}

plot_pca_eigenvectors_by_line_independent <- function(sample_info, quant_expr, eigenvectors_output_file, actual_pc_num) {
    #  Compute singular value decomposition
    svd1 <- svd(as.matrix(quant_expr))

    #  Scores of first 2 pc's across all samples
    pc <- svd1$v[,(actual_pc_num)] 
   

    #  Get unique cell lines
    cell_lines <- sort(unique(sample_info$cell_line))
    num_cell_lines <- length(cell_lines)


    
    #  Make pc plot for each cell line seperately (not automated yet..
     
    p1 <- make_eigenvector_plot(1, cell_lines, pc, sample_info) + labs(x="time step",y=paste0("PC",actual_pc_num)) # cell line 1
    p2 <- make_eigenvector_plot(2, cell_lines, pc, sample_info) + labs(x="time step",y=paste0("PC",actual_pc_num))+ theme(legend.position="none") # cell line 1
    p3 <- make_eigenvector_plot(3, cell_lines, pc, sample_info) + labs(x="time step",y=paste0("PC",actual_pc_num))+ theme(legend.position="none") # cell line 1
    p4 <- make_eigenvector_plot(4, cell_lines, pc, sample_info) + labs(x="time step",y=paste0("PC",actual_pc_num))+ theme(legend.position="none") # cell line 1
    p5 <- make_eigenvector_plot(5, cell_lines, pc, sample_info) + labs(x="time step",y=paste0("PC",actual_pc_num))+ theme(legend.position="none") # cell line 1
    p6 <- make_eigenvector_plot(6, cell_lines, pc, sample_info) + labs(x="time step",y=paste0("PC",actual_pc_num))+ theme(legend.position="none") # cell line 1
    p7 <- make_eigenvector_plot(7, cell_lines, pc, sample_info) + labs(x="time step",y=paste0("PC",actual_pc_num))+ theme(legend.position="none") # cell line 1
    p8 <- make_eigenvector_plot(8, cell_lines, pc, sample_info) + labs(x="time step",y=paste0("PC",actual_pc_num))+ theme(legend.position="none") # cell line 1
    p9 <- make_eigenvector_plot(9, cell_lines, pc, sample_info) + labs(x="time step",y=paste0("PC",actual_pc_num))+ theme(legend.position="none") # cell line 1
    p10 <- make_eigenvector_plot(10, cell_lines, pc, sample_info) + labs(x="time step",y=paste0("PC",actual_pc_num))+ theme(legend.position="none") # cell line 1
    p11 <- make_eigenvector_plot(11, cell_lines, pc, sample_info) + labs(x="time step",y=paste0("PC",actual_pc_num))+ theme(legend.position="none") # cell line 1
    p12 <- make_eigenvector_plot(12, cell_lines, pc, sample_info) + labs(x="time step",y=paste0("PC",actual_pc_num))+ theme(legend.position="none") # cell line 1
    p13 <- make_eigenvector_plot(13, cell_lines, pc, sample_info) + labs(x="time step",y=paste0("PC",actual_pc_num))+ theme(legend.position="none") # cell line 1
    p14 <- make_eigenvector_plot(14, cell_lines, pc, sample_info) + labs(x="time step",y=paste0("PC",actual_pc_num))+ theme(legend.position="none") # cell line 1



    legend <- get_legend(p1)

    # Merge all cell lines into one plot
    pdf(eigenvectors_output_file)
    gg <- plot_grid(p1 + theme(legend.position="none"),p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12,p13,p14,nrow=4,ncol=4)
    combined_gg <- ggdraw() + draw_plot(gg,0,0,1,1) + draw_plot(legend,.83,0,1,.3)
    print(combined_gg)
    dev.off()

}




plot_library_size <- function(sample_info, library_size_output_file) {
    #  Get unique cell lines
    cell_lines <- sort(unique(sample_info$cell_line))
    num_cell_lines <- length(cell_lines)

    #  Initialize arrays to store information
    ordered_library_sizes <- c()
    ordered_cell_lines <- c()
    ordered_time_steps <- c()
    ordered_sample_name <- c()
    
    # Loop through cell lines and samples
    for (cell_line_iter in 1:num_cell_lines) {
        for (time_step in 0:15) {
            # Extract cell line of ith iteration
            cell_line <- cell_lines[cell_line_iter]
            # Find which row this corresponds to
            correct_row <- which(sample_info$cell_line == cell_line & sample_info$time == time_step)
            # Check to make sure cellLine_timeStep exists in our data
            if (length(correct_row) > 0) {
                # If it does, append arrays
                ordered_library_sizes <- c(ordered_library_sizes,sample_info$lib.size[correct_row])
                ordered_cell_lines <- c(ordered_cell_lines, sample_info$cell_line[correct_row])
                ordered_time_steps <- c(ordered_time_steps, sample_info$time[correct_row])
                ordered_sample_name <- c(ordered_sample_name, sample_info$Sample_name[correct_row])
            }
        }
    }
    df <- data.frame(library_size = ordered_library_sizes, cell_line = factor(paste0("NA",ordered_cell_lines)), time_step = ordered_time_steps, sample_name = ordered_sample_name)

    colourCount = length(unique(df$cell_line))
    getPalette = colorRampPalette(brewer.pal(9, "Set1"))

    bar_plot <- ggplot(df, aes(x=sample_name, y=library_size, fill=cell_line)) + geom_bar(stat="identity")
    bar_plot <- bar_plot + scale_fill_manual(values = getPalette(colourCount))

    bar_plot <- bar_plot + theme(text = element_text(size=8),axis.text=element_text(size=7), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.text = element_text(size=7), legend.title = element_text(size=8)) 
    bar_plot <- bar_plot + labs(fill="Cell Line",x = "Sample", y = "Library Size") + theme(legend.key.size =  unit(0.16, "in"))

    ggsave(bar_plot, file=library_size_output_file, width=7.2, height=4.0, units="in")


}


plot_pca_real_valued_gene_filled <- function(sample_info, quant_expr, ensamble_id, gene_name, pc_num1, pc_num2, output_file) {
    #  Compute singular value decomposition
    svd1 <- svd(as.matrix(quant_expr))

    #  Scores of first 2 pc's across all samples
    pc1 <- svd1$v[,pc_num1]
    pc2 <- svd1$v[,pc_num2]

    row_label <- which(rownames(quant_expr) == ensamble_id)
    quant_expr <- as.matrix(quant_expr)

    # Put all information into data structure
    df <- data.frame(pc1 = pc1, pc2 = pc2, time_step = as.vector(quant_expr[row_label,]))


    #PLOT!
    pca_scatter <-  ggplot(df,aes(pc1,pc2)) + geom_point(aes(colour=time_step)) + theme(text = element_text(size=18), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black")) 
    pca_scatter <- pca_scatter + scale_color_gradient(low="pink",high="blue") + labs(colour="Expression",x = paste0("PC",pc_num1), title = gene_name,y = paste0("PC",pc_num2))


    ggsave(pca_scatter, file=output_file,width = 15,height=14,units="cm")
}

plot_cell_line_pca_real_valued_gene_filled <- function(cell_line_names, cell_line_expr, sample_info, quant_expr, ensamble_id,gene_name,pc_num1,pc_num2, output_file, time_step) {
    #  Compute singular value decomposition
    svd1 <- svd(as.matrix(cell_line_expr))

    #  Scores of first 2 pc's across all samples
    pc1 <- svd1$v[,pc_num1]
    pc2 <- svd1$v[,pc_num2]

    row_label <- which(rownames(quant_expr) == ensamble_id)
    quant_expr <- as.matrix(quant_expr)
    genes_specific_expr <- as.vector(quant_expr[row_label,])

    num_cell_lines <- length(cell_line_names)
    cell_line_expr_values <- c()
    for (counter in 1:num_cell_lines) {
        curr_cell_line <- cell_line_names[counter]
        index <- which(sample_info$cell_line == curr_cell_line & sample_info$time == time_step)
        cell_line_expr_values <- c(cell_line_expr_values, genes_specific_expr[index])
    }

    # Put all information into data structure
    df <- data.frame(pc1 = pc1, pc2 = pc2, time_step = cell_line_expr_values)


    #PLOT!
    pca_scatter <-  ggplot(df,aes(pc1,pc2)) + geom_point(aes(colour=time_step)) + theme(text = element_text(size=18), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black")) 
    pca_scatter <- pca_scatter + scale_color_gradient(low="pink",high="blue") + labs(colour="Expression",x = paste0("PC",pc_num1), title = paste0("t=", time_step, " ", gene_name),y = paste0("PC",pc_num2))


    ggsave(pca_scatter, file=output_file,width = 15,height=10.5,units="cm")


}



plot_ica_real_valued_gene_filled <- function(sample_info, quant_expr, ensamble_id, gene_name, pc_num1, pc_num2, output_file) {

    #  Compute singular value decomposition
    imod <- icafast(as.matrix(quant_expr),nc=5)

    #  Scores of first 2 pc's across all samples
    pc1 <- imod$M[,pc_num1]
    pc2 <- imod$M[,pc_num2]

    row_label <- which(rownames(quant_expr) == ensamble_id)
    quant_expr <- as.matrix(quant_expr)

    # Put all information into data structure
    df <- data.frame(pc1 = pc1, pc2 = pc2, time_step = as.vector(quant_expr[row_label,]))


    #PLOT!
    pca_scatter <-  ggplot(df,aes(pc1,pc2)) + geom_point(aes(colour=time_step)) + theme(text = element_text(size=18), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black")) 
    pca_scatter <- pca_scatter + scale_color_gradient(low="pink",high="blue") + labs(colour="Expression",x = paste0("PC",pc_num1), title = gene_name,y = paste0("PC",pc_num2))


    ggsave(pca_scatter, file=output_file,width = 15,height=10.5,units="cm")
}


gene_time_course_line_plot_grouped_by_cell_line <- function(sample_info, quant_expr, ensamble_id, gene_name, line_plot_file) {
    row_label <- which(rownames(quant_expr) == ensamble_id)
    quant_expr <- as.matrix(quant_expr)

    df <- data.frame(time = sample_info$time, expression = as.vector(quant_expr[row_label,]), cell_line = factor(paste0("NA",sample_info$cell_line)))

    colourCount = length(unique(df$cell_line))
    getPalette = colorRampPalette(brewer.pal(9, "Set1"))

    line_plot <- ggplot(df, aes(x=time, y=expression, group=cell_line)) + geom_line(aes(color=cell_line)) +
                geom_point(aes(color=cell_line)) +
                scale_colour_manual(values = getPalette(colourCount)) +
                theme(text = element_text(size=18), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black")) +
                labs(colour="Cell Line",x = "Time Step", y = "Normalized Expression", title = gene_name)
    ggsave(line_plot, file=line_plot_file,width = 15.9,height=11.5,units="cm")

}

gene_time_course_line_plot_grouped_by_cell_line_modular <- function(sample_info, quant_expr, ensamble_id, gene_name, line_plot_file) {
    row_label <- which(rownames(quant_expr) == ensamble_id)
    quant_expr <- as.matrix(quant_expr)

    df <- data.frame(time = sample_info$time, expression = as.vector(quant_expr[row_label,]), cell_line = factor(paste0("NA",sample_info$cell_line)))

    colourCount = length(unique(df$cell_line))
    getPalette = colorRampPalette(brewer.pal(9, "Set1"))

    line_plot <- ggplot(df, aes(x=time, y=expression, group=cell_line)) + geom_line(aes(color=cell_line)) +
                geom_point(aes(color=cell_line)) +
                scale_colour_manual(values = getPalette(colourCount)) + theme(plot.title = element_text(size=12)) + 
                theme(text = element_text(size=12), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black")) +
                labs(colour="Cell Line",x = "Day", y = "Expression", title = gene_name) +
                theme(plot.title=element_text(size=8,face="plain"),text = element_text(size=8),axis.text=element_text(size=7), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.text = element_text(size=7), legend.title = element_text(size=8)) 
    return(line_plot)

}



gene_time_course_line_plot_grouped_by_cell_line_old_vs_new <- function(sample_info, quant_expr, ensamble_id, gene_name, line_plot_file) {
    row_label <- which(rownames(quant_expr) == ensamble_id)
    quant_expr <- as.matrix(quant_expr)

    df <- data.frame(time = sample_info$time, expression = as.vector(quant_expr[row_label,]),old_vs_new=factor(sample_info$old_vs_new), cell_line = factor(paste0("NA",sample_info$cell_line)))
    line_plot <- ggplot(df, aes(x=time, y=expression, group=cell_line)) + geom_line(aes(color=old_vs_new)) +
                geom_point(aes(color=old_vs_new)) +
                theme(text = element_text(size=18), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black")) +
                labs(colour="old_vs_new",x = "Time Step", y = "Normalized Expression", title = gene_name)
    ggsave(line_plot, file=line_plot_file,width = 15.9,height=11.5,units="cm")


}




pc_gene_scatter <- function(sample_info, quant_expr, ensamble_id, gene_name, time_step, pc_num, pc_gene_scatter_output_file) {
    #  Compute singular value decomposition
    svd1 <- svd(as.matrix(quant_expr))

    #  Scores of pcs of interest
    pc_scores <- svd1$v[,pc_num]

    # Row corresponding to gene of interest
    row_label <- which(rownames(quant_expr) == ensamble_id)
    # Get matrix into correct format
    quant_expr <- as.matrix(quant_expr)


    df <- data.frame(time = sample_info$time, expression = as.vector(quant_expr[row_label,]), pc_scores = pc_scores, cell_line = factor(sample_info$cell_line))
    
    df <- df[df$time == time_step,]

    #PLOT!
    pca_scatter <- ggplot(df, aes(x = expression, y = pc_scores, colour = cell_line)) + geom_point() 
    pca_scatter <- pca_scatter + theme(text = element_text(size=18), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black")) 
    pca_scatter <- pca_scatter + labs(colour="Cell Line",x = paste0(gene_name, " expression"), y = paste0("PC",pc_num), title = paste0("Time step ", time_step))
    ggsave(pca_scatter, file=pc_gene_scatter_output_file,width = 15,height=10.5,units="cm")

}


pc_avg_troponin_scatter_colored_by_cell_line <- function(sample_info, quant_expr, pc_num, covariate_file, output_file) {
    #  Compute singular value decomposition
    svd1 <- svd(as.matrix(quant_expr))

    #  Scores of pcs of interest
    pc_scores <- svd1$v[,pc_num]
    
    covs <- read.table(covariate_file, header=TRUE)

    # Remove unimportant columns
    covs <- data.frame(as.matrix(covs[,3:dim(covs)[2]]))

    covs$avg_10_15_troponin_mRNA_expr <- abs(as.numeric(as.character(covs$avg_10_15_troponin_mRNA_expr)))
    covs$avg_10_15_sox2_mRNA_expr <- abs(as.numeric(as.character(covs$avg_10_15_sox2_mRNA_expr)))

    df <- data.frame(pc_scores = pc_scores, avg_troponin=covs$avg_10_15_troponin_mRNA_expr, cell_line = factor(sample_info$cell_line),time_step = sample_info$time)

    pca_scatter <- ggplot(df, aes(y = pc_scores, x = avg_troponin, colour = cell_line)) + geom_point() 
    pca_scatter <- pca_scatter + theme(text = element_text(size=18), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black")) 
    pca_scatter <- pca_scatter + labs(colour="Cell Line",x = "avg (t=10-15) troponin mRNA expr", y = paste0("PC",pc_num))
    ggsave(pca_scatter, file=output_file,width = 15,height=10.5,units="cm")

}

pc_avg_troponin_scatter_colored_by_time_step <- function(sample_info, quant_expr, pc_num, covariate_file, output_file) {
    #  Compute singular value decomposition
    svd1 <- svd(as.matrix(quant_expr))

    #  Scores of pcs of interest
    pc_scores <- svd1$v[,pc_num]
    
    covs <- read.table(covariate_file, header=TRUE)

    # Remove unimportant columns
    covs <- data.frame(as.matrix(covs[,3:dim(covs)[2]]))

    covs$avg_10_15_troponin_mRNA_expr <- abs(as.numeric(as.character(covs$avg_10_15_troponin_mRNA_expr)))
    covs$avg_10_15_sox2_mRNA_expr <- abs(as.numeric(as.character(covs$avg_10_15_sox2_mRNA_expr)))

    df <- data.frame(pc_scores = pc_scores, avg_troponin=covs$avg_10_15_troponin_mRNA_expr, cell_line = factor(sample_info$cell_line),time_step = sample_info$time)


    pca_scatter <- ggplot(df, aes(y = pc_scores, x = avg_troponin, colour = time_step)) + geom_point() 
    pca_scatter <- pca_scatter + theme(text = element_text(size=18), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black")) 
    pca_scatter <- pca_scatter + labs(colour="time step",x = "avg (t=10-15) troponin mRNA expr", y = paste0("PC",pc_num)) + scale_color_gradient(low="pink",high="blue")
    ggsave(pca_scatter, file=output_file,width = 15,height=10.5,units="cm")

}



pc_gene_scatter_all_time_colored_by_cell_line <- function(sample_info, quant_expr, ensamble_id, gene_name, pc_num, pc_gene_scatter_output_file) {
    #  Compute singular value decomposition
    svd1 <- svd(as.matrix(quant_expr))

    #  Scores of pcs of interest
    pc_scores <- svd1$v[,pc_num]

    # Row corresponding to gene of interest
    row_label <- which(rownames(quant_expr) == ensamble_id)
    # Get matrix into correct format
    quant_expr <- as.matrix(quant_expr)


    df <- data.frame(time = sample_info$time, expression = as.vector(quant_expr[row_label,]), pc_scores = pc_scores, cell_line = factor(sample_info$cell_line))
    

    #PLOT!
    pca_scatter <- ggplot(df, aes(x = expression, y = pc_scores, colour = cell_line)) + geom_point() 
    pca_scatter <- pca_scatter + theme(text = element_text(size=18), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black")) 
    pca_scatter <- pca_scatter + labs(colour="Cell Line",x = paste0(gene_name, " expression"), y = paste0("PC",pc_num))
    ggsave(pca_scatter, file=pc_gene_scatter_output_file,width = 15,height=10.5,units="cm")

}

pc_gene_scatter_all_time_colored_by_time <- function(sample_info, quant_expr, ensamble_id, gene_name, pc_num, pc_gene_scatter_output_file) {
    #  Compute singular value decomposition
    svd1 <- svd(as.matrix(quant_expr))

    #  Scores of pcs of interest
    pc_scores <- svd1$v[,pc_num]

    # Row corresponding to gene of interest
    row_label <- which(rownames(quant_expr) == ensamble_id)
    # Get matrix into correct format
    quant_expr <- as.matrix(quant_expr)


    df <- data.frame(time = sample_info$time, expression = as.vector(quant_expr[row_label,]), pc_scores = pc_scores, cell_line = factor(sample_info$cell_line))
    

    #PLOT!
    pca_scatter <- ggplot(df, aes(x = expression, y = pc_scores, colour = time)) + geom_point() 
    pca_scatter <- pca_scatter + theme(text = element_text(size=18), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black")) 
    pca_scatter <- pca_scatter + labs(colour="Time",x = paste0(gene_name, " expression"), y = paste0("PC",pc_num))+ scale_color_gradient(low="pink",high="blue")
    ggsave(pca_scatter, file=pc_gene_scatter_output_file,width = 15,height=10.5,units="cm")

}


#  Plot first two PC's. Color points by cell_line
plot_pca_categorical_covariate <- function(sample_info, quant_expr, output_file, covariate, covariate_name,pc_num1,pc_num2) {

    #  Compute singular value decomposition
    svd1 <- svd(as.matrix(quant_expr))

    #  Scores of first 2 pc's across all samples
    pc1 <- svd1$v[,pc_num1]
    pc2 <- svd1$v[,pc_num2]

    # Put all information into data structure
    df <- data.frame(pc1 = pc1, pc2 = pc2, covariate = covariate)

    colourCount = length(unique(df$covariate))
    getPalette = colorRampPalette(brewer.pal(9, "Set1"))

    #PLOT!
    pca_scatter <- ggplot(df, aes(x = pc1, y = pc2, colour = covariate)) + geom_point() 
    pca_scatter <- pca_scatter + scale_colour_manual(values = getPalette(colourCount))
    pca_scatter <- pca_scatter + theme(text = element_text(size=12), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black")) 
    pca_scatter <- pca_scatter + labs(colour=covariate_name,x = paste0("PC",pc_num1), y = paste0("PC",pc_num2))
    pca_scatter <- pca_scatter + theme(legend.key.size =  unit(0.16, "in"))
    pca_scatter <- pca_scatter + theme(plot.title=element_text(size=8,face="plain"),text = element_text(size=8),axis.text=element_text(size=7), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.text = element_text(size=7), legend.title = element_text(size=8)) 

    ggsave(pca_scatter, file=output_file, width=7.2, height=4.5,units="in")

}

#  Plot first two PC's. Color points by cell_line
plot_ica_categorical_covariate <- function(sample_info, quant_expr, output_file, covariate, covariate_name,pc_num1,pc_num2) {

    imod <- icafast(as.matrix(quant_expr),nc=5)

    #  Scores of first 2 pc's across all samples
    pc1 <- imod$M[,pc_num1]
    pc2 <- imod$M[,pc_num2]

    # Put all information into data structure
    df <- data.frame(pc1 = pc1, pc2 = pc2, covariate = covariate)

    #PLOT!
    pca_scatter <- ggplot(df, aes(x = pc1, y = pc2, colour = covariate)) + geom_point() 
    pca_scatter <- pca_scatter + theme(text = element_text(size=18), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black")) 
    pca_scatter <- pca_scatter + labs(colour=covariate_name,x = paste0("IC",pc_num1), y = paste0("IC",pc_num2))
    ggsave(pca_scatter, file=output_file,width = 15,height=10.5,units="cm")

}


# Make plot showing variance explained of first n pcs
plot_pca_variance_explained <- function(sample_info, quant_expr, n, output_file, x_axis_label) {
    sv <- svd(as.matrix(quant_expr))


    variance_explained <- (sv$d^2/sum(sv$d^2))[1:n]

    # Merge global vectors into a data frame
    df <- data.frame(variance_explained = variance_explained, pc_num = 1:n)

    # PLOT AWAY
    line_plot <- ggplot(data=df, aes(x=pc_num, y=variance_explained)) +
                geom_line() +
                geom_point() +
                ylim(0,max(variance_explained) + .01) + 
                scale_x_continuous(breaks=1:n) +
                labs(x = x_axis_label, y = "Variance Explained") + 
                 theme(plot.title=element_text(size=8,face="plain"),text = element_text(size=8),axis.text=element_text(size=7), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.text = element_text(size=7), legend.title = element_text(size=8), axis.text.x = element_text(vjust=.5)) 

    # SAVE PLOT
    ggsave(line_plot, file=output_file,width = 19,height=13.5,units="cm")
    return(line_plot)
}





# Helper method to make pca plot
make_pca_plot <- function(i, cell_lines, pc1, pc2, sample_info) {
    # Get name of ith cell line
    i_cell_line <- cell_lines[i]

    # Get indices  of all samples that are from ith cell line
    i_indices <- sample_info$cell_line == i_cell_line

    #  Extract 1st two pcs of all samples that belong to the ith cell lein
    i_pc1 <- pc1[i_indices]
    i_pc2 <- pc2[i_indices]

    # Get time steps of these samples
    i_time <- sample_info$time[i_indices]

    # Put into compact data frame for plotting
    df <- data.frame(pc1 = i_pc1, pc2 = i_pc2, time_step = i_time)

    #PLOT!
    pca_scatter <-  ggplot(df,aes(pc1,pc2)) + geom_point(aes(colour=time_step)) + theme(text = element_text(size=12), axis.text = element_text(size=7), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black")) 
    pca_scatter <- pca_scatter + scale_color_gradient(low="pink",high="blue") + ggtitle(i_cell_line) + xlim(min(pc1)-.01,max(pc1)+.01) + ylim(min(pc2)-.01,max(pc2) + .01)
    pca_scatter <- pca_scatter + labs(colour = " ")

    return(pca_scatter)
}


#  Perform PCA. Make one plot for each cell line of first two pcs.
plot_pca_seperate_cell_lines <- function(sample_info, quant_expr, output_file,pc_num1,pc_num2) {

    #  Compute singular value decomposition
    svd1 <- svd(as.matrix(quant_expr))

    #  Scores of first 2 pc's across all samples
    pc1 <- svd1$v[,pc_num1]
    pc2 <- svd1$v[,pc_num2]


    #  Get unique cell lines
    cell_lines <- sort(unique(sample_info$cell_line))
    num_cell_lines <- length(cell_lines)


    
    #  Make pc plot for each cell line seperately (not automated yet..
     
    p1 <- make_pca_plot(1, cell_lines, pc1, pc2, sample_info) + labs(x=paste0("PC",pc_num1),y=paste0("PC",pc_num2)) # cell line 1
    p2 <- make_pca_plot(2, cell_lines, pc1, pc2, sample_info)+ labs(x=paste0("PC",pc_num1),y=paste0("PC",pc_num2)) + theme(legend.position="none") # cell line 2..
    p3 <- make_pca_plot(3, cell_lines, pc1, pc2, sample_info)+ labs(x=paste0("PC",pc_num1),y=paste0("PC",pc_num2)) +theme(legend.position="none")
    p4 <- make_pca_plot(4, cell_lines, pc1, pc2, sample_info)+ labs(x=paste0("PC",pc_num1),y=paste0("PC",pc_num2)) +theme(legend.position="none")
    p5 <- make_pca_plot(5, cell_lines, pc1, pc2, sample_info)+ labs(x=paste0("PC",pc_num1),y=paste0("PC",pc_num2)) +theme(legend.position="none")
    p6 <- make_pca_plot(6, cell_lines, pc1, pc2, sample_info)+ labs(x=paste0("PC",pc_num1),y=paste0("PC",pc_num2)) +theme(legend.position="none")
    p7 <- make_pca_plot(7, cell_lines, pc1, pc2, sample_info)+labs(x=paste0("PC",pc_num1),y=paste0("PC",pc_num2)) + theme(legend.position="none")
    p8 <- make_pca_plot(8, cell_lines, pc1, pc2, sample_info)+ labs(x=paste0("PC",pc_num1),y=paste0("PC",pc_num2)) +theme(legend.position="none")
    p9 <- make_pca_plot(9, cell_lines, pc1, pc2, sample_info)+ labs(x=paste0("PC",pc_num1),y=paste0("PC",pc_num2)) +theme(legend.position="none")
    p10 <- make_pca_plot(10, cell_lines, pc1, pc2, sample_info)+ labs(x=paste0("PC",pc_num1),y=paste0("PC",pc_num2)) +theme(legend.position="none")
    p11 <- make_pca_plot(11, cell_lines, pc1, pc2, sample_info)+ labs(x=paste0("PC",pc_num1),y=paste0("PC",pc_num2)) +theme(legend.position="none")
    p12 <- make_pca_plot(12, cell_lines, pc1, pc2, sample_info)+ labs(x=paste0("PC",pc_num1),y=paste0("PC",pc_num2)) +theme(legend.position="none")
    p13 <- make_pca_plot(13, cell_lines, pc1, pc2, sample_info)+ labs(x=paste0("PC",pc_num1),y=paste0("PC",pc_num2)) +theme(legend.position="none")
    p14 <- make_pca_plot(14, cell_lines, pc1, pc2, sample_info)+ labs(x=paste0("PC",pc_num1),y=paste0("PC",pc_num2)) +theme(legend.position="none")
    p15 <- make_pca_plot(15, cell_lines, pc1, pc2, sample_info)+ labs(x=paste0("PC",pc_num1),y=paste0("PC",pc_num2)) +theme(legend.position="none")
    p16 <- make_pca_plot(16, cell_lines, pc1, pc2, sample_info)+ labs(x=paste0("PC",pc_num1),y=paste0("PC",pc_num2)) +theme(legend.position="none")
    p17 <- make_pca_plot(17, cell_lines, pc1, pc2, sample_info)+ labs(x=paste0("PC",pc_num1),y=paste0("PC",pc_num2)) +theme(legend.position="none")
    p18 <- make_pca_plot(18, cell_lines, pc1, pc2, sample_info)+ labs(x=paste0("PC",pc_num1),y=paste0("PC",pc_num2)) +theme(legend.position="none")
    p19 <- make_pca_plot(19, cell_lines, pc1, pc2, sample_info)+ labs(x=paste0("PC",pc_num1),y=paste0("PC",pc_num2)) +theme(legend.position="none")


    legend <- get_legend(p1)

    # Merge all cell lines into one plot
    pdf(output_file)
    gg <- plot_grid(p1 + theme(legend.position="none"),p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12,p13,p14,p15,p16,p17,p18,p19,nrow=5,ncol=4)
    combined_gg <- ggdraw() + draw_plot(gg,0,0,1,1) + draw_plot(legend,.83,0,1,.3)
    print(combined_gg)
    dev.off()
}



# Fit the current model using training data
fit_model <- function(x_train, y_train, lambda, model_type) {
    if (model_type == "linear_regression") {
        #  Fit lasso linear regression with glmnet
        fit <- glmnet(x=x_train, y=y_train, alpha=1, lambda=lambda)
        #  Extract coefficient matrix of dim (num_genes + 1, 1)
        beta <- as.matrix(fit$beta)
        #  Include intercept..
        beta <- rbind(fit$a0,beta)
    } else if (model_type == "beta_binomial_regression") {
        #  Organize data for stan
        concShape <- 1.0001  # Diffuse prior hyperparameters for BB
        concRate <- 1e-4  # Diffuse prior hyperparameters for BB
        ns <- numeric(length(y_train)) + 15  # n in beta binomial is max of all time steps
        dat <- list(N=length(y_train), P=ncol(x_train), ys=y_train, x=x_train, concShape=concShape, concRate=concRate, lambda=lambda, ns=ns)

        ############################################
        # Propper initialization
        ###########################################
        rat <- dat$ys/dat$ns # ratios
        # moment estimator of concentration parameter
        conc <- pmin( 1/var(rat, na.rm=T), 1000 )
        m <- pmin( pmax( mean(rat, na.rm=T), 1/1000 ), 1-1/1000)
        betaInit <- numeric(ncol(x_train))
        betaInit[1] <- log(m/(1.0-m))
        init=list(conc=conc, beta=as.array(betaInit),alpha=0)


        #  Run stan second order optimization
        fit <- optimizing(BETABINOMIAL_GLM, data=dat,iter=10000, init=init, algorithm="LBFGS", hessian=T, as_vector=F, verbose=TRUE)

        #  Extract Parameters
        beta <- as.matrix(fit$par$beta)
        #  Include intercept..
        beta <- rbind(fit$par$alpha, beta)
    }
    return(beta)
}

#  Use fitted model to make predictions on testing data
model_prediction <- function(x_test, model_type, beta) {
    #  Add intercept to model
    column_of_ones <- (numeric(nrow(x_test)) +1)
    x_test<- cbind(column_of_ones, x_test)
    if (model_type == "linear_regression") {
        # Compute linear predictions
        y_predicted <- x_test %*% beta
    } else if (model_type == "beta_binomial_regression") {
        y_predicted <- (1/(1 + exp(-(x_test %*% beta))))*15
    }

    # Return in vector form
    return(y_predicted[,1])
}

# Determine optimal choice of lambda through k-fold cross validation
# x is design matrix
# y is response vector
# lambdas is 1d grid of possible lambdas
# k puts the k in k-fold
# model_type is type of glm
select_lambda_through_k_fold_cv <- function(x, y, lambdas, k, model_type) {
    #  Shuffle data
    shuffle_indices <- sample(nrow(x))
    x_shuffled <- x[shuffle_indices,]
    y_shuffled <- y[shuffle_indices]

    # Assign samples to folds
    folds <- cut(seq(1, nrow(x)), breaks=k, labels=FALSE)

    # Keep track of MSE of every lambda in every fold (matrix of dimension (number of lambdas X number of folds))
    mse_mat <- matrix(0, length(lambdas), k)

    #  Loop through each lambda
    for (i in 1:length(lambdas)) {
        #  Value of lambda:
        i_lambda <- lambdas[i]

        #  Loop through each fold
        for (fold_number in 1:k) {
            #Split data into training data and validation data (based on fold)
            x_train <- x_shuffled[folds != fold_number,]
            y_train <- y_shuffled[folds != fold_number]

            x_test <- x_shuffled[folds == fold_number,]
            y_test <- y_shuffled[folds == fold_number]

            #  Fit the current model using training data
            beta <- fit_model(x_train, y_train, i_lambda, model_type)
            #  Use fitted model to make predictions on testing data
            y_predicted <- model_prediction(x_test, model_type, beta)

            #  Summarize accuracy
            mse <- mean((y_test - y_predicted)^2)
            mse_mat[i, fold_number] <- mse
        }
    }

    #  Compute average(mse) across all folds (lambda_errors is of length(lambdas))
    lambda_errors <- rowMeans(mse_mat)

    #  Compute optimal lambda (lambda with smallest error)
    optimal_lambda <- lambdas[which.min(lambda_errors)]
    return(optimal_lambda)
}


#  Driver to train the glm.
#  x is the expr matrix of dim (samples X genes)
#  y is the corresponding vector of time step of dim (samples)
#  We want to first learn lambda via K-fold CV
#  Then, using the optimal lambda, we will optimize beta
train_glm <- function(x, y, model_type) {

    # Grid of optional lambdas
    lambdas <- c(0, .001, .01, .1, 1)
    
    # Number of folds in CV
    k <- 3
    
    # Search grid of lambdas with k-fold cross validation for optimal lambda
    optimal_lambda <- select_lambda_through_k_fold_cv(x, y, lambdas, k, model_type)

    #  Train model on ALL of training data using optimal lambda
    beta <- fit_model(x, y, optimal_lambda, model_type)
    return(beta)
}





# Compute variance of each row of a matrix (x)
row_var <- function(x) {
  rowSums((x - rowMeans(x))^2)/(dim(x)[2] - 1)
}

# Return the n top indices of a vector x
get_top_indices <- function(x, n) {
  nx <- length(x)
  p <- nx-n
  xp <- sort(x, partial=p)[p]
  which(x > xp)
}


# Extract indices of the top num_genes with the largest variance.
select_genes_with_largest_variance <- function(rpkm_expr, num_genes) {
    variances <- row_var(rpkm_expr)

    indices <- get_top_indices(variances, num_genes)
    return(indices)
}

#  Visualize our time predictions in the form of a heatmap
plot_heatmap_of_time_predictions <- function(mat, file_name, model_type) {
    #  Convert from matrix form to data frame format
    melted_mat <- melt(mat)
    colnames(melted_mat) <- c("cellLine", "Time","predictedTime")

    #  Use factors to represent cell line and time step
    melted_mat$cellLine <- factor(melted_mat$cellLine)
    melted_mat$trueTime <- factor(melted_mat$Time)

    #  PLOT!
    heatmap <- ggplot(data=melted_mat, aes(x=cellLine, y=Time)) + geom_tile(aes(fill=predictedTime)) + scale_fill_gradient2(midpoint=8, guide="colorbar")
    heatmap <- heatmap + theme(text = element_text(size=18), panel.background = element_blank(), axis.text.x = element_text(angle = 90, vjust=.5))
    heatmap <- heatmap + ggtitle(paste0("    ",model_type))

    ggsave(heatmap, file=file_name,width = 15,height=13.5,units="cm")

}

time_step_prediction <- function(sample_info, quant_expr, rpkm_expr, model_type, num_genes) {
    # Extract indices of the top num_genes with the largest variance.
    gene_indices <- select_genes_with_largest_variance(rpkm_expr, num_genes)

    #  Filter expression data to only include those top num_genes genes
    quant_expr_filt <- t(quant_expr[gene_indices,])
    rpkm_expr_filt <- t(rpkm_expr[gene_indices,])

    #  Extract unique cell lines
    cell_lines <- sort(unique(sample_info$cell_line))

    #  Extract unique time steps
    time_steps <- sort(unique(sample_info$time))

    # Keep track of time predictions for each (cell_line,time_step)
    time_prediction_mat <- matrix(NA, length(cell_lines), length(time_steps))
    rownames(time_prediction_mat) <- cell_lines
    colnames(time_prediction_mat) <- time_steps

    #  Loop through each cell line
    for (i in 1:length(cell_lines)) {
        #  ith cell line
        i_cell_line <- cell_lines[i]

        #  Extract indices of samples used for training
        training_indices <- sample_info$cell_line != i_cell_line
        #  Extract indices of samples used for testing
        testing_indices <- sample_info$cell_line == i_cell_line

        #  Extract training data (expr data and corresponding time step)
        training_expr <- quant_expr_filt[training_indices,]
        training_time <- sample_info$time[training_indices]
        #  Extract testing data (expr data and corresponding time step)
        testing_expr <- quant_expr_filt[testing_indices,]
        testing_time <- sample_info$time[testing_indices]

        #  Train glm. This include learning:
        #    1. Coefficient vector (beta) through mle (of dimension (num_genes +1,1))
        #    2. Lasso regression parameter (lambda) through cross validation
        beta <- train_glm(training_expr, training_time, model_type)

        #  Use glm to make predictions
        predicted_times <- model_prediction(testing_expr, model_type, beta)

        #  Add time predictions to matrix (to keep track of everything)
        for (t_index in 1:length(testing_time)) {
            actual_time <- testing_time[t_index]
            predicted_time <- predicted_times[t_index]

            time_prediction_mat[i, actual_time + 1] <- predicted_time
        }
    }
    return(time_prediction_mat)
}




#  Main driver for time step prediction
time_step_prediction_heatmap_plot <- function(sample_info, quant_expr, rpkm_expr, model_type, num_genes, output_file) {
    # Make predictions
    time_prediction_mat <- time_step_prediction(sample_info, quant_expr, rpkm_expr, model_type, num_genes)
    #  Visualize our time predictions in the form of a heatmap
    plot_heatmap_of_time_predictions(time_prediction_mat, output_file, model_type)
}


compare_variance_of_actual_vs_predicted_time <- function(sample_info, quant_expr, rpkm_expr, model_type, num_genes, variance_metric, output_file) {
    #  Extract unique cell lines
    cell_lines <- sort(unique(sample_info$cell_line))
    #  Extract unique time steps
    time_steps <- sort(unique(sample_info$time))
    # Make predictions
    time_prediction_mat <- time_step_prediction(sample_info, quant_expr, rpkm_expr, model_type, num_genes)
    
    # Update sample_info matrix to include time predictions
    sample_info_predicted <- sample_info
    #  Loop through all samples
    for (sample_num in 1:dim(sample_info)[1]) {
        # Find row index of time_prediction_mat corresponding to this sample
        t_cell_line_index <- which(sample_info$cell_line[sample_num] == cell_lines)
        # Find column index of time_prediction_mat corresponding to this sample
        t_time_index <- which(sample_info$time[sample_num] == time_steps)
        # Extract predicted time corresponding to this sample
        sample_info_predicted$time[sample_num] <- round(time_prediction_mat[t_cell_line_index, t_time_index])
        if (sample_info_predicted$time[sample_num] > 15) {
            sample_info_predicted$time[sample_num] <- 15
        }
        if (sample_info_predicted$time[sample_num] < 0) {
            sample_info_predicted$time[sample_num] <- 0
        }
    }

    rpkm_expr <- filter_expr_to_non_zero_genes(rpkm_expr, sample_info)
    rpkm_expr <- filter_expr_to_non_zero_genes(rpkm_expr,sample_info_predicted)

    # Initialize vector of distances 
    distances <- c()
    # Initialize vector of time steps
    t_steps <- c()
    # Initialize predicted_vs_actual
    time_type <- c()

    # Need to compute  sqare distance from mean for each time step seperately..
    # so Loop through each time step
    for (t_step in 0:15) {
        #  Get indices of all samples that are from t_step
        t_sample_indices <- sample_info$time == t_step

        #  Subset rpkm_expr matrix so it only contains samples from t_step
        t_rpkm_expr <- rpkm_expr[,t_sample_indices]

        #  compute square distance from the mean for all (gene, sample) pairs from t_step
        if (variance_metric == "square_distance_from_mean") {
            t_disty <- compute_square_distance_from_mean(t_rpkm_expr)
        }
        if (variance_metric == "sdev") {
            t_disty <- compute_sdev(t_rpkm_expr)
        }
        if (variance_metric == "log_square_distance_from_mean") {
            t_disty <- log(compute_square_distance_from_mean(t_rpkm_expr))
        }
        if (variance_metric == "log_sdev") {
            t_disty <- log(compute_sdev(t_rpkm_expr))
        }
        if (variance_metric == "avg_square_distance_from_mean") {
            t_disty <- compute_avg_square_distance_from_mean(t_rpkm_expr)
        }

        distances <- c(distances, t_disty)
        t_steps <- c(t_steps, numeric(length(t_disty)) + t_step)
        time_type <- c(time_type,rep("actual", length(t_disty)))
    }
    for (t_step in 0:15) {
        #  Get indices of all samples that are from t_step
        t_sample_indices <- sample_info_predicted$time == t_step

        #  Subset rpkm_expr matrix so it only contains samples from t_step
        t_rpkm_expr <- rpkm_expr[,t_sample_indices]

        #  compute square distance from the mean for all (gene, sample) pairs from t_step
        if (variance_metric == "square_distance_from_mean") {
            t_disty <- compute_square_distance_from_mean(t_rpkm_expr)
        }
        if (variance_metric == "sdev") {
            t_disty <- compute_sdev(t_rpkm_expr)
        }
        if (variance_metric == "log_square_distance_from_mean") {
            t_disty <- log(compute_square_distance_from_mean(t_rpkm_expr))
        }
        if (variance_metric == "log_sdev") {
            t_disty <- log(compute_sdev(t_rpkm_expr))
        }
        if (variance_metric == "avg_square_distance_from_mean") {
            t_disty <- compute_avg_square_distance_from_mean(t_rpkm_expr)
        }

        distances <- c(distances, t_disty)
        t_steps <- c(t_steps, numeric(length(t_disty)) + t_step)
        time_type <- c(time_type,rep("predicted", length(t_disty)))
    }

    #  Convert data into data frame format
    df <- data.frame(distance=distances, time=factor(t_steps),time_type=factor(time_type))
    
    testy <- wilcox.test(distances ~ time_type,data=df)

    print(variance_metric)
    print(mean(df[df$time_type == "predicted",]$distance))
    print(mean(df[df$time_type == "actual",]$distance))
    print(median(df[df$time_type == "predicted",]$distance))
    print(median(df[df$time_type == "actual",]$distance))

    # PLOT
    boxplot <- ggplot(df, aes(x=time,y=distance,fill=time_type)) + geom_boxplot(alpha=.7) + labs(x = "Time Step", y = variance_metric) + scale_x_discrete(name="Time Step") + scale_fill_brewer(palette = "Accent")
    boxplot <- boxplot + theme(text = element_text(size=18)) + ggtitle(paste0("Wilcoxon pvalue = ",testy$p.value))
    ggsave(boxplot, file=output_file,width = 20,height=10.5,units="cm")

}





#  compute average square distance from the mean for all samples from t_step
#  Returns vector of length number of (sample,gene) pairs belonging to t_step
compute_square_distance_from_mean <- function(t_rpkm_expr) {

    #  Number of genes/loci
    L <- dim(t_rpkm_expr)[1]

    #  Number of samples belonging to time step t
    N <- dim(t_rpkm_expr)[2]

    # Initialize vector corresponding to the average_square distance from the mean (each element is one (sample, gene) pair)
    disty <- numeric(N*L)

    # compute mean of each gene/locus
    loci_means <- rowMeans(t_rpkm_expr)

    # Compute distance for each sample,gene pair
    counter <- 1
    for (n in 1:N) {
        for (l in 1:L) {
            disty[counter] <- ((t_rpkm_expr[l,n] - loci_means[l])**2)/(loci_means[l]**2)
            counter <- counter + 1
        }
    }
    return(disty)
}

#  compute average square distance from the mean for all samples from t_step
#  Returns vector of length number of (sample,gene) pairs belonging to t_step
compute_avg_square_distance_from_mean <- function(t_rpkm_expr) {

    #  Number of genes/loci
    L <- dim(t_rpkm_expr)[1]

    #  Number of samples belonging to time step t
    N <- dim(t_rpkm_expr)[2]

    # Initialize vector corresponding to the average_square distance from the mean (each element is one (sample, gene) pair)
    disty <- numeric(N)

    # compute mean of each gene/locus
    loci_means <- rowMeans(t_rpkm_expr)

    # Compute distance for each sample,gene pair
    counter <- 1
    for (n in 1:N) {
        disty[counter] <- (N/(L*(N-1)))*(sum(((t_rpkm_expr[,n] - loci_means)**2)/(loci_means**2)))
        counter <- counter + 1
    }
    return(disty)
}

compute_sdev <- function(t_rpkm_expr) {

    #  Number of genes/loci
    L <- dim(t_rpkm_expr)[1]

    #  Number of samples belonging to time step t
    N <- dim(t_rpkm_expr)[2]

    # Initialize vector corresponding to the average_square distance from the mean (each element is one (sample, gene) pair)
    disty <- numeric(L)

    # Compute variance for each gene pair
    counter <- 1
    for (l in 1:L) {
        disty[counter] <- sd(t_rpkm_expr[l,])
        counter <- counter + 1
    }
    return(disty)
}


# Filter to only genes that do not have rpkm == 0 for all sample subsets
# A sample subset is all the samples at time step t
filter_expr_to_non_zero_genes <- function(rpkm_expr_all, sample_info) {
    #  Number of samples
    N <- dim(rpkm_expr_all)[2]
    #  Number of genes
    D <- dim(rpkm_expr_all)[1]

    #  Initialize TRUE/FALSE vector to all TRUE
    #  TRUE represents sample is non-zero expression in each expr matrix subset
    binary_vec <- rep(TRUE, D)

    #  Need to loop through each expr matrix subset (so loop through time steps)
    for (t_step in 0:15) {
        # Sample indices of tth time step
        t_sample_indices <- sample_info$time == t_step

        #  Subseted expression matrix for t time step
        t_rpkm_subset <- rpkm_expr_all[,t_sample_indices]

        # Determine which genes have mapped reads for this sample subset
        t_binary_vec <- !apply(t_rpkm_subset, 1, function(x){all(x==0)})

        binary_vec <- binary_vec*t_binary_vec
    }
    boolean_vec <- binary_vec == 1.0
    return(rpkm_expr_all[boolean_vec,])
}


# Compute the  square distance from mean across all samples (in each time step seperately). Then plot distribution for each time
square_distance_from_mean_driver <- function(sample_info, rpkm_expr_all, metric, square_distance_from_mean_output_file) {
    # Filter to only genes that do not have rpkm == 0 for all sample subsets
    # A sample subset is all the samples at time step t
    rpkm_expr <- filter_expr_to_non_zero_genes(rpkm_expr_all, sample_info)

    # Initialize vector of distances 
    distances <- c()
    # Initialize vector of time steps
    t_steps <- c()

    # Need to compute  sqare distance from mean for each time step seperately..
    # so Loop through each time step
    for (t_step in 0:15) {
        #  Get indices of all samples that are from t_step
        t_sample_indices <- sample_info$time == t_step

        #  Subset rpkm_expr matrix so it only contains samples from t_step
        t_rpkm_expr <- rpkm_expr[,t_sample_indices]

        #  compute square distance from the mean for all (gene, sample) pairs from t_step
        if (metric == "square_distance_from_mean") {
            t_disty <- compute_square_distance_from_mean(t_rpkm_expr)
        }
        if (metric == "sdev") {
            t_disty <- compute_sdev(t_rpkm_expr)
        }
        if (metric == "log_square_distance_from_mean") {
            t_disty <- log(compute_square_distance_from_mean(t_rpkm_expr))
        }
        if (metric == "log_sdev") {
            t_disty <- log(compute_sdev(t_rpkm_expr))
        }
        if (metric == "avg_square_distance_from_mean") {
            t_disty <- compute_avg_square_distance_from_mean(t_rpkm_expr)
        }

        distances <- c(distances, t_disty)
        t_steps <- c(t_steps, numeric(length(t_disty)) + t_step)
    }

    #  Convert data into data frame format
    df <- data.frame(distance=distances, time=factor(t_steps))

    # Plot!

    violin <- ggplot(df, aes(x=time,y=distance)) + geom_boxplot() + labs(x = "Time Step", y = metric)
    violin <- violin + theme(text = element_text(size=18))
    ggsave(violin, file=square_distance_from_mean_output_file,width = 15,height=10.5,units="cm")

}

time_step_specific_covariate_pc_pve_heatmap <- function(pc_file, covariate_file, output_file, time_step) {
    pcs <- read.table(pc_file, header=TRUE)
    all_covs <- read.table(covariate_file, header=TRUE)



    covs <- all_covs[all_covs$time==time_step,]
    if (toString(pcs$Sample_id) != toString(covs$sample_names)) {
        print("ASSUMPTION ERROR IN TIME STEP SPECIFIC COVARIATE PC PVE HEATMAP")
    }

    # Remove unimportant columns
    pcs <- as.matrix(pcs[,2:dim(pcs)[2]])
    covs <- data.frame(as.matrix(covs[,3:dim(covs)[2]]))
    # Get covariates into propper class (only necessary for numeric)
    covs$time <- as.numeric(as.character(covs$time))
    covs$library_size <- as.numeric(as.character(covs$library_size))
    covs$feeder_passage <- as.numeric(as.character(covs$feeder_passage))
    covs$feeder_free_passage <- as.numeric(as.character(covs$feeder_free_passage))
    covs$rna_extraction_conc <- as.numeric(as.character(covs$rna_extraction_conc))
    covs$RIN <- as.numeric(as.character(covs$RIN))
    covs$beating <- as.numeric(as.character(covs$beating))
    covs$line_beating <- as.numeric(as.character(covs$line_beating))
    covs$flash_freezing <- as.numeric(as.character(covs$flash_freezing))
    covs$fastq_percent_duplicates <- as.numeric(as.character(covs$fastq_percent_duplicates))
    covs$percent_gc <- as.numeric(as.character(covs$percent_gc))
    covs$average_sequence_length <- as.numeric(as.character(covs$average_sequence_length))
    covs$total_sequences <- as.numeric(as.character(covs$total_sequences))
    covs$percent_fails <- as.numeric(as.character(covs$percent_fails))
    covs$avg_10_15_troponin_mRNA_expr <- abs(as.numeric(as.character(covs$avg_10_15_troponin_mRNA_expr)))
    covs$avg_10_15_nanog_mRNA_expr <- abs(as.numeric(as.character(covs$avg_10_15_nanog_mRNA_expr)))

    keep_cols <- c()

    num_samples <- dim(covs)[1]
    for (col_num in 1:(dim(covs)[2])) {
        if (sd(covs[,col_num]) == 0 | nlevels(covs[,col_num]) == num_samples) {
            keep_cols <- c(keep_cols, FALSE)
        } else {
            keep_cols <- c(keep_cols, TRUE)
        }
    }
    covs <- covs[,keep_cols]

    # Initialize PVE heatmap
    pve_map <- matrix(0, dim(covs)[2], dim(pcs)[2])
    colnames(pve_map) <- colnames(pcs)
    rownames(pve_map) <- colnames(covs)

    # Loop through each PC, COV Pair and take correlation
    num_pcs <- dim(pcs)[2]
    num_covs <- dim(covs)[2]
    for (num_pc in 1:num_pcs) {
        for (num_cov in 1:num_covs) {
            pc_vec <- pcs[,num_pc]
            cov_vec <- covs[,num_cov]
            lin_model <- lm(pc_vec ~ cov_vec)
            pve_map[num_cov, num_pc] <- summary(lin_model)$adj.r.squared
            if (summary(lin_model)$adj.r.squared < 0.0) {
                pve_map[num_cov, num_pc] <- 0.0
            }
        }
    }
    

    ord <- hclust( dist(scale(pve_map), method = "euclidean"), method = "ward.D" )$order

    melted_mat <- melt(pve_map)
    colnames(melted_mat) <- c("Covariate", "PC","PVE")

    #  Use factors to represent covariate and pc name
    melted_mat$Covariate <- factor(chartr("_", " ",melted_mat$Covariate), levels = chartr("_", " ",rownames(pve_map)[ord]))
    melted_mat$PC <- factor(melted_mat$PC)
    if (dim(pcs)[2] == 10) {
        levels(melted_mat$PC) <- c(levels(melted_mat$PC)[1],levels(melted_mat$PC)[3:10],levels(melted_mat$PC)[2])
    }
    if (dim(pcs)[2] == 21) {
        levels(melted_mat$PC) <- c(levels(melted_mat$PC)[1],levels(melted_mat$PC)[12],levels(melted_mat$PC)[15:21],levels(melted_mat$PC)[2:11], levels(melted_mat$PC)[13:14])
    }

    #  PLOT!
    title <- paste0("time_step_", time_step)
    heatmap <- ggplot(data=melted_mat, aes(x=Covariate, y=PC)) + geom_tile(aes(fill=PVE)) + scale_fill_gradient2(midpoint=-.05, guide="colorbar")
    heatmap <- heatmap + theme(text = element_text(size=14), panel.background = element_blank(), axis.text.x = element_text(angle = 90, vjust=.5))
    heatmap <- heatmap + labs(y="latent factor", title=title)

    # Save File
    ggsave(heatmap, file=output_file,width = 19,height=13.5,units="cm")


}

covariate_pc_pve_heatmap <- function(pc_file, covariate_file, output_file, title) {
    # Load in data
    pcs <- read.table(pc_file, header=TRUE)
    covs <- read.table(covariate_file, header=TRUE)


    # Remove unimportant columns
    pcs <- as.matrix(pcs[,2:dim(pcs)[2]])
    covs <- data.frame(as.matrix(covs[,3:dim(covs)[2]]))



    # Get covariates into propper class (only necessary for numeric)
    covs$time <- as.numeric(as.character(covs$time))
    covs$library_size <- as.numeric(as.character(covs$library_size))
    covs$feeder_passage <- as.numeric(as.character(covs$feeder_passage))
    covs$feeder_free_passage <- as.numeric(as.character(covs$feeder_free_passage))
    covs$rna_extraction_conc <- as.numeric(as.character(covs$rna_extraction_conc))
    covs$RIN <- as.numeric(as.character(covs$RIN))
    covs$beating <- as.numeric(as.character(covs$beating))
    covs$line_beating <- as.numeric(as.character(covs$line_beating))
    covs$flash_freezing <- as.numeric(as.character(covs$flash_freezing))
    covs$percent_duplicates <- as.numeric(as.character(covs$percent_duplicates))
    covs$percent_gc <- as.numeric(as.character(covs$percent_gc))
    covs$average_sequence_length <- as.numeric(as.character(covs$average_sequence_length))
    covs$total_sequences <- as.numeric(as.character(covs$total_sequences))
    covs$percent_fails <- as.numeric(as.character(covs$percent_fails))
    covs$avg_10_15_troponin_mRNA_expr <- (as.numeric(as.character(covs$avg_10_15_troponin_mRNA_expr)))
    covs$avg_10_15_nanog_mRNA_expr <- (as.numeric(as.character(covs$avg_10_15_nanog_mRNA_expr)))

    # Initialize PVE heatmap
    pve_map <- matrix(0, dim(covs)[2], dim(pcs)[2])
    colnames(pve_map) <- colnames(pcs)
    rownames(pve_map) <- colnames(covs)

    # Loop through each PC, COV Pair and take correlation
    num_pcs <- dim(pcs)[2]
    num_covs <- dim(covs)[2]
    for (num_pc in 1:num_pcs) {
        for (num_cov in 1:num_covs) {
            pc_vec <- pcs[,num_pc]
            cov_vec <- covs[,num_cov]
            lin_model <- lm(pc_vec ~ cov_vec)
            pve_map[num_cov, num_pc] <- summary(lin_model)$adj.r.squared
        }
    }
    
    ord <- hclust( dist(scale(pve_map), method = "euclidean"), method = "ward.D" )$order

    melted_mat <- melt(pve_map)
    colnames(melted_mat) <- c("Covariate", "PC","PVE")

    #  Use factors to represent covariate and pc name
    melted_mat$Covariate <- factor(chartr("_", " ",melted_mat$Covariate), levels = chartr("_", " ",rownames(pve_map)[ord]))
    # melted_mat$Covariate <- factor(melted_mat$Covariate, levels = rownames(pve_map)[ord])
    melted_mat$PC <- factor(melted_mat$PC)
    if (dim(pcs)[2] == 10) {
        levels(melted_mat$PC) <- c(levels(melted_mat$PC)[1],levels(melted_mat$PC)[3:10],levels(melted_mat$PC)[2])
    }
    if (dim(pcs)[2] == 21) {
        levels(melted_mat$PC) <- c(levels(melted_mat$PC)[1],levels(melted_mat$PC)[12],levels(melted_mat$PC)[15:21],levels(melted_mat$PC)[2:11], levels(melted_mat$PC)[13:14])
    }

    #  PLOT!
    heatmap <- ggplot(data=melted_mat, aes(x=Covariate, y=substr(PC, 3,5))) + geom_tile(aes(fill=PVE)) + scale_fill_gradient2(midpoint=-.05, guide="colorbar")
    heatmap <- heatmap + labs(y="PC number",fill="VE")
    heatmap <- heatmap + theme(text = element_text(size=8),axis.text=element_text(size=7), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.text = element_text(size=7), legend.title = element_text(size=8),  axis.text.x = element_text(angle = 90,hjust=1, vjust=.5)) 
    # Save File
    ggsave(heatmap, file=output_file, width=7.2, height=5.0,units="in")
    return(heatmap)
}

covariate_pc_specific_genes_pve_heatmap <- function(pc_file, quant_expr, covariate_file, output_file) {
    # Load in data
    pcs <- read.table(pc_file, header=TRUE)
    covs <- read.table(covariate_file, header=TRUE)

    # Remove unimportant columns
    pcs <- as.matrix(pcs[,2:dim(pcs)[2]])
    covs <- data.frame(as.matrix(covs[,3:dim(covs)[2]]))
    # Get covariates into propper class (only necessary for numeric)
    covs$time <- as.numeric(as.character(covs$time))
    covs$library_size <- as.numeric(as.character(covs$library_size))
    covs$feeder_passage <- as.numeric(as.character(covs$feeder_passage))
    covs$feeder_free_passage <- as.numeric(as.character(covs$feeder_free_passage))
    covs$rna_extraction_conc <- as.numeric(as.character(covs$rna_extraction_conc))
    covs$RIN <- as.numeric(as.character(covs$RIN))
    covs$volume_needed <- as.numeric(as.character(covs$volume_needed))
    covs$water_added <- as.numeric(as.character(covs$water_added))
    covs$beating <- as.numeric(as.character(covs$beating))
    covs$line_beating <- as.numeric(as.character(covs$line_beating))
    covs$flash_freezing <- as.numeric(as.character(covs$flash_freezing))
    covs$fastq_percent_duplicates <- as.numeric(as.character(covs$fastq_percent_duplicates))
    covs$percent_gc <- as.numeric(as.character(covs$percent_gc))
    covs$average_sequence_length <- as.numeric(as.character(covs$average_sequence_length))
    covs$total_sequences <- as.numeric(as.character(covs$total_sequences))
    covs$percent_fails <- as.numeric(as.character(covs$percent_fails))
    covs$avg_10_15_troponin_mRNA_expr <- abs(as.numeric(as.character(covs$avg_10_15_troponin_mRNA_expr)))
    covs$avg_10_15_sox2_mRNA_expr <- abs(as.numeric(as.character(covs$avg_10_15_sox2_mRNA_expr)))

    # Add specific gene info
    row_label <- which(rownames(quant_expr) == "ENSG00000118194")
    quant_expr <- as.matrix(quant_expr)
    covs$Troponin <- as.vector(quant_expr[row_label,])

    row_label <- which(rownames(quant_expr) == "ENSG00000181449")
    covs$Sox2 <- as.vector(quant_expr[row_label,])

    # Initialize PVE heatmap
    pve_map <- matrix(0, dim(covs)[2], dim(pcs)[2])
    colnames(pve_map) <- colnames(pcs)
    rownames(pve_map) <- colnames(covs)

    # Loop through each PC, COV Pair and take correlation
    num_pcs <- dim(pcs)[2]
    num_covs <- dim(covs)[2]
    for (num_pc in 1:num_pcs) {
        for (num_cov in 1:num_covs) {
            pc_vec <- pcs[,num_pc]
            cov_vec <- covs[,num_cov]
            lin_model <- lm(pc_vec ~ cov_vec)
            pve_map[num_cov, num_pc] <- summary(lin_model)$adj.r.squared
        }
    }
    
    ord <- hclust( dist(scale(pve_map), method = "euclidean"), method = "ward.D" )$order

    melted_mat <- melt(pve_map)
    colnames(melted_mat) <- c("Covariate", "PC","PVE")

    #  Use factors to represent covariate and pc name
    melted_mat$Covariate <- factor(melted_mat$Covariate, levels = rownames(pve_map)[ord])
    melted_mat$PC <- factor(melted_mat$PC)
    levels(melted_mat$PC) <- c(levels(melted_mat$PC)[1],levels(melted_mat$PC)[3:10],levels(melted_mat$PC)[2])



    #  PLOT!
    heatmap <- ggplot(data=melted_mat, aes(x=Covariate, y=PC)) + geom_tile(aes(fill=PVE)) + scale_fill_gradient2(midpoint=-.05, guide="colorbar")
    heatmap <- heatmap + theme(text = element_text(size=14), panel.background = element_blank(), axis.text.x = element_text(angle = 90, vjust=.5))

    # Save File
    ggsave(heatmap, file=output_file,width = 19,height=13.5,units="cm")
}


# Perform PCA on time-step independent quantile normalized matrix. Plot variance explained of first n PCs:
plot_pca_time_independent_variance_explained <- function(sample_info, time_step_independent_quant_expr, n, pca_plot_time_step_output_file) {
    # Initialize global vectors to keep track of relevent quantities
    time_step_vec <- c()
    pc_num_vec <- c()
    variance_explained_vec <- c()

    #  Compute variance explained for each time step independently
    #  So loop through time steps
    for (time_step in 0:15) {
        # Compute indices of samples that belong to this time step
        time_step_indices <- sample_info$time == time_step
        # Filter expression matrix to only samples from this time step
        time_step_expr <- time_step_independent_quant_expr[,time_step_indices]

        # Run PCA (ie svc) on filtered expression matrix
        sv <- svd(as.matrix(time_step_expr))
        # Compute variance_explained
        variance_explained <- (sv$d^2/sum(sv$d^2))[1:n]
        
        # Add results to global vectors to keep track
        variance_explained_vec <- c(variance_explained_vec, variance_explained)
        time_step_vec <- c(time_step_vec, rep(time_step, length(variance_explained)))
        pc_num_vec <- c(pc_num_vec, 1:n)
    }

    # Merge global vectors into a data frame
    df <- data.frame(variance_explained = variance_explained_vec, time_step = time_step_vec, pc_num = pc_num_vec)

    # PLOT AWAY
    line_plot <- ggplot(data=df, aes(x=pc_num, y=variance_explained, group=time_step)) +
                geom_line(aes(colour=time_step)) +
                geom_point(aes(colour=time_step)) +
                scale_color_gradient(low="pink",high="blue") + 
                ylim(0,.35) +
                theme(text = element_text(size=14), panel.background = element_blank(), axis.text.x = element_text(vjust=.5)) + 
                scale_x_continuous(breaks=1:n) +
                labs(colour="Time Step", x = "PC Number", title = "Independent Time Step Variance Explained",y = "Variance Explained")

    # SAVE PLOT
    ggsave(line_plot, file=pca_plot_time_step_output_file,width = 19,height=13.5,units="cm")

}

pca_correlation_heatmap_time_independent_v_global_helper <- function(time_step, correlation_map,global_pca_loadings,sample_info,n_global,n_independent, time_step_independent_quant_expr) {
    # Initialize correlation heatmap
    correlation_map <- matrix(0, n_global, n_independent)
    rownames(correlation_map) <- paste0(1:n_global)
    colnames(correlation_map) <- paste0(1:n_independent)

    # Extract indices of samples corresponding to this time step
    time_step_indices <- sample_info$time == time_step
    # Extract expression matrix corresponding to this time step
    time_step_expr <- time_step_independent_quant_expr[,time_step_indices]
    # Run PCA (ie svc) on filtered expression matrix
    time_step_sv_independent <- svd(as.matrix(time_step_expr))
    # Extract loadings from SVA
    time_step_pca_loadings <- time_step_sv_independent$v[,1:n_independent]
    # Now loop through global pcs
    for (global_pc_num in 1:n_global) {
        # Loop through local pcs
        for (local_pc_num in 1:n_independent) {
            row_num <- global_pc_num
            col_num <- local_pc_num 
            local_loadings <- time_step_pca_loadings[,local_pc_num]
            global_loadings <- global_pca_loadings[time_step_indices,global_pc_num]
            correlation_map[row_num,col_num] <- abs(cor(local_loadings, global_loadings))
        }
    }

    melted_mat <- melt(correlation_map)
    colnames(melted_mat) <- c("global_pc", "local_pc","ABS_CORR")

    #  Use factors to represent covariate and pc name
    #melted_mat$GLOBAL_PC <- factor(melted_mat$GLOBAL_PC, levels = rownames(correlation_map)[ord])
    melted_mat$GLOBAL_PC <- factor(melted_mat$global_pc)
    melted_mat$LOCAL_PC <- factor(melted_mat$local_pc)
    #levels(melted_mat$PC) <- c(levels(melted_mat$PC)[1],levels(melted_mat$PC)[3:10],levels(melted_mat$PC)[2])



    #  PLOT!
    heatmap <- ggplot(data=melted_mat, aes(x=GLOBAL_PC, y=LOCAL_PC)) + geom_tile(aes(fill=ABS_CORR)) + scale_fill_gradient2(midpoint=-.05, guide="colorbar") + labs(title=paste0("time step: ", time_step), x = "Aggregated PCs", y = "Time step PCs", fill = "Abs(correlation)")

    heatmap <- heatmap + theme(text = element_text(size=10), panel.background = element_blank(), axis.text.x = element_text(angle = 0, vjust=.5)) + theme(legend.position = "top")
    return(heatmap)
}

pca_correlation_heatmap_time_independent_v_global <- function(sample_info, quant_expr, time_step_independent_quant_expr, n_global, n_independent, pca_correlation_heatmap_output_file) {
    # Run SVD (PCA) on global (across all time steps) expression matrix
    svd_global <- svd(as.matrix(quant_expr))
    #  Extract scores (loadings) from first n_global pcs
    global_pca_loadings <- svd_global$v[,1:n_global]


    time_step <- 0
    t0 <- pca_correlation_heatmap_time_independent_v_global_helper(time_step, correlation_map,global_pca_loadings,sample_info,n_global,n_independent,time_step_independent_quant_expr )

    #ggsave(t0, file=pca_correlation_heatmap_output_file,width = 19,height=13.5,units="cm")
    time_step <- 1
    t1 <- pca_correlation_heatmap_time_independent_v_global_helper(time_step, correlation_map,global_pca_loadings,sample_info,n_global,n_independent,time_step_independent_quant_expr )+ theme(legend.position="none")
    
    time_step <- 2
    t2 <- pca_correlation_heatmap_time_independent_v_global_helper(time_step, correlation_map,global_pca_loadings,sample_info,n_global,n_independent,time_step_independent_quant_expr )+ theme(legend.position="none")
    
    time_step <- 3
    t3 <- pca_correlation_heatmap_time_independent_v_global_helper(time_step, correlation_map,global_pca_loadings,sample_info,n_global,n_independent,time_step_independent_quant_expr )+ theme(legend.position="none")
    
    time_step <- 4
    t4 <- pca_correlation_heatmap_time_independent_v_global_helper(time_step, correlation_map,global_pca_loadings,sample_info,n_global,n_independent,time_step_independent_quant_expr )+ theme(legend.position="none")
    
    time_step <- 5
    t5 <- pca_correlation_heatmap_time_independent_v_global_helper(time_step, correlation_map,global_pca_loadings,sample_info,n_global,n_independent,time_step_independent_quant_expr )+ theme(legend.position="none")
    
    time_step <- 6
    t6 <- pca_correlation_heatmap_time_independent_v_global_helper(time_step, correlation_map,global_pca_loadings,sample_info,n_global,n_independent,time_step_independent_quant_expr )+ theme(legend.position="none")
    
    time_step <- 7
    t7 <- pca_correlation_heatmap_time_independent_v_global_helper(time_step, correlation_map,global_pca_loadings,sample_info,n_global,n_independent,time_step_independent_quant_expr )+ theme(legend.position="none")
    
    time_step <- 8
    t8 <- pca_correlation_heatmap_time_independent_v_global_helper(time_step, correlation_map,global_pca_loadings,sample_info,n_global,n_independent,time_step_independent_quant_expr )+ theme(legend.position="none")
    
    time_step <- 9
    t9 <- pca_correlation_heatmap_time_independent_v_global_helper(time_step, correlation_map,global_pca_loadings,sample_info,n_global,n_independent,time_step_independent_quant_expr )+ theme(legend.position="none")
    
    time_step <- 10
    t10 <- pca_correlation_heatmap_time_independent_v_global_helper(time_step, correlation_map,global_pca_loadings,sample_info,n_global,n_independent,time_step_independent_quant_expr )+ theme(legend.position="none")
    
    time_step <- 11
    t11 <- pca_correlation_heatmap_time_independent_v_global_helper(time_step, correlation_map,global_pca_loadings,sample_info,n_global,n_independent,time_step_independent_quant_expr )+ theme(legend.position="none")
    
    time_step <- 12
    t12 <- pca_correlation_heatmap_time_independent_v_global_helper(time_step, correlation_map,global_pca_loadings,sample_info,n_global,n_independent,time_step_independent_quant_expr )+ theme(legend.position="none")
    
    time_step <- 13
    t13 <- pca_correlation_heatmap_time_independent_v_global_helper(time_step, correlation_map,global_pca_loadings,sample_info,n_global,n_independent,time_step_independent_quant_expr )+ theme(legend.position="none")
    
    time_step <- 14
    t14 <- pca_correlation_heatmap_time_independent_v_global_helper(time_step, correlation_map,global_pca_loadings,sample_info,n_global,n_independent,time_step_independent_quant_expr )+ theme(legend.position="none")
    
    time_step <- 15
    t15 <- pca_correlation_heatmap_time_independent_v_global_helper(time_step, correlation_map,global_pca_loadings,sample_info,n_global,n_independent,time_step_independent_quant_expr )+ theme(legend.position="none")

    legend <- get_legend(t0)

    pdf(pca_correlation_heatmap_output_file)
    gg <- plot_grid(t0+ theme(legend.position="none"),t1,t2,t3,t4,t5,t6,t7,t8,t9,t10,t11,t12,t13,t14,t15,NULL,legend,nrow=5,ncol=4,label_size=10,rel_heights=c(1,1,1,1,.23))
    combined_gg <- ggdraw() + draw_plot(gg,0,0,1,1) 
    print(combined_gg)
    dev.off()

}


cell_line_pc_colored_by_state_model <- function(cell_pc_file, state_file, output_file) {

        # Load in data
        pcs <- read.table(cell_pc_file, header=TRUE)
        states <- read.table(state_file, header=FALSE,sep=",")

        cell_lines_states <- states[,1]
        cell_lines_pcs <- states[,1]
        if (sum(cell_lines_states != cell_lines_pcs) != 0) {
            print("FUNDAMENTAL ASSUMPTION ERROR")
        }


        pc1_vec <- as.numeric(pcs[,2])
        pc2_vec <- as.numeric(pcs[,3])
        state_vec <- states[,2]
        for (index in 1:length(state_vec)) {
            if (state_vec[index] == 0) {
                state_vec[index] = 2
            }
        }

        df <- data.frame(pc1 = pc1_vec, pc2 = pc2_vec, state = factor(paste0("Cell line cluster ",state_vec)), cell_line=factor(cell_lines_states))

        pca_scatter <- ggplot(df, aes(x = pc1, y=pc2, colour = state)) + geom_point() 
        pca_scatter <- pca_scatter + labs(colour="",x = "Cell line collapsed PC1", y = "Cell line collapsed PC2")
        pca_scatter <- pca_scatter + scale_colour_manual(values=c("dodgerblue3","chartreuse4"))
        pca_scatter <- pca_scatter + xlim(-.35, .45) + ylim(-.45,.35)
        pca_scatter <- pca_scatter + theme(text = element_text(size=8),axis.text=element_text(size=7), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.text = element_text(size=7), legend.title = element_text(size=8)) 
        ggsave(pca_scatter, file=output_file,width = 15,height=10.5,units="cm")
        return(pca_scatter)

}

cell_line_pc_colored_by_state_model_real_valued <- function(cell_pc_file, state_file, output_file) {

        # Load in data
        pcs <- read.table(cell_pc_file, header=TRUE)
        states <- read.table(state_file, header=FALSE,sep=",")

        cell_lines_states <- states[,1]
        cell_lines_pcs <- states[,1]
        if (sum(cell_lines_states != cell_lines_pcs) != 0) {
            print("FUNDAMENTAL ASSUMPTION ERROR")
        }


        pc1_vec <- as.numeric(pcs[,2])
        pc2_vec <- as.numeric(pcs[,3])
        state_vec <- states[,2]


        df <- data.frame(pc1 = pc1_vec, pc2 = pc2_vec, state = state_vec, cell_line=factor(cell_lines_states))

        valid_rows <- !is.na(df$state)


        print(cor.test(state_vec[valid_rows],pcs[,2][valid_rows]))
        print(cor.test(state_vec[valid_rows],pcs[,3][valid_rows]))
        print(cor.test(state_vec[valid_rows],pcs[,4][valid_rows]))
        print(cor.test(state_vec[valid_rows],pcs[,5][valid_rows]))

        df <- df[valid_rows,]

        pca_scatter <- ggplot(df, aes(x = pc1, y=pc2, colour = state)) + geom_point() 
        pca_scatter <- pca_scatter + theme(text = element_text(size=8),axis.text=element_text(size=7), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.text = element_text(size=7), legend.title = element_text(size=8)) 
        pca_scatter <- pca_scatter + labs(colour="% live cells expressing TNNT2",x = "Cell line collapsed PC1", y = "Cell line collapsed PC2")
        pca_scatter <- pca_scatter + scale_color_gradient(low="darkgoldenrod1",high="cyan4")
        pca_scatter <- pca_scatter + xlim(-.35, .45) + ylim(-.45,.35)
        ggsave(pca_scatter, file=output_file,width = 15,height=10.5,units="cm")
        return(pca_scatter)

}



cell_line_pc_colored_by_state_model_hmm <- function(cell_pc_file, state_file, output_file) {
        print(cell_pc_file)
        print(state_file)
        print(output_file)
        # Load in data
        pcs <- read.table(cell_pc_file, header=TRUE)
        states <- read.table(state_file, header=TRUE,sep=",")

        cell_lines_states <- states[,1]
        cell_lines_pcs <- states[,1]
        if (sum(cell_lines_states != cell_lines_pcs) != 0) {
            print("FUNDAMENTAL ASSUMPTION ERROR")
        }

        print(pcs)

        pc1_vec <- as.numeric(pcs[,2])
        pc2_vec <- as.numeric(pcs[,3])
        state_vec <- states[,2]

        df <- data.frame(pc1 = pc1_vec, pc2 = pc2_vec, state = factor(state_vec), cell_line=factor(cell_lines_states))

        pca_scatter <- ggplot(df, aes(x = pc1, y=pc2, colour = state)) + geom_point() 
        pca_scatter <- pca_scatter + theme(text = element_text(size=18), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black")) 
        pca_scatter <- pca_scatter + labs(colour="HMM grouping",x = "PC1", y = "PC2")
        ggsave(pca_scatter, file=output_file,width = 15,height=10.5,units="cm")

}


cell_line_pc_colored_by_avg_troponin <- function(cell_line_pc_file, covariates, output_file) {
    pcs <- read.table(cell_pc_file, header=TRUE)
    cell_lines <- pcs[,1]
    pc1_vec <- as.numeric(pcs[,2])
    pc2_vec <- as.numeric(pcs[,3])
    troponin_vec <- c()

    for (cell_line_iter in 1:length(cell_lines)) {
        cell_line <- cell_lines[cell_line_iter]
        indices <- (covariates$cell_line == cell_line)
        avg_expr <- (covariates$avg_10_15_troponin_mRNA_expr[indices])[1]
        troponin_vec <- c(troponin_vec, avg_expr)
    }

    df <- data.frame(pc1 = pc1_vec, pc2 = pc2_vec, avg_10_15_troponin=as.numeric(troponin_vec), cell_line=factor(cell_lines))

    #PLOT!
    pca_scatter <-  ggplot(df,aes(pc1,pc2)) + geom_point(aes(colour=avg_10_15_troponin)) + theme(text = element_text(size=18), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black")) 
    pca_scatter <- pca_scatter + scale_color_gradient(low="pink",high="blue")


    ggsave(pca_scatter, file=output_file,width = 15,height=10.5,units="cm",dpi=600)

}


add_old_vs_new <- function(sample_info) {
    num_entries <- length(sample_info$cell_line)
    string_id <- rep("old", num_entries)
    for (iter in 1:num_entries) {
        curr <- (sample_info$cell_line)[iter]
        if (curr == "19190" | curr == "18907" | curr == "18517" | curr == "19127" | curr == "19193") {
            string_id[iter] = "new"
        }
    }
    sample_info$old_vs_new <- string_id
    print(summary(sample_info))
    return(sample_info)

}

scatter_of_avg_troponin_and_flow <- function(flow_file, covariates, output_file) {
    flow_results <- read.table(flow_file, header=FALSE,sep=",")

    cell_lines <- flow_results[,1]
    flow_vec <- flow_results[,2]

    troponin_vec <- c()

    for (cell_line_iter in 1:length(cell_lines)) {
        cell_line <- cell_lines[cell_line_iter]
        indices <- (covariates$cell_line == cell_line)
        avg_expr <- (covariates$avg_10_15_troponin_mRNA_expr[indices])[1]
        troponin_vec <- c(troponin_vec, avg_expr)
    }


    df <- data.frame(flow_results=flow_vec,avg_troponin=as.numeric(troponin_vec))

    valid_rows <- !is.na(df$flow_results)

    df <- df[valid_rows,]

    corry <- cor(df$flow_results, df$avg_troponin)

    pca_scatter <- ggplot(df, aes(x = flow_results, y=avg_troponin)) + geom_point() 
    pca_scatter <- pca_scatter + theme(text = element_text(size=18), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black")) 
    pca_scatter <- pca_scatter + labs(x="Flow results",y = "avg troponin (t=10-15)", title=paste0("pearson correlation=",corry))
    ggsave(pca_scatter, file=output_file,width = 15,height=10.5,units="cm")

}

scatter_of_pc_and_flow <- function(flow_file, pc_file, pc_num, output_file) {
    flow_results <- read.table(flow_file, header=FALSE,sep=",")

    cell_lines <- flow_results[,1]
    flow_vec <- flow_results[,2]

    pcs <- read.table(pc_file, header=TRUE)
    cell_lines <- pcs[,1]
    pc_vec <- as.numeric(pcs[,(pc_num+1)])



    df <- data.frame(flow_results=flow_vec,pc_vec=pc_vec)

    valid_rows <- !is.na(df$flow_results)

    df <- df[valid_rows,]

    corry <- cor(df$pc_vec, df$flow_results)

    pca_scatter <- ggplot(df, aes(x = flow_results, y=pc_vec)) + geom_point() 
    pca_scatter <- pca_scatter + theme(text = element_text(size=18), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black")) 
    pca_scatter <- pca_scatter + labs(x="Flow results",y = paste0("pc",pc_num), title=paste0("pearson correlation=",corry))
    ggsave(pca_scatter, file=output_file,width = 15,height=10.5,units="cm")

}

gene_cluster_plotter = function(melted,gene_cluster, ci=FALSE){
    data = melted[melted$L == gene_cluster, ]
    data0 = data[data$CellLineCluster == 'Cell Line Cluster 1',]
    data1 = data[data$CellLineCluster == 'Cell Line Cluster 2',]


    p <- ggplot(data, aes(x=X, y=value, color=CellLineCluster)) + geom_line() + theme(text = element_text(size=12), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"))
    if (ci){
        p <- p + geom_ribbon(data=data, aes(ymin=lower,ymax=upper, fill=CellLineCluster),alpha=0.3, colour=NA)
    }
    p <- p + labs(x="Day",y = "Expression",fill="",color="") + theme(legend.position = "bottom") + scale_fill_manual(values=c("slateblue3","chartreuse4"))+ scale_color_manual(values=c("slateblue3","chartreuse4"))
    p <- p + theme(legend.text = element_text(size=10))+ theme(legend.title = element_text(size=10)) + ylim(-1.5,2.5)
    return(p)
}

plot_karl_gene_cluster <- function(mixutre_hmm_cell_line_grouping_dir, gene_cluster_name, CI_boolean) {
    fmean <- read.csv(paste0(mixutre_hmm_cell_line_grouping_dir, "mixsvgp_K2_L20_0_29526888_fmean"))
    fvar <- read.csv(paste0(mixutre_hmm_cell_line_grouping_dir, "mixsvgp_K2_L20_0_29526888_fvar"))
    fmean$X <- fmean$X * 15/100

    L  <- paste('Gene Cluster ', rep(c(paste('0', seq(0, 9), sep=""), seq(10, 19)), 2), sep="")
    K <- paste('Cell Line Cluster ', c(rep(2, 20), rep(1, 20)), sep="")

    KL <- paste(K, L, sep="")

    colnames(fmean) <- c('X', KL)
    melted <- melt(fmean, id.vars = 'X')
    melted$CellLineCluster <- substr(melted$variable, 1, 19)
    melted$L <- substr(melted$variable, 20, 100)
    melted$sd <-sqrt(melt(fvar, id.vars = 'X')$value)
    melted$lower <- melted$value - 1.96 * melted$sd
    melted$upper <- melted$value + 1.96 * melted$sd


    return(gene_cluster_plotter(melted,gene_cluster_name,CI_boolean))

}



plot_karl_gene_cluster_grid <- function(mixutre_hmm_cell_line_grouping_dir,ci=FALSE) {
    fmean <- read.csv(paste0(mixutre_hmm_cell_line_grouping_dir, "mixsvgp_K2_L20_0_29526888_fmean"))
    fvar <- read.csv(paste0(mixutre_hmm_cell_line_grouping_dir, "mixsvgp_K2_L20_0_29526888_fvar"))
    fmean$X <- fmean$X * 15/100



    L  <- paste('Gene Cluster ', rep(c(paste('0', seq(0, 9), sep=""), seq(10, 19)), 2), sep="")
    K <- paste('cell line cluster ', c(rep(2, 20), rep(1, 20)), sep="")

    KL <- paste(K, L, sep="")

    colnames(fmean) <- c('X', KL)
    melted <- melt(fmean, id.vars = 'X')
    melted$CellLineCluster <- paste0("Inferred ",substr(melted$variable, 1, 19))
    melted$L <- substr(melted$variable, 20, 100)
    melted$sd <-sqrt(melt(fvar, id.vars = 'X')$value)
    melted$lower <- melted$value - 1.96 * melted$sd
    melted$upper <- melted$value + 1.96 * melted$sd


    p <- ggplot(melted, aes(x=X, y=value, color=CellLineCluster)) + geom_line() + facet_wrap(~ L, ncol=10) 
    if (ci){
        p <- p + geom_ribbon(data=data, aes(ymin=lower,ymax=upper, fill=CellLineCluster),alpha=0.3, colour=NA)
    }
    p <- p + labs(x="Day",y = "Expression",fill="",color="") + theme(legend.position = "bottom") + scale_fill_manual(values=c("dodgerblue3","chartreuse4"))+ scale_color_manual(values=c("dodgerblue3","chartreuse4"))
    p <- p + theme(text = element_text(size=8),axis.text=element_text(size=7), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.text = element_text(size=7), legend.title = element_text(size=8)) 
    return(p)


}


plot_karl_gene_cluster_grid_v2 <- function(mixutre_hmm_cell_line_grouping_dir, row_name,ci=FALSE) {
    fmean <- read.csv(paste0(mixutre_hmm_cell_line_grouping_dir, "mixsvgp_K2_L20_0_29526888_fmean"))
    fvar <- read.csv(paste0(mixutre_hmm_cell_line_grouping_dir, "mixsvgp_K2_L20_0_29526888_fvar"))
    fmean$X <- fmean$X * 15/100

    L  <- paste('Gene Cluster ', rep(c(paste('0', seq(0, 9), sep=""), seq(10, 19)), 2), sep="")
    if (row_name == "row1") {
        L  <- rep(paste0("Gene Cluster 0" , seq(0,9)), 2)
    } 
    if (row_name == "row2") {
        L  <- rep(paste0("Gene Cluster " , seq(10,19)), 2)
    }
    K <- paste('Cell Line Cluster ', c(rep(2, 20), rep(1, 20)), sep="")

    KL <- paste(K, L, sep="")

    colnames(fmean) <- c('X', KL)
    melted <- melt(fmean, id.vars = 'X')
    melted$CellLineCluster <- substr(melted$variable, 1, 19)
    melted$L <- substr(melted$variable, 20, 100)
    melted$sd <-sqrt(melt(fvar, id.vars = 'X')$value)
    melted$lower <- melted$value - 1.96 * melted$sd
    melted$upper <- melted$value + 1.96 * melted$sd


    p <- ggplot(melted, aes(x=X, y=value, color=CellLineCluster)) + geom_line() + facet_wrap(~ L, ncol=10) + theme(text = element_text(size=12),axis.text=element_text(size=7), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"))
    if (ci){
        p <- p + geom_ribbon(data=data, aes(ymin=lower,ymax=upper, fill=CellLineCluster),alpha=0.3, colour=NA)
    }
    p <- p + labs(x="Day",y = "Expression",fill="",color="") + theme(legend.position = "bottom") + scale_fill_manual(values=c("dodgerblue3","chartreuse4"))+ scale_color_manual(values=c("dodgerblue3","chartreuse4"))
    p <- p + theme(legend.text = element_text(size=10))+ theme(legend.title = element_text(size=10)) + ylim(-1.5,2.6)
    return(p)


}

make_edf_troponin_nanog_time_course <- function(sample_info, quant_expr, output_file) {

    ensamble_id <- "ENSG00000111704"
    gene_name <- "Nanog"
    nanog_line_plot <- gene_time_course_line_plot_grouped_by_cell_line_modular(sample_info, quant_expr, ensamble_id, gene_name)



    ensamble_id <- "ENSG00000118194"
    gene_name <- "Troponin T2"
    troponin_line_plot <- gene_time_course_line_plot_grouped_by_cell_line_modular(sample_info, quant_expr, ensamble_id, gene_name)

    legend = get_legend(troponin_line_plot + theme(legend.key.size =  unit(0.2, "in")) + theme(legend.text = element_text(size=8))+ theme(legend.title = element_text(size=8)))

    combined <- ggdraw() + 
                draw_plot(nanog_line_plot + theme(legend.position='none'), 0,.5,.8,.5) + 
                draw_plot(troponin_line_plot + theme(legend.position='none'),0,0,.8,.5) + 
                draw_plot(legend,.8,0,1,1) + 
                draw_plot_label(c("A","B"),c(0,0),c(1,.5),size=12)

    ggsave(combined, file=output_file, width=7.2, height=5.5,units="in")


}

make_gene_cluster_color_bar <- function(gene_clusters) {
    ordered_cluster_ids <- unique(gene_clusters)
    fraction_arr <- c()
    cluster_arr <- c()
    pos_arr <- c()
    for (counter in 1:length(ordered_cluster_ids)) {
        cluster_id <- ordered_cluster_ids[counter]
        cluster_count <- sum(gene_clusters==cluster_id)
        fraction_arr <- c(fraction_arr, cluster_count)
        cluster_arr <- c(cluster_arr, cluster_id)
        pos_arr <- c(pos_arr, 1)

    }

    colourCount = length(ordered_cluster_ids)
    getPalette = colorRampPalette(brewer.pal(9, "Set1"))

    df <- data.frame(position=pos_arr,cluster=factor(cluster_arr,levels=ordered_cluster_ids), values=fraction_arr)
    bar <- ggplot(df, aes(x=position, y=values, fill=cluster)) + geom_bar(stat="identity") + theme_void() #+ scale_fill_manual(values=getPalette(colourCount))
    bar <- bar + scale_fill_manual(values=c("goldenrod3","chartreuse4"))
    return(bar)
}

make_cell_line_cluster_color_bar_cross_all_genes <- function(cell_line_names, gene_cluster_names) {
    ordered_cluster_ids <- c()
    fraction_arr <- c()
    cluster_arr <- c()
    pos_arr <- c()
    for (gene_cluster in 1:20) {
        gene_cluster_id <- as.character(gene_cluster)
        cell_line_cluster_id <- as.character("1")
        joint_cluster_id <- paste0(gene_cluster_id,"_",cell_line_cluster_id)
        count <- sum(gene_cluster_names == gene_cluster_id & cell_line_names == cell_line_cluster_id)
        fraction_arr <- c(fraction_arr,count)
        cluster_arr <- c(cluster_arr, joint_cluster_id)
        pos_arr <- c(pos_arr,1)

        cell_line_cluster_id <- as.character("2")
        joint_cluster_id <- paste0(gene_cluster_id,"_",cell_line_cluster_id)
        count <- sum(gene_cluster_names == gene_cluster_id & cell_line_names == cell_line_cluster_id)
        fraction_arr <- c(fraction_arr,count)
        cluster_arr <- c(cluster_arr, joint_cluster_id)
        pos_arr <- c(pos_arr,1)
    }

    df <- data.frame(position=pos_arr,cluster=factor(cluster_arr,levels=cluster_arr), values=fraction_arr)

    bar <- ggplot(df, aes(x=position, y=values, fill=cluster)) + geom_bar(stat="identity") + theme_void() 
    bar <- bar + scale_fill_manual(values=rep(c("goldenrod3","chartreuse4"),20))
    return(bar)
    
}

heatmap_showing_average_expression_for_cell_line_cluster <- function(cell_line_cluster_average_expression_file, cell_line_cluster) {
    avg_expr <- read.table(cell_line_cluster_average_expression_file, header=TRUE)
    gene_names <- avg_expr[,1]
    gene_clusters <- avg_expr[,2]
    expr <- as.matrix(avg_expr[,3:(dim(avg_expr)[2])])
    colnames(expr) <- as.character(0:15)

    melted_mat <- melt(expr)
    colnames(melted_mat) <- c("Gene", "Day","Expression")

    gene_cluster_color_bar <- make_gene_cluster_color_bar(gene_clusters)


    melted_mat$Gene <- factor(gene_names, levels=rev(gene_names))
    melted_mat$Day <- factor(melted_mat$Day, levels=as.character(0:15))

    heatmap <- ggplot(data=melted_mat, aes(x=Day, y=Gene)) + geom_tile(aes(fill=Expression)) + scale_fill_gradient2(midpoint=0, guide="colorbar")
    heatmap <- heatmap + theme(text = element_text(size=12), panel.background = element_blank())
    heatmap <- heatmap + labs(x="Day",y="Gene",fill="Expression") + theme(legend.text = element_text(size=10))+ theme(legend.title = element_text(size=10))
    # Remove marks on y-axis
    heatmap <- heatmap +  theme(axis.text.y=element_blank(),axis.ticks.y=element_blank())
    heatmap <- heatmap + scale_x_discrete(breaks=c("0","1","2","3","4","5","6","7","8","9","10","11","12","13","14","15"),labels=c("0","","","","","5","","","","","10","","","","","15"))
        
    combined <- ggdraw() + 
        draw_plot(heatmap + theme(legend.position='none'),.02,0,.98,1) + 
        draw_plot(gene_cluster_color_bar + theme(legend.position='none'),.095,.18,.03,.82) 


    return(combined)


}

split_gpm_concordance_heatmap <- function(input_file, l) {
    incidence_mat <- read.table(input_file, header=TRUE, sep=",")
    cell_lines <- paste0("NA",as.character(incidence_mat[,1]))
    incidence <- as.matrix(incidence_mat[,2:(dim(incidence_mat)[2])])
    colnames(incidence) <- cell_lines
    rownames(incidence) <- cell_lines

    melted_incidence <- melt(incidence)
    colnames(melted_incidence) <- c("Cell_line_1","Cell_line_2","overlap")

    melted_incidence$Cell_line_1 <- factor(melted_incidence$Cell_line_1, levels=paste0("NA",c("18489", "18907", "19127", "18508","19108", "19093","18505","19159","18855","18912","18520","18511","19190","18858","18517","18870","19193","18499","19209")))
    melted_incidence$Cell_line_2 <- factor(melted_incidence$Cell_line_2, levels=rev(paste0("NA",c("18489", "18907", "19127", "18508","19108", "19093","18505","19159","18855","18912","18520","18511","19190","18858","18517","18870","19193","18499","19209"))))
    
    heatmap <- ggplot(data=melted_incidence, aes(x=Cell_line_1, y=Cell_line_2)) + geom_tile(aes(fill=overlap)) #+ scale_fill_gradient(low="grey",high="plum2")
    #heatmap <- heatmap + scale_fill_distiller()
    #heatmap <- heatmap + scale_fill_brewer(values = brewer.pal(3,"RdPu"))
    heatmap <- heatmap + scale_fill_distiller(palette = "Blues", direction=1) + theme(plot.title = element_text(face="plain",size=8))
    heatmap <- heatmap + theme(text = element_text(size=8),axis.text=element_text(size=7), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.text = element_text(size=7), legend.title = element_text(size=8), axis.text.x = element_text(angle = 90,hjust=1, vjust=.5)) 
    heatmap <- heatmap + labs(x = "Cell line", y = "Cell line", fill= "Overlap", title=paste0(l, " gene clusters"))
    return(heatmap)

}


make_heatmap_for_gene_cluster <- function(gene_cluster, expr_cluster1, expr_cluster2, gene_names_cluster1, gene_names_cluster2,gene_clusters_cluster1, gene_clusters_cluster2) {
    gene_cluster_indices <- gene_clusters_cluster1 == gene_cluster
    gene_cluster_expr_cluster1 <- expr_cluster1[gene_cluster_indices,]
    gene_cluster_expr_cluster2 <- expr_cluster2[gene_cluster_indices,]
    expr <- rbind(gene_cluster_expr_cluster1, gene_cluster_expr_cluster2)
    gene_names <- c(paste0(gene_names_cluster1[gene_cluster_indices], "_cluster1"),paste0(gene_names_cluster2[gene_cluster_indices], "_cluster2"))
    cell_line_names <- c(rep("1",sum(gene_cluster_indices)), rep("2",sum(gene_cluster_indices)))



    cell_line_color_bar <- make_gene_cluster_color_bar(cell_line_names)

    colnames(expr) <- as.character(0:15)
    melted_mat <- melt(expr)
    colnames(melted_mat) <- c("Gene", "Day","Expression")

    melted_mat$Gene <- factor(gene_names, levels=gene_names)
    melted_mat$Day <- factor(melted_mat$Day, levels=as.character(0:15))

    heatmap <- ggplot(data=melted_mat, aes(x=Day, y=Gene)) + geom_tile(aes(fill=Expression)) + scale_fill_gradient2(midpoint=0, guide="colorbar")
    heatmap <- heatmap + theme(text = element_text(size=12), panel.background = element_blank())
    heatmap <- heatmap + labs(x="",y="",fill="") + theme(legend.text = element_text(size=10))+ theme(legend.title = element_text(size=10))
    # Remove marks on y-axis
    heatmap <- heatmap +  theme(axis.text.y=element_blank(),axis.ticks.y=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank()) + theme_void()
    heatmap <- heatmap + scale_x_discrete(breaks=c("0","1","2","3","4","5","6","7","8","9","10","11","12","13","14","15"),labels=c("0","","","","","5","","","","","10","","","","","15"))
    
  
    heatmap <- heatmap #+ theme(panel.border = element_rect(colour = "black", fill=NA, size=5))
    combined <- ggdraw() + 
        draw_plot(heatmap + theme(legend.position='none') ,.02,.045,.98,.91) + 
        draw_plot(cell_line_color_bar + theme(legend.position='none'),0,0,.03,1) 
    return(combined)
}

heatmap_showing_average_expression_for_gene_cluster <- function(cell_line_cluster_1_average_expression_file, cell_line_cluster_2_average_expression_file) {
    avg_expr_cluster1 <- read.table(cell_line_cluster_1_average_expression_file, header=TRUE)
    gene_names_cluster1 <- avg_expr_cluster1[,1]
    gene_clusters_cluster1 <- avg_expr_cluster1[,2]
    expr_cluster1 <- as.matrix(avg_expr_cluster1[,3:(dim(avg_expr_cluster1)[2])])
    

    avg_expr_cluster2 <- read.table(cell_line_cluster_2_average_expression_file, header=TRUE)
    gene_names_cluster2 <- avg_expr_cluster2[,1]
    gene_clusters_cluster2 <- avg_expr_cluster2[,2]
    expr_cluster2 <- as.matrix(avg_expr_cluster2[,3:(dim(avg_expr_cluster2)[2])])


    li = list() 
    for (i in 1:20) {
        li[[i]] = make_heatmap_for_gene_cluster(as.character(i), expr_cluster1, expr_cluster2, gene_names_cluster1, gene_names_cluster2,gene_clusters_cluster1, gene_clusters_cluster2)
    }



    # Should return a vector of heatmaps. 19 and 20 should have x-axis labels

    return(li)


}

heatmap_showing_average_expression_for_all_clusters <- function(cell_line_cluster_1_average_expression_file, cell_line_cluster_average_expression_file) {
    avg_expr_cluster1 <- read.table(cell_line_cluster_1_average_expression_file, header=TRUE)
    gene_names_cluster1 <- avg_expr_cluster1[,1]
    gene_clusters_cluster1 <- avg_expr_cluster1[,2]
    expr_cluster1 <- as.matrix(avg_expr_cluster1[,3:(dim(avg_expr_cluster1)[2])])
    

    avg_expr_cluster2 <- read.table(cell_line_cluster_2_average_expression_file, header=TRUE)
    gene_names_cluster2 <- avg_expr_cluster2[,1]
    gene_clusters_cluster2 <- avg_expr_cluster2[,2]
    expr_cluster2 <- as.matrix(avg_expr_cluster2[,3:(dim(avg_expr_cluster2)[2])])

    fraction_from_gene_cluster <- c()

    for (gene_cluster in 1:20) {
        gene_cluster_indices <- gene_clusters_cluster1 == gene_cluster
        gene_cluster_expr_cluster1 <- expr_cluster1[gene_cluster_indices,]
        gene_cluster_expr_cluster2 <- expr_cluster2[gene_cluster_indices,]
        fraction_from_gene_cluster <- c(fraction_from_gene_cluster, sum(gene_clusters_cluster1 == gene_cluster)/length(gene_clusters_cluster1))
        if (gene_cluster == 1) {
            expr <- rbind(gene_cluster_expr_cluster1, gene_cluster_expr_cluster2)
            gene_names <- c(paste0(gene_names_cluster1[gene_cluster_indices], "_cluster1"),paste0(gene_names_cluster2[gene_cluster_indices], "_cluster2"))
            cell_line_names <- c(rep("1",sum(gene_cluster_indices)), rep("2",sum(gene_cluster_indices)))
            gene_cluster_names <- c(rep(as.character(gene_cluster), 2*sum(gene_cluster_indices)))
        }
        else {
            expr <- rbind(expr, gene_cluster_expr_cluster1, gene_cluster_expr_cluster2)
            gene_names <- c(gene_names, paste0(gene_names_cluster1[gene_cluster_indices], "_cluster1"),paste0(gene_names_cluster2[gene_cluster_indices], "_cluster2"))
            cell_line_names <- c(cell_line_names, rep("1",sum(gene_cluster_indices)), rep("2",sum(gene_cluster_indices)))
            gene_cluster_names <- c(gene_cluster_names, rep(as.character(gene_cluster), 2*sum(gene_cluster_indices)))
        }
    }

    color_bar <- make_cell_line_cluster_color_bar_cross_all_genes(cell_line_names, gene_cluster_names)


    colnames(expr) <- as.character(0:15)
    melted_mat <- melt(expr)
    colnames(melted_mat) <- c("Gene", "Day","Expression")

    melted_mat$Gene <- factor(gene_names, levels=rev(gene_names))
    melted_mat$Day <- factor(melted_mat$Day, levels=as.character(0:15))

    heatmap <- ggplot(data=melted_mat, aes(x=Day, y=Gene)) + geom_tile(aes(fill=Expression))+ scale_fill_gradient2(midpoint=0, guide="colorbar")#+ scale_fill_distiller(palette = "PRGn") 

    heatmap <- heatmap + theme(text = element_text(size=12), panel.background = element_blank())
    heatmap <- heatmap + labs(x="Day",y="",fill="") + theme(legend.text = element_text(size=10))+ theme(legend.title = element_text(size=10))
    # Remove marks on y-axis
    heatmap <- heatmap +  theme(axis.text.y=element_blank(),axis.ticks.y=element_blank())
    heatmap <- heatmap + scale_x_discrete(breaks=c("0","1","2","3","4","5","6","7","8","9","10","11","12","13","14","15"),labels=c("0","","","","","5","","","","","10","","","","","15"))
    

    heatmap_legend <- get_legend(heatmap)
    heatmap <- heatmap + theme(legend.position="none")

    combined <- ggdraw() + 
        draw_plot(heatmap + theme(legend.position='none') ,0,0,1,.99) + 
        draw_plot(color_bar + theme(legend.position='none'),0,.14,.05,.8572) 
        #draw_line(x=c(.0602,.981),y=c(.9,.9),color="black",size=.2)

    start_y <- .959
    scaling_factor <- .765
    for (i in 1:19) {
        new_y <- start_y - ((sum(fraction_from_gene_cluster[1:i]))*scaling_factor)
        combined <- combined + draw_line(x=c(.0602,.981),y=c(new_y,new_y),color="black",size=.2)
    }

  
    return(combined)

}

make_figure_1_alt_a <- function(sample_info, quant_expr, mixutre_hmm_cell_line_grouping_dir,cell_line_cluster_1_average_expression_file, cell_line_cluster_2_average_expression_file, output_file,fig1a_image_file) {
    fig1b <- plot_pca_time_step_modular(sample_info, quant_expr)

    #fig1c <- plot_karl_gene_cluster(mixutre_hmm_cell_line_grouping_dir, 'Gene Cluster 01')

    #fig1d <- plot_karl_gene_cluster(mixutre_hmm_cell_line_grouping_dir, 'Gene Cluster 03')

    #cell_line_legend <- get_legend(fig1d)

    fig1c <- heatmap_showing_average_expression_for_gene_cluster(cell_line_cluster_1_average_expression_file, cell_line_cluster_2_average_expression_file)

    #fig1d <- heatmap_showing_average_expression_for_gene_cluster(cell_line_cluster_average_expression_file, "2")

    pca_legend <- get_legend(fig1b)

    combined <- ggdraw() + 
                draw_image(fig1a_image_file,0,.56,.4,.4) +
                draw_plot(fig1b + theme(legend.position='none'), .505,.65,.4,.35) + 
                draw_plot(pca_legend,.9,.43,.9,.9) + 
                draw_plot_label(c("a","b"),c(0,.5),c(1,1),size=12)
    counter <- 1
    start_y <- .66
    step_size <- .05
    for (i in 1:10) {
        combined <- combined + draw_plot(fig1c[[counter]]+ theme(legend.position='none'), 0,start_y-(i*step_size),.5,.06) 
        counter <- counter + 1
        combined <- combined + draw_plot(fig1c[[counter]]+ theme(legend.position='none'), .5,start_y-(i*step_size),.5,.06) 
        counter <- counter + 1
    }
    #counter <- 1
    #start_y <- .74
    #for (i in 1:10) {
    #    combined <- combined + draw_plot_label(paste0("Gene cluster ", counter),c(.17),c(start_y-(i*step_size)),size=6)
    #    counter <- counter + 1
    #    combined <- combined + draw_plot_label(paste0("Gene cluster ", counter),c(.67),c(start_y-(i*step_size)),size=6)
    #    counter <- counter + 1
   # }


    ggsave(combined, file=output_file, width=7.2, height=5.1, units="in")
}

make_figure_1_alt_b <- function(sample_info, quant_expr, mixutre_hmm_cell_line_grouping_dir,cell_line_cluster_1_average_expression_file, cell_line_cluster_2_average_expression_file, output_file,fig1a_image_file) {
    fig1b <- plot_pca_time_step_modular(sample_info, quant_expr)

    #fig1c <- plot_karl_gene_cluster(mixutre_hmm_cell_line_grouping_dir, 'Gene Cluster 01')

    #fig1d <- plot_karl_gene_cluster(mixutre_hmm_cell_line_grouping_dir, 'Gene Cluster 03')

    #cell_line_legend <- get_legend(fig1d)

    #fig1c <- heatmap_showing_average_expression_for_gene_cluster(cell_line_cluster_1_average_expression_file, cell_line_cluster_2_average_expression_file)

    fig1c <- heatmap_showing_average_expression_for_cell_line_cluster(cell_line_cluster_1_average_expression_file, "1")

    fig1d <- heatmap_showing_average_expression_for_cell_line_cluster(cell_line_cluster_2_average_expression_file, "2")

    pca_legend <- get_legend(fig1b)

    combined <- ggdraw() + 
                draw_image(fig1a_image_file,0,.6,.5,.5) +
                draw_plot(fig1b + theme(legend.position='none'), .505,.65,.4,.35) + 
                draw_plot(pca_legend,.9,.43,.9,.9) + 
                draw_plot_label(c("a","b"),c(0,.5),c(1,1),size=12) + 
                draw_plot_label("figure 1a not done yet",c(0.12),c(.98),size=7) +
                draw_plot(fig1c+ theme(legend.position='none'),0,0,.5,.5) +
                draw_plot(fig1d+ theme(legend.position='none'),.5,0,.5,.5) +
                draw_plot_label("Cell line cluster 1",c(.16),c(.52),size=8) + 
                draw_plot_label("Cell line cluster 2",c(.66),c(.52),size=8)

    ggsave(combined, file=output_file, width=7.2, height=5.1, units="in")
}

make_figure_1_alt_c <- function(sample_info, quant_expr, mixutre_hmm_cell_line_grouping_dir,cell_line_cluster_1_average_expression_file, cell_line_cluster_2_average_expression_file, output_file,fig1a_image_file) {
    fig1b <- plot_pca_time_step_modular(sample_info, quant_expr)


    fig1c <- heatmap_showing_average_expression_for_all_clusters(cell_line_cluster_1_average_expression_file, cell_line_cluster_average_expression_file)

    pca_legend <- get_legend(fig1b)
    print("hello")

    combined <- ggdraw() + 
                draw_image(fig1a_image_file,0,.6,.5,.5) +
                draw_plot(fig1b + theme(legend.position='none'), .505,.65,.4,.35) + 
                draw_plot(pca_legend,.9,.43,.9,.9) + 
                draw_plot_label(c("a","b"),c(0,.5),c(1,1),size=12) + 
                draw_plot_label("figure 1a not done yet",c(0.12),c(.98),size=7) + 
                draw_plot(fig1c,0,0,1,.6)


    ggsave(combined, file=output_file, width=7.2, height=5.1, units="in")
}



make_figure_1 <- function(sample_info, quant_expr, mixutre_hmm_cell_line_grouping_dir, output_file,fig1a_image_file) {
    fig1b <- plot_pca_time_step_modular(sample_info, quant_expr)
    pca_legend <- get_legend(fig1b)


    fig1c <- plot_karl_gene_cluster_grid(mixutre_hmm_cell_line_grouping_dir) + theme(strip.background = element_blank(),strip.text.x = element_blank())
    cell_line_cluster_legend <- get_legend(fig1c+ theme(legend.position="bottom")) 


    combined <- ggdraw() + 
                draw_image(fig1a_image_file,.003,.545,.47,.47) +
                draw_plot(fig1b + theme(legend.position='none'), .505,.54,.46,.45) + 
                draw_plot(pca_legend,.945,.37,.9,.9) + 
                draw_plot_label(c("A","B","C"),c(0,.5,0),c(1,1,.6),size=12) + 
                draw_plot(fig1c + theme(legend.position="none"),-.02,.05,1.02,.5) + 
                draw_plot(cell_line_cluster_legend,.32,-.47,1,1)
    ggsave(combined, file=output_file, width=7.2, height=4.3, units="in")
}


make_figure_1_alt_e <- function(sample_info, quant_expr, mixutre_hmm_cell_line_grouping_dir, output_file,fig1a_image_file) {
    fig1b <- plot_pca_time_step_modular(sample_info, quant_expr)
    pca_legend <- get_legend(fig1b)


    fig1c_row1 <- plot_karl_gene_cluster_grid_v2(mixutre_hmm_cell_line_grouping_dir, "row1") + theme(strip.background = element_blank(),strip.text.x = element_blank()) + theme(axis.text.x=element_blank()) + labs(x="",y="")
    fig1c_row2 <- plot_karl_gene_cluster_grid_v2(mixutre_hmm_cell_line_grouping_dir, "row2") + theme(strip.background = element_blank(),strip.text.x = element_blank()) + labs(y="")

    cell_line_cluster_legend <- get_legend(fig1c_row1+ theme(legend.position="bottom")) 


    combined <- ggdraw() + 
                draw_image(fig1a_image_file,0,.5,.5,.5) +
                draw_plot(fig1b + theme(legend.position='none'), .505,.54,.4,.45) + 
                draw_plot(pca_legend,.9,.39,.9,.9) + 
                draw_plot_label(c("a","b","c"),c(0,.5,0),c(1,1,.6),size=12) + 
                draw_plot_label("figure 1a not done yet",c(0.12),c(.98),size=7) + 
                draw_plot(fig1c_row1 + theme(legend.position="none"),-.02,.27,1.02,.29) + 
                draw_plot(fig1c_row2 + theme(legend.position="none"),-.02,.03,1.02,.32) + 
                draw_plot(cell_line_cluster_legend,.32,-.47,1,1) + 
                draw_plot_label("Expression", c(0),c(.17),size=12, angle=90,fontface="plain")

    #start_x <- -.05
    #scale <- .092
    #y1 <- .54 
    #y2 <- .35
    #for (gene_cluster in 1:10) {
    #    combined <- combined + draw_plot_label(paste0("Cluster ", gene_cluster), c(start_x + scale*gene_cluster), c(y1), size=8,fontface="plain")
    #    combined <- combined + draw_plot_label(paste0("Cluster ", 10+ gene_cluster), c(start_x + scale*gene_cluster), c(y2),size=8,fontface="plain")
    #}


    ggsave(combined, file=output_file, width=7.2, height=4.3, units="in")
}











#####################################################################################################
# Load Data
#####################################################################################################

options(warn=1)
preprocess_total_expression_dir = args[1]  # Where total expression processed data is_autosomal
visualize_total_expression_dir = args[2]  # Ouputdir to save images
covariate_dir = args[3]  # Input dir with covariate information
mixutre_hmm_cell_line_grouping_dir = args[4]  # Directory containing files assigning cell lines to groupings

#  Get sample information 
sample_info_file <- paste0(preprocess_total_expression_dir, "sample_info.txt")
sample_info <- read.table(sample_info_file, header=TRUE)

#  Get quantile normalized expression data
quantile_normalized_exp_file <- paste0(preprocess_total_expression_dir, "quantile_normalized.txt")
quant_expr <- read.csv(quantile_normalized_exp_file, header=TRUE, sep=" ")



##### USED to save time if already run
#saveRDS(sample_info, paste0(preprocess_total_expression_dir,"sample_info.rds"))
#saveRDS(quant_expr, paste0(preprocess_total_expression_dir,"quant_expr.rds"))

#sample_info <- readRDS(paste0(preprocess_total_expression_dir,"sample_info.rds"))
#quant_expr <- readRDS(paste0(preprocess_total_expression_dir,"quant_expr.rds"))


#  Get quantile normalized expression (done seperately for each time point) data
quantile_normalized_time_independent_expression_file <- paste0(preprocess_total_expression_dir, "time_step_independent_quantile_normalized.txt")
time_step_independent_quant_expr <- read.csv(quantile_normalized_time_independent_expression_file, header=TRUE, sep=" ")

#  Get rpkm expression_data
rpkm_exp_file <- paste0(preprocess_total_expression_dir, "rpkm.txt")
rpkm_expr <- read.csv(rpkm_exp_file, header=TRUE, sep=" ")

#  Get covariate file
covariate_file <- paste0(covariate_dir, "processed_covariates_categorical.txt")
covariates <- read.table(covariate_file,header=TRUE)

# Get standardized cell line specific expression
cell_line_expression_ignore_missing_file <- paste0(preprocess_total_expression_dir,"cell_line_expression_ignore_missing.txt")
cell_line_expression_ignore_missing <- read.csv(cell_line_expression_ignore_missing_file, header=TRUE, sep="\t")
cell_line_expression_ignore_missing <- cell_line_expression_ignore_missing[,2:(dim(cell_line_expression_ignore_missing)[2])]
colnames(cell_line_expression_ignore_missing) <- substr(colnames(cell_line_expression_ignore_missing), 2, 1000)


#####################################################################################################
# Run Analysis / Create Plots
#####################################################################################################



####################################################################
# Plot Figure 1
####################################################################
output_file <- paste0(visualize_total_expression_dir, "figure1.png")
fig1a_image_file <- "/project2/gilad/bstrober/ipsc_differentiation_19_lines/preprocess_input_data/fig1a_mock.png"
make_figure_1(sample_info, quant_expr, mixutre_hmm_cell_line_grouping_dir, output_file,fig1a_image_file)



####################################################################
# Plot EDF showing Nanog and troponin time courses for each cell line
####################################################################
output_file <- paste0(visualize_total_expression_dir, "edf_troponin_nanog_time_course.png")
make_edf_troponin_nanog_time_course(sample_info, quant_expr, output_file)




####################################################################
# Plot library size
####################################################################
################
# Make barplot showing library sizes of each sample
library_size_output_file <- paste0(visualize_total_expression_dir, "library_size.pdf")
plot_library_size(sample_info, library_size_output_file)



####################################################################
# PCA plots with samples labeled/colored by various covariates
####################################################################

##################
#  Perform PCA. Plot first 2 pcs as a function of time step 
pca_plot_time_step_output_file <- paste0(visualize_total_expression_dir, "pca_plot_1_2_time_step.png")
plot_pca_time_step(sample_info, quant_expr, pca_plot_time_step_output_file)




#################
# Perform PCA. Plot specified PCs with samples colored/labeled by their gene expression of:
# a. Troponin (gene expressed in cardiomyocytes)
# b. sox2 (gene expressed in ipscs)
# c. nanog (gene expressed in ipscs)
pc_num1 <- 1
pc_num2 <- 2
ensamble_id <- "ENSG00000118194"
gene_name <- "Troponin"
pca_plot_gene_filled_output_file <- paste0(visualize_total_expression_dir, "pca_plot_",pc_num1,"_",pc_num2,"_",gene_name,"_gene_filled.png")
plot_pca_real_valued_gene_filled(sample_info, quant_expr, ensamble_id,gene_name,pc_num1,pc_num2,pca_plot_gene_filled_output_file)


pc_num1 <- 1
pc_num2 <- 2
ensamble_id <- "ENSG00000181449"
gene_name <- "sox2"
pca_plot_gene_filled_output_file <- paste0(visualize_total_expression_dir, "pca_plot_",pc_num1,"_",pc_num2,"_",gene_name,"_gene_filled.png")
plot_pca_real_valued_gene_filled(sample_info, quant_expr, ensamble_id,gene_name,pc_num1,pc_num2,pca_plot_gene_filled_output_file)


pc_num1 <- 1
pc_num2 <- 2
ensamble_id <- "ENSG00000111704"
gene_name <- "nanog"
pca_plot_gene_filled_output_file <- paste0(visualize_total_expression_dir, "pca_plot_",pc_num1,"_",pc_num2,"_",gene_name,"_gene_filled.png")
plot_pca_real_valued_gene_filled(sample_info, quant_expr, ensamble_id,gene_name,pc_num1,pc_num2,pca_plot_gene_filled_output_file)



#################
# Perform Cell line pca PCA. Plot specified PCs as a function of gene expression of:
# a. Troponin (gene expressed in cardiomyocytes)
# b. sox2 (gene expressed in ipscs)
pc_num1 <- 1
pc_num2 <- 2
time_step <- 15
ensamble_id <- "ENSG00000118194"
gene_name <- "Troponin"
cell_line_pca_plot_gene_filled_output_file <- paste0(visualize_total_expression_dir, "cell_line_ignore_missing_pca_plot_",pc_num1,"_",pc_num2,"_time_",time_step,"_",gene_name,"_gene_filled.png")
plot_cell_line_pca_real_valued_gene_filled(colnames(cell_line_expression_ignore_missing), cell_line_expression_ignore_missing, sample_info, quant_expr, ensamble_id,gene_name,pc_num1,pc_num2,cell_line_pca_plot_gene_filled_output_file, time_step)




#################
# Perform PCA. Plot specified PCS as a function of various covariates:
# a. cell line
# b. rna extraction person
# c. rna extraction round
# d. differentiation batch

pc_num1 <- 1
pc_num2 <- 2

pca_plot_cell_line_output_file <- paste0(visualize_total_expression_dir, "pca_plot_",pc_num1,"_",pc_num2,"_cell_line.png")
plot_pca_categorical_covariate(sample_info, quant_expr, pca_plot_cell_line_output_file,factor(paste0("NA",sample_info$cell_line)), "Cell Line", pc_num1,pc_num2)

pca_plot_cell_line_output_file <- paste0(visualize_total_expression_dir, "pca_plot_",pc_num1,"_",pc_num2,"_rna_extraction_persion.png")
plot_pca_categorical_covariate(sample_info, quant_expr, pca_plot_cell_line_output_file,factor(covariates$RNA_extraction_person), "rna_extraction_person", pc_num1,pc_num2)

pca_plot_cell_line_output_file <- paste0(visualize_total_expression_dir, "pca_plot_",pc_num1,"_",pc_num2,"_rna_extraction_round.png")
plot_pca_categorical_covariate(sample_info, quant_expr, pca_plot_cell_line_output_file,factor(covariates$RNA_extraction_round), "rna_extraction_round", pc_num1,pc_num2)

pca_plot_cell_line_output_file <- paste0(visualize_total_expression_dir, "pca_plot_",pc_num1,"_",pc_num2,"_differentiation_batch.png")
plot_pca_categorical_covariate(sample_info, quant_expr, pca_plot_cell_line_output_file,factor(covariates$differentiation_batch), "differentiation_batch", pc_num1,pc_num2)




###############################
# Plot eigenvectors
################################
eigenvectors_output_file <- paste0(visualize_total_expression_dir, "pca_plot_eigenvector_viz.pdf")
plot_pca_eigenvectors(sample_info, quant_expr, eigenvectors_output_file)

eigenvectors_output_file <- paste0(visualize_total_expression_dir, "pca_plot_eigenvector_viz_by_line.pdf")
plot_pca_eigenvectors_by_line(sample_info, quant_expr, eigenvectors_output_file)


####################################################################
# Variance explained line plots
####################################################################

# Variance explained of first 7 pcs from cell PCA
n <- 7
cell_line_pca_plot_variance_explained_output_file <- paste0(visualize_total_expression_dir, "cell_line_ignore_missing_pca_plot_variance_explained", n, ".png")
cell_line_pc_pve <- plot_pca_variance_explained(colnames(cell_line_expression_ignore_missing), cell_line_expression_ignore_missing, n, cell_line_pca_plot_variance_explained_output_file, "Cell line collapsed PC number")


#################
# Perform PCA on full quantile normalized matrix. Plot variance explained of the first n PCs:
n <- 20
pca_plot_variance_explained_output_file <- paste0(visualize_total_expression_dir, "pca_plot_variance_explained", n, ".png")
plot_pca_variance_explained(sample_info, quant_expr, n, pca_plot_variance_explained_output_file)

n <- 10
pca_plot_variance_explained_output_file <- paste0(visualize_total_expression_dir, "pca_plot_variance_explained", n, ".png")
pc_pve <- plot_pca_variance_explained(sample_info, quant_expr, n, pca_plot_variance_explained_output_file, "PC number")





####################################################################
# Time course of expression of specific genes seperated by cell line
####################################################################
ensamble_id <- "ENSG00000118194"
gene_name <- "Troponin"
line_plot_file <- paste0(visualize_total_expression_dir, gene_name,"_time_course_grouped_by_cell_line.png")
gene_time_course_line_plot_grouped_by_cell_line(sample_info, quant_expr, ensamble_id, gene_name, line_plot_file)

ensamble_id <- "ENSG00000181449"
gene_name <- "sox2"
line_plot_file <- paste0(visualize_total_expression_dir, gene_name,"_time_course_grouped_by_cell_line.png")
gene_time_course_line_plot_grouped_by_cell_line(sample_info, quant_expr, ensamble_id, gene_name, line_plot_file)


ensamble_id <- "ENSG00000111704"
gene_name <- "nanog"
line_plot_file <- paste0(visualize_total_expression_dir, gene_name,"_time_course_grouped_by_cell_line.png")
gene_time_course_line_plot_grouped_by_cell_line(sample_info, quant_expr, ensamble_id, gene_name, line_plot_file)




####################################################################
# Covariates explaining variance in principle components
####################################################################


# Make heatmap showing PVE between pcs and covariates
pc_file <- paste0(covariate_dir,"principal_components_10.txt")
covariate_file <- paste0(covariate_dir, "processed_covariates_categorical.txt")
output_file <- paste0(visualize_total_expression_dir, "pc_covariate_pve_heatmap.png")
heatmap <- covariate_pc_pve_heatmap(pc_file, covariate_file,output_file, "PCA")

combined_output_file <- paste0(visualize_total_expression_dir, "pc_covariate_pve_heatmap_joint.png")
combined <- plot_grid(pc_pve, heatmap, labels = c("A", "B"), ncol=1,rel_heights = c(.7, 1.2))
ggsave(combined, file=combined_output_file, width=7.2, height=5.5,units="in")








####################################################################
# Cell line PC1 vs PC2 colored by:
####################################################################
######### 1. Troponin expression
######### 2. Karl's groupings
######### 3. Flow results

# File containing cell line PCs
cell_pc_file <- paste0(covariate_dir,"cell_line_ignore_missing_principal_components_9.txt")

# Results from Karl's clustering
bi_clustering_state_file <- paste0(mixutre_hmm_cell_line_grouping_dir, "mixsvgp_K2_L100_1_28542829_assignments")

# Flow cytometry results
flow_file <- paste0(mixutre_hmm_cell_line_grouping_dir, "flow_results.txt")



# avg10-15 troponin expression
output_file <- paste0(visualize_total_expression_dir, "cell_line_pc1_2_colored_by_troponin_expression.png")
cell_line_pc_colored_by_avg_troponin(cell_pc_file, covariates, output_file)

# bi-clustering model
output_file <- paste0(visualize_total_expression_dir, "cell_line_pc1_2_colored_by_gp_mixture_model.png")
cell_line_pc_scatter_model_xx <- cell_line_pc_colored_by_state_model(cell_pc_file, bi_clustering_state_file, output_file)

# FLow results
output_file <- paste0(visualize_total_expression_dir, "cell_line_pc1_2_colored_by_flow_results.png")
cell_line_pc_scatter_flow <- cell_line_pc_colored_by_state_model_real_valued(cell_pc_file, flow_file, output_file)



# Make cowplot combined plot of cell_line_pc_scatter_model_xx and cell_line_pc_scatter_flow
output_file <- paste0(visualize_total_expression_dir, "cell_line_pc1_2_model_xx_and_flow_and_pve_joint.png")
combined <- plot_grid(cell_line_pc_scatter_flow+ theme(legend.position='bottom'), cell_line_pc_scatter_model_xx+ theme(legend.position='bottom'), labels = c("B", "C"), align = "h")
combined2 <- plot_grid(cell_line_pc_pve, combined, labels=c("A", ""), ncol=1, rel_heights = c(.9, 1.2))
ggsave(combined2, file=output_file, width=7.2, height=5.5,units="in")



####################################################################
# Correlate flow results with PC results, as well as avg troponin
####################################################################
output_file <- paste0(visualize_total_expression_dir, "flow_avg_10_15_troponin_scatter.png")
scatter_of_avg_troponin_and_flow(flow_file, covariates, output_file)

pc_num <- 2
output_file <- paste0(visualize_total_expression_dir, "flow_pc",pc_num,"_troponin_scatter.png")
scatter_of_pc_and_flow(flow_file, cell_pc_file, pc_num, output_file)








####################################################################
# Look at concordance of karls results across different runs
####################################################################
output_file <- paste0(visualize_total_expression_dir,"split_gpm_concordance.png")

l <- 10
input_file <- paste0(mixutre_hmm_cell_line_grouping_dir,"K2L",l,"_incidence")
l10_concordance_heatmap <- split_gpm_concordance_heatmap(input_file, l)

l <- 20
input_file <- paste0(mixutre_hmm_cell_line_grouping_dir,"K2L",l,"_incidence")
l20_concordance_heatmap <- split_gpm_concordance_heatmap(input_file, l)

l <- 50
input_file <- paste0(mixutre_hmm_cell_line_grouping_dir,"K2L",l,"_incidence")
l50_concordance_heatmap <- split_gpm_concordance_heatmap(input_file, l)

l <- 100
input_file <- paste0(mixutre_hmm_cell_line_grouping_dir,"K2L",l,"_incidence")
l100_concordance_heatmap <- split_gpm_concordance_heatmap(input_file, l)


combined <- plot_grid(l10_concordance_heatmap, l20_concordance_heatmap, l50_concordance_heatmap, l100_concordance_heatmap, labels = c("A", "B", "C", "D"), ncol=2)

ggsave(combined, file=output_file,width=7.2,height=5.5,units="in")












































































#####################################
# Retired scripts
#####################################

#################
# Perform PCA. Make seperate plot for each cell line:
pc_num1<-1
pc_num2<-2
pca_plot_cell_line_output_file <- paste0(visualize_total_expression_dir, "pca_plot_",pc_num1,"_",pc_num2,"_seperate_cell_lines.pdf")
#plot_pca_seperate_cell_lines(sample_info, quant_expr, pca_plot_cell_line_output_file,pc_num1,pc_num2)



####################################################################
# Make plot showing average expression in each cell line cluster
####################################################################
cell_line_cluster <- "1"
cell_line_cluster_average_expression_file <- paste0(preprocess_total_expression_dir, "cell_line_cluster", cell_line_cluster, "_average_expression.txt")
output_file <- paste0(visualize_total_expression_dir, "cell_line_cluster", cell_line_cluster,"_heatmap.png")
#heatmap_cluster1 <- heatmap_showing_average_expression_for_cell_line_cluster(cell_line_cluster_average_expression_file, cell_line_cluster)
# ggsave(heatmap_cluster1, file=output_file, width=7.2, height=1.5,units="in")

cell_line_cluster <- "2"
cell_line_cluster_average_expression_file <- paste0(preprocess_total_expression_dir, "cell_line_cluster", cell_line_cluster, "_average_expression.txt")
output_file <- paste0(visualize_total_expression_dir, "cell_line_cluster", cell_line_cluster,"_heatmap.png")
#heatmap_cluster1 <- heatmap_showing_average_expression_for_cell_line_cluster(cell_line_cluster_average_expression_file, cell_line_cluster)
#ggsave(heatmap_cluster1, file=output_file, width=7.2, height=5.5,units="in")




####################################################################
# Correlation heatmap between time step independent latent factors and global latent factors
####################################################################


n_independent <- 5
n_global <- 7
pca_correlation_heatmap_output_file <- paste0(visualize_total_expression_dir, "pca_correlation_heatmap_time_independent_v_global.pdf")
# pca_correlation_heatmap_time_independent_v_global(sample_info, quant_expr, time_step_independent_quant_expr, n_global, n_independent, pca_correlation_heatmap_output_file)



# Perform PCA on time-step independent quantile normalized matrix. Plot variance explained of first n PCs:
n <- 6
pca_plot_time_independent_output_file <- paste0(visualize_total_expression_dir, "pca_plot_time_step_independent_variance_explained",n,".png")
#plot_pca_time_independent_variance_explained(sample_info, time_step_independent_quant_expr, n, pca_plot_time_independent_output_file)



pc_file <- paste0(covariate_dir,"cell_line_ignore_missing_principal_components_5_time_time_per_sample.txt")
covariate_file <- paste0(covariate_dir, "processed_covariates_categorical.txt")
output_file <- paste0(visualize_total_expression_dir, "cell_line_pcXtime_covariate_pve_heatmap.png")
#covariate_pc_pve_heatmap(pc_file, covariate_file,output_file, "cell-line PCA")




# File containing assignment of each cell line to states in the 12 state model
state_file_12 <- paste0(mixutre_hmm_cell_line_grouping_dir, "12state")

# File containing assignmnets of each cell line to states in the 8 state model
state_file_8 <- paste0(mixutre_hmm_cell_line_grouping_dir, "8state")

# 8 state model
output_file <- paste0(visualize_total_expression_dir, "cell_pc1_2_colored_by_8_state_model.png")
#cell_line_pc_colored_by_state_model_hmm(cell_pc_file, state_file_8, output_file)

# 12 state model
output_file <- paste0(visualize_total_expression_dir, "cell_pc1_2_colored_by_12_state_model.png")
#cell_line_pc_colored_by_state_model_hmm(cell_pc_file, state_file_12, output_file)


# 0:15
for (time_step in 0:15) {
    pc_file <- paste0(covariate_dir, "pca_loadings_time_step_independent_", time_step, ".txt")
    covariate_file <- paste0(covariate_dir, "processed_covariates_categorical.txt")
    output_file <- paste0(visualize_total_expression_dir, "time_step_", time_step, "_pc_covariate_pve_heatmap.png")
    #time_step_specific_covariate_pc_pve_heatmap(pc_file, covariate_file, output_file, time_step)
}




# Make heatmap showing PVE between pcs & (covariates and troponin/sox2 expression)
pc_file <- paste0(covariate_dir,"principal_components_10.txt")
covariate_file <- paste0(covariate_dir, "processed_covariates_categorical.txt")
output_file <- paste0(visualize_total_expression_dir, "pc_covariate_troponin_sox2_pve_heatmap.png")
#covariate_pc_specific_genes_pve_heatmap(pc_file, quant_expr, covariate_file,output_file)















pc_num <-1
eigenvectors_output_file <- paste0(visualize_total_expression_dir, "pca_plot_", pc_num, "_eigenvector_viz_by_line_independent.pdf")
# plot_pca_eigenvectors_by_line_independent(sample_info, quant_expr, eigenvectors_output_file, pc_num)


pc_num <-2
eigenvectors_output_file <- paste0(visualize_total_expression_dir, "pca_plot_", pc_num, "_eigenvector_viz_by_line_independent.pdf")
# plot_pca_eigenvectors_by_line_independent(sample_info, quant_expr, eigenvectors_output_file, pc_num)


pc_num <-3
eigenvectors_output_file <- paste0(visualize_total_expression_dir, "pca_plot_", pc_num, "_eigenvector_viz_by_line_independent.pdf")
# plot_pca_eigenvectors_by_line_independent(sample_info, quant_expr, eigenvectors_output_file, pc_num)


pc_num <-4
eigenvectors_output_file <- paste0(visualize_total_expression_dir, "pca_plot_", pc_num, "_eigenvector_viz_by_line_independent.pdf")
#plot_pca_eigenvectors_by_line_independent(sample_info, quant_expr, eigenvectors_output_file, pc_num)



###########################################################################################
# Scatter plot of 2nd (and others) PC loading vs a genes expression levels in t = 15 samples
###########################################################################################
ensamble_id <- "ENSG00000118194"
gene_name <- "Troponin"
time_step <- 15
pc_num <- 2

pc_gene_scatter_output_file <- paste0(visualize_total_expression_dir, "pc_num_", pc_num, "_",gene_name,"_time_step_",time_step,"_scatter.pdf")
#pc_gene_scatter(sample_info, quant_expr, ensamble_id, gene_name, time_step, pc_num, pc_gene_scatter_output_file)


ensamble_id <- "ENSG00000181449"
gene_name <- "sox2"
time_step <- 15
pc_num <- 2

pc_gene_scatter_output_file <- paste0(visualize_total_expression_dir, "pc_num_", pc_num, "_",gene_name,"_time_step_",time_step,"_scatter.pdf")
#pc_gene_scatter(sample_info, quant_expr, ensamble_id, gene_name, time_step, pc_num, pc_gene_scatter_output_file)


###########################################################################################
# Scatter plot of 2nd (and others) PC loading vs a genes expression levels at all time steps
###########################################################################################


ensamble_id <- "ENSG00000118194"
gene_name <- "Troponin"
pc_num <- 2

pc_gene_scatter_output_file <- paste0(visualize_total_expression_dir, "pc_num_", pc_num, "_",gene_name,"_colored_by_cell_line_scatter.pdf")
#pc_gene_scatter_all_time_colored_by_cell_line(sample_info, quant_expr, ensamble_id, gene_name, pc_num, pc_gene_scatter_output_file)


pc_gene_scatter_output_file <- paste0(visualize_total_expression_dir, "pc_num_", pc_num, "_",gene_name,"_colored_by_time_scatter.pdf")
#pc_gene_scatter_all_time_colored_by_time(sample_info, quant_expr, ensamble_id, gene_name, pc_num, pc_gene_scatter_output_file)

###########################################################################################
# Scatter plot of 2nd (and others) PC loading vs a genes AVERAGE expression
###########################################################################################
pc_num <- 2
covariate_file <- paste0(covariate_dir, "processed_covariates_categorical.txt")
output_file <- paste0(visualize_total_expression_dir,"pc_num_", pc_num,"_avg_10_15_troponin_colored_by_cell_line_scatter.pdf")
#pc_avg_troponin_scatter_colored_by_cell_line(sample_info, quant_expr, pc_num, covariate_file, output_file)
pc_num <- 2
covariate_file <- paste0(covariate_dir, "processed_covariates_categorical.txt")
output_file <- paste0(visualize_total_expression_dir,"pc_num_", pc_num,"_avg_10_15_troponin_colored_by_time_step_scatter.pdf")
#pc_avg_troponin_scatter_colored_by_time_step(sample_info, quant_expr, pc_num, covariate_file, output_file)



####################################################################
# Comparison of old vs new cell lines (will be removed)
####################################################################

# Add old vs new to sample info
revised_sample_info <- add_old_vs_new(sample_info)


##################
#  Perform PCA. Plot first 2 pcs as a function of time step 
pca_plot_time_step_output_file <- paste0(visualize_total_expression_dir, "pca_plot_1_2_old_vs_new.png")
#plot_pca_old_vs_new(revised_sample_info, quant_expr, pca_plot_time_step_output_file)


###################
# Time course of expression of specific genes seperated by cell line (and colored by old vs_new)
ensamble_id <- "ENSG00000118194"
gene_name <- "Troponin"
line_plot_file <- paste0(visualize_total_expression_dir, gene_name,"_time_course_grouped_by_cell_line_old_vs_new.png")
#gene_time_course_line_plot_grouped_by_cell_line_old_vs_new(revised_sample_info, quant_expr, ensamble_id, gene_name, line_plot_file)

ensamble_id <- "ENSG00000181449"
gene_name <- "sox2"
line_plot_file <- paste0(visualize_total_expression_dir, gene_name,"_time_course_grouped_by_cell_line_old_vs_new.png")
#gene_time_course_line_plot_grouped_by_cell_line_old_vs_new(revised_sample_info, quant_expr, ensamble_id, gene_name, line_plot_file)


ensamble_id <- "ENSG00000111704"
gene_name <- "nanog"
line_plot_file <- paste0(visualize_total_expression_dir, gene_name,"_time_course_grouped_by_cell_line_old_vs_new.png")
#gene_time_course_line_plot_grouped_by_cell_line_old_vs_new(revised_sample_info, quant_expr, ensamble_id, gene_name, line_plot_file)





####################################################################
# Plot Figure 1_alt
####################################################################
cell_line_cluster_1_average_expression_file <- paste0(preprocess_total_expression_dir, "cell_line_cluster1_average_expression.txt")
cell_line_cluster_2_average_expression_file <- paste0(preprocess_total_expression_dir, "cell_line_cluster2_average_expression.txt")


output_file <- paste0(visualize_total_expression_dir, "figure1_alt_a.png")
fig1a_image_file <- "/project2/gilad/bstrober/ipsc_differentiation_19_lines/preprocess_input_data/fig1a.png"
#make_figure_1_alt_a(sample_info, quant_expr, mixutre_hmm_cell_line_grouping_dir, cell_line_cluster_1_average_expression_file, cell_line_cluster_2_average_expression_file, output_file,fig1a_image_file)


output_file <- paste0(visualize_total_expression_dir, "figure1_alt_b.png")
fig1a_image_file <- "/project2/gilad/bstrober/ipsc_differentiation_19_lines/preprocess_input_data/fig1a.png"
#make_figure_1_alt_b(sample_info, quant_expr, mixutre_hmm_cell_line_grouping_dir, cell_line_cluster_1_average_expression_file, cell_line_cluster_2_average_expression_file, output_file,fig1a_image_file)

output_file <- paste0(visualize_total_expression_dir, "figure1_alt_a.png")
fig1a_image_file <- "/project2/gilad/bstrober/ipsc_differentiation_19_lines/preprocess_input_data/fig1a.png"
#make_figure_1_alt_a(sample_info, quant_expr, mixutre_hmm_cell_line_grouping_dir, cell_line_cluster_1_average_expression_file, cell_line_cluster_2_average_expression_file, output_file,fig1a_image_file)


output_file <- paste0(visualize_total_expression_dir, "figure1_alt_c.png")
fig1a_image_file <- "/project2/gilad/bstrober/ipsc_differentiation_19_lines/preprocess_input_data/fig1a.png"
#make_figure_1_alt_c(sample_info, quant_expr, mixutre_hmm_cell_line_grouping_dir, cell_line_cluster_1_average_expression_file, cell_line_cluster_2_average_expression_file, output_file,fig1a_image_file)


output_file <- paste0(visualize_total_expression_dir, "figure1_alt_e.png")
fig1a_image_file <- "/project2/gilad/bstrober/ipsc_differentiation_19_lines/preprocess_input_data/fig1a.png"
#make_figure_1_alt_e(sample_info, quant_expr, mixutre_hmm_cell_line_grouping_dir, output_file,fig1a_image_file)







