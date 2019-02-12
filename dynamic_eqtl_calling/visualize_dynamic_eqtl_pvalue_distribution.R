args = commandArgs(trailingOnly=TRUE)
library(ggplot2)
library(ggthemes)
library(cowplot)
library(reshape)




pvalue_histogram <- function(pvalues, output_file) {
    df <- data.frame(pvalue=pvalues)

    histo <- ggplot(df, aes(x = pvalues)) +
        geom_histogram(aes(y = ..density..), breaks = seq(0,1,by=.01)) +
        theme(text = element_text(size=18), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black")) +
        labs(x = "pvalue", y = "density", title="nominal pvalues")
    ggsave(histo, file=output_file,width = 20,height=10.5,units="cm")
}


qq_plot_cmp_to_uniform <- function(pvalues, perm_pvalue_arr, output_file) {
    real_pvalues <- sort(pvalues)

    perm_pvalues <- sort(perm_pvalue_arr)

    uniform_pvalues_1 <- sort(runif(length(real_pvalues)))
    uniform_pvalues_2 <- sort(runif(length(perm_pvalues)))

    cat_pvalues <- c(real_pvalues, perm_pvalues)
    cat_uniform_pvalues <- c(uniform_pvalues_1, uniform_pvalues_2)

    version <- c(rep("real", length(real_pvalues)), rep("permuted", length(perm_pvalues)))

    df <- data.frame(pvalues=-log10(cat_pvalues + .000000000001), uniform_pvalues=-log10(cat_uniform_pvalues + .000000000001), version=factor(version))
    
    # PLOT!
    max_val <-max(max(-log10(cat_uniform_pvalues + .000000000001)), max(-log10(cat_pvalues + .000000000001)))

    #PLOT!
    scatter <- ggplot(df, aes(x = uniform_pvalues, y = pvalues, colour=version)) + geom_point() 
    scatter <- scatter + theme(text = element_text(size=14), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black")) 
    scatter <- scatter + labs(x = "Uniform", y = "Real")
    scatter <- scatter + geom_abline() 
    scatter <- scatter + scale_x_continuous(limits = c(-.1, max_val + .1), breaks = round(seq(0, max_val, by = 5),1)) + scale_y_continuous(limits = c(-.1,max_val+.1), breaks = round(seq(0, max_val, by = 5),1))

    ggsave(scatter, file=output_file, width=20, height=10.5, units="cm")
}


qq_plot_cmp_true_to_uniform <- function(pvalues, output_file) {
    real_pvalues <- sort(pvalues)


    uniform_pvalues_1 <- sort(runif(length(real_pvalues)))

    cat_pvalues <- c(real_pvalues)
    cat_uniform_pvalues <- c(uniform_pvalues_1)


    df <- data.frame(pvalues=-log10(cat_pvalues + .000000000001), uniform_pvalues=-log10(cat_uniform_pvalues + .000000000001))
    
    # PLOT!
    max_val <-max(max(-log10(cat_uniform_pvalues + .000000000001)), max(-log10(cat_pvalues + .000000000001)))

    #PLOT!
    scatter <- ggplot(df, aes(x = uniform_pvalues, y = pvalues)) + geom_point() 
    scatter <- scatter + theme(text = element_text(size=14), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black")) 
    scatter <- scatter + labs(x = "Uniform", y = "Real")
    scatter <- scatter + geom_abline() 
    scatter <- scatter + scale_x_continuous(limits = c(-.1, max_val + .1), breaks = round(seq(0, max_val, by = 5),1)) + scale_y_continuous(limits = c(-.1,max_val+.1), breaks = round(seq(0, max_val, by = 5),1))

    ggsave(scatter, file=output_file, width=20, height=10.5, units="cm")
}



qq_plot_real_vs_permuted <- function(pvalues, perm_pvalue_arr, output_file) {
    real_pvalues <- sort(pvalues)

    perm_pvalues <- sort(perm_pvalue_arr)


    df <- data.frame(pvalues=-log10(real_pvalues + .000000000001), uniform_pvalues=-log10(perm_pvalues + .000000000001))
    
    # PLOT!
    max_val <-max(max(-log10(real_pvalues + .000000000001)), max(-log10(perm_pvalues + .000000000001)))

    #PLOT!
    scatter <- ggplot(df, aes(x = uniform_pvalues, y = pvalues)) + geom_point() 
    scatter <- scatter + theme(text = element_text(size=14), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black")) 
    scatter <- scatter + labs(x = "Permuted", y = "Real")
    scatter <- scatter + geom_abline() 
    scatter <- scatter + scale_x_continuous(limits = c(-.1, max_val + .1), breaks = round(seq(0, max_val, by = 5),1)) + scale_y_continuous(limits = c(-.1,max_val+.1), breaks = round(seq(0, max_val, by = 5),1))

    ggsave(scatter, file=output_file, width=20, height=10.5, units="cm")
}






# Command Line Args
real_merged_results = args[1]
perm_merged_results = args[2]
parameter_string = args[3]
qtl_visualization_dir = args[4]

# Load in real results
real_data <- read.table(real_merged_results,header=FALSE, sep = "\t")
perm_data <- read.table(perm_merged_results, header=FALSE, sep = "\t")





#############################
# QQ-plot showing real and permuted compared to uniform
#############################

output_file <- paste0(qtl_visualization_dir, parameter_string, "_qq_plot_cmp_true_to_uniform.png")
qq_plot_cmp_true_to_uniform(real_data$V4, output_file)


output_file <- paste0(qtl_visualization_dir, parameter_string, "_qq_plot_cmp_to_uniform.png")
qq_plot_cmp_to_uniform(real_data$V4, perm_data$V4, output_file)

#############################
# QQ-plot showing real and permuted compared to uniform
#############################
output_file <- paste0(qtl_visualization_dir, parameter_string, "_qq_plot_real_vs_permuted.png")
qq_plot_real_vs_permuted(real_data$V4, perm_data$V4, output_file)