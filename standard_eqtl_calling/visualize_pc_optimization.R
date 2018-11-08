args = commandArgs(trailingOnly=TRUE)
library(ggplot2)
library(ggthemes)
library(reshape)
library(cowplot)



number_egenes_greater_than_thresh_boxplot <- function(directory, parameter_string, pvalue_thresh, output_file) {
    num_gene_vec <- c()
    num_pc_vec <- c()
    time_step_vec <- c()
    for (temp_time_step in 0:15) {
        for (temp_pc_num in 0:3) {
            print(temp_time_step)
            print(temp_pc_num)
            input_file_name <- paste0(directory,"cht_results_",parameter_string,"_num_pc_",temp_pc_num,"_time_",temp_time_step,"_eqtl_results.txt")
            data <- read.table(input_file_name,header=TRUE)
            indexes <- data$pvalue <= pvalue_thresh
            num_genes <- length(unique(data[indexes,]$gene_id))

            num_gene_vec <- c(num_gene_vec, num_genes)
            num_pc_vec <- c(num_pc_vec, temp_pc_num)
            time_step_vec <- c(time_step_vec, temp_time_step)
        }
    }
    df <- data.frame(num_genes = as.numeric(num_gene_vec), num_pcs <- factor(num_pc_vec), time_step = factor(time_step_vec))

    box_plot <- ggplot(df, aes(x=num_pcs, y=num_genes)) + geom_boxplot(width=.54)
    box_plot <- box_plot + theme(text = element_text(size=18), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black")) 
    box_plot <- box_plot + labs(x = "Number of PCs", y = "Number of genes")
    ggsave(box_plot, file=output_file,width = 20,height=10.5,units="cm")
}




directory = args[1]
parameter_string = args[2]



pvalue_thresh <- .0001
output_file <- paste0(directory,"optimize_pcs_",pvalue_thresh,"_boxplot.png")
number_egenes_greater_than_thresh_boxplot(directory, parameter_string, pvalue_thresh, output_file)