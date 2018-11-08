args = commandArgs(trailingOnly=TRUE)
library(ggplot2)
library(ggthemes)
library(cowplot)
library(reshape)



load_in_odds_ratios <- function(file_name, adding_constant) {
    aa <- read.table(file_name,header=TRUE)

    real_overlaps <- as.numeric(aa$real_overlaps) + adding_constant
    real_misses <- as.numeric(aa$real_misses) + adding_constant
    perm_overlaps <- as.numeric(aa$perm_overlaps) + adding_constant
    perm_misses <- as.numeric(aa$perm_misses) + adding_constant
    odds_ratios <- (real_overlaps/real_misses)/(perm_overlaps/perm_misses)
    return(odds_ratios)
}

all_available_enrichments_boxplot <- function(cre, cell_lines, hits_versions, adding_constant, num_permutations, input_root, output_file, marker_type) {
    odds_ratios <- c()
    roadmap_cell_types <- c()
    dynamic_qtl_versions <- c()
    for (cell_line_counter in 1:length(cell_lines)) {
        for (hits_version_counter in 1:length(hits_versions)) {
            cell_line <- cell_lines[cell_line_counter]
            hits_version <- hits_versions[hits_version_counter]
            input_file <- paste0(input_root, cre, "_", cell_line, "_cell_lines_", hits_version,"_hits_", num_permutations, "_",cluster_assignment, ".txt" )
            or <- load_in_odds_ratios(input_file, adding_constant)
            odds_ratios <- c(odds_ratios, or)
            roadmap_cell_types <- c(roadmap_cell_types, rep(cell_line, length(or)))
            dynamic_qtl_versions <- c(dynamic_qtl_versions, rep(hits_version, length(or)))
        }
    }
    df <- data.frame(odds_ratio=odds_ratios,roadmap_cell_type=factor(roadmap_cell_types,levels=c("all", "ipsc", "heart", "ipsc_only", "heart_only", "heart_and_ipsc")), qtl_version=factor(dynamic_qtl_versions))
    # PLOT
    boxplot <- ggplot(df, aes(x=roadmap_cell_type,y=odds_ratio,fill=qtl_version)) + geom_boxplot() + labs(x = "Roadmap Cell Type", y = "Odds Ratio", title= marker_type)
    boxplot <- boxplot + theme(text = element_text(size=18))
    boxplot <- boxplot + geom_hline(yintercept = 1.0) 
    boxplot <- boxplot 
    boxplot <- boxplot
    ggsave(boxplot, file=output_file,width = 26,height=10.5,units="cm")

}
enrichment_boxplot_cross_time <- function(marker_type, cell_lines, time_steps, cht_enrichment_dir, parameter_string, pc_num, num_permutations, efdr, output_file) {
    odds_ratios <- c()
    roadmap_cell_types <- c()
    time_step_vec <- c()
    adding_constant <- 0
    for (cell_line_counter in 1:length(cell_lines)) {
        for (time_step_counter in 1:length(time_steps)) {
            cell_line <- cell_lines[cell_line_counter]
            time_step <- time_steps[time_step_counter]
            input_file <- paste0(cht_enrichment_dir, parameter_string, "_num_pc_", pc_num, "_time_", time_step, "_efdr_thresh_",efdr,"_", cell_line, "_cell_lines_", marker_type, "_", num_permutations,"_enrich.txt")
            or <- load_in_odds_ratios(input_file, adding_constant)
            odds_ratios <- c(odds_ratios, or)
            roadmap_cell_types <- c(roadmap_cell_types, rep(cell_line, length(or)))
            time_step_vec <- c(time_step_vec, rep(time_step, length(or)))
        }
    }
    df <- data.frame(odds_ratio=odds_ratios,roadmap_cell_type=factor(roadmap_cell_types,levels=c("all", "ipsc", "heart", "ipsc_only", "heart_only", "heart_and_ipsc")), time_step=factor(time_step_vec))

    boxplot <- ggplot(df, aes(x=time_step,y=odds_ratio,fill=roadmap_cell_types)) + geom_boxplot() + labs(x = "Time", y = "Odds Ratio", title= marker_type)
    boxplot <- boxplot + theme(text = element_text(size=18))
    boxplot <- boxplot + geom_hline(yintercept = 1.0) 
    boxplot <- boxplot 
    boxplot <- boxplot
    ggsave(boxplot, file=output_file,width = 26,height=10.5,units="cm")

}


parameter_string = args[1]
cht_visualization_dir = args[2]
pc_num = args[3]
cht_enrichment_dir = args[4]
num_permutations = args[5]
efdr = args[6]


marker_type <- "promotor"
cell_lines <- c("all", "ipsc", "heart", "ipsc_only", "heart_only", "heart_and_ipsc")
time_steps <- 0:15
output_file <- paste0(cht_visualization_dir, parameter_string, "_num_pc_", pc_num, "_efdr_thresh_",efdr,"_", marker_type, "_", num_permutations, "_all_cell_line_options_enrich.png")
print(output_file)
enrichment_boxplot_cross_time(marker_type, cell_lines, time_steps, cht_enrichment_dir, parameter_string, pc_num, num_permutations, efdr, output_file)

marker_type <- "enhancer"
cell_lines <- c("all", "ipsc", "heart", "ipsc_only", "heart_only", "heart_and_ipsc")
time_steps <- 0:15
output_file <- paste0(cht_visualization_dir, parameter_string, "_num_pc_", pc_num, "_efdr_thresh_",efdr,"_", marker_type, "_", num_permutations, "_all_cell_line_options_enrich.png")
print(output_file)
enrichment_boxplot_cross_time(marker_type, cell_lines, time_steps, cht_enrichment_dir, parameter_string, pc_num, num_permutations, efdr, output_file)


marker_type <- "promotor"
cell_lines <- c("ipsc", "heart")
time_steps <- 0:15
output_file <- paste0(cht_visualization_dir, parameter_string, "_num_pc_", pc_num, "_efdr_thresh_",efdr,"_", marker_type, "_", num_permutations, "_cell_line_specific_enrich.png")
print(output_file)
enrichment_boxplot_cross_time(marker_type, cell_lines, time_steps, cht_enrichment_dir, parameter_string, pc_num, num_permutations, efdr, output_file)

marker_type <- "enhancer"
cell_lines <- c("ipsc", "heart")
time_steps <- 0:15
output_file <- paste0(cht_visualization_dir, parameter_string, "_num_pc_", pc_num, "_efdr_thresh_",efdr,"_", marker_type, "_", num_permutations, "_cell_line_specific_enrich.png")
print(output_file)
enrichment_boxplot_cross_time(marker_type, cell_lines, time_steps, cht_enrichment_dir, parameter_string, pc_num, num_permutations, efdr, output_file)


marker_type <- "promotor"
cell_lines <- c("ipsc_only", "heart_only")
time_steps <- 0:15
output_file <- paste0(cht_visualization_dir, parameter_string, "_num_pc_", pc_num, "_efdr_thresh_",efdr,"_", marker_type, "_", num_permutations, "_cell_line_only_enrich.png")
print(output_file)
enrichment_boxplot_cross_time(marker_type, cell_lines, time_steps, cht_enrichment_dir, parameter_string, pc_num, num_permutations, efdr, output_file)

marker_type <- "enhancer"
cell_lines <- c("ipsc_only", "heart_only")
time_steps <- 0:15
output_file <- paste0(cht_visualization_dir, parameter_string, "_num_pc_", pc_num, "_efdr_thresh_",efdr,"_", marker_type, "_", num_permutations, "_cell_line_only_enrich.png")
print(output_file)
enrichment_boxplot_cross_time(marker_type, cell_lines, time_steps, cht_enrichment_dir, parameter_string, pc_num, num_permutations, efdr, output_file)

# file_name <- paste0(cht_enrichment_dir, parameter_string, "_num_pc_", pc_num, "_time_", time_step, "_efdr_thresh_",efdr,"_", cell_line_name, "_", marker_type, "_", num_permutations,"_enrich.txt")

