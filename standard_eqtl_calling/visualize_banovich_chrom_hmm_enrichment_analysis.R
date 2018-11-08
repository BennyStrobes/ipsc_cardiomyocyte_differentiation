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


odds_ratio_cell_line_specific_boxplot <- function(ipsc_early_or, ipsc_late_or, cardio_early_or, cardio_late_or, ipsc_change_or, cardio_change_or, output_file, marker_type) {
    odds_ratios <- c()
    roadmap_cell_types <- c()
    dynamic_qtl_versions <- c()


    odds_ratios <- c(odds_ratios, ipsc_early_or)
    roadmap_cell_types <- c(roadmap_cell_types, rep("ipsc", length(ipsc_early_or)))
    dynamic_qtl_versions <- c(dynamic_qtl_versions, rep("early_qtl", length(ipsc_early_or)))


    odds_ratios <- c(odds_ratios, ipsc_late_or)
    roadmap_cell_types <- c(roadmap_cell_types, rep("ipsc", length(ipsc_late_or)))
    dynamic_qtl_versions <- c(dynamic_qtl_versions, rep("late_qtl", length(ipsc_late_or)))

    odds_ratios <- c(odds_ratios, cardio_early_or)
    roadmap_cell_types <- c(roadmap_cell_types, rep("heart", length(cardio_early_or)))
    dynamic_qtl_versions <- c(dynamic_qtl_versions, rep("early_qtl", length(cardio_early_or)))

    odds_ratios <- c(odds_ratios, cardio_late_or)
    roadmap_cell_types <- c(roadmap_cell_types, rep("heart", length(cardio_early_or)))
    dynamic_qtl_versions <- c(dynamic_qtl_versions, rep("late_qtl", length(cardio_late_or)))


    odds_ratios <- c(odds_ratios, ipsc_change_or)
    roadmap_cell_types <- c(roadmap_cell_types, rep("ipsc", length(ipsc_change_or)))
    dynamic_qtl_versions <- c(dynamic_qtl_versions, rep("change_qtl", length(ipsc_change_or)))

    odds_ratios <- c(odds_ratios, cardio_change_or)
    roadmap_cell_types <- c(roadmap_cell_types, rep("heart", length(cardio_change_or)))
    dynamic_qtl_versions <- c(dynamic_qtl_versions, rep("change_qtl", length(cardio_change_or)))

    df <- data.frame(odds_ratio=odds_ratios,roadmap_cell_type=factor(roadmap_cell_types,levels=c("ipsc","heart")), qtl_version=factor(dynamic_qtl_versions))
    # PLOT
    boxplot <- ggplot(df, aes(x=roadmap_cell_type,y=odds_ratio,fill=qtl_version)) + geom_boxplot() + labs(x = "Roadmap Cell Type", y = "Odds Ratio", title= marker_type)
    boxplot <- boxplot + theme(text = element_text(size=18))
    boxplot <- boxplot + geom_hline(yintercept = 1.0) 
    boxplot <- boxplot 
    boxplot <- boxplot
    ggsave(boxplot, file=output_file,width = 20,height=10.5,units="cm")


}


all_available_enrichments_boxplot <- function(cre, cell_lines, adding_constant, num_permutations, input_root, output_file, title) {
    odds_ratios <- c()
    roadmap_cell_types <- c()
    for (cell_line_counter in 1:length(cell_lines)) {
        cell_line <- cell_lines[cell_line_counter]
        input_file <- paste0(input_root, cell_line, "_cell_lines_", cre, "_", num_permutations, "_enrich.txt")
        or <- load_in_odds_ratios(input_file, adding_constant)
        odds_ratios <- c(odds_ratios, or)
        roadmap_cell_types <- c(roadmap_cell_types, rep(cell_line, length(or)))
    }
    df <- data.frame(odds_ratio=odds_ratios,roadmap_cell_type=factor(roadmap_cell_types,levels=c("all", "ipsc", "heart", "ipsc_only", "heart_only", "heart_and_ipsc")))
    # PLOT
    boxplot <- ggplot(df, aes(x=roadmap_cell_type,y=odds_ratio)) + geom_boxplot() + labs(x = "Roadmap Cell Type", y = "Odds Ratio", title= title)
    boxplot <- boxplot + theme(text = element_text(size=18))
    boxplot <- boxplot + geom_hline(yintercept = 1.0) 
    boxplot <- boxplot 
    boxplot <- boxplot
    ggsave(boxplot, file=output_file,width = 26,height=10.5,units="cm")

}

all_hits_enrichments_boxplot <- function(cre, cell_lines, hits_versions, adding_constant, num_permutations, input_root, output_file, marker_type) {
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
    boxplot <- boxplot + theme(legend.position="none")
    boxplot <- boxplot
    ggsave(boxplot, file=output_file,width = 26,height=10.5,units="cm")

}

only_all_enrichment_boxplot <- function(cell_line, adding_constant, num_permutations, input_root, output_file, title) {
    odds_ratios <- c()
    cre_version <- c()

    cre <- "promotor"    
    input_file <- paste0(input_root, cell_line, "_cell_lines_", cre, "_", num_permutations, "_enrich.txt")
    or <- load_in_odds_ratios(input_file, adding_constant)
    odds_ratios <- c(odds_ratios, or)
    cre_version <- c(cre_version, rep(cre, length(or)))

    cre <- "enhancer"
    input_file <- paste0(input_root, cell_line, "_cell_lines_", cre, "_", num_permutations, "_enrich.txt")
    or <- load_in_odds_ratios(input_file, adding_constant)
    odds_ratios <- c(odds_ratios, or)
    cre_version <- c(cre_version, rep(cre, length(or)))


    df <- data.frame(odds_ratio=odds_ratios, cre_type=factor(cre_version))
    # PLOT
    boxplot <- ggplot(df, aes(x=cre_type,y=odds_ratio,fill=cre_type)) + geom_boxplot() + labs(x = "CRE Type", y = "Odds Ratio",title=title)
    boxplot <- boxplot + theme(text = element_text(size=18))
    boxplot <- boxplot + geom_hline(yintercept = 1.0) 
    boxplot <- boxplot + theme(legend.position="none")
    boxplot <- boxplot
    ggsave(boxplot, file=output_file,width = 15,height=10.5,units="cm") 
}




input_directory <- args[1]
visualization_directory <- args[2]
num_permutations <- args[3]

cell_line <- "all"
adding_constant <- 0
output_file <- paste0(visualization_directory, "banovich_ipsc_prom_enh_all_enrichments_num_perm_", num_permutations, "_odds_ratios.png")
only_all_enrichment_boxplot(cell_line, adding_constant, num_permutations, paste0(input_directory, "banovich_ipsc_"), output_file, "banovich_ipsc")

cell_line <- "all"
adding_constant <- 0
output_file <- paste0(visualization_directory, "banovich_cm_prom_enh_all_enrichments_num_perm_", num_permutations, "_odds_ratios.png")
only_all_enrichment_boxplot(cell_line, adding_constant, num_permutations, paste0(input_directory, "banovich_cm_"), output_file, "banovich_cm")




cell_lines <- c("all", "ipsc", "heart", "ipsc_only", "heart_only", "heart_and_ipsc")
adding_constant <- 1


########################
# Promoter
########################
cre <- "promotor"
output_file <- paste0(visualization_directory, "banovich_ipsc_",cre,"_num_perm_",num_permutations,"_odds_ratios.png")
all_available_enrichments_boxplot(cre, cell_lines, adding_constant, num_permutations, paste0(input_directory, "banovich_ipsc_"), output_file, "banovich_ipsc promotor")

########################
# Enhancer
########################
cre <- "enhancer"
output_file <- paste0(visualization_directory, "banovich_ipsc_",cre,"_num_perm_",num_permutations,"_odds_ratios.png")
all_available_enrichments_boxplot(cre, cell_lines, adding_constant, num_permutations, paste0(input_directory, "banovich_ipsc_"), output_file, "banovich_ipsc enhancer")

########################
# Promoter
########################
cre <- "promotor"
output_file <- paste0(visualization_directory, "banovich_cm_",cre,"_num_perm_",num_permutations,"_odds_ratios.png")
all_available_enrichments_boxplot(cre, cell_lines, adding_constant, num_permutations, paste0(input_directory, "banovich_cm_"), output_file, "banovich_cm promotor")

########################
# Enhancer
########################
cre <- "enhancer"
output_file <- paste0(visualization_directory, "banovich_cm_",cre,"_num_perm_",num_permutations,"_odds_ratios.png")
all_available_enrichments_boxplot(cre, cell_lines, adding_constant, num_permutations, paste0(input_directory, "banovich_cm_"), output_file, "banovich_cm enhancer")















