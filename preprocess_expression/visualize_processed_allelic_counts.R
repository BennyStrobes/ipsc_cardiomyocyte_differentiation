args = commandArgs(trailingOnly=TRUE)
library(Rsubread)
library(ggplot2)
library(ggthemes)
library(reshape)
library(mvtnorm)
library(cowplot)


# Visualize the percent of het snps that show biallelic expression. Color points by cell line
percent_biallelic_het_snps_scatter_cell_line <- function(ref_counts, total_counts, sample_info, percent_biallelic_ouptut_file, num_read_threshold) {
    N <- dim(sample_info)[1]  # Num samples
    #  Add new column to sample_info to keep track of the percent bi-allelic in each sample
    sample_info$percent_biallelic <- numeric(N)

    # Loop through each sample (compute percent biallleic in each)
    for (n in 1:N) {
        n_ref <- ref_counts[,n]  # reference allele counts for the nth sample
        n_total <- total_counts[,n]  # total counts for the nth sample
        # Consider only heterozygous sites (the nan check) and sites that have greater than num_read_threshold reads mapped totally
        n_observed_indices <- !is.nan(n_total) & (n_total > num_read_threshold)

        # Compute whether each sites is biallelic in terms of reads mapped
        biallelic_sites <- (n_ref[n_observed_indices] != n_total[n_observed_indices]) & (n_ref[n_observed_indices] != 0)

        # Compute percent of het snps that show biallelic expression
        sample_info$percent_biallelic[n] <- sum(biallelic_sites)/length(biallelic_sites)
    }

    # Now reorder percent_biallelic so they are ordered by cell line, and also time
    ordered_biallelic <- c()
    ordered_cell_lines <- c()
    unique_cell_lines <- unique(sample_info$cell_line)
    # Loop through each cell line
    for (i in 1:length(unique_cell_lines)) {
        i_cell_line <- unique_cell_lines[i]  # current cell line
        i_sample_info <- sample_info[sample_info$cell_line == i_cell_line,]  # Subset sample_info to only include samples belonging to this cell line

        i_times <- sort(i_sample_info$time)  # Now sort times within this cell line in numerical order
        #  Loop through sorted times
        for (j in 1:length(i_times)) {
            ij_time <- i_times[j]  # current time
            # calculate percent biallelic in this cell line at this time step
            ordered_biallelic <- c(ordered_biallelic, i_sample_info[i_sample_info$time == ij_time,]$percent_biallelic) 
            ordered_cell_lines <- c(ordered_cell_lines, i_cell_line)
        }

    }


    # Put data into data.frame for plotting
    df <- data.frame(sample_num = as.numeric(rownames(sample_info)), percent_biallelic = ordered_biallelic, cell_line = factor(ordered_cell_lines))

    #PLOT!
    scatter <- ggplot(df, aes(x = sample_num, y = percent_biallelic, colour = cell_line)) + geom_point() 
    scatter <- scatter + theme(text = element_text(size=18), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black")) 
    scatter <- scatter  + labs(x = "Sample", y = "Percent Biallelic", colour = "Cell Line")
    ggsave(scatter, file=percent_biallelic_ouptut_file,width = 15,height=10.5,units="cm")
}



# Visualize the percent of het snps that show biallelic expression. Color points by cell line
percent_biallelic_het_snps_scatter_library_size <- function(ref_counts, total_counts, sample_info, percent_biallelic_ouptut_file, num_read_threshold) {
    N <- dim(sample_info)[1]  # Num samples
    #  Add new column to sample_info to keep track of the percent bi-allelic in each sample
    sample_info$percent_biallelic <- numeric(N)

    # Loop through each sample (compute percent biallleic in each)
    for (n in 1:N) {
        n_ref <- ref_counts[,n]  # reference allele counts for the nth sample
        n_total <- total_counts[,n]  # total counts for the nth sample
        # Consider only heterozygous sites (the nan check) and sites that have greater than num_read_threshold reads mapped totally
        n_observed_indices <- !is.nan(n_total) & (n_total > num_read_threshold)

        # Compute whether each sites is biallelic in terms of reads mapped
        biallelic_sites <- (n_ref[n_observed_indices] != n_total[n_observed_indices]) & (n_ref[n_observed_indices] != 0)

        # Compute percent of het snps that show biallelic expression
        sample_info$percent_biallelic[n] <- sum(biallelic_sites)/length(biallelic_sites)
    }

    # Now reorder percent_biallelic so they are ordered by cell line, and also time
    ordered_biallelic <- c()
    ordered_lib_sizes <- c()
    unique_cell_lines <- unique(sample_info$cell_line)
    # Loop through each cell line
    for (i in 1:length(unique_cell_lines)) {
        i_cell_line <- unique_cell_lines[i]  # current cell line
        i_sample_info <- sample_info[sample_info$cell_line == i_cell_line,]  # Subset sample_info to only include samples belonging to this cell line

        i_times <- sort(i_sample_info$time)  # Now sort times within this cell line in numerical order
        #  Loop through sorted times
        for (j in 1:length(i_times)) {
            ij_time <- i_times[j]  # current time
            # calculate percent biallelic in this cell line at this time step
            ordered_biallelic <- c(ordered_biallelic, i_sample_info[i_sample_info$time == ij_time,]$percent_biallelic) 
            ordered_lib_sizes <- c(ordered_lib_sizes, i_sample_info[i_sample_info$time == ij_time,]$lib.size)
        }

    }


    # Put data into data.frame for plotting
    df <- data.frame(sample_num = as.numeric(rownames(sample_info)), percent_biallelic = ordered_biallelic, library_size = ordered_lib_sizes)
    #PLOT!
    scatter <-  ggplot(df,aes(sample_num,percent_biallelic)) + geom_point(aes(colour=library_size)) + theme(text = element_text(size=18), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black")) 
    scatter <- scatter + scale_color_gradient(low="pink",high="blue") 
    scatter <- scatter + labs(x = "Sample", y = "Percent Biallelic", colour = "Lib. Size")
    ggsave(scatter, file=percent_biallelic_ouptut_file,width = 15,height=10.5,units="cm")

}

# Compute fraction of imputated genotype sites for each sample
fraction_of_hard_coded_genotype_sites <- function(het_prob_file, sample_info, output_file) {
    # Stream het prob file (each line is a site)
    stop = FALSE
    f = file(het_prob_file, "r")
    total_sites <- 0
    while(!stop) {
        next_line = readLines(f, n = 1)
        if(length(next_line) == 0) {  # Stopping criteria
            stop = TRUE
            close(f)
            break
        }
        data <- unlist(strsplit(next_line,"\t"))
        if (startsWith(data[1],'#CHROM')) {  # HEADER
            cell_lines <- data[10:length(data)]  # Keep track of names of cell lines
            hard_coded <- numeric(length(cell_lines))  # Keep track of how many of the genotype calls are hard-coded
        }
        if (startsWith(data[1],'#') == FALSE) {  # Standard line
            total_sites <- total_sites + 1 # Keep track of total sites
            probs <- as.numeric(data[10:length(data)])
            hard_coded_1 <- 1.0*(probs==1.0)  # Binary variable whether each cell line is hard coded 1 for this genotype
            hard_coded_0 <- 1.0*(probs==0.0)  # Binary variable whether each cell line is hard coded 0 for this genotype
            hard_coded <- hard_coded + hard_coded_1 + hard_coded_0  # Keep track
        }
    }

    df <- data.frame(cell_line = factor(cell_lines), fraction_hard_coded = hard_coded/total_sites)

    bar_plot <- ggplot(df, aes(x=cell_line, y=fraction_hard_coded, fill=cell_line)) + geom_bar(stat="identity")

    bar_plot <- bar_plot + theme(text = element_text(size=18), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"),axis.text.x = element_text(angle = 90, vjust=.43,hjust = -1)) 
    bar_plot <- bar_plot + labs(colour="Cell Line",x = "Cell Line", y = "% hard-coded sites") + theme(legend.position="none")

    #pdf(output_file)
    #print(bar_plot)
    #dev.off()
    ggsave(bar_plot, file=output_file,width = 15,height=10.5,units="cm")
}

# Make boxplot of number of expressed het-snps per individual with one box for every read cutoff (that defines what is an expressed het-snp)
number_of_expressed_het_snps_per_individual <- function(total_counts, output_filer, het_thresh) {
    # Keep track of variables (initialization)
    num_expressed_sites <- c()
    thresholds <- c()
    cell_line <- c()
    
    #Get order of cel_line ids
    ordered_cell_lines <- substr(colnames(total_counts),2,6)

    # Compute number of expressed het_snps for each of these thresholds
    read_threshs <- c(2,6,10,15,30)
    for (iter in 1:length(read_threshs)) {
        read_thresh <- read_threshs[iter]
        # Compute which (site,samples) have expression
        binary_matrix <- total_counts > read_thresh
        # Create array of length number of samples that is the number of expressed sites in each sample
        num_expressed <- colSums(binary_matrix,na.rm=TRUE)

        # keep track of data
        num_expressed_sites <- c(num_expressed_sites, num_expressed)
        cell_line <- c(cell_line, ordered_cell_lines)
        thresholds <- c(thresholds, numeric(length(num_expressed)) + read_thresh)
    }
    # PLOT
    df <- data.frame(cell_line = factor(cell_line), num_expressed_sites = num_expressed_sites, thresholds=factor(thresholds))
    box_plot <- ggplot(df, aes(x=thresholds, y=num_expressed_sites)) + geom_boxplot(notch=TRUE)
    box_plot <- box_plot + theme(text = element_text(size=18), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black")) 
    box_plot <- box_plot + labs(x = "Read Threshold", y = "Expressed het sites / Sample", title=paste0("Het threshold = .",het_thresh))
    ggsave(box_plot, file=output_filer,width = 15,height=10.5,units="cm")
}

# Make boxplot of number of expressed het-snps per individual with one box for every read cutoff (that defines what is an expressed het-snp) and also per cell_line
number_of_expressed_het_snps_per_individual_cell_line_binned <- function(total_counts, output_filer, het_thresh) {
    # Keep track of variables (initialization)
    num_expressed_sites <- c()
    thresholds <- c()
    cell_line <- c()
    
    #Get order of cel_line ids
    ordered_cell_lines <- substr(colnames(total_counts),2,6)

    # Compute number of expressed het_snps for each of these thresholds
    read_threshs <- c(2,6,10)
    for (iter in 1:length(read_threshs)) {
        read_thresh <- read_threshs[iter]
        # Compute which (site,samples) have expression
        binary_matrix <- total_counts > read_thresh
        # Create array of length number of samples that is the number of expressed sites in each sample
        num_expressed <- colSums(binary_matrix,na.rm=TRUE)

        # keep track of data
        num_expressed_sites <- c(num_expressed_sites, num_expressed)
        cell_line <- c(cell_line, ordered_cell_lines)
        thresholds <- c(thresholds, numeric(length(num_expressed)) + read_thresh)
    }
    # PLOT
    df <- data.frame(cell_line = factor(cell_line), num_expressed_sites = num_expressed_sites, thresholds=factor(thresholds))
    box_plot <- ggplot(df, aes(x=thresholds, y=num_expressed_sites, fill=cell_line)) + geom_boxplot()
    box_plot <- box_plot + theme(text = element_text(size=18), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black")) 
    box_plot <- box_plot + labs(fill= "Cell Line",x = "Read Threshold", y = "Expressed het sites / Sample", title=paste0("Het threshold = .",het_thresh))
    ggsave(box_plot, file=output_filer,width = 15,height=10.5,units="cm")
}


compute_number_of_sites_that_pass_filter <- function(ref_counts,total_counts, sample_info, number_of_heterozygous_lines_i, fraction_of_biallelic_samples_i, min_reads) {
    ordered_cell_lines <- factor(sample_info$cell_line)
    # Keep track of number of sites that pass
    num_sites_that_pass <- 0



    aa <- rowSums(total_counts >= min_reads & ref_counts/total_counts >= .01 & (total_counts-ref_counts)/total_counts >= .01,na.rm=TRUE)
    bb <- rowSums(!is.na(total_counts))
    biallelic_fractions <- aa/bb

    num_cell_lines_with_het_sites <- ceiling(bb/16)

    num_sites_that_pass <- sum(biallelic_fractions >= fraction_of_biallelic_samples_i & num_cell_lines_with_het_sites >= number_of_heterozygous_lines_i)

    return(num_sites_that_pass)
}

compute_number_of_sites_that_pass_independent_filter_helper <- function(ref_counts,total_counts, sample_info, number_of_heterozygous_lines_i, fraction_of_biallelic_samples_i, min_reads) {
    aa <- rowSums(total_counts >= min_reads & ref_counts/total_counts >= .01 & (total_counts-ref_counts)/total_counts >= .01,na.rm=TRUE)
    bb <- rowSums(!is.na(total_counts))
    biallelic_fractions <- aa/bb

    num_cell_lines_with_het_sites <- bb


    time_step_truth_vec <- biallelic_fractions >= fraction_of_biallelic_samples_i & num_cell_lines_with_het_sites >= number_of_heterozygous_lines_i
    return(time_step_truth_vec)
}

compute_number_of_sites_that_pass_independent_filter <- function(ref_counts,total_counts, sample_info, number_of_heterozygous_lines_i, fraction_of_biallelic_samples_i, min_reads) {
    ordered_cell_lines <- factor(sample_info$cell_line)
    time_steps <- sample_info$time
    # Keep track of number of sites that pass
    num_sites_that_pass <- 0
    truth_vector <- rep(1.0,dim(ref_counts)[1])

    for (time_step in 0:15) {

        observed_samples <- time_steps == time_step
        time_step_ref_counts <- ref_counts[,observed_samples]
        time_step_total_counts <- total_counts[,observed_samples]
        time_step_truth_vec <- compute_number_of_sites_that_pass_independent_filter_helper(time_step_ref_counts, time_step_total_counts, sample_info, number_of_heterozygous_lines_i, fraction_of_biallelic_samples_i, min_reads)
        truth_vector <- truth_vector*time_step_truth_vec

    }
    num_sites_that_pass <- sum(truth_vector)
    return(num_sites_that_pass)
}

compute_number_of_unique_genes_that_pass_filter <- function(ref_counts,total_counts, sample_info, number_of_heterozygous_lines_i, fraction_of_biallelic_samples_i, min_reads, gene_names) {
    ordered_cell_lines <- factor(sample_info$cell_line)
    # Keep track of number of sites that pass
    num_sites_that_pass <- 0



    aa <- rowSums(total_counts >= min_reads & ref_counts/total_counts >= .01 & (total_counts-ref_counts)/total_counts >= .01,na.rm=TRUE)
    bb <- rowSums(!is.na(total_counts))
    biallelic_fractions <- aa/bb

    num_cell_lines_with_het_sites <- ceiling(bb/16)

    sites_that_pass <- biallelic_fractions >= fraction_of_biallelic_samples_i & num_cell_lines_with_het_sites >= number_of_heterozygous_lines_i


    # Simple error checking
    if (length(sites_that_pass) != length(gene_names)) {
        print('ERROR!!!')
    }

    # Initialize vector of genes to keep track of those that pass
    genes_that_pass <- c()

    # Loop through sites
    for (site_index in 1:length(sites_that_pass)) {
        # Check to see if site passes filter
        if (sites_that_pass[site_index] == TRUE) {
            # Full gene string is name of gene(s) at site. If more than 1 gene, then they are comma seperated
            full_gene_string <- gene_names[site_index]
            # Add each element of the ',' seperated array to  genes_that_pass
            full_gene_info <- strsplit(full_gene_string,",")[[1]]
            for (i in 1:length(full_gene_info)) {
                genes_that_pass <- c(genes_that_pass, full_gene_info[i])
            }

        }
    }
    return(length(unique(genes_that_pass)))
}

# Experiment with number of heterozygous sites that remain when we use various filters.
# Specifically wish to vary:
### a. Minimum number of cell lines heterozygous
### b. Minimum fraction of remaining samples that have bi-allelic expression (fewer than 3 reads mapping, or less than 1% of reads mapping to one allele)
number_of_heterozygous_sites_at_various_filters_lineplot <- function(ref_counts, total_counts, sample_info, output_file, het_thresh) {
    # Keep track of quantities of interest
    number_of_sites <- c()
    number_of_heterozygous_lines <- c()
    fraction_of_biallelic_samples <- c()
    min_num_reads <- c()

    #num_het_lines <- c(6,7,8,9,10)
    #fraction_of_samples <- c(.6,.7,.8,.9)
    num_het_lines <- c(4,5,6,7,8,9)
    fraction_of_samples <- c(.5,.6,.7,.8,.9)
    min_reads <- c(2,3)
    for (i in 1:length(num_het_lines)) {
        for (j in 1:length(fraction_of_samples)) {
            for (k in 1:length(min_reads)) {
                # Threshold for min number of heterozygous cell lines
                number_of_heterozygous_lines_i <- num_het_lines[i]
                # Threshold for min fraction of heterozygous samples that show biallelic expression
                fraction_of_biallelic_samples_i <- fraction_of_samples[j]
                # Threshold for min reads
                min_reads_i <- min_reads[k]

                # Compute the number of sites that pass this threshold
                num_sites_that_pass <- compute_number_of_sites_that_pass_filter(ref_counts,total_counts, sample_info, number_of_heterozygous_lines_i, fraction_of_biallelic_samples_i, min_reads_i)
            
                # Store results
                number_of_sites <- c(number_of_sites,num_sites_that_pass)
                number_of_heterozygous_lines <- c(number_of_heterozygous_lines, number_of_heterozygous_lines_i)
                fraction_of_biallelic_samples <- c(fraction_of_biallelic_samples, fraction_of_biallelic_samples_i)
                min_num_reads <- c(min_num_reads,min_reads_i)
            }
        }
    }


    df <- data.frame(number_of_sites = number_of_sites,number_of_reads = factor(min_num_reads), number_of_heterozygous_lines = factor(number_of_heterozygous_lines), fraction_of_biallelic_samples = factor(fraction_of_biallelic_samples))

    #PLOT
    line_plot <- ggplot(df, aes(x=number_of_heterozygous_lines, y=number_of_sites, colour=fraction_of_biallelic_samples, shape=number_of_reads, group=interaction(fraction_of_biallelic_samples,number_of_reads))) + geom_line() +
                geom_point(aes(color=fraction_of_biallelic_samples)) +
                theme(text = element_text(size=16), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black")) +
                labs(colour="Min % biallelic", shape = "Min reads",x = "Min # of het cell lines", y = "Number of sites") + ylim(0,14000)
    ggsave(line_plot, file=output_file,width = 15,height=10.5,units="cm")
}

# Experiment with number of heterozygous sites that remain when we use various filters.
# Specifically wish to vary:
### a. Minimum number of cell lines heterozygous
### b. Minimum fraction of remaining samples that have bi-allelic expression (fewer than 3 reads mapping, or less than 1% of reads mapping to one allele)
number_of_heterozygous_sites_at_various_filters_independent_time_step_lineplot <- function(ref_counts, total_counts, sample_info, output_file, het_thresh) {
    # Keep track of quantities of interest
    number_of_sites <- c()
    number_of_heterozygous_lines <- c()
    fraction_of_biallelic_samples <- c()
    min_num_reads <- c()

    #num_het_lines <- c(6,7,8,9,10)
    #fraction_of_samples <- c(.6,.7,.8,.9)
    num_het_lines <- c(4,5,6,7,8,9)
    fraction_of_samples <- c(.5,.6,.7,.8,.9)
    min_reads <- c(2,3)
    for (i in 1:length(num_het_lines)) {
        for (j in 1:length(fraction_of_samples)) {
            for (k in 1:length(min_reads)) {
                # Threshold for min number of heterozygous cell lines
                number_of_heterozygous_lines_i <- num_het_lines[i]
                # Threshold for min fraction of heterozygous samples that show biallelic expression
                fraction_of_biallelic_samples_i <- fraction_of_samples[j]
                # Threshold for min reads
                min_reads_i <- min_reads[k]

                # Compute the number of sites that pass this threshold
                num_sites_that_pass <- compute_number_of_sites_that_pass_independent_filter(ref_counts,total_counts, sample_info, number_of_heterozygous_lines_i, fraction_of_biallelic_samples_i, min_reads_i)
            
                # Store results
                number_of_sites <- c(number_of_sites,num_sites_that_pass)
                number_of_heterozygous_lines <- c(number_of_heterozygous_lines, number_of_heterozygous_lines_i)
                fraction_of_biallelic_samples <- c(fraction_of_biallelic_samples, fraction_of_biallelic_samples_i)
                min_num_reads <- c(min_num_reads,min_reads_i)
            }
        }
    }


    df <- data.frame(number_of_sites = number_of_sites,number_of_reads = factor(min_num_reads), number_of_heterozygous_lines = factor(number_of_heterozygous_lines), fraction_of_biallelic_samples = factor(fraction_of_biallelic_samples))

    #PLOT
    line_plot <- ggplot(df, aes(x=number_of_heterozygous_lines, y=number_of_sites, colour=fraction_of_biallelic_samples, shape=number_of_reads, group=interaction(fraction_of_biallelic_samples,number_of_reads))) + geom_line() +
                geom_point(aes(color=fraction_of_biallelic_samples)) +
                theme(text = element_text(size=16), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black")) +
                labs(colour="Min % biallelic", shape = "Min reads",x = "Min # of het cell lines", y = "Number of sites") + ylim(0,4000)
    ggsave(line_plot, file=output_file,width = 15,height=10.5,units="cm")
}



# Experiment with number of heterozygous sites (IN TERMS OF NUMBER OF UNIQUE GENES) that remain when we use various filters.
# Specifically wish to vary:
### a. Minimum number of cell lines heterozygous
### b. Minimum fraction of remaining samples that have bi-allelic expression (fewer than 3 reads mapping, or less than 1% of reads mapping to one allele)
number_of_genes_at_various_filters_lineplot <- function(ref_counts, total_counts, sample_info, output_file, het_thresh, gene_names) {
    # Keep track of quantities of interest
    number_of_sites <- c()
    number_of_heterozygous_lines <- c()
    fraction_of_biallelic_samples <- c()
    min_num_reads <- c()

    #num_het_lines <- c(6,7,8,9,10)
    #fraction_of_samples <- c(.6,.7,.8,.9)
    num_het_lines <- c(4,5,6,7,8,9)
    fraction_of_samples <- c(.5,.6,.7,.8,.9)
    min_reads <- c(2,3)
    for (i in 1:length(num_het_lines)) {
        for (j in 1:length(fraction_of_samples)) {
            for (k in 1:length(min_reads)) {
                # Threshold for min number of heterozygous cell lines
                number_of_heterozygous_lines_i <- num_het_lines[i]
                # Threshold for min fraction of heterozygous samples that show biallelic expression
                fraction_of_biallelic_samples_i <- fraction_of_samples[j]
                # Threshold for min reads
                min_reads_i <- min_reads[k]

                # Compute the number of unique genes that pass this threshold
                num_unique_genes <- compute_number_of_unique_genes_that_pass_filter(ref_counts,total_counts, sample_info, number_of_heterozygous_lines_i, fraction_of_biallelic_samples_i, min_reads_i, gene_names)
                # Store results
                number_of_sites <- c(number_of_sites,num_unique_genes)
                number_of_heterozygous_lines <- c(number_of_heterozygous_lines, number_of_heterozygous_lines_i)
                fraction_of_biallelic_samples <- c(fraction_of_biallelic_samples, fraction_of_biallelic_samples_i)
                min_num_reads <- c(min_num_reads,min_reads_i)
            }
        }
    }


    df <- data.frame(number_of_genes = number_of_sites,number_of_reads = factor(min_num_reads), number_of_heterozygous_lines = factor(number_of_heterozygous_lines), fraction_of_biallelic_samples = factor(fraction_of_biallelic_samples))

    #PLOT
    line_plot <- ggplot(df, aes(x=number_of_heterozygous_lines, y=number_of_genes, colour=fraction_of_biallelic_samples, shape=number_of_reads, group=interaction(fraction_of_biallelic_samples,number_of_reads))) + geom_line() +
                geom_point(aes(color=fraction_of_biallelic_samples)) +
                theme(text = element_text(size=16), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black")) +
                labs(colour="Min % biallelic", shape = "Min reads",x = "Min # of het cell lines", y = "Number of unique genes")
    ggsave(line_plot, file=output_file,width = 15,height=10.5,units="cm")
}


extract_gene_names_from_site_ids <- function(site_ids) {
    num_sites <- length(site_ids)
    # Initialize output vector
    gene_names <- rep("init",num_sites)
    # split up each element of site_ids with the delimiter "_" (the genes will contain distinct elements now)
    split_names <- strsplit(site_ids,"_")

    # Loop through each site seperately
    for (index in 1:num_sites) {
        # Find out which indices contain gene names (potentially more than 1 index)
        gene_indices <- grep("ENSG",split_names[[index]])

        # Append gene string to gene_names vector
        gene_names[index] <- paste(split_names[[index]][gene_indices], collapse=",")

        # Simple error checking
        if (length(gene_indices) == 0) {
            print("ERRROR in extracting gene names")
        }
    }
    return(gene_names)
}

# Histogram of number of mapped genes per heterozygous site
number_of_mapped_genes_per_site_histogram <- function(gene_names, output_file) {
    num_sites <- length(gene_names)
    # Initialize vector to keep track of number of genes per site
    num_genes_per_site <- c()
    
    # Looop through ach of sites
    for (index in 1:num_sites) {
        # Gene name(s). If more than one gene, they are seperated by ','
        full_gene_string <- gene_names[index]
        # Convert ',' seperated string to array
        full_gene_info <- strsplit(full_gene_string,",")[[1]]
        # Compute number of genes per site and add to num_genes_per_site array
        num_genes_per_site <- c(num_genes_per_site, length(full_gene_info))
    }

    df <- data.frame(num_genes_per_site = num_genes_per_site)

    # PLOT
    p <- ggplot(df, aes(x=num_genes_per_site)) + geom_histogram() +
        theme(text = element_text(size=16), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black")) +
        labs(x = "Number of genes per site", y = "Number of sites")
    ggsave(p, file=output_file,width = 15,height=10.5,units="cm")
}





###########################################################################
# Command line arguments
###########################################################################
processed_allelic_counts_dir = args[1]  # Contains allelic counts
visualize_allelic_counts_dir = args[2]  # Directory to save output files to
genotype_dir = args[3]  # genotype file
preprocess_total_expression_dir = args[4]  # contains sample_info file


het_thresh <- 999 #  Heterozygous probability threshold




###########################################################################
# Load in data
###########################################################################

#  Get sample information 
sample_info_file <- paste0(preprocess_total_expression_dir, "sample_info.txt")
sample_info <- read.table(sample_info_file, header=TRUE)


# total_counts_file contains total number of reads mapping to each site in each sample
# ref_counts_file contains number of reads mapping to reference allele at each site in each sample
total_counts_file <- paste0(processed_allelic_counts_dir, "allelic_counts_gene_mapped_het_prob_", het_thresh, "_total_counts.txt")
ref_counts_file <- paste0(processed_allelic_counts_dir, "allelic_counts_gene_mapped_het_prob_", het_thresh, "_ref_counts.txt")
total_counts <- read.table(total_counts_file, header=TRUE, sep="\t")
total_counts <- total_counts[,2:dim(total_counts)[2]]  # Skip first column due to transition from python (no information lost)
ref_counts <- read.table(ref_counts_file, header=TRUE, sep="\t")
ref_counts <- ref_counts[,2:dim(ref_counts)[2]] # Skip first column due to transition from python (no information lost)


#  Get the gene corresponding to each heterozygous site (row label). Note there can be multiple genes per site_id
# These can be found in:
full_allelic_counts_file <- paste0(processed_allelic_counts_dir, "allelic_counts_gene_mapped_het_prob_", het_thresh,".txt")
full_allelic_count_matrix <- read.table(full_allelic_counts_file, header=TRUE, sep="\t")
# Get row labels (site_ids) of ref_counts and total_counts matrices
site_ids <- as.character(full_allelic_count_matrix$siteID)
gene_names <- extract_gene_names_from_site_ids(site_ids)





###########################################################################
# Plotting functions
###########################################################################

# Compute fraction of imputated genotype sites for each sample
het_prob_file <- paste0(genotype_dir, "YRI_het_prob_genotype.vcf")
fraction_hard_coded_output_file <- paste0(visualize_allelic_counts_dir, "percent_hard_coded_genotypes_by_cell_line.pdf")
fraction_of_hard_coded_genotype_sites(het_prob_file, sample_info, fraction_hard_coded_output_file)


# Visualize the percent of het snps that show biallelic expression. Color points by cell line
num_read_threshold <- 5  # Only consider sites that have at least num_read_threshold reads mapping to both alleles
percent_biallelic_ouptut_file <- paste0(visualize_allelic_counts_dir, "percent_biallelic_het_snps_by_cell_line_",het_thresh, ".png")
percent_biallelic_het_snps_scatter_cell_line(ref_counts, total_counts, sample_info, percent_biallelic_ouptut_file, num_read_threshold)

# Visualize the percent of het snps that show biallelic expression. Color points by total read-depth
num_read_threshold <- 5  # Only consider sites that have at least num_read_threshold reads mapping to both alleles
percent_biallelic_ouptut_file <- paste0(visualize_allelic_counts_dir, "percent_biallelic_het_snps_by_library_size_",het_thresh, ".png")
percent_biallelic_het_snps_scatter_library_size(ref_counts, total_counts, sample_info, percent_biallelic_ouptut_file, num_read_threshold)



# Make boxplot of number of expressed het-snps per individual with one box for every read cutoff (that defines what is an expressed het-snp)
number_of_expressed_het_snps_output_file <- paste0(visualize_allelic_counts_dir, "number_of_expressed_het_snps_per_individual_boxplot_",het_thresh,".png")
number_of_expressed_het_snps_per_individual(total_counts, number_of_expressed_het_snps_output_file, het_thresh)

# Make boxplot of number of expressed het-snps per individual with one box for every read cutoff (that defines what is an expressed het-snp) and also per cell_line
number_of_expressed_het_snps_output_file <- paste0(visualize_allelic_counts_dir, "number_of_expressed_het_snps_per_individual_cell_line_boxplot_",het_thresh,".png")
number_of_expressed_het_snps_per_individual_cell_line_binned(total_counts, number_of_expressed_het_snps_output_file, het_thresh)


# Experiment with number of heterozygous sites that remain when we use various filters.
# Specifically wish to vary:
### a. Minimum number of cell lines heterozygous
### b. Minimum fraction of remaining samples that have bi-allelic expression (fewer than n reads mapping, or less than 1% of reads mapping to one allele)
### c. The n reads in step b
number_of_heterozygous_sites_after_filters_output_file <- paste0(visualize_allelic_counts_dir, "number_of_heterozygous_sites_at_various_filters_lineplot_min_reads_",het_thresh,".png")
number_of_heterozygous_sites_at_various_filters_lineplot(ref_counts, total_counts, sample_info, number_of_heterozygous_sites_after_filters_output_file, het_thresh)

# Experiment with number of heterozygous sites (IN TERMS OF UNIQUE GENES) that remain when we use various filters.
# Specifically wish to vary:
### a. Minimum number of cell lines heterozygous
### b. Minimum fraction of remaining samples that have bi-allelic expression (fewer than n reads mapping, or less than 1% of reads mapping to one allele)
### c. The n reads in step b
number_of_genes_after_filters_output_file <- paste0(visualize_allelic_counts_dir, "number_of_genes_at_various_filters_lineplot_min_reads_",het_thresh,".png")
number_of_genes_at_various_filters_lineplot(ref_counts, total_counts, sample_info, number_of_genes_after_filters_output_file, het_thresh, gene_names)


# Experiment with number of heterozygous sites that remain when we use various filters.
# Specifically wish to vary:
### a. Minimum number of cell lines heterozygous
### b. Minimum fraction of remaining samples that have bi-allelic expression (fewer than n reads mapping, or less than 1% of reads mapping to one allele)
### c. The n reads in step b
number_of_heterozygous_sites_after_filters_output_file <- paste0(visualize_allelic_counts_dir, "number_of_heterozygous_sites_at_various_independent_filters_lineplot_min_reads_",het_thresh,".png")
number_of_heterozygous_sites_at_various_filters_independent_time_step_lineplot(ref_counts, total_counts, sample_info, number_of_heterozygous_sites_after_filters_output_file, het_thresh)



# Histogram of number of mapped genes per heterozygous site
number_of_mapped_genes_per_site_output_file <- paste0(visualize_allelic_counts_dir, "number_of_mapped_genes_per_site_", het_thresh, ".png")
number_of_mapped_genes_per_site_histogram(gene_names, number_of_mapped_genes_per_site_output_file)
