args = commandArgs(trailingOnly=TRUE)
library(edgeR)
library(Rsubread)
library(Biobase)
library(preprocessCore)

# Helper method to save DGE data structure to tab-deliminated text file
save_DGE_matrix <- function(counts, output_file) {
    #  Convert DGE data structure to matrix
    temp_mat <- as.matrix(counts)

    #  Edit colnames to include a header over the row-labels.
    revised_column_names <- colnames(temp_mat)
    revised_column_names[1] <- paste0("Gene_id\t",revised_column_names[1])

    write.table(temp_mat, output_file, quote=FALSE,col.names = revised_column_names)

}

# Helper method to save DGE data structure to tab-deliminated text file
save_DGE_matrix2 <- function(counts, output_file) {
    #  Convert DGE data structure to matrix
    temp_mat <- as.matrix(counts)

    #  Edit colnames to include a header over the row-labels.
    revised_column_names <- colnames(temp_mat)
    revised_column_names[1] <- paste0("Sample_name\t",revised_column_names[1])

    write.table(temp_mat, output_file, quote=FALSE,col.names = revised_column_names)

}


# Make several organizational changes to count data including:
#   1. Remove "Undetermined samples"
#   2. Remove replicated samples (19238 was sequenced twice. We are excluding the lane 8 version for downstream analysis)
#   3. Remove samples with less than 5 million mapped reads
#   4. Rename sample ids from sequencing name to name of form cell_line"_"time_step
organize_count_data <- function(counts, preprocess_total_expression_dir) {
    # Remove uninformative prefix from column names
    full_names <- substring(colnames(counts), 72, nchar(colnames(counts))-11)
    colnames(counts) <- full_names

    #  Assign cell line id and time-step to each sample
    num_samples <- nrow(counts$samples)
    counts$samples$cell_line <- substr(colnames(counts),1,5)  # cell line
    counts$samples$time <- substr(colnames(counts),7,8)  # time
    #  If time is a 1 digit number, it will have a "_" at the end. Remove the "_" if it exists
    for (i in 1:num_samples) {if (substring(counts$samples$time[i],2,2) == "_") {counts$samples$time[i] <- substring(counts$samples$time[i],1,1)}}

    #  Save count matrix containing all possible genes to text file
    output_file <- paste0(preprocess_total_expression_dir, "raw_counts_all_genes.txt")
    save_DGE_matrix(counts, output_file)

    return(counts)
}



#  Perform RPKM transformation
rpkm_transformation <- function(counts, preprocess_total_expresssion_dir) {
    #  Perform RPKM transformation using edgeR
    rpkm_data <- rpkm(counts, counts$genes$Length)

    #  Save RPKM transformed data to text file
    output_file <- paste0(preprocess_total_expression_dir,"rpkm_all_genes.txt")
    save_DGE_matrix(rpkm_data, output_file)

    return(rpkm_data)
}

#  Add biotype category to counts data structure (ie. is gene protein coding)
add_biotype_to_count_data_structure <- function(counts, exon_table) {
    #  Loop through all genes
    biotype <- NULL
    for (i in 1:nrow(counts)) {
        #  GeneID
        i_gene_id <- counts$genes$GeneID[i]

        #  Table where each row is one of the exons composing i_gene_id
        i_exons <- exon_table[exon_table$GeneID == i_gene_id,]

        # If all exons are classified as "protein coding", we have protein coding gene
        if (sum(i_exons$Biotype != "protein_coding") == 0) {
            biotype <- c(biotype,"protein_coding")
        } else {  # Not all exons are classified as protein coding
            biotype <- c(biotype, "not_protein_coding")
        }

    }
    counts$genes$biotype <- biotype
    return(counts)
}

#  Check if ith gene is autosomal
#  Chromosome string is a ";" seperated string (a vector). 
#  Where each element in the vector the chromosome number of the jth exon for this gene
check_if_gene_is_autosomal <- function(chromosome_string) {
    # Initialize check to be true
    autosomal <- TRUE

    #  Convery chromosome string into a vector, and return only unique elements (chromosomes) in the vector
    chromosome_array <- unique(strsplit(chromosome_string,";")[[1]])

    #  Just a check to ensure that every exon composing a gene belongs to the same chromosome
    if (length(chromosome_array) != 1) {
        print("ASSUMPTION ERROR")
        print(chromosome_string)
    }

    #  The chromosome number of the ith gene
    chromer <- chromosome_array[1]

    #  Check if the chromosome is non-autosomal
    if (chromer == "X" || chromer ==  "Y" || chromer == "MT") {
        autosomal <- FALSE
    }

    return(autosomal)
}


#  Check if ith gene passes min-read thresholds (criteria number three) (TRUE if passes)
#  count_vector is vector of read counts for ith gene
#  rpkm_vector is vecotr of rpkms for ith gene
#  Pass filter if greater than or equal to min_samples that have both:
#      a. rpkm greater than or equal to  rpkm_threshold
#      b. read count greater than or equal to count_threshold
check_if_gene_passes_min_read_threshold <- function(count_vector, rpkm_vector, min_samples, rpkm_threshold, count_threshold) {
    # Initialize boolean output variable to be true
    pass_threshold <- TRUE

    #  Create vector of length samples, where element is 1 if sample has count and rpkm greater than or equal to thresholds
    binary_samples <- (count_vector >= count_threshold) * (rpkm_vector >= rpkm_threshold)

    #  Check to make sure enough samples pass threshold
    if (sum(binary_samples) < min_samples) {
        pass_threshold <- FALSE
    }
    return(pass_threshold)
}

# Check if ith gene is protein_coding
check_if_gene_is_protein_coding <- function(biotype) {
    #  Initialize boolean to true
    boolean <- TRUE

    if (biotype != "protein_coding") {
        boolean <- FALSE
    }

    return(boolean)
}

# Peform filtering on genes. We include genes such that:
#     1. Genes are protein-coding
#     2. genes are autosomal 
#     3. genes have at least 10 samples such that RPKM >= .1 and counts >= 6
filter_genes <- function(counts, rpkm_data, preprocess_total_expression_dir, normalization_method) {

    #  Kepp list of indices of genes that pass our three filters
    pass_filters <- NULL

    #  Put count data temporarily into matrix form for computational ease.
    count_matrix <- as.matrix(counts)

    # Loop through each gene
    for (i in 1:nrow(counts)) {

        # Check if ith gene is protein coding (Criteria number 1)
        is_protein_coding <- check_if_gene_is_protein_coding(counts$genes$biotype[i])

        #  Check if ith gene is autosomal (criteria number 2) (TRUE if autosomal)
        is_autosomal <- check_if_gene_is_autosomal(counts$genes$Chr[i])

        #  Check if ith gene passes min-read thresholds (criteria number three) (TRUE if passes)
        min_samples <- 10 
        rpkm_threshold <- .1
        count_threshold <- 6
        pass_min_read_threshold <- check_if_gene_passes_min_read_threshold(count_matrix[i,], rpkm_data[i,], min_samples, rpkm_threshold, count_threshold)


        #  Check if ith gene passes all three criteria
        if (is_protein_coding == TRUE && is_autosomal == TRUE && pass_min_read_threshold == TRUE) {
            pass_filters <- c(pass_filters, i)
        }
    }

    #  Filter genes in both count data and rpkm data
    counts <- counts[pass_filters,]
    rpkm_data <- rpkm_data[pass_filters,]

    if (normalization_method == "log_rpkm") {
        rpkm_data <- log2(rpkm_data + 1)
    }

    #  Write filtered matrices to output files
    count_output_file <- paste0(preprocess_total_expression_dir,"raw_counts.txt")
    save_DGE_matrix(counts, count_output_file)

    rpkm_output_file <- paste0(preprocess_total_expression_dir,"rpkm.txt")
    save_DGE_matrix(rpkm_data, rpkm_output_file)

    return(list(counts,rpkm_data))
}

#  Peform quantile normalization, and then standardize each row 
quantile_normalize_and_standardize <- function(rpkm_data, preprocess_total_expression_dir, normalization_method) {

    if (normalization_method == "inverse_gaussian") {
        # Quantile normalize (so the samples have equivalent variance)
        quantile_normalized_samples <- normalize.quantiles(as.matrix(rpkm_data))

        # Quantile Normalized without standardization or projection onto gaussian
        colnames(quantile_normalized_samples) <- colnames(rpkm_data)
        rownames(quantile_normalized_samples) <- rownames(rpkm_data)

        #  Write results to output file
        output_file <- paste0(preprocess_total_expression_dir,"quantile_normalized_no_projection.txt")
        save_DGE_matrix(quantile_normalized_samples, output_file) 


        temp_mat <- t(apply(quantile_normalized_samples, 1, rank, ties.method = "average"))
        quantile_normalized <- qnorm(temp_mat / (ncol(temp_mat)+1));


        # Quantile normalization with gaussian projection
        colnames(quantile_normalized) <- colnames(rpkm_data)
        rownames(quantile_normalized) <- rownames(rpkm_data)

        #  Write results to output file
        output_file <- paste0(preprocess_total_expression_dir,"quantile_normalized.txt")
        save_DGE_matrix(quantile_normalized, output_file)
    # Applies to "rpkm" and "log_rpkm"
    } else {
        # Quantile normalize (so the samples have equivalent variance)
        quantile_normalized_samples <- normalize.quantiles(as.matrix(rpkm_data))
    
        # Quantile Normalized without standardization or projection onto gaussian
        colnames(quantile_normalized_samples) <- colnames(rpkm_data)
        rownames(quantile_normalized_samples) <- rownames(rpkm_data)

        #  Write results to output file
        output_file <- paste0(preprocess_total_expression_dir,"quantile_normalized_no_projection.txt")
        save_DGE_matrix(quantile_normalized_samples, output_file) 


        quantile_normalized <- quantile_normalized_samples
        # Standardize row by row (ie gene by gene)
        for (i in 1:nrow(quantile_normalized)) {
            quantile_normalized[i,] <- (quantile_normalized[i,] - mean(quantile_normalized[i,]))/sd(quantile_normalized[i,])
        }

        # Quantile normalization with gaussian projection
        colnames(quantile_normalized) <- colnames(rpkm_data)
        rownames(quantile_normalized) <- rownames(rpkm_data)

        #  Write results to output file
        output_file <- paste0(preprocess_total_expression_dir,"quantile_normalized.txt")
        save_DGE_matrix(quantile_normalized, output_file)
    }


    return(quantile_normalized)
}


#  Peform quantile normalization, and then standardize each row 
#  Do this for each time step independently
quantile_normalize_and_standardize_time_step_independent <- function(rpkm_data, preprocess_total_expression_dir, sample_info) {

    # Initialize matrix that will contain quantile normalized expression
    quantile_normalized <- matrix(0, dim(rpkm_data)[1], dim(rpkm_data)[2])

    # Loop through time steps..
    for (temp_time_step in 0:15) {
        # Get matrix of samples belonging to only the current time step
        time_step_indices <- sample_info$time == as.character(temp_time_step)
        time_step_rpkm_matrix <- as.matrix(rpkm_data[,time_step_indices])
        # Quantile normalize (so the samples have equivalent variance)
        time_step_quantile_normalized_samples <- normalize.quantiles(time_step_rpkm_matrix)

        # Project each gene onto a gaussian
        time_step_temp_mat <- t(apply(time_step_quantile_normalized_samples, 1, rank, ties.method = "average"))
        time_step_quantile_normalized <- qnorm(time_step_temp_mat / (ncol(time_step_temp_mat)+1));

        # Now add time step independent results back into the full matrix for storage purposes
        quantile_normalized[,time_step_indices] <- time_step_quantile_normalized
    }

    colnames(quantile_normalized) <- colnames(rpkm_data)
    rownames(quantile_normalized) <- rownames(rpkm_data)

    #  Write results to output file
    output_file <- paste0(preprocess_total_expression_dir,"time_step_independent_quantile_normalized.txt")
    save_DGE_matrix(quantile_normalized, output_file)

    return(quantile_normalized)
}


preprocess_total_expression_dir = args[1]
exon_file = args[2]
bam_dir = args[3]
normalization_method = args[4]


# Load in table regarding exon information
exon_table <- read.table(exon_file,header=TRUE)




##############################################################################################################
#  Convert from BAMs to count based data (Using edgeR)
fc <- featureCounts(Sys.glob(paste0(bam_dir,"*bam")),annot.ext=exon_file)
counts <- DGEList(counts=fc$counts, genes=fc$annotation[,c("GeneID","Length","Chr","Start","End","Strand")])

#  Add biotype information to count data structure (ie. whether gene is protein_coding or not_protein_coding)
counts <- add_biotype_to_count_data_structure(counts, exon_table)
saveRDS(counts, paste0(preprocess_total_expression_dir,"raw_counts.rds"))
saveRDS(fc, paste0(preprocess_total_expression_dir,"fc.rds"))
################################################################################################################
# Use this part if already converted from bams once!
#fc <- readRDS(paste0(preprocess_total_expression_dir,"fc.rds"))
#counts <- readRDS(paste0(preprocess_total_expression_dir,"raw_counts.rds"))
################################################################################################################



# Make several organizational changes to count data including:
#   1. change sample names from full path of bam file to just cellLine_timeStep
#   2. Add cell_line and time_step columns to counts$samples data table
counts <- organize_count_data(counts, preprocess_total_expression_dir)

#  Perform RPKM transformation
rpkm_data <- rpkm_transformation(counts, preprocess_total_expression_dir)


# Peform filtering on genes. We include genes such that:
#     1. Genes are protein-coding
#     2. genes are autosomal 
#     3. genes have at least 10 samples such that RPKM >= .1 and counts >= 6
temp_data_struct <- filter_genes(counts, rpkm_data, preprocess_total_expression_dir, normalization_method)
counts <- temp_data_struct[[1]]
rpkm_data <- temp_data_struct[[2]]


#  Peform quantile normalization, and then standardize each row 
#  Do this for each time step independently
quantile_normalized_time_step_independent_data <- quantile_normalize_and_standardize_time_step_independent(rpkm_data, preprocess_total_expression_dir, counts$samples)


#  Peform quantile normalization, and then standardize each row 
quantile_normalized_data <- quantile_normalize_and_standardize(rpkm_data, preprocess_total_expression_dir, normalization_method)




#  Save sample information to tab-delimited output file
save_DGE_matrix2(counts$samples, paste0(preprocess_total_expression_dir, "sample_info.txt"))



