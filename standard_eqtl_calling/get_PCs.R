args = commandArgs(trailingOnly=TRUE)


samples_file = args[1]
cht_input_file_dir = args[2]
parameter_string = args[3]
time_step = args[4]

# read names of individuals
indv.tab <- read.table(samples_file, header=F)

# build matrix of read counts for all target regions
# and individuals. Number of rows is same as number
# of rows in CHT input files. One column per individual
read.count.matrix <- NULL
for(i in 1:nrow(indv.tab)) {
    indv.id <- indv.tab$V1[i]
    cat("  ", indv.id, "\n", file=stderr())
    
    # read combined haplotype test input file for this individual
    filename <- paste(cht_input_file_dir,"haplotype_read_counts_",parameter_string,".", indv.id, "_", time_step,
                      ".txt.gz", sep="")
    cht.tab <- read.table(filename, header=T)
    if(is.null(read.count.matrix)) {
        # initialize matrix if first individual
        read.count.matrix <- matrix(nrow=nrow(cht.tab), ncol=nrow(indv.tab))
    }
    read.count.matrix[,i] <- cht.tab$REGION.READ.COUNT
}

# perform principal component analysis
pca <- prcomp(read.count.matrix)

# write PC loadings to stdout
write.table(pca$rotation, file="", sep=" ",
            append=FALSE, quote=FALSE,
            col.names=FALSE, row.names=FALSE)

