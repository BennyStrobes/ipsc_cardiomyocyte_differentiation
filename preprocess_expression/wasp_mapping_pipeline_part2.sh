#!/bin/bash
#SBATCH --mem=30G --time=15:30:00 --partition=broadwl
standard_id="$1"
genotype_dir="$2"
fastq_dir="$3"
wasp_intermediate_dir="$4"
genome_dir="$5"
vcf_file="$6"
preprocess_allelic_counts_dir="$7"
chrom_info_file="$8"


date
module load samtools/1.4.1


#########################################################
# A. Map the fastq files using your favorite mapper/options and filter for quality using a cutoff of your choice (SUBREAD)
#########################################################
# This part takes around 25 minutes
echo "WASP STEP 2: Initial Mapping of fastq file (with subread)"

fq_file=$fastq_dir$standard_id".fastq.gz"
echo $fq_file

Rscript run-subread.R $fq_file $wasp_intermediate_dir $genome_dir


#  Sort and index the bams
samtools sort -o $wasp_intermediate_dir$standard_id.sort.bam $wasp_intermediate_dir$standard_id.bam
samtools index $wasp_intermediate_dir$standard_id.sort.bam



#########################################################
# B. Use find_intersecting_snps.py (WASP script) to identify reads that may have mapping bias
#########################################################
# This part takes around 10 minutes
echo "WASP STEP 3: find_intersecting_snps.py"

python find_intersecting_snps.py \
    --is_sorted \
    --output_dir $wasp_intermediate_dir \
    --snp_tab $genotype_dir"snp_tab.h5" \
    --snp_index $genotype_dir"snp_index.h5" \
    --haplotype $genotype_dir"haps.h5" \
    --samples $genotype_dir"used_samples.txt" \
    $wasp_intermediate_dir$standard_id.sort.bam


#########################################################
# C. Map the filtered fastq file a second time using the same arguments as the previous mapping
#########################################################
# This part takes around 10 minutes
echo "WASP STEP 4: Second mapping"

fq_file=$wasp_intermediate_dir$standard_id".sort.remap.fastq.gz"

Rscript run-subread.R $fq_file $wasp_intermediate_dir $genome_dir

#  Sort and index the bams
samtools sort -o $wasp_intermediate_dir$standard_id.sort.remap.sort.bam $wasp_intermediate_dir$standard_id.sort.remap.bam
samtools index $wasp_intermediate_dir$standard_id.sort.remap.sort.bam



#########################################################
# D. Use filter_remapped_reads.py (WASP script) to filter out reads where one or more of the allelic versions of the reads fail to map back to the same location as the original read
#########################################################
#  This part takes around 3 minutes
echo "WASP STEP 5: filter_remapped_reads.py"
python filter_remapped_reads.py \
        $wasp_intermediate_dir$standard_id.sort.to.remap.bam \
        $wasp_intermediate_dir$standard_id.sort.remap.sort.bam \
        $wasp_intermediate_dir$standard_id.keep.bam



#########################################################
# E. Merge the bams we plan to use. Then sort and index
#########################################################
#  This part takes around 5 minutes
echo "WASP STEP 6: Merge bams"

samtools merge $wasp_intermediate_dir$standard_id.keep.merge.bam \
                $wasp_intermediate_dir$standard_id.keep.bam \
                $wasp_intermediate_dir$standard_id.sort.keep.bam

samtools sort -o $wasp_intermediate_dir$standard_id.keep.merge.sort.bam $wasp_intermediate_dir$standard_id.keep.merge.bam
samtools index $wasp_intermediate_dir$standard_id.keep.merge.sort.bam



#########################################################
# F. Filter Duplicate Reads using rmdup.py (WASP script)
#########################################################
#  This part takes around 5 minutes
echo "WASP STEP 7: Filter duplicate reads using rmdup.py"
python rmdup.py $wasp_intermediate_dir$standard_id.keep.merge.sort.bam $wasp_intermediate_dir$standard_id.wasp_corrected.bam




#########################################################
# G. Run GATK ASE read counter (using output from wasp mapping pipeline)
#########################################################

echo "Adding uninformative groupnames to bam file (so GATK accepts the bams)"
#  This is only done so GATK accepts the bams. The group names mean nothing!
java -jar picard.jar AddOrReplaceReadGroups \
      I=$wasp_intermediate_dir$standard_id.wasp_corrected.bam  \
      O=$wasp_intermediate_dir$standard_id.wasp_corrected2.bam \
      RGID=4 \
      RGLB=lib1 \
      RGPL=illumina \
      RGPU=unit1 \
      RGSM=20

echo "Sorting one more time before running GATK"
samtools sort -o $wasp_intermediate_dir$standard_id.wasp_corrected3.bam $wasp_intermediate_dir$standard_id.wasp_corrected2.bam
samtools index $wasp_intermediate_dir$standard_id.wasp_corrected3.bam


echo "Reordering bam file using picard: reorder sam (ensures correct format for GATK)"
java -jar picard.jar ReorderSam \
    I=$wasp_intermediate_dir$standard_id.wasp_corrected3.bam \
    O=$wasp_intermediate_dir$standard_id.wasp_corrected_gatk_ready.bam \
    R=$genome_dir"Homo_sapiens.GRCh37.dna_sm.fa" \
    CREATE_INDEX=TRUE


echo "Running GATK"
 java -jar GenomeAnalysisTK.jar \
   -R $genome_dir"Homo_sapiens.GRCh37.dna_sm.fa" \
   -T ASEReadCounter \
   -o $preprocess_allelic_counts_dir$standard_id"_wasp_corrected.txt" \
   -I $wasp_intermediate_dir$standard_id.wasp_corrected_gatk_ready.bam \
   -sites $vcf_file \
   --outputFormat "TABLE"






#########################################################
# H. Convert bam file from this individual into h5 format using bam2h5_tables_update.py (WASP script)
#########################################################

# Extact cell line name (INDIVIDUAL) from $standard_id
INDIVIDUAL=`sed 's/_/\n/g' <<< $standard_id | head -1`

time_step=`sed 's/_/\n/g' <<< $standard_id | head -2 | tail -1`

ALL_SAMPLES_FILE=$genotype_dir"all_genotyped_samples.txt"


python bam2h5_tables_update.py --chrom $chrom_info_file \
    --snp_index $genotype_dir"snp_index.h5" \
    --snp_tab $genotype_dir"snp_tab.h5" \
    --haplotype $genotype_dir"haps.h5" \
    --samples $ALL_SAMPLES_FILE \
    --individual $INDIVIDUAL \
    --ref_as_counts $preprocess_allelic_counts_dir"ref_as_counts."$INDIVIDUAL"_"$time_step".h5" \
    --alt_as_counts $preprocess_allelic_counts_dir"alt_as_counts."$INDIVIDUAL"_"$time_step".h5" \
    --other_as_counts $preprocess_allelic_counts_dir"other_as_counts."$INDIVIDUAL"_"$time_step".h5" \
    --read_counts $preprocess_allelic_counts_dir"read_counts."$INDIVIDUAL"_"$time_step".h5" \
    $wasp_intermediate_dir$INDIVIDUAL"_"$time_step"_merged.wasp_corrected_gatk_ready.bam"


echo "DONE!"
date

