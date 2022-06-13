#!/bin/bash
 

#UPPMAX commands (Uppsala Multidisciplinary Center for Advanced Computational Science)
#SBATCH -A proj_num
#SBATCH -p core
#SBATCH -n 8
#SBATCH -t 01:00:00
#SBATCH -J 210723_FC_SE
#SBATCH --output=210723_FC_SE.out
#SBATCH --error=210723_FC_SE.err

#load packages. bioinfo-tools is loaded on uppmax in order to load all other packages used.
module load bioinfo-tools
module load subread/1.5.2

#file paths
BAM_PATH=/proj/alignment/SE
OUTPUT_PATH=/proj/featurecount
REF_PATH=/proj/reference/genecode_GTF_annotation

#run featureCounts
featureCounts \
	-t exon \
	-g gene_id \
	-s 2 \
	-T 8 \
	-M \
	-O \
	--fraction \
	-C \
	-a ${REF_PATH}/gencode.v32.annotation.gtf \
	-o ${OUTPUT_PATH}/GSE103405_GSE148175_210723.readCount \
	${BAM_PATH}/*.bam \
	&> ${OUTPUT_PATH}/GSE103405_GSE148175_210723.readCount.log



#Readme
#-t: Specify feature type in GTF annotation.
#-g: Specify attribute type in GTF annotation. Here we could chose e.g. transcript ID or gene ID. I chose gene ID, because I want to do DE analysis on gene level.
#-s: use '-s 2' if reversely stranded 
#-T: Number of computational cores/threads used for the analysis
#-p: The experiment is paired end, Keyi note: --countReadPairs is to count the fragment rather than reads, I think it is more appropiate for pair-end dataset.
#-M: Multi-mapping reads will also be counted. Each alignment will have 1 count or a fractional count if --fraction is specified
#-O: Allow reads that overlaps multiple features to be counted
#-C: If specified, the chimeric fragments (those fragments that have their two ends aligned to different chromosomes) will NOT be counted.
#-a: Name of the annotation file.
#-o: output file name

