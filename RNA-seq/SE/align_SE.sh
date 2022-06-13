#!/bin/bash
#Made on 181215

#UPPMAX commands (Uppsala Multidisciplinary Center for Advanced Computational Science)
#SBATCH -A proj_num
#SBATCH -p core
#SBATCH -n 8
#SBATCH -t 8:00:00
#SBATCH -J 210723_HiSAT2_ali_SE
#SBATCH --output=210723_HiSAT2_ali_SE.out
#SBATCH --error=210723_HiSAT2_ali_SE.err

#load packages. bioinfo-tools is loaded on uppmax in order to load all other packages used.
module load bioinfo-tools
module load HISAT2/2.1.0 
module load samtools/1.12

#file paths
FQ_PATH=/proj/rmrRNA/SE
OUTPUT_PATH=/proj/alignment/SE
REF_PATH=/proj/reference


cd ${FQ_PATH}
for i in $(ls *_tc_rmrRNA.fastq.gz | sed 's/_tc_rmrRNA.fastq.gz//g'); do	
	hisat2 \
		-p 8 \
		--rna-strandness R \
		-k 5 \
		--summary-file ${OUTPUT_PATH}/${i}_tc_align.txt \
		-x ${REF_PATH}/hg38_HISAT \
		-U ${i}_tc_rmrRNA.fastq.gz \
		-S ${OUTPUT_PATH}/${i}_tc_rmrRNA.sam
	samtools sort -@ 8 -o ${OUTPUT_PATH}/${i}_tc_rmrRNA.bam ${OUTPUT_PATH}/${i}_tc_rmrRNA.sam
	rm ${OUTPUT_PATH}/${i}_tc_rmrRNA.sam
done


#Readme:
# -k Default: 5 (linear index) or 10 (graph index).
#-p specifies the number of computational cores/threads that will be used by the program
#--rna-strandness: strand-specific information. R or RF is the choice for dUTP based method which is used by most of the stranded lib prep kit. R for SE and RF for PE.
#--summary-file: Print alignment summary to this file
#-x path to the pre-built genome index. Note that the index consists of multiple files ending in .ht2 , and only the shared part of the filename should be indicated (e.g. genome if the files are called genome.1.ht2 , genome.2.ht2 , etc).
#-1 the first-read mate FASTQ file
#-2 the second-read mate FASTQ file
#-S name of the result file that will be created
#>> send all messages from HISAT2 (including errors and warnings) into the specified file
