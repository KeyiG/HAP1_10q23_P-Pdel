#!/bin/bash
#Script used to remove ribosomal RNA still left in the samples after library preparation
#Made on 181215

#UPPMAX commands (Uppsala Multidisciplinary Center for Advanced Computational Science)
#SBATCH -A proj_num
#SBATCH -p core
#SBATCH -n 8
#SBATCH -t 6:00:00
#SBATCH -J 210721_HiSAT2_rRNA_SE
#SBATCH --output=210721_HiSAT2_rRNA.out
#SBATCH --error=210721_HiSAT2_rRNA.err

#load packages. bioinfo-tools is loaded on uppmax in order to load all other packages used.
module load bioinfo-tools
module load HISAT2/2.1.0 

#file paths
FQ_PATH=/proj/trim_fastq/SE
OUTPUT_PATH=/proj/rmrRNA/SE
REF_PATH=/proj/reference/rRNA.dir

#loop to run HiSAT2 alignment to rRNA, and keeping unaligned reads.

cd ${FQ_PATH}
for i in $(ls *_tc.fq.gz | sed 's/_tc.fq.gz//g'); do	
	hisat2 \
		-p 8 \
		--rna-strandness R \
		--un-gz ${OUTPUT_PATH}/${i}_tc_rmrRNA.fastq.gz \
        	--summary-file ${OUTPUT_PATH}/${i}_tc_rRNAcontam.txt \
		-x ${REF_PATH}/hrRNA_combined \
		-U ${i}_tc.fq.gz \
		-S ${OUTPUT_PATH}/${i}_tc_rRNA.sam \
		>> ${OUTPUT_PATH}/${i}_rRNAalign_stdout.stderr.txt 2>&1
done


#Readme:
#-p specifies the number of computational cores/threads that will be used by the program
#--rna-strandness: strand-specific information. R or RF is the choice for dUTP based method which is used by most of the stranded lib prep kit. R for SE and RF for PE.
#--un-conc-gz: Write paired-end reads that fail to align concordantly to file(s) at <path>. Useful for rRNA removal. for SE, use --un-gz
#--summary-file: Print alignment summary to this file
#-x path to the pre-built genome index. Note that the index consists of multiple files ending in .ht2 , and only the shared part of the filename should be indicated (e.g. genome if the files are called genome.1.ht2 , genome.2.ht2 , etc).
#hrRNA_combined is made by combining all RNA fasta files from the databases HGNC, ENA, and SILVA on 171009.
#-1 the first-read mate FASTQ file
#-2 the second-read mate FASTQ file
#-S name of the result file that will be created
#>> send all messages from HISAT2 (including errors and warnings) into the specified file
