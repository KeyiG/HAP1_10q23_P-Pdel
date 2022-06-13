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

#loop to run HiSAT2 alignment to rRNA, and keeping unaligned reads.

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
#--rna-strandness: strand-specific information. Needs to be RF if using Illumina Truseq library preparation. R or RF is the choice for dUTP based method which is used by most of the stranded lib prep kit. R for SE and RF for PE.
#--un-conc-gz: Write paired-end reads that fail to align concordantly to file(s) at <path>. Useful for rRNA removal
#--summary-file: Print alignment summary to this file
#-x path to the pre-built genome index. Note that the index consists of multiple files ending in .ht2 , and only the shared part of the filename should be indicated (e.g. genome if the files are called genome.1.ht2 , genome.2.ht2 , etc).
#--dta: this option is present in JOnas script. it is for With this option, HISAT2 requires longer anchor lengths for de novo discovery of splice sites. This leads to fewer alignments with short-anchors, which helps transcript assemblers improve significantly in computation and memory usage. for transcript assembly. I don't think it is very necessary for our purpose. so I skip it.
#-1 the first-read mate FASTQ file
#-2 the second-read mate FASTQ file
#-S name of the result file that will be created
#>> send all messages from HISAT2 (including errors and warnings) into the specified file
