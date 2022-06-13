#!/bin/bash
#Made by keyi Geng


#UPPMAX commands
#SBATCH -A proj_num
#SBATCH -p core
#SBATCH -n 8
#SBATCH -t 12:00:00
#SBATCH -J BWA_MEM_211001
#SBATCH --output=bwa_mem_211001.out
#SBATCH --error=bwa_mem_211001.err

#module 
module load bioinfo-tools
module load bwa/0.7.17
module load samtools/1.9

FQ_PATH=/proj/fastq
OUTPUT_PATH=/proj/alignment

cd ${FQ_PATH}
for i in $(ls *.fq.gz | sed 's/.fq.gz//g'); do
	bwa mem -t 6 hg38.fa ${i}.fq.gz | samtools sort -@ 2 -o ${OUTPUT_PATH}/${i}.bam
done


module unload bwa/0.7.17
module unload samtools/1.9
