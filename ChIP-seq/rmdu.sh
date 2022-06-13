#!/bin/bash
#made by Keyi

#UPPMAX commands
#SBATCH -A proj_num
#SBATCH -p core
#SBATCH -n 8
#SBATCH -t 5:00:00
#SBATCH -J rm_dubl
#SBATCH --output=rm_dubl.out
#SBATCH --error=rm_dubl.err

#module
module load bioinfo-tools
module load samtools/1.9
module load NGSUtils/0.5.9

BAM_PATH=/proj/alignment
OUTPUT_PATH=/proj/rmdubl

cd ${BAM_PATH}
for i in $(ls *.bam | sed 's/.bam//g'); do
	samtools sort -@ 2 -n -o ${OUTPUT_PATH}/${i}_nsort.bam ${i}.bam
	samtools fixmate -@ 2 -m ${OUTPUT_PATH}/${i}_nsort.bam ${OUTPUT_PATH}/${i}_fixmate.bam
	samtools sort -@ 2 -o ${OUTPUT_PATH}/${i}_positionsort.bam ${OUTPUT_PATH}/${i}_fixmate.bam
        samtools markdup -r -@ 2 ${OUTPUT_PATH}/${i}_positionsort.bam ${OUTPUT_PATH}/${i}_markdup.bam 
	bamutils filter ${OUTPUT_PATH}/${i}_markdup.bam \
		${OUTPUT_PATH}/${i}_rmdubl.bam \
		-excludebed hg38.blacklist.bed nostrand
done

module unload bioinfo-tools
module unload samtools/1.9
