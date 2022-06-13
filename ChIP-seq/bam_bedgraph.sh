#!/bin/bash
#Made by keyi Geng, converting processed bam files to bedgraph.


#UPPMAX commands
#SBATCH -A proj_num
#SBATCH -p core
#SBATCH -n 8
#SBATCH -t 06:00:00
#SBATCH -J bedgraph_220223
#SBATCH --output=bedgraph_220223.out
#SBATCH --error=bedgraph_220223.err


#modules
module load bioinfo-tools
module load deepTools/3.1.0


BAM_PATH=/proj/rmdubl
OUTPUT=/proj/bedgraph

cd ${BAM_PATH}
for i in $(ls *rmdubl.bam | sed 's/.bam//g'); do \
	bamCoverage \
		-p 8 \
		--bam ${i}.bam \
		--outFileName ${OUTPUT}/${i}_MQ20.bedgraph \
		--effectiveGenomeSize 2913022398 \
		--normalizeUsing RPGC \
		--minMappingQuality 20 \
		--binSize 50 \
		--outFileFormat bedgraph \

done

#samtools used for remove soft-clipped read. Then use deep tools to generate bedgraph.
