#!/bin/bash
#Made by keyi Geng


#UPPMAX commands
#SBATCH -A proj_num
#SBATCH -p core
#SBATCH -n 8
#SBATCH -t 15:00:00
#SBATCH -J bedgraph_220210
#SBATCH --output=bedgraph_220210.out
#SBATCH --error=bedgraph_220210.err


#modules
module load bioinfo-tools
module load deepTools/3.1.0
module load samtools/1.9
#follow this example: https://deeptools.readthedocs.io/en/develop/content/tools/bamCoverage.html#usage-examples-for-rna-seq

INPUT_PATH=/proj/alignment/PE
OUTPUT_PATH=/proj/bedgraph
cd ${INPUT_PATH}
for i in $(ls *.bam | sed 's/.bam//g'); do 
	samtools view -@ 6 -h ${i}.bam | awk '$6 !~ /H|S/{print}' | samtools view -bS > ${i}_rmsc.bam
	samtools index -@ 8 ${i}_rmsc.bam
	bamCoverage \
		-p 8 \
		--bam ${i}_rmsc.bam \
		--outFileName ${OUTPUT_PATH}/${i}_rmsc_fwd.bedgraph \
		--effectiveGenomeSize 2913022398 \
		--normalizeUsing BPM \
		--binSize 50 \
		--filterRNAstrand forward \
		--outFileFormat bedgraph
	bamCoverage \
		-p 8 \
                --bam ${i}_rmsc.bam \
                --outFileName ${OUTPUT_PATH}/${i}_rmsc_rvs.bedgraph \
                --effectiveGenomeSize 2913022398 \
                --normalizeUsing BPM \
                --binSize 50 \
                --filterRNAstrand reverse \
                --outFileFormat bedgraph

done

module unload bioinfo-tools
module unload deepTools/3.1.0

#in older version, argument normalizeTo1x are used followed by effective size, but in lastest version 3.1.3, two argument effectiveGenomeSize and normalizeUsing are used. and RPGC is 1x normalized method also recommended in ChIP-seq sample of their manual.

