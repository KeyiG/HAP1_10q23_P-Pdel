#!/bin/bash
#Updated on 20201006 by Keyi Geng

#UPPMAX commands
#SBATCH -A proj_num
#SBATCH -p core
#SBATCH -n 8
#SBATCH -t 4:00:00

#modules
module load bioinfo-tools
module load trimmomatic/0.36
#end of modules

FQ_PATH=/proj/publish_fastq/PE
OUTPUT_PATH=/proj/trim_fastq/PE

cd ${FQ_PATH}
for i in $(ls *_1.fastq.gz | sed 's/_1.fastq.gz//g'); do
	java -jar $TRIMMOMATIC_HOME/trimmomatic.jar \
		PE \
		-threads 8\
		-phred33 \
		${i}_1.fastq.gz \
		${i}_2.fastq.gz \
		${OUTPUT_PATH}/${i}_tc_R1.fastq.gz \
		${OUTPUT_PATH}/${i}_tc_unparied_R1.fastq.gz \
		${OUTPUT_PATH}/${i}_tc_R2.fastq.gz \
		${OUTPUT_PATH}/${i}_tc_unparied_R2.fastq.gz \
		ILLUMINACLIP:$TRIMMOMATIC_HOME/adapters/TruSeq3-PE.fa:2:30:10 \
		HEADCROP:10 \
		TRAILING:10 \
		MINLEN:36
done

