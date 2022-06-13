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

FQ_PATH=/proj/publish_fastq/SE
OUTPUT_PATH=/proj/trim_fastq/SE

cd ${FQ_PATH}
for i in $(ls *.fq.gz | sed 's/.fq.gz//g'); do
	java -jar $TRIMMOMATIC_HOME/trimmomatic.jar \
		SE \
		-threads 8\
		-phred33 \
		${i}.fq.gz \
		${OUTPUT_PATH}/${i}_tc.fq.gz \
		ILLUMINACLIP:$TRIMMOMATIC_HOME/adapters/TruSeq3-SE.fa:2:30:10 \
		HEADCROP:10 \
		MINLEN:36 \
		>>${OUTPUT_PATH}/${i}.trimmomatic.stdout.stderr.txt 2>&1
done
