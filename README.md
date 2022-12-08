This repository contains the script for the work: “Cas9-mediated genome engineering exaggerates collateral deletion at the PTEN gene locus mimicking cancer profiles”.

Three folders ChIP-seq, RNA-seq and patient_sample_analysis contain the script used to analyze the data.

ChIP-seq folder: BWA_align_MEM.sh is for alignment; rmdu.sh is for remove the PCR duplicates and the read aligned to exclusion list; bam_bedgraph is for generating bedgraph for visualization; MACS_HAP1.sh is for calling peaks; ChIP-seq_analysis.Rmd is for differential acetylation of H3K27 (DAc) analysis and visualization; files folder contains the files obtained from RNA-seq analysis and used here, as well as the sample sheet for DAc analysis.

RNA-seq folder: PE folder and SE folder contain the script for analysing pair-end (PE) and single-end (SE) data: trim_PE(SE).sh is for trim the fastq files; rm_rRNA_PE(SE).sh is for remove the read aligned to rRNA; align_PE(SE).sh is for alignment; bam_bedgraph_PE(SE)_rmsc is for generating bedgraph files to visualize; featurecount_genecode_PE(SE).sh is for read counting, and its output counting table is used directly for downstream differential expression analysis in R. RNA-seq_DEanalysis.Rmd is for pooling all the dataset together. RNA-seq_individual_DEanalysis.Rmd is for performing differential expression analysis to each dataset separately.

Patient_sample_analysis folder: patient_sample_analysis.Rmd contains the script analysing the patient data and visualization. sampleID_used_info.txt contains the sample ID used for generating figure to make it reproducible.


