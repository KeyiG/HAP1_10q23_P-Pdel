#!/bin/bash
#made on 20190131

#UPPMAX commands
#SBATCH -A snic2020-15-292
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 8:00:00
#SBATCH -J MACS2_211004
#SBATCH --output=MACS2_211004.out
#SBATCH --error=MACS2_211004.err

#modules
module load bioinfo-tools
module load MACS/2.1.2

BAM_PATH=/proj/rmdubl

macs2 callpeak -t ${BAM_PATH}/CK0356_HAP1_t1_H3K27ac_rep2_rmdubl.bam -c ${BAM_PATH}/CK0361_rmbldu.bam -f BAM -g 2.7e9 -n CK0356_HAP1_t1_H3K27ac_rep2_CK361_deleinput.macs2 -q 0.01
macs2 callpeak -t ${BAM_PATH}/CK0358_HAP1_ctrl_H3K27ac_rep1_rmdubl.bam -c ${BAM_PATH}/CK0218_rmbldu.bam -f BAM -g 2.7e9 -n CK0358_HAP1_ctrl_H3K27ac_rep1.macs2 -q 0.01
macs2 callpeak -t ${BAM_PATH}/CK0359_HAP1_ctrl_H3K27ac_rep2_rmdubl.bam -c ${BAM_PATH}/CK0218_rmbldu.bam -f BAM -g 2.7e9 -n CK0359_HAP1_ctrl_H3K27ac_rep2.macs2 -q 0.01
macs2 callpeak -t ${BAM_PATH}/CK0360_HAP1_t1_H3K27ac_rep1_rmdubl.bam -c ${BAM_PATH}/CK0361_rmbldu.bam -f BAM -g 2.7e9 -n CK0360_HAP1_t1_H3K27ac_rep1_CK361_deleinput.macs2 -q 0.01

macs2 callpeak -t ${BAM_PATH}/SMARCA4KO_SMARCA2dTAG_DMSO_H3K27ac_rmdubl.bam -c ${BAM_PATH}/SMARCA4KO_SMARCA2dTAG_DMSO_IgG_rmdubl.bam -f BAM -g 2.7e9 -n SMARCA4KO_SMARCA2dTAG_DMSO_H3K27ac.macs2 -q 0.01 
macs2 callpeak -t ${BAM_PATH}/SMARCA4KO_SMARCA2dTAG_dTAG47_72h_H3K27ac_rmdubl.bam -c ${BAM_PATH}/SMARCA4KO_SMARCA2dTAG_dTAG47_72h_IgG_rmdubl.bam -f BAM -g 2.7e9 -n SMARCA4KO_SMARCA2dTAG_dTAG47_72h_H3K27ac.macs2 -q 0.01
macs2 callpeak -t ${BAM_PATH}/SMARCC1KO_SMARCC2dTAG_DMSO_H3K27ac_rmdubl.bam -c ${BAM_PATH}/SMARCC1KO_SMARCC2dTAG_DMSO_IgG_rmdubl.bam -f BAM -g 2.7e9 -n SMARCC1KO_SMARCC2dTAG_DMSO_H3K27ac.macs2 -q 0.01
macs2 callpeak -t ${BAM_PATH}/SMARCC1KO_SMARCC2dTAG_dTAG47_72h_H3K27ac_rmdubl.bam -c ${BAM_PATH}/SMARCC1KO_SMARCC2dTAG_dTAG47_72h_IgG_rmdubl.bam -f BAM -g 2.7e9 -n SMARCC1KO_SMARCC2dTAG_dTAG47_72h_H3K27ac.macs2 -q 0.01 --scale-to large
#the last one, control sequencing less than enrichment, so I need to ask macs2 to scale to large to avoid downsampling IP data. defaul is to small, downsizing.







