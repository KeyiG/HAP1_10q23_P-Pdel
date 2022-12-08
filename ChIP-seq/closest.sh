module load bioinfo-tools
module load BEDTools/2.29.2

tail -n +2 pooling_peakanno_delvsctrl.txt > pooling_peakanno_delvsctrl.bed
bedtools sort -i allfour_HAP1_RNAseq_DE_loc.bed > allfour_HAP1_RNAseq_DE_loc_sorted.bed
bedtools sort -i pooling_peakanno_delvsctrl.bed > pooling_peakanno_delvsctrl_sorted.bed
bedtools closest -d -a allfour_HAP1_RNAseq_DE_loc_sorted.bed -b pooling_peakanno_delvsctrl_sorted.bed > DE_DB_dis.txt
bedtools closest -d -a pooling_peakanno_delvsctrl_sorted.bed -b allfour_HAP1_RNAseq_DE_loc_sorted.bed > DB_DE_dis.txt











