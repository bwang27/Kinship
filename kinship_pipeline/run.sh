
############## For the kinship analysis, it needs to check the chip sample names and cell lines, change it formally, if there's new sample sequenced each time, check whether that cell line is in WGS or not, add each bam size and cell line information to files in feature folder; the input file for kinship is using count file from dipseek, so if there is new sample, need to run dipseek first to get the txt count file ##################

DATA="/research/groups/xmagrp/projects/SensitiveDetection/common/bwang/kinship/20250825/feature";

CODE="/research/groups/xmagrp/projects/SensitiveDetection/common/bwang/kinship/20250825/pipeline/tmp";

perl creat_job.pl CHIPSeq.lst | sh;

cat out* >merge_v_CHIPWGS.txt;

#alternatively, you can clean the files names by sed,this is very variable due to the dataset, e.g. 
#sed -i "s/ai_check_PedDepAcc_//g" merge_v_CHIPWGS.txt;
#sed -i "s/.txt//g" merge_v_CHIPWGS.txt;
#sed -i "s/ai_check_CHIP_//g" merge_v_CHIPWGS.txt;
#sed -i "s/ai_check_osteo_//g" merge_v_CHIPWGS.txt;
#sed -i "s/PedDep_CHIPSEQ_//g" merge_v_CHIPWGS.txt;
#sed -i "s/epig_//g" merge_v_CHIPWGS.txt;
#sed -i "s/ai_check_WGS_//g" merge_v_CHIPWGS.txt;

mkdir output;

mv out.* output;

perl convert_2_matrix.pl merge_v_CHIPWGS.txt thematrix.txt;

Rscript clustering_bw.R;

perl get_pair_value.pl thematrix_clustered.txt;

perl pair_uniq.pl;

perl faster_add_more_infotokinship_chip.pl merge_v_CHIPWGS.txt final_pairs_4plot.txt $DATA/wgs_chip.bam.size $DATA/wgs.chip.cellline.txt

perl better.get_overlapped.chip.wgs.kinship.pl CHIP_WGS_pair.txt overlapped.wgs.id;

perl get_realid.pl $DATA/namekey/namekey.txt v2_ovlap_WGS_kinship_pair.txt; 

perl get_each_cluster.pl;

perl find_2ndkinship.pl;

################ merged version: in this part I merged some codes into one as merged.pl which has ~400 lines, alternatively you can prepair the input file and feed into the merged code as required, the output is sufficient to tell which chip sample closest to which wgs #############################

DATA="/research/groups/xmagrp/projects/SensitiveDetection/common/bwang/kinship/20250825/feature";

perl creat_job.pl CHIPSeq.lst | sh;

cat out* >merge_v_CHIPWGS.txt;

perl convert_2_matrix.pl merge_v_CHIPWGS.txt thematrix.txt;

Rscript clustering_bw.R;

perl main_kinship.pl thematrix_clustered.txt merge_v_CHIPWGS.txt final_pairs_4plot.txt $DATA/wgs_chip.bam.size $DATA/wgs.chip.cellline.txt $DATA/namekey.txt;

###the main output from here is newfinal_pair.txt which includes all the most important information for each kinship pairs, and final.corrected_samples_with_closest_kinship_in_diffcellline.txt which is the summary file of those correctly labeled CHIP samples and their 2nd ranked wgs samples to further confirm;
