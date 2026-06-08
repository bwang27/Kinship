module unload sjcb/python/3.7.4
module load python/3.9.9
module load gcc/6.3.0

code_dir='/research/groups/xmagrp/projects/SensitiveDetection/xmagrp/code/dipseek'
cp $code_dir/run_dipseek_v2.py .
cp $code_dir/ngs_lib.py .
cp $code_dir/dipseek_lib.py .

python run_dipseek_v2.py blat test_indel_meta 75 test_indel_meta_blatresult.txt hg38
python run_dipseek_v2.py omni indel B500_Indel_step0_1.txt NameMap_WGS_344.lst SJLIFE_WGS_bam.3.lst B500_Indel hg38 WGS

