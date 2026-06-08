import sys
import glob
import os
import re
import subprocess
if not (len(sys.argv)==3):
	print("please input SNVIndel_fn,hg_version")
	sys.exit()
SNVIndel_fn=sys.argv[1]
genome_version=sys.argv[2]

#if not os.path.isdir(out_folder):
#	os.mkdir(out_folder)
left_flanking_length=75
right_flanking_length=75
if genome_version =="hg38":
	bit_fn="/research_jude/rgs01_jude/groups/xmagrp/projects/AMLRelapse/common/NeoSplice/code/yliu1/reference/hg38.gz"
	bit_fa_fn="/research_jude/rgs01_jude/groups/xmagrp/projects/AMLRelapse/common/NeoSplice/code/yliu1/reference/hg38.2bit"
if genome_version =="hg19":
	bit_fn="/research_jude/rgs01_jude/groups/xmagrp/projects/AMLRelapse/common/NeoSplice/code/yliu1/reference/hg19.gz"
	bit_fa_fn="/research_jude/rgs01_jude/groups/xmagrp/projects/AMLRelapse/common/NeoSplice/code/yliu1/reference/hg19.2bit"

def get_revcmp(strs):
	"""return reverse complement of seq"""
	cpm={"A":"T","T":"A","C":"G","G":"C","N":"N"}
	rev=strs[::-1]
	revcmp=""
	for i in rev:
		if i in cpm:
			revcmp+=cpm[i]
		else:
			revcmp+="N"
	return revcmp

def extract_seq(chro,start,end,strand):
	"""return seq in the defined region (1 based)"""
	if "chr" not in chro:
		chro="chr"+chro
	seq_command="samtools faidx "+bit_fn+" "+chro+":"+str(start)+"-"+str(end)
	seq=("").join(str(subprocess.check_output(seq_command,shell=True),"utf-8").strip().split("\n")[1:]).upper()
	pos_range=[chro+"."+str(i) for i in range(int(start),int(end)+1)]
	if strand=="-":
		seq=get_revcmp(seq)
		pos_range=pos_range[::-1]
	return (seq,pos_range)

def make_contig_4SNVIndel(SNVIndel):
	chro,pos,ref_allele,mut_allele=SNVIndel.split(".")
	left_end_pos=int(pos)-1	
	right_start_pos=int(pos)+len(ref_allele)
	left_seq,left_pos_range=extract_seq(chro,left_end_pos-left_flanking_length+1,left_end_pos,"+")
		
	right_seq,right_pos_range=extract_seq(chro,right_start_pos,right_start_pos+right_flanking_length-1,"+")
	SNVIndel_seq=left_seq+mut_allele+right_seq
	SNVIndel_seq_2=left_seq+"_"+mut_allele+"_"+right_seq
	SNVIndel_seq_pos_range=left_pos_range+["ins."+str(i) for i in range(len(mut_allele))]+right_pos_range

	#ref_allel_check=extract_seq(chro,int(pos),int(pos)+len(ref_allele)-1,"+")[0]
	#if ref_allel_check== ref_allele:
	#	ref_allel_check_result="good_ref_allele"
	#else:
	#	ref_allel_check_result="bad_ref_allele"
	#return (SNVIndel_seq,SNVIndel_seq_pos_range,ref_allel_check_result)
	return (SNVIndel_seq,SNVIndel_seq_2,SNVIndel_seq_pos_range)

SNVIndel_candidates={}
with open (SNVIndel_fn,"r")as in_fh:
	for line in in_fh:
		mut=line.strip()
		chro,pos,ref_allele,mut_allele=mut.split(".")
		if chro not in (["chr"+str(i) for i in range(1,23)]+["chrX","chrY"]):
			continue
		if len(ref_allele)==len(mut_allele)==1:
			mut_type="SNV" #without padding
		else:
			mut_type="Indel" #with padding including MNV
		mut_seq,mut_seq_2,mut_seq_pos_range=make_contig_4SNVIndel(mut)	
		SNVIndel_candidates[mut]=[mut_seq,mut_seq_2,mut_seq_pos_range,mut_type]

SNVIndel_fa_fn=SNVIndel_fn+".fa"
with open(SNVIndel_fa_fn,"w")as mut_fa_fh:
	for mut in SNVIndel_candidates:
		mut_seq=SNVIndel_candidates[mut][0]
		mut_fa_fh.write(">"+mut+"\n")
		mut_fa_fh.write(mut_seq+"\n")

blat_fn=SNVIndel_fa_fn+".blast8"
blat_command="/research/groups/xmagrp/projects/SensitiveDetection/xmagrp/code/dipseek/apps/blat "+bit_fa_fn+" "+SNVIndel_fa_fn+" -out=blast8 "+blat_fn
os.system(blat_command)

blat_data={}
with open (blat_fn,"r")as in_fh:
	for line in in_fh:
		lst=line.strip().split("\t")
		mut=lst[0]
		hit_chro=lst[1]
		if hit_chro not in (["chr"+str(i) for i in range(1,23)]+["chrX","chrY"]):
			continue
		blat_score=float(lst[-1])
		blat_identity=float(lst[2])
		blat_exact_match_len=int(lst[3])
		mut_match_start=int(lst[6])
		mut_match_end=int(lst[7])
		mut_hit_len=abs(mut_match_end-mut_match_start)+1
		genome_match_start=int(lst[8])
		genome_match_end=int(lst[9])
		hit_center_pos=0.5*(int(genome_match_start)+int(genome_match_end))
		if mut not in SNVIndel_candidates:
			continue
		if mut not in blat_data:
			blat_data[mut]={
				"hit_max_score":0,
				"nonhit_max_score":0,
				"hit_score_lst":[],
				"nonhit_score_lst":[]
			}	
		mut_seq,mut_seq_2,mut_seq_pos_range,mut_type=SNVIndel_candidates[mut]
		chro,pos,ref_allele,mut_allele=mut.split(".")
		if hit_chro !=chro:
			blat_data[mut]["nonhit_score_lst"].append(blat_score)
		else:	
			if mut_type=="SNV":
				if hit_center_pos ==int(pos):
					blat_data[mut]["hit_score_lst"].append(blat_score)
				else:
					blat_data[mut]["nonhit_score_lst"].append(blat_score)
			if mut_type=="Indel":
				if (int(pos)-30)<hit_center_pos<(int(pos)+30):
					blat_data[mut]["hit_score_lst"].append(blat_score)
				else:
					blat_data[mut]["nonhit_score_lst"].append(blat_score)
for mut in blat_data:
	hit_score_lst=blat_data[mut]["hit_score_lst"]
	nonhit_score_lst=blat_data[mut]["nonhit_score_lst"]
	if len(hit_score_lst)>0:
		hit_score_lst=sorted(hit_score_lst)
		blat_data[mut]["hit_max_score"]=hit_score_lst[-1]
	if len(nonhit_score_lst)>0:
		nonhit_score_lst=sorted(nonhit_score_lst)
		blat_data[mut]["nonhit_max_score"]=nonhit_score_lst[-1]

blat_out_fn=blat_fn+".blat_check.txt"
blat_out_fh=open(blat_out_fn,"w")
blat_out_fh.write("SNVIndel\tSNVIndel_contig\tmut_type\tblat_hit_bitscore\tblat_nonhit_max_bitscore\n")

with open (SNVIndel_fn,"r") as in_fh:
	for line in in_fh:
		lst=line.strip().split("\t")
		mut=lst[0]
		mut_seq,mut_seq_2,mut_seq_pos_range,mut_type=SNVIndel_candidates[mut]
		hit_score="0"
		nonit_score="0"
		if mut in blat_data:
			hit_score=blat_data[mut]["hit_max_score"]
			nonhit_score=blat_data[mut]["nonhit_max_score"]
		blat_out_fh.write(mut+"\t"+mut_seq_2+"\t"+mut_type+"\t"+str(hit_score)+"\t"+str(nonhit_score)+"\n")
blat_out_fh.close()
fn_to_delete=[SNVIndel_fa_fn,blat_fn]
for fn in fn_to_delete:
	rm_command="rm "+fn
	# os.system(rm_command)
