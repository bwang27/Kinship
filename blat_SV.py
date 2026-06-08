import sys
import glob
import os
import re
import subprocess

if not (len(sys.argv)==5):
	print("please input SV_fn,hg_version,flanking_length,out_folder")
	sys.exit()

SV_fn=sys.argv[1]
genome_version=sys.argv[2]
flanking_length=sys.argv[3]
out_folder=sys.argv[4]

if not os.path.isdir(out_folder):
	os.mkdir(out_folder)

left_flanking_length=int(flanking_length)
right_flanking_length=int(flanking_length)

if genome_version =="hg38":
	bit_fn="/research_jude/rgs01_jude/groups/xmagrp/projects/AMLRelapse/common/NeoSplice/code/yliu1/reference/hg38.gz"
	bit_fa_fn="/research_jude/rgs01_jude/groups/xmagrp/projects/AMLRelapse/common/NeoSplice/code/yliu1/reference/hg38.2bit"
elif genome_version =="hg19":
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

SV_candidates={}
with open (SV_fn,"r")as in_fh:
	for line in in_fh:
		SV=line.strip().split("\t")[0] # modified
		chrA,posA,strandA,chrB,posB,strandB,insertion=SV.split(".")
		if strandA=="+":
			startA=str(int(posA)-left_flanking_length)
			endA=posA
		else:
			startA=posA
			endA=str(int(posA)+left_flanking_length)
		seqA,posA_range=extract_seq(chrA,startA,endA,strandA)

		if strandB=="+":
			startB=posB
			endB=str(int(posB)+right_flanking_length)
		else:
			startB=str(int(posB)-right_flanking_length)
			endB=posB
		seqB,posB_range=extract_seq(chrB,startB,endB,strandB)
		SV_seq=seqA+insertion+seqB
		SV_seq_pos_range=posA_range+["ins."+str(i) for i in range(len(insertion))]+posB_range
		SV_candidates[SV]=[SV_seq,SV_seq_pos_range]

SV_fa_fn=out_folder+"/"+SV_fn+".fa"
with open(SV_fa_fn,"w")as SV_fa_fh:
	for SV in SV_candidates:
		SV_seq=SV_candidates[SV][0]
		SV_fa_fh.write(">"+SV+"\n")
		SV_fa_fh.write(SV_seq+"\n")

blat_fn=SV_fa_fn+".blast8"
blat_command="blat "+bit_fa_fn+" "+SV_fa_fn+" -out=blast8 "+blat_fn
os.system(blat_command)

SVs_blat={}
with open (blat_fn,"r")as in_fh:
	for line in in_fh:
		lst=line.strip().split("\t")
		SV=lst[0]
		chro=lst[1]
		if float(lst[2])<95:
			continue
		blat_exact_match_len=int(lst[3])
		if chro not in (["chr"+str(i) for i in range(1,23)]+["chrX","chrY"]):
			continue
		SV_match_start=int(lst[6])
		SV_match_end=int(lst[7])
		SV_hit_len=abs(SV_match_end-SV_match_start)+1
		genome_match_start=int(lst[8])
		genome_match_end=int(lst[9])
		if genome_match_start<=genome_match_end:
			strand="+"
		else:
			strand="-"
		if SV not in SVs_blat:
			SVs_blat[SV]=[]
		SVs_blat[SV].append([chro,genome_match_start,genome_match_end,strand,SV_match_start,SV_match_end,SV_hit_len,blat_exact_match_len])

blat_out_fn=blat_fn+".blat_check.txt"

blat_out_fh=open(blat_out_fn,"w")
blat_out_fh.write("SV\tSV_contig\tblat_result\n")

with open (SV_fn,"r") as in_fh:
	for line in in_fh:
		lst=line.strip().split("\t")
		SV=lst[0] # modified
		ini_chrA,ini_posA,ini_strandA,ini_chrB,ini_posB,ini_strandB,insertion=SV.split(".")
		ini_posA=int(ini_posA)
		ini_posB=int(ini_posB)

		SV_seq=SV_candidates[SV][0]
		SV_seq_len=len(SV_seq)
		blated_SV="no_hit"
		SV_blat_check="bad_blat"
		
		if SV not in SVs_blat:
			blat_out_fh.write(SV+"\t"+SV_seq+"\t"+blated_SV+"\t"+SV_blat_check+"\n")
			continue
		sorted_blat_segs=sorted(SVs_blat[SV],key=lambda x:x[4])
		if (len(sorted_blat_segs)==1):
			chro,genome_match_start,genome_match_end,strand,SV_match_start,SV_match_end,SV_hit_len,blat_exact_match_len=sorted_blat_segs[0]
			if SV_hit_len<(SV_seq_len+30):
				blated_SV="single_hit"
				SV_blat_check="bad_covered"
			else:
				blated_SV="single_hit"
				SV_blat_check="focal_indel"
		else:
			clean_sorted_blat_segs=[]
			for i in range(len(sorted_blat_segs)):
				chro,genome_match_start,genome_match_end,strand,SV_match_start,SV_match_end,SV_hit_len,blat_exact_match_len=sorted_blat_segs[i]
				if SV_hit_len<left_flanking_length-20:
					continue
				elif (SV_hit_len>(SV_seq_len-8)):
					clean_sorted_blat_segs.append(sorted_blat_segs[i])
					break
				else:
					if len(clean_sorted_blat_segs)==0:
						clean_sorted_blat_segs.append(sorted_blat_segs[i])
					else:
						if (SV_match_start>clean_sorted_blat_segs[-1][4]) and (SV_match_end<=clean_sorted_blat_segs[-1][5]):
							continue
						if (SV_match_start>=clean_sorted_blat_segs[-1][4]) and (SV_match_end<clean_sorted_blat_segs[-1][5]):
							continue
						else:
							clean_sorted_blat_segs.append(sorted_blat_segs[i])
			#print(clean_sorted_blat_segs)	
			if len(clean_sorted_blat_segs)>0:
				if len(clean_sorted_blat_segs)==2:
					chrA=clean_sorted_blat_segs[0][0]
					posA=clean_sorted_blat_segs[0][2]
					strandA=clean_sorted_blat_segs[0][3]
					posA_hit_len=clean_sorted_blat_segs[0][6]
					posA_blat_exact_match_len=clean_sorted_blat_segs[0][7]
					
					chrB=clean_sorted_blat_segs[1][0]
					posB=clean_sorted_blat_segs[1][1]
					strandB=clean_sorted_blat_segs[1][3]
					posB_hit_len=clean_sorted_blat_segs[1][6]
					posB_blat_exact_match_len=clean_sorted_blat_segs[1][7]
					insert_seq=""
					posA_end_in_SV=clean_sorted_blat_segs[0][5]
					posB_start_in_SV=clean_sorted_blat_segs[1][4]
					overlapped_len=posA_end_in_SV+1-posB_start_in_SV
					if overlapped_len==0: 
						if chrA==ini_chrA and chrB==ini_chrB and strandA==ini_strandA and strandB==ini_strandB and abs(posA-ini_posA)<15 and abs(posB-ini_posB)<15: 
							if (clean_sorted_blat_segs[-1][5]-clean_sorted_blat_segs[-2][4])>0.9*SV_seq_len:
								SV_blat_check="perfect_blat"
								blated_SV=(".").join([chrA,str(posA),strandA,chrB,str(posB),strandB,insert_seq])
							else:
								SV_blat_check="perfect_blat;poor_covered"
								blated_SV=(".").join([chrA,str(posA),strandA,chrB,str(posB),strandB,insert_seq])
					else:
						if overlapped_len<0:
							insert_seq=SV_seq[posA_end_in_SV:(posB_start_in_SV-1)]
						if overlapped_len>0:
							if posA_hit_len==posA_blat_exact_match_len:
								if strandB=="+":
									posB=posB+overlapped_len
								else:
									posB=posB-overlapped_len
						if chrA==ini_chrA and chrB==ini_chrB and strandA==ini_strandA and strandB==ini_strandB and abs(posA-ini_posA)<15 and abs(posB-ini_posB)<15: 
							SV_blat_check="perfect_blat"
							blated_SV=(".").join([chrA,str(posA),strandA,chrB,str(posB),strandB,insert_seq])
				else:
					if len(clean_sorted_blat_segs)==1:
						chro,genome_match_start,genome_match_end,strand,SV_match_start,SV_match_end,SV_hit_len,blat_exact_match_len=sorted_blat_segs[0]
						if SV_hit_len<(SV_seq_len+30):
							blated_SV="single_hit"
							SV_blat_check="bad_covered"
						else:
							blated_SV="single_hit"
							SV_blat_check="focal_indel"
					else:
						blated_SV="multiple_hits"
						SV_blat_check="bad_blat"
		blat_out_fh.write(SV+"\t"+SV_seq+"\t"+blated_SV+"\t"+SV_blat_check+"\n")
blat_out_fh.close()
fn_to_delete=[SV_fa_fn,blat_fn]
#for fn in fn_to_delete:
#	rm_command="rm "+fn
#	os.system(rm_command)
