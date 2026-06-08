import os
import sys
import re
import subprocess
import json
import math
import random,string

PATH="/research/groups/xmagrp/projects/SensitiveDetection/xmagrp/code/dipseek/"

SAMBAMBA      = PATH + 'apps/sambamba'
SAMTOOLS      = PATH + 'apps/samtools'
BWA           = PATH + 'apps/bwa'
BLAT          = PATH + "apps/blat"

if len(sys.argv)!=2:
	print ("please input indel list fn")
	sys.exit()

indel_fn=sys.argv[1]

indel_padding=40 #final indel allele with 60 bp flanking seq chrA.posA.starndA.chrB.posB.starndB.flanking_A_seq.insertion_seq.flanking_B_seq  
blast_trim=8 # trim blast midlines at breakpoints to 8bp exact match '||||||||' 
flanking_length=200
qual_cutoff=20
qual_perc=0.8
MNV_len_cutoff=[2,10]

reference_folder="/research_jude/rgs01_jude/groups/xmagrp/projects/AMLRelapse/common/NeoSplice/code/yliu1/reference/"
blast_exe_folder="/research_jude/rgs01_jude/groups/xmagrp/projects/AMLRelapse/common/NeoSplice/code/yliu1/tools/"
centromere_regions={
	"hg19":{
		"chr1":[120000000,145000000],
		"chr2":[87000000,98000000],
		"chr3":[88000000,96000000],
		#"chr4":[48000000,54000000],
		"chr4":[48000000,52900000],
		"chr5":[44000000,52000000],
		"chr6":[57000000,64000000],
		"chr7":[57000000,64000000],
		"chr8":[42000000,48000000],
		"chr9":[44000000,70000000],
		"chr10":[38000000,44000000],
		"chr11":[50000000,56000000],
		"chr12":[34000000,40000000],
		"chr13":[0,20000000],
		"chr14":[0,20000000],
		"chr15":[0,20000000],
		#"chr16":[32000000,50000000],
		"chr16":[32000000,49500000],
		"chr17":[22000000,26000000],
		"chr18":[14000000,20000000],
		"chr19":[23000000,29000000],
		"chr20":[25000000,30000000],
		"chr21":[0,15000000],
		"chr22":[0,17000000],
		#"chrX":[57000000,63000000],
		"chrX":[57000000,62500000],
		"chrY":[9000000,14000000]
	},
	"hg38":{
		"chr1":[120000000,145000000],
		"chr2":[87000000,98000000],
		"chr3":[88000000,96000000],
		#"chr4":[48000000,54000000],
		"chr4":[48600000,54000000],
		"chr5":[44000000,52000000],
		"chr6":[57000000,62000000],
		"chr7":[57000000,64000000],
		#"chr8":[42000000,48000000],
		"chr8":[43000000,46400000],
		#"chr9":[40000000,70000000],
		"chr9":[42000000,70000000],
		"chr10":[38000000,43100000], # yf edited 44000000
		"chr11":[50000000,56000000],
		"chr12":[32000000,40000000],
		"chr13":[0,20000000],
		"chr14":[0,20000000],
		"chr15":[0,20000000],
		"chr16":[34000000,48000000],
		"chr17":[22000000,28000000],
		#"chr18":[14000000,23000000],
		"chr18":[14000000,22000000],
		"chr19":[23000000,29000000],
		"chr20":[25000000,32400000], # yf edited 33000000
		"chr21":[0,14000000],
		"chr22":[0,17000000],
		"chrX":[57000000,63000000],
		"chrY":[6000000,11000000]
	}
}


def get_bit_fn(bam):
	bit_fn=reference_folder+"hg19.gz"
	genome_version="hg19"
	chr_in_bam="yes"
	genome_version_command=f"{SAMTOOLS} idxstats "+bam+"|head -1"
	# genome_version_command="samtools idxstats "+bam+"|head -1"
	chr_checkresult =str(subprocess.check_output(genome_version_command,shell=True),"utf-8")
	if "chr1\t248956422" in chr_checkresult or "1\t248956422" in chr_checkresult:
		bit_fn=reference_folder+"hg38.gz"
		genome_version="hg38"

	if chr_checkresult.count("chr")==0:
		chr_in_bam="no"

	return (bit_fn,genome_version,chr_in_bam)

def get_revcmp(strs):
	cpm={"A":"T","T":"A","C":"G","G":"C","N":"N","":"","-":"-"}
	rev=strs[::-1]
	revcmp=""
	for i in rev:
		revcmp+=cpm[i]
	return revcmp
		
def make_WT_contig(indel,bam,bit_fn):
	chro,pos=indel.split(".")[0:2]
	sample=os.path.basename(bam)
	WT_contig_fn=("").join(random.choice(string.ascii_lowercase) for x in range(7))+"."+sample+"."+indel+".fa"
	start=int(pos)-200
	end=int(pos)+200
	WT_contig_command=f"{SAMTOOLS} faidx "+bit_fn+" "+chro+":"+str(start)+"-"+str(end)
	# WT_contig_command="samtools faidx "+bit_fn+" "+chro+":"+str(start)+"-"+str(end)
	WT_contig=("").join(str(subprocess.check_output(WT_contig_command,shell=True),"utf-8").strip().split("\n")[1:]).upper()
	genome_pos=[chro+"."+str(i) for i in range(start,end+1)]
	with open(WT_contig_fn,"w") as out_fh:
		out_fh.write(">"+indel+"\n")
		out_fh.write(WT_contig)
	return (WT_contig_fn,genome_pos)
	
def extract_sam(indel,bam):
	indel_lst=indel.split(".")
	sample=os.path.basename(bam)
	WT_contig_fn=sample+"."+indel+".fa"
	chro=indel_lst[0]
	pos=indel_lst[1]
	start=int(pos)-200
	end=int(pos)+200
	get_chro_command=f"{SAMTOOLS} idxstats "+bam+"|cut -f1"
	# get_chro_command="samtools idxstats "+bam+"|cut -f1"
	chr_checkresult =str(subprocess.check_output(get_chro_command,shell=True),"utf-8")
	if ("chr" in chr_checkresult) and ("chr" not in indel):
		chro="chr"+chro
	elif ("chr" not in chr_checkresult) and ("chr" in indel):
		chro=chro.replace("chr","")
	else:
		chro=chro
	sam_fn=sample+"."+indel+".sam"
	sam_command=f"{SAMTOOLS} view "+bam+" "+str(chro)+":"+str(int(pos)-200)+"-"+str(int(pos)+200)+">"+sam_fn
	# sam_command="samtools view "+bam+" "+str(chro)+":"+str(int(pos)-200)+"-"+str(int(pos)+200)+">"+sam_fn
	os.system(sam_command)
	return sam_fn


def extract_fasta(sam_fn):
	seq_data={}
	seqid_data={}
	candidate_fa_fn=sam_fn.replace(".sam",".candidate_indel.fa")
	candidate_fa_fh=open(candidate_fa_fn,'w')
	with open (sam_fn,'r') as sam_fh:
		for line in sam_fh:
			lst=line.strip().split('\t')
			seqid=lst[0]+"_"+lst[1]
			chro=lst[2]
			pos=lst[3]
			map_qual=lst[4]
			cigar=lst[5]
			paired_read_chr=lst[6]
			paired_read_pos=lst[7]
			segment_length=lst[8]
			seq=lst[9]
			qual_str=lst[10]
			ngood=0
			#if re.findall("^\d+M$",cigar):
			#	continue
			for nn_qual in qual_str:
				qscore=ord(nn_qual)-33
				if qscore >qual_cutoff:
					ngood+=1
			if ngood<qual_perc*len(qual_str):
				continue

			if seq not in seq_data:
				seq_data[seq]=[seqid]
				candidate_fa_fh.write(">"+seqid+"\n")
				candidate_fa_fh.write(seq+"\n")
				seqid_data[seqid]=seq
			else:
				seq_data[seq].append(seqid)

	candidate_fa_fh.close()
	return (seq_data,seqid_data,candidate_fa_fn)

def do_blast(WT_contig_fn,candidate_fa_fn):
	blast_fn=WT_contig_fn.replace(".fa",".blast")
	blast_command=blast_exe_folder+"blastn -task blastn -word_size 4 -evalue 0.0001 -penalty -3 -reward 2 -gapopen 6 -gapextend 2 -dust no -soft_masking false -query "+candidate_fa_fn+" -subject "+WT_contig_fn+" -outfmt 15 > "+blast_fn
	os.system(blast_command)
	return (blast_fn)

def adjust_indel_pos(qseq,hseq,midline):
	hseq_lst=list(hseq)
	qseq_lst=list(qseq)
	midline_lst=list(midline)
	for i in range(len(qseq_lst)):
		if list(qseq_lst)[i]==list(hseq_lst)[i]:
			continue
			
		if hseq_lst[i]=="-" and i <len(hseq_lst)-1:
			j=1
			while hseq_lst[i+j]=="-" and (i+j)<len(hseq_lst)-1:
				j+=1
			if (i+j)<=(len(hseq_lst)-1) and (hseq_lst[i+j]==qseq_lst[i]):
				hseq_lst[i]=hseq_lst[i+j]
				hseq_lst[i+j]="-"
				midline_lst[i]="|"
				midline_lst[i+j]=" "
			else:
				break

		if qseq_lst[i]=="-" and (i <len(qseq_lst)-1):
			j=1
			while (qseq_lst[i+j]=="-") and (i+j)<len(qseq_lst)-1:
				j+=1
			if (i+j)<=(len(qseq_lst)-1) and (qseq_lst[i+j]==hseq_lst[i]):
				qseq_lst[i]=qseq_lst[i+j]
				qseq_lst[i+j]="-"
				midline_lst[i]="|"
				midline_lst[i+j]=" "
			else:
				break
	adjust_qseq=("").join(qseq_lst)
	adjust_hseq=("").join(hseq_lst)
	adjust_midline=("").join(midline_lst)

	return (adjust_qseq,adjust_hseq,adjust_midline)

def adjust_1seg_indel_pos(seqid,qseq,hseq,midline):
	hseq_lst=list(hseq)
	qseq_lst=list(qseq)
	midline_lst=list(midline)

	for i in range(len(qseq_lst)):
		if list(qseq_lst)[i]==list(hseq_lst)[i]:
			continue
			
		if hseq_lst[i]=="-" and i <len(hseq_lst)-1:
			j=1
			while hseq_lst[i+j]=="-" and (i+j)<len(hseq_lst)-1:
				j+=1
			if (i+j)<=(len(hseq_lst)-1) and (hseq_lst[i+j]==qseq_lst[i]):
				hseq_lst[i]=hseq_lst[i+j]
				hseq_lst[i+j]="-"
				midline_lst[i]="|"
				midline_lst[i+j]=" "
			else:
				break

		if qseq_lst[i]=="-" and (i <len(qseq_lst)-1):
			j=1
			while (qseq_lst[i+j]=="-") and (i+j)<len(qseq_lst)-1:
				j+=1
			if (i+j)<=(len(qseq_lst)-1) and (qseq_lst[i+j]==hseq_lst[i]):
				qseq_lst[i]=qseq_lst[i+j]
				qseq_lst[i+j]="-"
				midline_lst[i]="|"
				midline_lst[i+j]=" "
			else:
				break
	
	qseq=("").join(qseq_lst)
	hseq=("").join(hseq_lst)
	midline=("").join(midline_lst)
	left_index=0
	for ii in range(len(qseq_lst)):
		if list(qseq_lst)[ii]==list(hseq_lst)[ii]:
			continue
		if not(hseq_lst[ii]=="-" or qseq_lst[ii]=="-"):
			continue
		left_index=ii-1
		break
	left_qseq=("").join(qseq_lst[0:(left_index+1)])
	left_hseq=("").join(hseq_lst[0:(left_index+1)])
	left_midline=("").join(midline_lst[0:(left_index+1)])
	while (left_index-blast_trim+1>=0) and (left_midline[(left_index-blast_trim+1):(left_index+1)]!="|"*blast_trim):
		left_index=left_index-1
	left_qseq=("").join(qseq_lst[0:(left_index+1)])
	left_hseq=("").join(hseq_lst[0:(left_index+1)])
	left_midline=("").join(midline_lst[0:(left_index+1)])
		
	right_index=left_index+1
	right_qseq_rev,right_hseq_rev,right_midline_rev=adjust_indel_pos(qseq[right_index:][::-1],hseq[right_index:][::-1],midline[right_index:][::-1])
	right_qseq,right_hseq,right_midline=right_qseq_rev[::-1],right_hseq_rev[::-1],right_midline_rev[::-1]

	adjust_qseq=left_qseq+right_qseq
	adjust_hseq=left_hseq+right_hseq
	adjust_midline=left_midline+right_midline

	return (adjust_qseq,adjust_hseq,adjust_midline)

def locate_1seg_indel(initial_indel,seqid,seqid_data,seq_data,pos_blast_data,genome_pos):
	read_seq=seqid_data[seqid]
	all_read_ids=seq_data[read_seq]
	initial_chro,initial_pos=initial_indel.split(".")[0:2]

	query_from=int(pos_blast_data[seqid]["query_from"])
	query_to=int(pos_blast_data[seqid]["query_to"])
	query_strand=pos_blast_data[seqid]["query_strand"]
	hit_from=int(pos_blast_data[seqid]["hit_from"])
	hit_to=int(pos_blast_data[seqid]["hit_to"])
	hit_strand=pos_blast_data[seqid]["hit_strand"]
	blast_midline=pos_blast_data[seqid]["midline"]
	blast_gaps=pos_blast_data[seqid]["gaps"]
	blast_qseq=pos_blast_data[seqid]["qseq"]
	blast_hseq=pos_blast_data[seqid]["hseq"]
	align_len=pos_blast_data[seqid]["align_len"]
	identity=pos_blast_data[seqid]["identity"]

	curated_indel_lst=[]
	all_read_ids_lst=[]

	#deal with indel in the alignment
	if "-" in blast_qseq or "-" in blast_hseq:
		if hit_strand=="+": 
			blast_qseq,blast_hseq,blast_midline=adjust_1seg_indel_pos(seqid,blast_qseq,blast_hseq,blast_midline)
		if hit_strand=="-": 
			blast_qseq_rev,blast_hseq_rev,blast_midline_rev=adjust_1seg_indel_pos(seqid,blast_qseq[::-1],blast_hseq[::-1],blast_midline[::-1])
			blast_qseq,blast_hseq,blast_midline=blast_qseq_rev[::-1],blast_hseq_rev[::-1],blast_midline_rev[::-1]

	#trim breakpoint gaps	
	indels_index_in_align=[]
	if "-" in blast_qseq or "-" in blast_hseq:
		indel_left_index="ND"
		indel_right_index="ND"
		for i in range(blast_trim,(len(blast_midline)-blast_trim)):
			indel_left_index=i-1
			if blast_midline[i]!="|" and (blast_qseq[i]=="-" or blast_hseq[i]=="-"):
				break
		while (indel_left_index-blast_trim+1)>=0 and (blast_midline[(indel_left_index-blast_trim+1):(indel_left_index+1)]!="|"*blast_trim):
			indel_left_index=indel_left_index-1

		for j in range(blast_trim,(len(blast_midline)-blast_trim)):
			indel_right_index=len(blast_midline)-j
			if blast_midline[::-1][j]!="|" and (blast_qseq[::-1][j]=="-" or blast_hseq[::-1][j]=="-"):
				break
		while (indel_right_index+blast_trim)<=len(blast_midline) and (blast_midline[indel_right_index:(indel_right_index+blast_trim)]!="|"*blast_trim):
			indel_right_index=indel_right_index+1
		if (blast_midline[(indel_left_index-blast_trim+1):(indel_left_index+1)]!="|"*blast_trim) or (blast_midline[indel_right_index:(indel_right_index+blast_trim)]!="|"*blast_trim):
			indel_left_index="ND"
			indel_right_index="ND"
		else:
			indels_index_in_align.append([indel_left_index,indel_right_index])
			if hit_strand=="+":
				for indel in indels_index_in_align:
					align_left=indel[0]
					align_right=indel[1]
					ref_allele=blast_hseq[align_left:(align_right+1)].replace("-","").upper()
					mut_allele=blast_qseq[align_left:(align_right+1)].replace("-","").upper()
					ref_pos=genome_pos[int(hit_from)+len(blast_hseq[0:(align_left+1)].replace("-",""))-1-1]
					curated_indel=ref_pos+"."+ref_allele+"."+mut_allele
					if len(ref_allele)>0 and len(mut_allele)>0:
						curated_indel_lst.append(curated_indel)
						all_read_ids_lst.append(all_read_ids)
			if hit_strand=="-":
				for indel in indels_index_in_align:
					align_left=indel[0]
					align_right=indel[1]
					ref_allele=get_revcmp(blast_hseq[align_left:(align_right+1)].replace("-","")).upper()
					mut_allele=get_revcmp(blast_qseq[align_left:(align_right+1)].replace("-","")).upper()
					ref_pos=genome_pos[int(hit_to)+len(blast_hseq[align_right:].replace("-",""))-1-1]
					curated_indel=ref_pos+"."+ref_allele+"."+mut_allele
					if len(ref_allele)>0 and len(mut_allele)>0:
						curated_indel_lst.append(curated_indel)
						all_read_ids_lst.append(all_read_ids)
	else:
		#updated line #365-406 in Nov15_2023 
		indel_left_index="ND"
		indel_right_index="ND"
		if int(align_len)-int(identity)<MNV_len_cutoff[0]:
			pass
		else:
			for i in range(blast_trim,(len(blast_midline)-blast_trim)):
				indel_left_index=i-1
				if blast_midline[i]!="|":
					break
			for j in range(blast_trim,(len(blast_midline)-blast_trim)):
				indel_right_index=len(blast_midline)-j
				if blast_midline[::-1][j]!="|":
					break
			if (blast_midline[(indel_left_index-blast_trim+1):(indel_left_index+1)]!="|"*blast_trim) or (blast_midline[indel_right_index:(indel_right_index+blast_trim)]!="|"*blast_trim):
				indel_left_index="ND"
				indel_right_index="ND"
				pass
			else:
				indels_index_in_align.append([indel_left_index,indel_right_index])
				if hit_strand=="+":
					for indel in indels_index_in_align:
						align_left=indel[0]
						align_right=indel[1]
						ref_allele=blast_hseq[align_left:(align_right+1)].upper()
						mut_allele=blast_qseq[align_left:(align_right+1)].upper()
						ref_pos=genome_pos[int(hit_from)+len(blast_hseq[0:(align_left+1)])-1-1]
						curated_indel=ref_pos+"."+ref_allele+"."+mut_allele
						if len(ref_allele)>0 and len(mut_allele)>0:
							curated_indel_lst.append(curated_indel)
							all_read_ids_lst.append(all_read_ids)
				if hit_strand=="-":
					for indel in indels_index_in_align:
						align_left=indel[0]
						align_right=indel[1]
						ref_allele=get_revcmp(blast_hseq[align_left:(align_right+1)]).upper()
						mut_allele=get_revcmp(blast_qseq[align_left:(align_right+1)]).upper()
						ref_pos=genome_pos[int(hit_to)+len(blast_hseq[align_right:])-1-1]
						curated_indel=ref_pos+"."+ref_allele+"."+mut_allele
						if len(ref_allele)>0 and len(mut_allele)>0:
							curated_indel_lst.append(curated_indel)
							all_read_ids_lst.append(all_read_ids)
			
	return(curated_indel_lst,all_read_ids_lst)

def locate_2seg_indel(initial_indel,seqid,seqid_data,seq_data,posA_blast_data,posB_blast_data,genome_pos):
	read_seq=seqid_data[seqid]
	all_read_ids=seq_data[read_seq]
	initial_chro,initial_pos=initial_indel.split(".")[0:2]

	A_query_from=int(posA_blast_data[seqid]["query_from"])
	A_query_to=int(posA_blast_data[seqid]["query_to"])
	A_query_strand=posA_blast_data[seqid]["query_strand"]
	A_hit_from=int(posA_blast_data[seqid]["hit_from"])
	A_hit_to=int(posA_blast_data[seqid]["hit_to"])
	A_hit_strand=posA_blast_data[seqid]["hit_strand"]
	A_blast_midline=posA_blast_data[seqid]["midline"]
	A_blast_gaps=posA_blast_data[seqid]["gaps"]
	A_blast_qseq=posA_blast_data[seqid]["qseq"]
	A_blast_hseq=posA_blast_data[seqid]["hseq"]

	B_query_from=int(posB_blast_data[seqid]["query_from"])
	B_query_to=int(posB_blast_data[seqid]["query_to"])
	B_query_strand=posB_blast_data[seqid]["query_strand"]
	B_hit_from=int(posB_blast_data[seqid]["hit_from"])
	B_hit_to=int(posB_blast_data[seqid]["hit_to"])
	B_hit_strand=posB_blast_data[seqid]["hit_strand"]
	B_blast_midline=posB_blast_data[seqid]["midline"]
	B_blast_gaps=posB_blast_data[seqid]["gaps"]
	B_blast_qseq=posB_blast_data[seqid]["qseq"]
	B_blast_hseq=posB_blast_data[seqid]["hseq"]
	curated_chr="ND"
	curated_posA="ND"
	curated_posB="ND"
	ref_allele="ND"
	mut_allele="ND"
	breakpoint_A="ND"
	breakpoint_B="ND"

	#deal with overlapped seq mapping issue from blast result
	if A_hit_strand==B_hit_strand=="+":
		if A_query_to-B_query_from>=0:
			overlapped_bp=A_query_to-B_query_from+1
			try:
				compare_len=min([(overlapped_bp+5),len(A_blast_midline),len(B_blast_midline)])
				if A_blast_midline[-compare_len:].count("|")<B_blast_midline[0:compare_len].count("|"):
					A_query_to=A_query_to-overlapped_bp
					i=0
					overlapped_align_bp=0
					for nn in A_blast_qseq[::-1]:
						if nn !="-":
							i+=1
						overlapped_align_bp+=1
						if i==overlapped_bp:
							break
					A_hit_to=A_hit_to-overlapped_align_bp+A_blast_hseq[-overlapped_align_bp:].count("-")
					A_blast_midline=A_blast_midline[:-overlapped_align_bp]
					A_blast_qseq=A_blast_qseq[:-overlapped_align_bp]
					A_blast_hseq=A_blast_hseq[:-overlapped_align_bp]
				else:
					B_query_from=B_query_from+overlapped_bp
					i=0
					overlapped_align_bp=0
					for nn in B_blast_qseq:
						if nn !="-":
							i+=1
						overlapped_align_bp+=1
						if i==overlapped_bp:
							break
					B_hit_from=B_hit_from+overlapped_align_bp-B_blast_hseq[0:overlapped_align_bp].count("-")
					B_blast_midline=B_blast_midline[overlapped_align_bp:]
					B_blast_qseq=B_blast_qseq[overlapped_align_bp:]
					B_blast_hseq=B_blast_hseq[overlapped_align_bp:]
			except:
				pass
				#print("error in dealing with overlapped seq mapping issue",initial_indel,read_seq)
			

	if A_hit_strand==B_hit_strand=="-":
		if B_query_to-A_query_from>=0:
			overlapped_bp=B_query_to-A_query_from+1
			try:
				compare_len=min([(overlapped_bp+5),len(A_blast_midline),len(B_blast_midline)])
				if A_blast_midline[0:compare_len].count("|")<B_blast_midline[-compare_len:].count("|"):
					A_query_from=A_query_from+overlapped_bp
					i=0
					overlapped_align_bp=0
					for nn in A_blast_qseq:
						if nn !="-":
							i+=1
						overlapped_align_bp+=1
						if i==overlapped_bp:
							break
					#A_hit_from=A_hit_from+overlapped_align_bp-A_blast_hseq[:overlapped_align_bp].count("-")
					A_hit_from=A_hit_from-overlapped_align_bp+A_blast_hseq[:overlapped_align_bp].count("-")
					A_blast_midline=A_blast_midline[overlapped_align_bp:]
					A_blast_qseq=A_blast_qseq[overlapped_align_bp:]
					A_blast_hseq=A_blast_hseq[overlapped_align_bp:]
				else:
					B_query_to=B_query_to-overlapped_bp
					i=0
					overlapped_align_bp=0
					for nn in B_blast_qseq[::-1]:
						if nn !="-":
							i+=1
						overlapped_align_bp+=1
						if i==overlapped_bp:
							break
					B_hit_to=B_hit_to+overlapped_align_bp-B_blast_hseq[-overlapped_align_bp:].count("-")
					B_blast_midline=B_blast_midline[0:-overlapped_align_bp]
					B_blast_qseq=B_blast_qseq[0:-overlapped_align_bp]
					B_blast_hseq=B_blast_hseq[0:-overlapped_align_bp]
			except:
				pass
				#print("error in dealing with overlapped seq mapping issue",initial_indel)
			
	#deal with indel in the alignment
	if "-" in A_blast_qseq or "-" in A_blast_hseq:
		if A_hit_strand==B_hit_strand=="+": 
			A_blast_qseq,A_blast_hseq,A_blast_midline=adjust_indel_pos(A_blast_qseq,A_blast_hseq,A_blast_midline)
		if A_hit_strand==B_hit_strand=="-": 
			A_blast_qseq_rev,A_blast_hseq_rev,A_blast_midline_rev=adjust_indel_pos(A_blast_qseq[::-1],A_blast_hseq[::-1],A_blast_midline[::-1])
			A_blast_qseq,A_blast_hseq,A_blast_midline=A_blast_qseq_rev[::-1],A_blast_hseq_rev[::-1],A_blast_midline_rev[::-1]

	if "-" in B_blast_qseq or "-" in B_blast_hseq:
		if A_hit_strand==B_hit_strand=="+": 
			B_blast_qseq_rev,B_blast_hseq_rev,B_blast_midline_rev=adjust_indel_pos(B_blast_qseq[::-1],B_blast_hseq[::-1],B_blast_midline[::-1])
			B_blast_qseq,B_blast_hseq,B_blast_midline=B_blast_qseq_rev[::-1],B_blast_hseq_rev[::-1],B_blast_midline_rev[::-1]
		if A_hit_strand== B_hit_strand=="-": 
			B_blast_qseq,B_blast_hseq,B_blast_midline=adjust_indel_pos(B_blast_qseq,B_blast_hseq,B_blast_midline)

	#trim breakpoint gaps	
	if A_hit_strand==B_hit_strand=="+":
		try:
			if A_blast_midline[-blast_trim:]!=blast_trim*"|":
				drop_A_len=A_blast_midline[::-1].index(blast_trim*"|")
				A_query_to=A_query_to-drop_A_len+A_blast_qseq[-drop_A_len:].count("-")
				A_hit_to=A_hit_to-drop_A_len+A_blast_hseq[-drop_A_len:].count("-")
				A_blast_midline=A_blast_midline[:-drop_A_len]
				A_blast_qseq=A_blast_qseq[:-drop_A_len]
				A_blast_hseq=A_blast_hseq[:-drop_A_len]
		except:
			pass
			#print("error in trimming breakpoint gaps",seqid,initial_indel,read_seq)
		

		try:
			if B_blast_midline[0:blast_trim]!=blast_trim*"|":
				drop_B_len=B_blast_midline.index(blast_trim*"|")
				B_query_from=B_query_from+drop_B_len-B_blast_qseq[0:drop_B_len].count("-")
				B_hit_from=B_hit_from+drop_B_len-B_blast_hseq[0:drop_B_len].count("-")
				B_blast_midline=B_blast_midline[drop_B_len:]
				B_blast_qseq=B_blast_qseq[drop_B_len:]
				B_blast_hseq=B_blast_hseq[drop_B_len:]
		except:
			pass
			#print("error in trimming breakpoint gaps",initial_indel,read_seq,)

		if A_blast_midline.count(" ")<=5 or B_blast_midline.count(" ")<=5:
			try:
				if A_query_to<B_query_from:
					breakpoint_A=genome_pos[A_hit_to-1]
					breakpoint_B=genome_pos[B_hit_from-1]
					mut_allele=read_seq[(A_query_to-1):B_query_from]
					curated_posA=breakpoint_A.split(".")[1]
					curated_posB=breakpoint_B.split(".")[1]
					ref_allele_command=f"{SAMTOOLS} faidx "+bit_fn+" "+initial_chro+":"+str(curated_posA)+"-"+str(curated_posB)
					# ref_allele_command="samtools faidx "+bit_fn+" "+initial_chro+":"+str(curated_posA)+"-"+str(curated_posB)
					ref_allele=("").join(str(subprocess.check_output(ref_allele_command,shell=True),"utf-8").strip().split("\n")[1:]).upper()

			except:
				pass
				#print("error in locating breakpoint",initial_indel)
		if A_blast_midline[-blast_trim:]!=blast_trim*"|" or B_blast_midline[0:blast_trim]!=blast_trim*"|":
			ref_allele="ND"
	if A_hit_strand== B_hit_strand=="-":
		try:
			if A_blast_midline[0:blast_trim]!=blast_trim*"|":
				drop_A_len=A_blast_midline.index(blast_trim*"|")
				A_query_from=A_query_from+drop_A_len-A_blast_qseq[0:drop_A_len].count("-")
				#A_hit_from=A_hit_from+drop_A_len-A_blast_hseq[0:drop_A_len].count("-")
				A_hit_from=A_hit_from-drop_A_len+A_blast_hseq[0:drop_A_len].count("-")
				A_blast_midline=A_blast_midline[drop_A_len:]
				A_blast_qseq=A_blast_qseq[drop_A_len:]
				A_blast_hseq=A_blast_hseq[drop_A_len:]
		except:
			pass
			#print("error in trimming breakpoint gaps",initial_indel,read_seq)

		try:
			if B_blast_midline[-blast_trim:]!=blast_trim*"|":
				drop_B_len=B_blast_midline[::-1].index(blast_trim*"|")
				B_query_to=B_query_to-drop_B_len+B_blast_qseq[-drop_B_len:].count("-")
				B_hit_to=B_hit_to-drop_B_len+B_blast_hseq[-drop_B_len:].count("-")
				B_blast_midline=B_blast_midline[:-drop_B_len]
				B_blast_qseq=B_blast_qseq[:-drop_B_len]
				B_blast_hseq=B_blast_hseq[:-drop_B_len]
		except:
			pass
			#print("error in trimming breakpoint gaps",initial_indel,read_seq)

		if A_blast_midline.count(" ")<=5 and B_blast_midline.count(" ")<=5:
			try:
				if B_query_to<A_query_from:
					breakpoint_A=genome_pos[A_hit_from-1]
					breakpoint_B=genome_pos[B_hit_to-1]
					mut_allele=get_revcmp(read_seq[(B_query_to-1):A_query_from])
					curated_posA=breakpoint_A.split(".")[1]
					curated_posB=breakpoint_B.split(".")[1]
					ref_allele_command=f"{SAMTOOLS} faidx "+bit_fn+" "+initial_chro+":"+str(curated_posA)+"-"+str(curated_posB)
					# ref_allele_command="samtools faidx "+bit_fn+" "+initial_chro+":"+str(curated_posA)+"-"+str(curated_posB)
					ref_allele=("").join(str(subprocess.check_output(ref_allele_command,shell=True),"utf-8").strip().split("\n")[1:]).upper()

			except:	
				pass
				#print("error in locating breakpoint",initial_indel)
		if A_blast_midline[0:blast_trim]!=blast_trim*"|" or B_blast_midline[-blast_trim:]!=blast_trim*"|":
			ref_allele="ND"

	curated_indel=initial_chro+"."+curated_posA+"."+ref_allele+"."+mut_allele
	if "ND" in curated_indel:
		all_read_ids="ND"
	return(curated_indel,all_read_ids)

def parse_blast(blast_fn,seqid_data,seq_data,initial_indel,genome_pos):
	strand_transform={"Plus":"+","Minus":"-"}
	curated_indel="ND"
	curated_indels_data={}
	initial_chr,initial_pos=initial_indel.split(".")[0:2]
	curated_indel_lst=[]
	all_read_ids_lst=[]
	blast_data={}
	try:
		blast_json=json.load(open(blast_fn))["BlastOutput2"]
		for i in range(0,len(blast_json)):
			blast_report= blast_json[i]["report"]
			blast_result= blast_report["results"]
			bl2seq=blast_result["bl2seq"][0]
			seqid=bl2seq["query_title"]
			if seqid not in blast_data:
				blast_data[seqid]={}
			query_len=bl2seq["query_len"]

			hsps="ND"
			align_len="ND"
			identity="ND"
			gaps="ND"
			query_from="ND"
			query_to="ND"
			query_strand="ND"
			hit_from="ND"
			hit_to="ND"
			hit_strand="ND"
			qseq="ND"
			hseq="ND"
			midline="ND"
			bit_score="ND"
			score="ND"
			evalue="ND"

			if ("hits" in bl2seq) and len(bl2seq["hits"])>0:
				hits=bl2seq["hits"]
				hsps=hits[0]["hsps"]
				if len(hsps)>0:
					for j in range(len(hsps)):
						align_len=hsps[j]["align_len"]
						identity=hsps[j]["identity"]
						gaps=hsps[j]["gaps"]
						query_from=hsps[j]["query_from"]
						query_to=hsps[j]["query_to"]
						query_strand=hsps[j]["query_strand"]
						hit_from=hsps[j]["hit_from"]
						hit_to=hsps[j]["hit_to"]
						hit_strand=hsps[j]["hit_strand"]
						qseq=hsps[j]["qseq"]
						hseq=hsps[j]["hseq"]
						midline=hsps[j]["midline"]
						bit_score=hsps[j]["bit_score"]
						score=hsps[j]["score"]
						evalue=hsps[j]["evalue"]
						blast_data[seqid][j]={
							"initial_indel":initial_indel,
							"align_len":align_len,
							"identity":identity,
							"gaps":gaps,
							"query_from":query_from,
							"query_to":query_to,
							"query_strand":strand_transform[query_strand],
							"hit_from":hit_from,
							"hit_to":hit_to,
							"hit_strand":strand_transform[hit_strand],
							"qseq":qseq,
							"hseq":hseq,
							"midline":midline,
							"bit_score":bit_score,
							"score":score,
							"evalue":evalue
						}
						#updated by adding 2 lines below
						if (int(query_from)<5) and (query_to>len(seqid_data[seqid])-5):
							break
	except:
		print(initial_indel,bam)

	for seqid in blast_data:
		hsp_nums=blast_data[seqid].keys()
		seq=seqid_data[seqid]
		hsp_range_lst=[]
		if len(blast_data[seqid])==1:
			if int(blast_data[seqid][0]["query_from"])>5 or int(blast_data[seqid][0]["query_to"])<len(seq)-5:
				continue
			hsp_range_lst.append(0)
		if len(blast_data[seqid])>1:
			for ii in blast_data[seqid]:
				ii_query_from=int(blast_data[seqid][ii]["query_from"])
				ii_query_to=int(blast_data[seqid][ii]["query_to"])
				if ii_query_from<5 and ii_query_to>len(seq)-5:
					hsp_range_lst.append(ii)
					continue
				ii_hit_from=int(blast_data[seqid][ii]["hit_from"])
				ii_hit_to=int(blast_data[seqid][ii]["hit_to"])
				ii_query_strand=blast_data[seqid][ii]["query_strand"]
				ii_hit_strand=blast_data[seqid][ii]["hit_strand"]
				for jj in blast_data[seqid]:
					if int(ii)>=int(jj):
						continue
					jj_query_from=int(blast_data[seqid][jj]["query_from"])
					jj_query_to=int(blast_data[seqid][jj]["query_to"])
					jj_hit_from=int(blast_data[seqid][jj]["hit_from"])
					jj_hit_to=int(blast_data[seqid][jj]["hit_to"])
					jj_query_strand=blast_data[seqid][jj]["query_strand"]
					jj_hit_strand=blast_data[seqid][jj]["hit_strand"]
					if ii_query_from<=jj_query_from and ii_query_to>=jj_query_to:
						continue
					if ii_hit_strand!=jj_hit_strand:
						continue
					if ii_query_to>len(seq)-10 and jj_query_from<10:
						hsp_range_lst+=[jj,ii]
						continue
					if ii_query_from<10 and jj_query_to>len(seq)-10:
						hsp_range_lst+=[ii,jj]
						continue
		#if "A00214:50:H5NKVDSXX:1:2367:2193:3474_99" in seqid:
		#	print ("A00214:50:H5NKVDSXX:1:2367:2193:3474_99",hsp_range_lst)
		if len(hsp_range_lst)==1:
			hsp_order=hsp_range_lst[0]
			pos_blast_data={
				seqid:{
				"initial_indel":initial_indel,
				"align_len":blast_data[seqid][hsp_order]["align_len"],
				"identity":blast_data[seqid][hsp_order]["identity"],
				"gaps":blast_data[seqid][hsp_order]["gaps"],
				"query_from":blast_data[seqid][hsp_order]["query_from"],
				"query_to":blast_data[seqid][hsp_order]["query_to"],
				"query_strand":blast_data[seqid][hsp_order]["query_strand"],
				"hit_from":blast_data[seqid][hsp_order]["hit_from"],
				"hit_to":blast_data[seqid][hsp_order]["hit_to"],
				"hit_strand":blast_data[seqid][hsp_order]["hit_strand"],
				"qseq":blast_data[seqid][hsp_order]["qseq"],
				"hseq":blast_data[seqid][hsp_order]["hseq"],
				"midline":blast_data[seqid][hsp_order]["midline"],
				"bit_score":blast_data[seqid][hsp_order]["bit_score"],
				"score":blast_data[seqid][hsp_order]["score"],
				"evalue":blast_data[seqid][hsp_order]["evalue"]
				}
			}
			
			curated_indel_lst,all_read_ids_lst=locate_1seg_indel(initial_indel,seqid,seqid_data,seq_data,pos_blast_data,genome_pos)
			for i in range(len(curated_indel_lst)):
				curated_indel=curated_indel_lst[i]
				all_read_ids=all_read_ids_lst[i]
				if curated_indel not in curated_indels_data:
					curated_indels_data[curated_indel]=all_read_ids
				else:
					curated_indels_data[curated_indel]+=all_read_ids

		if len(hsp_range_lst)==2:
			hsp_left_order=hsp_range_lst[0]
			hsp_right_order=hsp_range_lst[1]
			if blast_data[seqid][hsp_left_order]["hit_strand"]=="+":
				posA_blast_data={
					seqid:{
						"initial_indel":initial_indel,
						"align_len":blast_data[seqid][hsp_left_order]["align_len"],
						"identity":blast_data[seqid][hsp_left_order]["identity"],
						"gaps":blast_data[seqid][hsp_left_order]["gaps"],
						"query_from":blast_data[seqid][hsp_left_order]["query_from"],
						"query_to":blast_data[seqid][hsp_left_order]["query_to"],
						"query_strand":blast_data[seqid][hsp_left_order]["query_strand"],
						"hit_from":blast_data[seqid][hsp_left_order]["hit_from"],
						"hit_to":blast_data[seqid][hsp_left_order]["hit_to"],
						"hit_strand":blast_data[seqid][hsp_left_order]["hit_strand"],
						"qseq":blast_data[seqid][hsp_left_order]["qseq"],
						"hseq":blast_data[seqid][hsp_left_order]["hseq"],
						"midline":blast_data[seqid][hsp_left_order]["midline"],
						"bit_score":blast_data[seqid][hsp_left_order]["bit_score"],
						"score":blast_data[seqid][hsp_left_order]["score"],
						"evalue":blast_data[seqid][hsp_left_order]["evalue"]
					}
				}
				posB_blast_data={
					seqid:{
						"initial_indel":initial_indel,
						"align_len":blast_data[seqid][hsp_right_order]["align_len"],
						"identity":blast_data[seqid][hsp_right_order]["identity"],
						"gaps":blast_data[seqid][hsp_right_order]["gaps"],
						"query_from":blast_data[seqid][hsp_right_order]["query_from"],
						"query_to":blast_data[seqid][hsp_right_order]["query_to"],
						"query_strand":blast_data[seqid][hsp_right_order]["query_strand"],
						"hit_from":blast_data[seqid][hsp_right_order]["hit_from"],
						"hit_to":blast_data[seqid][hsp_right_order]["hit_to"],
						"hit_strand":blast_data[seqid][hsp_right_order]["hit_strand"],
						"qseq":blast_data[seqid][hsp_right_order]["qseq"],
						"hseq":blast_data[seqid][hsp_right_order]["hseq"],
						"midline":blast_data[seqid][hsp_right_order]["midline"],
						"bit_score":blast_data[seqid][hsp_right_order]["bit_score"],
						"score":blast_data[seqid][hsp_right_order]["score"],
						"evalue":blast_data[seqid][hsp_right_order]["evalue"]
					}
				}
				curated_indel,all_read_ids=locate_2seg_indel(initial_indel,seqid,seqid_data,seq_data,posA_blast_data,posB_blast_data,genome_pos)
				if curated_indel not in curated_indels_data:
					curated_indels_data[curated_indel]=all_read_ids
				else:
					curated_indels_data[curated_indel]+=all_read_ids
			if blast_data[seqid][hsp_left_order]["hit_strand"]==blast_data[seqid][hsp_right_order]["hit_strand"]=="-":
				posB_blast_data={
					seqid:{
						"initial_indel":initial_indel,
						"align_len":blast_data[seqid][hsp_left_order]["align_len"],
						"identity":blast_data[seqid][hsp_left_order]["identity"],
						"gaps":blast_data[seqid][hsp_left_order]["gaps"],
						"query_from":blast_data[seqid][hsp_left_order]["query_from"],
						"query_to":blast_data[seqid][hsp_left_order]["query_to"],
						"query_strand":blast_data[seqid][hsp_left_order]["query_strand"],
						"hit_from":blast_data[seqid][hsp_left_order]["hit_from"],
						"hit_to":blast_data[seqid][hsp_left_order]["hit_to"],
						"hit_strand":blast_data[seqid][hsp_left_order]["hit_strand"],
						"qseq":blast_data[seqid][hsp_left_order]["qseq"],
						"hseq":blast_data[seqid][hsp_left_order]["hseq"],
						"midline":blast_data[seqid][hsp_left_order]["midline"],
						"bit_score":blast_data[seqid][hsp_left_order]["bit_score"],
						"score":blast_data[seqid][hsp_left_order]["score"],
						"evalue":blast_data[seqid][hsp_left_order]["evalue"]
					}
				}
				posA_blast_data={
					seqid:{
						"initial_indel":initial_indel,
						"align_len":blast_data[seqid][hsp_right_order]["align_len"],
						"identity":blast_data[seqid][hsp_right_order]["identity"],
						"gaps":blast_data[seqid][hsp_right_order]["gaps"],
						"query_from":blast_data[seqid][hsp_right_order]["query_from"],
						"query_to":blast_data[seqid][hsp_right_order]["query_to"],
						"query_strand":blast_data[seqid][hsp_right_order]["query_strand"],
						"hit_from":blast_data[seqid][hsp_right_order]["hit_from"],
						"hit_to":blast_data[seqid][hsp_right_order]["hit_to"],
						"hit_strand":blast_data[seqid][hsp_right_order]["hit_strand"],
						"qseq":blast_data[seqid][hsp_right_order]["qseq"],
						"hseq":blast_data[seqid][hsp_right_order]["hseq"],
						"midline":blast_data[seqid][hsp_right_order]["midline"],
						"bit_score":blast_data[seqid][hsp_right_order]["bit_score"],
						"score":blast_data[seqid][hsp_right_order]["score"],
						"evalue":blast_data[seqid][hsp_right_order]["evalue"]
					}
				}
				curated_indel,all_read_ids=locate_2seg_indel(initial_indel,seqid,seqid_data,seq_data,posA_blast_data,posB_blast_data,genome_pos)
				if curated_indel not in curated_indels_data:
					curated_indels_data[curated_indel]=all_read_ids
				else:
					curated_indels_data[curated_indel]+=all_read_ids
	return (curated_indels_data)


def get_indel_padding(indel):
	chro,pos,ref_allele,mut_allele=indel.split(".")
	ref_allele_len=len(ref_allele)
	mut_allele_len=len(mut_allele)
	left_padding_command=f"{SAMTOOLS} faidx "+bit_fn+" "+chro+":"+str(int(pos)-indel_padding)+"-"+str(int(pos)-1)
	# left_padding_command="samtools faidx "+bit_fn+" "+chro+":"+str(int(pos)-indel_padding)+"-"+str(int(pos)-1)
	left_padding_seq=("").join(str(subprocess.check_output(left_padding_command,shell=True),"utf-8").strip().split("\n")[1:]).upper()
	right_padding_command=f"{SAMTOOLS} faidx "+bit_fn+" "+chro+":"+str(int(pos)+ref_allele_len)+"-"+str(int(pos)+ref_allele_len+indel_padding-1)
	# right_padding_command="samtools faidx "+bit_fn+" "+chro+":"+str(int(pos)+ref_allele_len)+"-"+str(int(pos)+ref_allele_len+indel_padding-1)
	right_padding_seq=("").join(str(subprocess.check_output(right_padding_command,shell=True),"utf-8").strip().split("\n")[1:]).upper()
	indel_with_padding=chro+"."+pos+"."+left_padding_seq+"."+ref_allele+"."+mut_allele+"."+right_padding_seq
	return indel_with_padding
		
def get_depth(indel):
	chro,pos=indel.split(".")[0:2]
	if chr_in_bam =="no":
		depth_command=f"{SAMTOOLS} depth -r "+chro.replace("chr","")+":"+pos+"-"+pos+" "+bam+"|cut -f3"
		# depth_command="samtools depth -r "+chro.replace("chr","")+":"+pos+"-"+pos+" "+bam+"|cut -f3"
	else:
		depth_command=f"{SAMTOOLS} depth -r "+chro+":"+pos+"-"+pos+" "+bam+"|cut -f3"
		# depth_command="samtools depth -r "+chro+":"+pos+"-"+pos+" "+bam+"|cut -f3"
	try:
		depth_indel=str(subprocess.check_output(depth_command,shell=True),"utf-8").strip()
		if depth_indel=="":
			depth_indel="0"
	except:
		depth_indel="0"
	return (depth_indel)


out_indel_fn=indel_fn+"_curated_indel"
out_indel_fh=open(out_indel_fn,"w")
out_indel_fh.write("sample\tinitial_indel\tcurated_indel\tgenome_version\tcurated_indel_reads\tdepth\tcurated_indel_with_padding\n")

out_indel_reads_fn=indel_fn+"_curated_indel_reads.txt"
out_indel_reads_fh=open(out_indel_reads_fn,"w")
out_indel_reads_fh.write("curated_indel_reads\tread_id\tbam\n")
#with open ("indel_test.txt","r") as indel_fh:
with open (indel_fn,"r") as indel_fh:
	for line in indel_fh:
		lst=line.strip().split("\t")
		initial_indel=lst[0]
		initial_chro,initial_pos=initial_indel.strip().split(".")[0:2]
		bam=lst[1]
		sample=os.path.basename(bam).replace(".bam", "")
		bit_fn,genome_version,chr_in_bam=get_bit_fn(bam)
		if not (initial_chro in centromere_regions[genome_version]):
			out_indel="%s\t%s\t%s\t%s\t%s\t%s\t%s\n"%(sample,initial_indel,"badchr","ND","ND","ND","ND")
			out_indel_fh.write(out_indel)
			continue
		if centromere_regions[genome_version][initial_chro][0]<int(initial_pos)<centromere_regions[genome_version][initial_chro][1]:
			out_indel="%s\t%s\t%s\t%s\t%s\t%s\t%s\n"%(sample,initial_indel,"centromere","ND","ND","ND","ND")
			out_indel_fh.write(out_indel)
			continue
		WT_contig_fn,genome_pos=make_WT_contig(initial_indel,bam,bit_fn)
		sam_fn=extract_sam(initial_indel,bam)
		seq_data,seqid_data,candidate_fa_fn=extract_fasta(sam_fn)
		blast_fn=do_blast(WT_contig_fn,candidate_fa_fn)
		curated_indels_data=parse_blast(blast_fn,seqid_data,seq_data,initial_indel,genome_pos)
		initial_indel_depth=get_depth(initial_indel)
		initial_indel_padding=get_indel_padding(initial_indel)
		files_to_delete=[sam_fn,candidate_fa_fn,WT_contig_fn,blast_fn]
		if len(curated_indels_data)==0:
			out_indel="%s\t%s\t%s\t%s\t%s\t%s\t%s\n"%(sample,initial_indel,"ND","ND","ND","ND","ND")
			out_indel_fh.write(out_indel)
			for fn_to_delete in set(files_to_delete):
				os.remove(fn_to_delete)
			continue
		for indel in curated_indels_data:
			ref_allel,mut_allel=indel.split(".")[-2:]
			#if ("N" in ref_allel) or ("N" in mut_allel) or ("ND" in indel) or (len(ref_allel)>20) or (len(mut_allel)>20):
			if ("N" in ref_allel) or ("N" in mut_allel) or ("ND" in indel):
				out_indel="%s\t%s\t%s\t%s\t%s\t%s\t%s\n"%(sample,initial_indel,"ND","ND","ND","ND","ND")
			elif len(ref_allel)==len(mut_allel) and len(ref_allel)>MNV_len_cutoff[1]:
				out_indel="%s\t%s\t%s\t%s\t%s\t%s\t%s\n"%(sample,initial_indel,"ND","ND","ND","ND","ND")
			else:
				curated_indel_depth=get_depth(indel)
				curated_indel_padding=get_indel_padding(indel)
				curated_indel_reads=str(len(set(curated_indels_data[indel])))
				if int(curated_indel_reads)>=0:
					out_indel="%s\t%s\t%s\t%s\t%s\t%s\t%s\n"%(sample,initial_indel,indel,genome_version,curated_indel_reads,curated_indel_depth,curated_indel_padding)
				else:
					out_indel="%s\t%s\t%s\t%s\t%s\t%s\t%s\n"%(sample,initial_indel,"ND","ND","ND","ND","ND")
			if "ND" in out_indel:
				continue
			out_indel_fh.write(out_indel)
			for read_id in set(curated_indels_data[indel]):
				out_indel_reads_fh.write(indel+"\t"+read_id+"\t"+bam+"\n")

		for fn_to_delete in set(files_to_delete):
			os.remove(fn_to_delete)

out_indel_fh.close()
out_indel_reads_fh.close()


