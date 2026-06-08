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
	print ("please input SV list fn")
	sys.exit()

SV_fn=sys.argv[1]

SV_padding=60 #final SV allele with 60 bp flanking seq chrA.posA.starndA.chrB.posB.starndB.flanking_A_seq.insertion_seq.flanking_B_seq  
blast_trim=8 # trim blast midlines at breakpoints to 8bp exact match '||||||||' 
flanking_length=200
qual_cutoff=20
qual_perc=0.9
#sam_reads_limit=50000
sam_reads_limit=500000

reference_folder="/research_jude/rgs01_jude/groups/xmagrp/projects/AMLRelapse/common/NeoSplice/code/yliu1/reference/"
blast_exe_folder="/research_jude/rgs01_jude/groups/xmagrp/projects/AMLRelapse/common/NeoSplice/code/yliu1/tools/"
centromere_regions={
	"hg19":{
		"chr1":[120000000,145000000],
		"chr2":[87000000,98000000],
		"chr3":[88000000,96000000],
		#"chr4":[48000000,54000000],
		"chr4":[48000000,52900000],
		#"chr5":[44000000,52000000],
		"chr5":[44000000,51600000],
		"chr6":[57000000,64000000],
		"chr7":[57000000,64000000],
		"chr8":[42500000,48000000],
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
		"chrY":[9000000,14000000],
		"chrM":[0,0]	#added Oct2,2024
	},
	"hg38":{
		"chr1":[120000000,145000000],
		"chr2":[87000000,98000000],
		"chr3":[88000000,96000000],
		#"chr4":[48000000,54000000],
		"chr4":[48600000,54000000],
		#"chr5":[44000000,52000000],
		"chr5":[44000000,51600000],
		"chr6":[57000000,62000000],
		"chr7":[57000000,64000000],
		#"chr8":[42000000,48000000],
		"chr8":[43000000,46400000],
		#"chr9":[40000000,70000000],
		"chr9":[42000000,70000000],
		"chr10":[38000000,44000000],
		"chr11":[50000000,56000000],
		"chr12":[32000000,40000000],
		"chr13":[0,20000000],
		"chr14":[0,20000000],
		"chr15":[0,20000000],
		"chr16":[34000000,48000000],
		"chr17":[22000000,28000000],
		#"chr18":[14000000,23000000],
		"chr18":[14000000,19000000],
		"chr19":[23000000,29000000],
		#"chr20":[25000000,33000000],
		"chr20":[26000000,32000000],
		"chr21":[0,14000000],
		"chr22":[0,17000000],
		"chrX":[57000000,63000000],
		"chrY":[6000000,11000000],
		"chrM":[0,0]	#added Oct2,2024
	}
}


def get_bit_fn(bam):
	bit_fn=reference_folder+"hg19.gz"
	genome_version="hg19"
	chr_in_bam="yes"
	genome_version_command=f"{SAMTOOLS} idxstats "+bam+"|head -1"
	#genome_version_command="samtools idxstats "+bam+"|head -1"
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
def adjust_SV(SV):
	if SV.count(".")==5:
		SV=SV+"."
	cmp_strand={"+":"-","-":"+"}
	SV_lst=SV.split(".")
	chrA=SV_lst[0]
	if "chr" not in chrA:
		chrA="chr"+chrA
	if chrA=="chrX":
		chrA="chr23"
	if chrA=="chrY":
		chrA="chr24"
	if chrA=="chrM": #added Oct2,2024
		chrA="chr25"
	posA=SV_lst[1]
	strandA=SV_lst[2]
	chrB=SV_lst[3]
	if "chr" not in chrB:
		chrB="chr"+chrB
	if chrB=="chrX":
		chrB="chr23"
	if chrB=="chrY":
		chrB="chr24"
	if chrB=="chrM": #added Oct2,2024
		chrB="chr25"
	posB=SV_lst[4]
	strandB=SV_lst[5]
	insertion=SV_lst[6]

	adjusted_chrA=chrA
	adjusted_posA=posA
	adjusted_strandA=strandA
	adjusted_chrB=chrB
	adjusted_posB=posB
	adjusted_strandB=strandB
	adjusted_insertion=insertion

	if (int(chrA.replace("chr",""))>int(chrB.replace("chr",""))) or (int(chrA.replace("chr",""))==int(chrB.replace("chr","")) and int(posA)>int(posB)):
		adjusted_chrA=chrB
		adjusted_posA=posB
		adjusted_strandA=cmp_strand[strandB]
		adjusted_chrB=chrA
		adjusted_posB=posA
		adjusted_strandB=cmp_strand[strandA]
		adjusted_insertion=get_revcmp(insertion)
	if adjusted_chrA=="chr23":
		adjusted_chrA="chrX"
	if adjusted_chrA=="chr24":
		adjusted_chrA="chrY"
	if adjusted_chrA=="chr25":
		adjusted_chrA="chrM"	#added Oct2,2024
	if adjusted_chrB=="chr23":
		adjusted_chrB="chrX"
	if adjusted_chrB=="chr24":
		adjusted_chrB="chrY"
	if adjusted_chrB=="chr25":
		adjusted_chrB="chrM"	#added Oct2,2024

	adjusted_SV="%s.%s.%s.%s.%s.%s.%s"%(adjusted_chrA,adjusted_posA,adjusted_strandA,adjusted_chrB,adjusted_posB,adjusted_strandB,adjusted_insertion)
	return adjusted_SV
		
def make_focalSV_ref(SV,bam,bit_fn):
	SV_lst=SV.split(".")
	chrA=SV_lst[0]
	posA=SV_lst[1]
	strandA=SV_lst[2]
	chrB=SV_lst[3]
	posB=SV_lst[4]
	strandB=SV_lst[5]
	insertion=SV_lst[6]
	sample=os.path.basename(bam)
	if "chr" not in SV:
		chrA="chr"+chrA
		chrB="chr"+chrB
	startA=int(posA)-200
	endA=int(posA)+200
	startB=int(posB)-200
	endB=int(posB)+200
	posA_ref_command=f"{SAMTOOLS} faidx "+bit_fn+" "+chrA+":"+str(startA)+"-"+str(endA)
	posB_ref_command=f"{SAMTOOLS} faidx "+bit_fn+" "+chrB+":"+str(startB)+"-"+str(endB)
	# posA_ref_command="samtools faidx "+bit_fn+" "+chrA+":"+str(startA)+"-"+str(endA)
	# posB_ref_command="samtools faidx "+bit_fn+" "+chrB+":"+str(startB)+"-"+str(endB)
	posA_ref_seq=("").join(str(subprocess.check_output(posA_ref_command,shell=True),"utf-8").strip().split("\n")[1:]).upper()
	posB_ref_seq=("").join(str(subprocess.check_output(posB_ref_command,shell=True),"utf-8").strip().split("\n")[1:]).upper()
	if strandA=="+" and strandB=="+":
		SV_ref_seq=posA_ref_seq+posB_ref_seq
		genome_pos=[chrA+"."+str(i) for i in range(startA,endA+1)]+[chrB+"."+str(i) for i in range(startB,endB+1)]
	if strandA=="+" and strandB=="-":
		posB_ref_seq=get_revcmp(posB_ref_seq)
		SV_ref_seq=posA_ref_seq+posB_ref_seq
		genome_pos=[chrA+"."+str(i) for i in range(startA,endA+1)]+[chrB+"."+str(i) for i in range(startB,endB+1)][::-1]
	if strandA=="-" and strandB=="+":
		posA_ref_seq=get_revcmp(posA_ref_seq)
		SV_ref_seq=posA_ref_seq+posB_ref_seq
		genome_pos=[chrA+"."+str(i) for i in range(startA,endA+1)][::-1]+[chrB+"."+str(i) for i in range(startB,endB+1)]
	if strandA=="-" and strandB=="-":
		posA_ref_seq=get_revcmp(posA_ref_seq)
		posB_ref_seq=get_revcmp(posB_ref_seq)
		SV_ref_seq=posA_ref_seq+posB_ref_seq
		genome_pos=[chrA+"."+str(i) for i in range(startA,endA+1)][::-1]+[chrB+"."+str(i) for i in range(startB,endB+1)][::-1]

	SV_ref_fn=('').join(random.choice(string.ascii_lowercase) for x in range(7))+"."+sample+"."+SV+"ref.fa"
	with open(SV_ref_fn,"w")as SV_ref_fh:
		SV_ref_fh.write(">"+sample+"."+SV+"\n")
		SV_ref_fh.write(SV_ref_seq+"\n")
	return (SV_ref_fn,genome_pos)

	
		
def make_ref(SV,bam,bit_fn):
	SV_lst=SV.split(".")
	chrA=SV_lst[0]
	posA=SV_lst[1]
	strandA=SV_lst[2]
	chrB=SV_lst[3]
	posB=SV_lst[4]
	strandB=SV_lst[5]
	insertion=SV_lst[6]
	sample=os.path.basename(bam)
	if "chr" not in SV:
		chrA="chr"+chrA
		chrB="chr"+chrB
	startA=int(posA)-200
	endA=int(posA)+200
	startB=int(posB)-200
	endB=int(posB)+200
	distance=abs(int(posB)-int(posA))
	if chrA==chrB and 300<distance<400:
		if strandA==strandB:
			#deal with small dels
			if strandA=="+" and (int(posB)-int(posA))>0:
				endA=int(posA)+math.floor(0.5*distance)
				startB=int(posB)-math.floor(0.5*distance)
			if strandA=="-" and (int(posA)-int(posB))>0:
				startA=int(posA)-math.floor(0.5*distance)
				endB=int(posB)+math.floor(0.5*distance)
			#deal with small ins
			if strandA=="+" and (int(posA)-int(posB))>0:
				startA=int(posA)-math.floor(0.5*distance)
				endB=int(posB)+math.floor(0.5*distance)
				#print(startA,endB)
			if strandA=="-" and (int(posB)-int(posA))>0:
				endA=int(posA)+math.floor(0.5*distance)
				startB=int(posB)-math.floor(0.5*distance)
		else:
			#deal with inv
			if strandA=="+" and (int(posB)-int(posA))>0:
				endA=int(posA)+math.floor(0.5*distance)
				startB=int(posB)-math.floor(0.5*distance)
			if strandA=="+" and (int(posA)-int(posB))>0:
				startA=int(posA)-math.floor(0.5*distance)
				endB=int(posB)+math.floor(0.5*distance)
			if strandA=="-" and (int(posB)-int(posA))>0:
				startA=int(posA)-math.floor(0.5*distance)
				endB=int(posB)+math.floor(0.5*distance)
			if strandA=="-" and (int(posA)-int(posB))>0:
				endA=int(posA)+math.floor(0.5*distance)
				startB=int(posB)-math.floor(0.5*distance)

	genomeA_pos=[chrA+"."+str(i) for i in range(startA,endA+1)]
	genomeB_pos=[chrB+"."+str(i) for i in range(startB,endB+1)]
	posA_ref_command=f"{SAMTOOLS} faidx "+bit_fn+" "+chrA+":"+str(startA)+"-"+str(endA)
	posB_ref_command=f"{SAMTOOLS} faidx "+bit_fn+" "+chrB+":"+str(startB)+"-"+str(endB)
	# posA_ref_command="samtools faidx "+bit_fn+" "+chrA+":"+str(startA)+"-"+str(endA)
	# posB_ref_command="samtools faidx "+bit_fn+" "+chrB+":"+str(startB)+"-"+str(endB)
	posA_ref_seq=("").join(str(subprocess.check_output(posA_ref_command,shell=True),"utf-8").strip().split("\n")[1:]).upper()
	posB_ref_seq=("").join(str(subprocess.check_output(posB_ref_command,shell=True),"utf-8").strip().split("\n")[1:]).upper()
	posA_ref_fn=('').join(random.choice(string.ascii_lowercase) for x in range(7))+"."+sample+"."+chrA+"."+posA+"."+strandA+".ref.fa"
	posB_ref_fn=('').join(random.choice(string.ascii_lowercase) for x in range(7))+"."+sample+"."+chrB+"."+posB+"."+strandB+".ref.fa"
	with open(posA_ref_fn,"w")as posA_ref_fh:
		posA_ref_fh.write(">"+sample+"."+chrA+"."+posA+"."+strandA+"\n")
		posA_ref_fh.write(posA_ref_seq+"\n")

	with open(posB_ref_fn,"w")as posB_ref_fh:
		posB_ref_fh.write(">"+sample+"."+chrB+"."+posB+"."+strandB+"\n")
		posB_ref_fh.write(posB_ref_seq+"\n")

	return (posA_ref_fn,posB_ref_fn,genomeA_pos,genomeB_pos)

def extract_sam(SV,bam,chr_in_bam):
	SV_lst=SV.split(".")
	chrA=SV_lst[0]
	posA=SV_lst[1]
	strandA=SV_lst[2]
	chrB=SV_lst[3]
	posB=SV_lst[4]
	strandB=SV_lst[5]
	sample=os.path.basename(bam).replace(".bam","")
	if chr_in_bam=="yes" and ("chr" not in SV):
		chrA="chr"+chrA
		chrB="chr"+chrB
	if (chr_in_bam=="no") and ("chr" in SV):
		chrA=chrA.replace("chr","")
		chrB=chrB.replace("chr","")
	sam_fn=('').join(random.choice(string.ascii_lowercase) for x in range(7))+"."+SV+"."+os.path.basename(bam).replace(".bam",".sam")
	command_A=f"{SAMTOOLS} view "+bam+" "+str(chrA)+":"+str(int(posA)-100)+"-"+str(int(posA)+100)+"|sort -k1,1 >"+sam_fn
	command_B=f"{SAMTOOLS} view "+bam+" "+str(chrB)+":"+str(int(posB)-100)+"-"+str(int(posB)+100)+"|sort -k1,1 >>"+sam_fn
	# command_A="samtools view "+bam+" "+str(chrA)+":"+str(int(posA)-100)+"-"+str(int(posA)+100)+"|sort -k1,1 >"+sam_fn
	# command_B="samtools view "+bam+" "+str(chrB)+":"+str(int(posB)-100)+"-"+str(int(posB)+100)+"|sort -k1,1 >>"+sam_fn
	os.system(command_A)
	os.system(command_B)
	return sam_fn


def extract_fasta(sam_fn):
	seq_data={}
	seqid_data={}
	candidate_reads_fn=sam_fn.replace(".sam",".candidate_SV.reads")
	candidate_SV_reads_fh=open(candidate_reads_fn,'w')
	candidate_fa_fn=sam_fn.replace(".sam",".candidate_SV.fa")
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
			if re.findall("^\d+M$",cigar):
				continue
			for nn_qual in qual_str:
				qscore=ord(nn_qual)-33
				if qscore >qual_cutoff:
					ngood+=1
			if ngood<qual_perc*len(qual_str):
				continue

			candidate_SV_reads_fh.write(line.strip()+"\n")
			if seq not in seq_data:
				seq_data[seq]=[seqid]
				candidate_fa_fh.write(">"+seqid+"\n")
				candidate_fa_fh.write(seq+"\n")
				seqid_data[seqid]=seq
			else:
				seq_data[seq].append(seqid)

	candidate_SV_reads_fh.close()
	candidate_fa_fh.close()
	return (seq_data,seqid_data,candidate_reads_fn,candidate_fa_fn)

def do_blast(posA_ref_fn,posB_ref_fn,candidate_fa_fn):
	posA_blast_fn=posA_ref_fn.replace(".fa",".blast")
	posB_blast_fn=posB_ref_fn.replace(".fa",".blast")
	command_posA=blast_exe_folder+"blastn -task blastn -word_size 4 -evalue 0.0001 -penalty -3 -reward 2 -gapopen 6 -gapextend 4 -dust no -soft_masking false -query "+candidate_fa_fn+" -subject "+posA_ref_fn+" -outfmt 15 > "+posA_blast_fn
	command_posB=blast_exe_folder+"blastn -task blastn -word_size 4 -evalue 0.0001 -penalty -3 -reward 2 -gapopen 6 -gapextend 4 -dust no -soft_masking false -query "+candidate_fa_fn+" -subject "+posB_ref_fn+" -outfmt 15 > "+posB_blast_fn
	os.system(command_posA)
	os.system(command_posB)
	return (posA_blast_fn,posB_blast_fn)

def do_focalSV_blast(SV_ref_fn,candidate_fa_fn):
	SV_blast_fn=SV_ref_fn.replace(".fa",".blast")
	command_SV=blast_exe_folder+"blastn -task blastn -word_size 4 -evalue 0.0001 -penalty -3 -reward 2 -gapopen 6 -gapextend 4 -dust no -soft_masking false -query "+candidate_fa_fn+" -subject "+SV_ref_fn+" -outfmt 15 > "+SV_blast_fn
	os.system(command_SV)
	return SV_blast_fn

def adjust_indel_pos(qseq,hseq,midline):
	hseq_lst=list(hseq)
	qseq_lst=list(qseq)
	midline_lst=list(midline)
	for i in range(len(qseq)):
		if list(qseq)[i]==list(hseq)[i]:
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

def locate_SV(initial_SV,seqid,seqid_data,seq_data,posA_blast_data,posB_blast_data):
	read_seq=seqid_data[seqid]
	all_read_ids=seq_data[read_seq]
	initial_A_chr,initial_A_pos,initial_A_strand,initial_B_chr,initial_B_pos,initial_B_strand,insertion=initial_SV.split(".")

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
	breakpoint_A="ND"
	breakpoint_B="ND"
	insertion="ND"
	#deal with overlapped seq mapping issue from blast result
	if A_hit_strand==initial_A_strand and B_hit_strand==initial_B_strand and A_query_to>=B_query_from and A_query_from<B_query_to:
		overlapped_bp=A_query_to-B_query_from+1
		try:
			if A_blast_midline[-(overlapped_bp+5):].count("|")<B_blast_midline[0:(overlapped_bp+5)].count("|"):
				A_query_to=A_query_to-overlapped_bp
				i=0
				overlapped_align_bp=0
				for nn in A_blast_qseq[::-1]:
					if nn !="-":
						i+=1
					overlapped_align_bp+=1
					if i==overlapped_bp:
						break
				if (A_hit_from<A_hit_to):
					A_hit_to=A_hit_to-overlapped_align_bp+A_blast_hseq[-overlapped_align_bp:].count("-")
				else:
					A_hit_to=A_hit_to+overlapped_align_bp-A_blast_hseq[-overlapped_align_bp:].count("-")

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
				if (B_hit_from<B_hit_to):
					B_hit_from=B_hit_from+overlapped_align_bp-B_blast_hseq[0:overlapped_align_bp].count("-")
				else:
					B_hit_from=B_hit_from-overlapped_align_bp+B_blast_hseq[0:overlapped_align_bp].count("-")

				B_blast_midline=B_blast_midline[overlapped_align_bp:]
				B_blast_qseq=B_blast_qseq[overlapped_align_bp:]
				B_blast_hseq=B_blast_hseq[overlapped_align_bp:]
		except:
			breakpoint_A="ND"
			breakpoint_B="ND"
			insertion="ND"
			
	if A_hit_strand!=initial_A_strand and B_hit_strand!=initial_B_strand and B_query_to>=A_query_from and B_query_from<A_query_to:
		overlapped_bp=B_query_to-A_query_from+1
		try:
			if A_blast_midline[0:(overlapped_bp+5)].count("|")<B_blast_midline[-(overlapped_bp+5):].count("|"):
				A_query_from=A_query_from+overlapped_bp
				i=0
				overlapped_align_bp=0
				for nn in A_blast_qseq:
					if nn !="-":
						i+=1
					overlapped_align_bp+=1
					if i==overlapped_bp:
						break
				if (A_hit_from<A_hit_to):
					A_hit_from=A_hit_from+overlapped_align_bp-A_blast_hseq[0:overlapped_align_bp].count("-")
				else:
					A_hit_from=A_hit_from-overlapped_align_bp+A_blast_hseq[0:overlapped_align_bp].count("-")
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
				if (B_hit_from<B_hit_to):
					B_hit_to=B_hit_to-overlapped_align_bp+B_blast_hseq[-overlapped_align_bp:].count("-")
				else:
					B_hit_to=B_hit_to+overlapped_align_bp-B_blast_hseq[-overlapped_align_bp:].count("-")
				B_blast_midline=B_blast_midline[:-overlapped_align_bp]
				B_blast_qseq=B_blast_qseq[:-overlapped_align_bp]
				B_blast_hseq=B_blast_hseq[:-overlapped_align_bp]
		except:
			breakpoint_A="ND"
			breakpoint_B="ND"
			insertion="ND"

	#deal with indel in the alignment
	if "-" in A_blast_qseq or "-" in A_blast_hseq:
		if A_hit_strand==initial_A_strand and B_hit_strand==initial_B_strand: 
			A_blast_qseq,A_blast_hseq,A_blast_midline=adjust_indel_pos(A_blast_qseq,A_blast_hseq,A_blast_midline)
		if A_hit_strand!=initial_A_strand and B_hit_strand!=initial_B_strand: 
			A_blast_qseq_rev,A_blast_hseq_rev,A_blast_midline_rev=adjust_indel_pos(A_blast_qseq[::-1],A_blast_hseq[::-1],A_blast_midline[::-1])
			A_blast_qseq,A_blast_hseq,A_blast_midline=A_blast_qseq_rev[::-1],A_blast_hseq_rev[::-1],A_blast_midline_rev[::-1]

	if "-" in B_blast_qseq or "-" in B_blast_hseq:
		if A_hit_strand==initial_A_strand and B_hit_strand==initial_B_strand: 
			B_blast_qseq_rev,B_blast_hseq_rev,B_blast_midline_rev=adjust_indel_pos(B_blast_qseq[::-1],B_blast_hseq[::-1],B_blast_midline[::-1])
			B_blast_qseq,B_blast_hseq,B_blast_midline=B_blast_qseq_rev[::-1],B_blast_hseq_rev[::-1],B_blast_midline_rev[::-1]
		if A_hit_strand!=initial_B_strand and B_hit_strand!=initial_B_strand: 
			B_blast_qseq,B_blast_hseq,B_blast_midline=adjust_indel_pos(B_blast_qseq,B_blast_hseq,B_blast_midline)

	#trim breakpoint gaps	
	if A_hit_strand==initial_A_strand and B_hit_strand==initial_B_strand:
		try:
			if A_blast_midline[-blast_trim:]!=blast_trim*"|":
				drop_A_len=A_blast_midline[::-1].index(blast_trim*"|")
				A_query_to=A_query_to-drop_A_len+A_blast_qseq[-drop_A_len:].count("-")
				if (A_hit_from<A_hit_to):
					A_hit_to=A_hit_to-drop_A_len+A_blast_hseq[-drop_A_len:].count("-")
				else:
					A_hit_to=A_hit_to+drop_A_len-A_blast_hseq[-drop_A_len:].count("-")
				A_blast_midline=A_blast_midline[:-drop_A_len]
				A_blast_qseq=A_blast_qseq[:-drop_A_len]
				A_blast_hseq=A_blast_hseq[:-drop_A_len]
		except:
				breakpoint_A="ND"
				breakpoint_B="ND"
				insertion="ND"

		try:
			if B_blast_midline[0:blast_trim]!=blast_trim*"|":
				drop_B_len=B_blast_midline.index(blast_trim*"|")
				B_query_from=B_query_from+drop_B_len-B_blast_qseq[0:drop_B_len].count("-")
				if (B_hit_from<B_hit_to):
					B_hit_from=B_hit_from+drop_B_len-B_blast_hseq[0:drop_B_len].count("-")
				else:
					B_hit_from=B_hit_from-drop_B_len+B_blast_hseq[0:drop_B_len].count("-")
				B_blast_midline=B_blast_midline[drop_B_len:]
				B_blast_qseq=B_blast_qseq[drop_B_len:]
				B_blast_hseq=B_blast_hseq[drop_B_len:]
		except:
				breakpoint_A="ND"
				breakpoint_B="ND"
				insertion="ND"

		if A_blast_midline.count(" ")>5 or B_blast_midline.count(" ")>5:
			breakpoint_A="ND"
			breakpoint_B="ND"
			insertion="ND"
		else:
			try:
				if A_query_to<B_query_from:
					breakpoint_A_used_seq=read_seq[:A_query_to]
					breakpoint_B_used_seq=read_seq[(B_query_from-1):]
					breakpoint_A=genomeA_pos[A_hit_to-1]+"."+initial_A_strand
					breakpoint_B=genomeB_pos[B_hit_from-1]+"."+initial_B_strand
					insertion=read_seq[A_query_to:(B_query_from-1)]
				else:
					breakpoint_A="ND"
					breakpoint_A="ND"
					insertion="ND"
			except:
				breakpoint_A="ND"
				breakpoint_B="ND"
				insertion="ND"

	if A_hit_strand!=initial_A_strand and B_hit_strand!=initial_B_strand:
		try:
			if A_blast_midline[0:blast_trim]!=blast_trim*"|":
				drop_A_len=A_blast_midline.index(blast_trim*"|")
				A_query_from=A_query_from+drop_A_len-A_blast_qseq[0:drop_A_len].count("-")
				if (A_hit_from<A_hit_to):
					A_hit_from=A_hit_from+drop_A_len-A_blast_hseq[0:drop_A_len].count("-")
				else:
					A_hit_from=A_hit_from-drop_A_len+A_blast_hseq[0:drop_A_len].count("-")

				A_blast_midline=A_blast_midline[drop_A_len:]
				A_blast_qseq=A_blast_qseq[drop_A_len:]
				A_blast_hseq=A_blast_hseq[drop_A_len:]
		except:
			breakpoint_A="ND"
			breakpoint_B="ND"
			insertion="ND"

		try:
			if B_blast_midline[-blast_trim:]!=blast_trim*"|":
				drop_B_len=B_blast_midline[::-1].index(blast_trim*"|")
				B_query_to=B_query_to-drop_B_len+B_blast_qseq[-drop_B_len:].count("-")
				if (B_hit_from<B_hit_to):
					B_hit_to=B_hit_to-drop_B_len+B_blast_hseq[-drop_B_len:].count("-")
				else:
					B_hit_to=B_hit_to+drop_B_len-B_blast_hseq[-drop_B_len:].count("-")
				B_blast_midline=B_blast_midline[:-drop_B_len]
				B_blast_qseq=B_blast_qseq[:-drop_B_len]
				B_blast_hseq=B_blast_hseq[:-drop_B_len]
		except:
			breakpoint_A="ND"
			breakpoint_B="ND"
			insertion="ND"

		if A_blast_midline.count(" ")>5 or B_blast_midline.count(" ")>5:
			breakpoint_A="ND"
			breakpoint_B="ND"
			insertion="ND"
		else:
			try:
				if B_query_to<A_query_from:
					breakpoint_A_used_seq=read_seq[(A_query_from-1):]
					breakpoint_B_used_seq=read_seq[:B_query_to]
					breakpoint_A=genomeA_pos[A_hit_from-1]+"."+initial_A_strand
					breakpoint_B=genomeB_pos[B_hit_to-1]+"."+initial_B_strand
					insertion=get_revcmp(read_seq[B_query_to:(A_query_from-1)])
				else:
					breakpoint_A="ND"
					breakpoint_B="ND"
					insertion="ND"
			except:
				breakpoint_A="ND"
				breakpoint_B="ND"
				insertion="ND"

	curated_SV=breakpoint_A+"."+breakpoint_B+"."+insertion
	if "ND" in curated_SV:
		all_read_ids="ND"
	return(curated_SV,all_read_ids)

def locate_focalSV(initial_SV,seqid,seqid_data,seq_data,posA_blast_data,posB_blast_data,genome_pos):
	read_seq=seqid_data[seqid]
	all_read_ids=seq_data[read_seq]
	initial_A_chr,initial_A_pos,initial_A_strand,initial_B_chr,initial_B_pos,initial_B_strand,insertion=initial_SV.split(".")

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
	breakpoint_A="ND"
	breakpoint_B="ND"
	insertion="ND"
	#if seqid=="A00641:249:HVM2FDMXX:1:1118:11406:17425_99":
	#	print (A_query_from,B_query_from,A_query_to,B_query_to)
	#deal with overlapped seq mapping issue from blast result
	if A_hit_strand==B_hit_strand=="+":
		if A_query_to-B_query_from>=0:
			overlapped_bp=A_query_to-B_query_from+1
			try:
				if A_blast_midline[-(overlapped_bp+5):].count("|")<B_blast_midline[0:(overlapped_bp+5)].count("|"):
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
				breakpoint_A="ND"
				breakpoint_B="ND"
				insertion="ND"
			

	if A_hit_strand==B_hit_strand=="-":
		if B_query_to-A_query_from>=0:
			overlapped_bp=B_query_to-A_query_from+1
			try:
				if A_blast_midline[0:(overlapped_bp+5)].count("|")<B_blast_midline[-(overlapped_bp+5):].count("|"):
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
				breakpoint_A="ND"
				breakpoint_B="ND"
				insertion="ND"

	#if seqid=="A00641:249:HVM2FDMXX:2:1123:1036:25833_1107":
	#	print("minus",seqid)
	#	print(seqid,A_query_from,A_query_to,B_query_from,B_query_to)
	#	print(overlapped_bp)
			
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
				breakpoint_A="ND"
				breakpoint_B="ND"
				insertion="ND"

		try:
			if B_blast_midline[0:blast_trim]!=blast_trim*"|":
				drop_B_len=B_blast_midline.index(blast_trim*"|")
				B_query_from=B_query_from+drop_B_len-B_blast_qseq[0:drop_B_len].count("-")
				B_hit_from=B_hit_from+drop_B_len-B_blast_hseq[0:drop_B_len].count("-")
				B_blast_midline=B_blast_midline[drop_B_len:]
				B_blast_qseq=B_blast_qseq[drop_B_len:]
				B_blast_hseq=B_blast_hseq[drop_B_len:]
		except:
				breakpoint_A="ND"
				breakpoint_B="ND"
				insertion="ND"

		if A_blast_midline.count(" ")>5 or B_blast_midline.count(" ")>5:
			breakpoint_A="ND"
			breakpoint_B="ND"
			insertion="ND"
		else:
			try:
				if A_query_to<B_query_from:
					breakpoint_A=genome_pos[A_hit_to-1]+"."+initial_A_strand
					breakpoint_B=genome_pos[B_hit_from-1]+"."+initial_B_strand
					insertion=read_seq[A_query_to:(B_query_from-1)]
				else:
					breakpoint_A="ND"
					breakpoint_A="ND"
					insertion="ND"
			except:
				breakpoint_A="ND"
				breakpoint_B="ND"
				insertion="ND"

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
			breakpoint_A="ND"
			breakpoint_B="ND"
			insertion="ND"

		try:
			if B_blast_midline[-blast_trim:]!=blast_trim*"|":
				drop_B_len=B_blast_midline[::-1].index(blast_trim*"|")
				B_query_to=B_query_to-drop_B_len+B_blast_qseq[-drop_B_len:].count("-")
				B_hit_to=B_hit_to-drop_B_len+B_blast_hseq[-drop_B_len:].count("-")
				B_blast_midline=B_blast_midline[:-drop_B_len]
				B_blast_qseq=B_blast_qseq[:-drop_B_len]
				B_blast_hseq=B_blast_hseq[:-drop_B_len]
		except:
			breakpoint_A="ND"
			breakpoint_B="ND"
			insertion="ND"

		if A_blast_midline.count(" ")>5 or B_blast_midline.count(" ")>5:
			breakpoint_A="ND"
			breakpoint_B="ND"
			insertion="ND"
		else:
			try:
				if B_query_to<A_query_from:
					breakpoint_A=genome_pos[A_hit_from-1]+"."+initial_A_strand
					breakpoint_B=genome_pos[B_hit_to-1]+"."+initial_B_strand
					insertion=get_revcmp(read_seq[B_query_to:(A_query_from-1)])
				else:
					breakpoint_A="ND"
					breakpoint_B="ND"
					insertion="ND"
			except:
				breakpoint_A="ND"
				breakpoint_B="ND"
				insertion="ND"
			#print(curated_SV)

	curated_SV=breakpoint_A+"."+breakpoint_B+"."+insertion
	if "ND" in curated_SV:
		all_read_ids="ND"
	return(curated_SV,all_read_ids)

def parse_focalSV_blast(blast_fn,seqid_data,seq_data,SV,genome_pos):
	strand_transform={"Plus":"+","Minus":"-"}
	curated_SV="ND"
	curated_SVs_data={}
	target_pos=blast_fn.replace(".blast","")
	initial_A_chr,initial_A_pos,initial_A_strand,initial_B_chr,initial_B_pos,initial_B_strand,insertion=SV.split(".")
	#anno_blast_fn=blast_fn+".anno"
	#anno_blast_fh=open(anno_blast_fn,"w")
	#anno_blast_fh.write("SV\ttarget_pos\tseqid\tseq\talign_len\tidentity\tgaps\tquery_from\tquery_to\tquery_strand\thit_from\thit_to\thit_strand\tqseq\thseq\tmidline\tbit_score\tscore\tevalue\n")
	blast_json=json.load(open(blast_fn))["BlastOutput2"]
	blast_data={}
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
			if len(hsps)>1:
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
						"initial_SV":SV,
						"target_pos":target_pos,
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
	
	for seqid in blast_data:
		hsp_nums=blast_data[seqid].keys()
		seq=seqid_data[seqid]
		drop_hsps=[]
		for ii in blast_data[seqid]:
			ii_query_from=int(blast_data[seqid][ii]["query_from"])
			ii_query_to=int(blast_data[seqid][ii]["query_to"])
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
				if ii_query_from==jj_query_from and ii_query_to==jj_query_to and ii_hit_strand==jj_hit_strand:
					if ii_query_strand==ii_hit_strand:
						if ii_query_to>len(seq)-10:
							if 0.5*(ii_hit_from+ii_hit_to)<0.5*(jj_hit_from+jj_hit_to):
								drop_hsps.append(ii)
							else:
								drop_hsps.append(jj)
						else:
							if 0.5*(ii_hit_from+ii_hit_to)>0.5*(jj_hit_from+jj_hit_to):
								drop_hsps.append(ii)
							else:
								drop_hsps.append(jj)
					else:
						if ii_query_to>len(seq)-10:
							if 0.5*(ii_hit_from+ii_hit_to)>0.5*(jj_hit_from+jj_hit_to):
								drop_hsps.append(ii)
							else:
								drop_hsps.append(jj)
						else:
							if 0.5*(ii_hit_from+ii_hit_to)<0.5*(jj_hit_from+jj_hit_to):
								drop_hsps.append(ii)
							else:
								drop_hsps.append(jj)
						
				if ii_query_from<jj_query_from and ii_query_to>jj_query_to:
					drop_hsps.append(jj)
		for drop_hsp in set(drop_hsps):
			del blast_data[seqid][drop_hsp]
		if len (blast_data[seqid])>1:
			hsp_range_lst=[]
			for hsp in blast_data[seqid]:
				hsp_range_lst.append((hsp,blast_data[seqid][hsp]["query_from"],blast_data[seqid][hsp]["query_to"]))
			hsp_range_lst=sorted(hsp_range_lst,key=lambda x:x[1])[:2]
			#print(SV,hsp_range_lst)
			hsp_left_order=hsp_range_lst[0][0]
			hsp_right_order=hsp_range_lst[1][0]
			if blast_data[seqid][hsp_left_order]["hit_strand"]==blast_data[seqid][hsp_right_order]["hit_strand"]=="+":
				posA_blast_data={
				seqid:{
						"initial_SV":SV,
						"target_pos":target_pos,
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
						"initial_SV":SV,
						"target_pos":target_pos,
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
				curated_SV,all_read_ids=locate_focalSV(SV,seqid,seqid_data,seq_data,posA_blast_data,posB_blast_data,genome_pos)
				if curated_SV not in curated_SVs_data:
					curated_SVs_data[curated_SV]=all_read_ids
				else:
					curated_SVs_data[curated_SV]+=all_read_ids
			if blast_data[seqid][hsp_left_order]["hit_strand"]==blast_data[seqid][hsp_right_order]["hit_strand"]=="-":
				posB_blast_data={
					seqid:{
						"initial_SV":SV,
						"target_pos":target_pos,
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
						"initial_SV":SV,
						"target_pos":target_pos,
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
				curated_SV,all_read_ids=locate_focalSV(SV,seqid,seqid_data,seq_data,posA_blast_data,posB_blast_data,genome_pos)
				if curated_SV not in curated_SVs_data:
					curated_SVs_data[curated_SV]=all_read_ids
				else:
					curated_SVs_data[curated_SV]+=all_read_ids
		#if seqid=="A00641:249:HVM2FDMXX:1:1219:10836:19038_1171":
		#	print(seqid,curated_SV,seq)
	return (curated_SVs_data)


def parse_blast(blast_fn,seqid_data,SV):
	strand_transform={"Plus":"+","Minus":"-"}
	target_pos=blast_fn.replace(".blast","")
	#anno_blast_fn=blast_fn+".anno"
	#anno_blast_fh=open(anno_blast_fn,"w")
	#anno_blast_fh.write("SV\ttarget_pos\tseqid\tseq\talign_len\tidentity\tgaps\tquery_from\tquery_to\tquery_strand\thit_from\thit_to\thit_strand\tqseq\thseq\tmidline\tbit_score\tscore\tevalue\n")
	blast_json=json.load(open(blast_fn))["BlastOutput2"]
	blast_data={}
	for i in range(0,len(blast_json)):
		blast_report= blast_json[i]["report"]
		blast_result= blast_report["results"]
		bl2seq=blast_result["bl2seq"][0]
		seqid=bl2seq["query_title"]
		seq=seqid_data[seqid]
		query_len=bl2seq["query_len"]


		if ("hits" in bl2seq) and len(bl2seq["hits"])>0:
			hits=bl2seq["hits"]
			hsps=hits[0]["hsps"][0]
			align_len=hsps["align_len"]
			identity=hsps["identity"]
			gaps=hsps["gaps"]
			query_from=hsps["query_from"]
			query_to=hsps["query_to"]
			
			query_strand=strand_transform[hsps["query_strand"]]
			hit_from=hsps["hit_from"]
			hit_to=hsps["hit_to"]
			hit_strand=strand_transform[hsps["hit_strand"]]
			qseq=hsps["qseq"]
			hseq=hsps["hseq"]
			midline=hsps["midline"]
			bit_score=hsps["bit_score"]
			score=hsps["score"]
			evalue=hsps["evalue"]
		else:
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
			
		blast_data[seqid]={
			"initial_SV":SV,
			"target_pos":target_pos,
			"hsps":hsps,
			"align_len":align_len,
			"identity":identity,
			"gaps":gaps,
			"query_from":query_from,
			"query_to":query_to,
			"query_strand":query_strand,
			"hit_from":hit_from,
			"hit_to":hit_to,
			"hit_strand":hit_strand,
			"qseq":qseq,
			"hseq":hseq,
			"midline":midline,
			"bit_score":bit_score,
			"score":score,
			"evalue":evalue
		}
		#out="%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n"%(SV,target_pos,seqid,seq,align_len,identity,gaps,query_from,query_to,query_strand,hit_from,hit_to,hit_strand,qseq,hseq,midline,bit_score,score,evalue)
		#anno_blast_fh.write(out)
	#anno_blast_fh.close()
	return blast_data
def make_reference_seq_4inv(chro,pos):
	if "chr" not in chro:
		chro="chr"+chro
	start=str(int(pos)-flanking_length)
	end=str(int(pos)+flanking_length)
	seq_command=f"{SAMTOOLS} faidx "+bit_fn+" "+chro+":"+start+"-"+end
	# seq_command="samtools faidx "+bit_fn+" "+chro+":"+start+"-"+end
	genome_seq=("").join(str(subprocess.check_output(seq_command,shell=True),"utf-8").strip().split("\n")[1:]).upper()
	genome_pos_range=range(int(start),int(end)+1)
	genome_revcmpseq=get_revcmp(genome_seq)
	genome_revpos_range=genome_pos_range[::-1]
	return (genome_seq,genome_pos_range,genome_revcmpseq,genome_revpos_range)

def collect_SV_reads_4inv(candidate_reads_fn,SV,posA_genome_seq,posA_genome_revcmpseq,posB_genome_seq,posB_genome_revcmpseq):
	chrA,posA,strandA,chrB,posB,strandB,insertion=SV.split(".")
	candidate_invSV_reads_fn=sample+"_"+SV+"_candidate_invSV_reads"
	candidate_invSV_reads_fh=open(candidate_invSV_reads_fn,"w")
	with open (candidate_reads_fn,"r") as in_fh:
		for line in in_fh:
			lst=line.strip().split("\t")
			read_seq=lst[9]
			if strandA=="+" and strandB=="-":
				if (read_seq[0:20] in posA_genome_seq) and (read_seq[-20:] in posB_genome_revcmpseq):
					candidate_invSV_reads_fh.write("initial_read\t"+line)
				elif (read_seq[0:20] in posB_genome_seq) and (read_seq[-20:] in posA_genome_revcmpseq):
					lst[9]=get_revcmp(read_seq)
					candidate_invSV_reads_fh.write("revcmp_read\t"+("\t").join(lst)+"\n")
				else:
					continue
			if strandA=="-" and strandB=="+":
				if (read_seq[0:20] in posA_genome_revcmpseq) and (read_seq[-20:] in posB_genome_seq):
					candidate_invSV_reads_fh.write("initial_read\t"+line)
				elif (read_seq[0:20] in posB_genome_revcmpseq) and (read_seq[-20:] in posA_genome_seq):
					lst[9]=get_revcmp(read_seq)
					candidate_invSV_reads_fh.write("revcmp_read\t"+("\t").join(lst)+"\n")
				else:
					continue
	candidate_invSV_reads_fh.close()
	return candidate_invSV_reads_fn
def curate_long_insertions(SV,left_seq_index,left_genome_seq,left_genome_pos_range,right_seq_index,right_genome_seq,right_genome_pos_range):
	raw_chrA,raw_posA,raw_strandA,raw_chrB,raw_posB,raw_strandB,insertion=SV.split(".")
	curated_left_genome_pos=raw_chrA+"."+raw_posA+"."+raw_strandA
	curated_right_genome_pos=raw_chrB+"."+raw_posB+"."+raw_strandB
	SNPs=[]
	tmp_insertion=insertion
	if tmp_insertion[1:8]==left_genome_seq[(left_seq_index+2):(left_seq_index+9)]:
		left_SNP_pos=raw_chrA+"."+str(left_genome_pos_range[left_seq_index+1])+"."+left_genome_seq[(left_seq_index+1)]+"."+tmp_insertion[0]
		if raw_strandA=="-":
			left_SNP_pos=raw_chrA+"."+str(left_genome_pos_range[left_seq_index+1])+"."+get_revcmp(left_genome_seq[(left_seq_index+1)])+"."+get_revcmp(tmp_insertion[0])
		SNPs.append(left_SNP_pos)
		for i in range(8,len(insertion)+1):
			if tmp_insertion[1:i] == left_genome_seq[(left_seq_index+2):(left_seq_index+1+i)]:
				if i<len(insertion):
					continue
				else:
					curated_left_seq=tmp_insertion
					curated_left_genome_pos=raw_chrA+"."+str(left_genome_pos_range[left_seq_index+i])+"."+raw_strandA
					tmp_insertion=tmp_insertion.replace(curated_left_seq,"")
					break
			else:
				try:
					curated_left_seq=tmp_insertion[0:(i-1)]
					curated_left_genome_pos=raw_chrA+"."+str(left_genome_pos_range[left_seq_index+i])+"."+raw_strandA
					tmp_insertion=tmp_insertion.replace(curated_left_seq,"")
					break
				except:
					curated_left_genome_pos="ND"
					tmp_insertion="ND"
	if len(tmp_insertion)>=8:
		if tmp_insertion[-8:-1]==right_genome_seq[(right_seq_index-8):(right_seq_index-1)]:
			right_SNP_pos=raw_chrB+"."+str(right_genome_pos_range[right_seq_index-1])+"."+right_genome_seq[(right_seq_index-1)]+"."+tmp_insertion[-1]
			if raw_strandB=="-":
				right_SNP_pos=raw_chrB+"."+str(right_genome_pos_range[right_seq_index-1])+"."+get_revcmp(right_genome_seq[(right_seq_index-1)])+"."+get_revcmp(tmp_insertion[-1])
			SNPs.append(right_SNP_pos)
			for i in range(8,len(tmp_insertion)+1):
				if tmp_insertion[-i:-1]==right_genome_seq[(right_seq_index-i):(right_seq_index-1)]:
					if i<len(tmp_insertion):
						continue
					else:
						curated_right_seq=tmp_insertion
						curated_right_genome_pos=raw_chrB+"."+str(right_genome_pos_range[right_seq_index-i])+"."+raw_strandB
						tmp_insertion=tmp_insertion.replace(curated_right_seq,"")
						break
				else:
					try:
						curated_right_seq=tmp_insertion[(-i+1):]
						curated_right_genome_pos=raw_chrB+"."+str(right_genome_pos_range[right_seq_index-i+1])+"."+raw_strandB
						tmp_insertion=tmp_insertion.replace(curated_right_seq,"")
						break
					except:
						curated_right_genome_pos="ND"
						tmp_insertion="ND"
	curated_insertion=tmp_insertion
	curated_SV=curated_left_genome_pos+"."+curated_right_genome_pos+"."+tmp_insertion
	return(curated_SV,curated_insertion,SNPs)




def locate_invSV(read_seq,left_chr,left_genome_seq,left_strand,left_genome_pos_range,right_chr,right_genome_seq,right_strand,right_genome_pos_range):
	tmp_read_seq=read_seq
	left_seq=""
	left_seq_index=[]
	left_genome_pos=""
	right_seq=""
	right_seq_index=[]
	right_genome_pos=""
	curated_SV=""
	SNPs=""

	for i in range(20,len(read_seq)):
		if read_seq[0:i] in left_genome_seq:
			continue
		else:
			left_seq=tmp_read_seq[0:(i-1)]
			left_seq_index=[(nn.end()-1) for nn in re.finditer(left_seq,left_genome_seq)]
			if len(left_seq_index)!=1:
				continue
			left_genome_pos=left_chr+"."+str(left_genome_pos_range[left_seq_index[0]])
			break

	tmp_read_seq=tmp_read_seq.replace(left_seq,"")
	for i in range(20,(len(tmp_read_seq)+1)):
		if tmp_read_seq[-i:] in right_genome_seq:
			if i<len(tmp_read_seq):
				continue
			else:
				right_seq=tmp_read_seq[-(len(tmp_read_seq)):]
				right_seq_index=[(nn.start()) for nn in re.finditer(right_seq,right_genome_seq)]
				if len(right_seq_index)!=1:
					continue
				right_genome_pos=right_chr+"."+str(right_genome_pos_range[right_seq_index[0]])
		else:
			right_seq=tmp_read_seq[-(i-1):]
			right_seq_index=[(i.start()) for i in re.finditer(right_seq,right_genome_seq)]
			if len(right_seq_index)!=1:
				continue
			right_genome_pos=right_chr+"."+str(right_genome_pos_range[right_seq_index[0]])
			break
	tmp_read_seq=tmp_read_seq.replace(right_seq,"")
	insertion_seq=tmp_read_seq
	curated_SV=left_genome_pos+"."+left_strand+"."+right_genome_pos+"."+right_strand+"."+insertion_seq
	if curated_SV.count("chr")!=2:
		return ["ND","ND"]
	else:
		SNPs=[]
		if len(insertion_seq)>=8:
			curated_SV,insertion_seq,SNPs=curate_long_insertions(curated_SV,left_seq_index[0],left_genome_seq,left_genome_pos_range,right_seq_index[0],right_genome_seq,right_genome_pos_range)
		if "ND" in curated_SV:
			return ["ND","ND"]
		else:
			left_genome_pos=int(curated_SV.split(".")[1])
			right_genome_pos=int(curated_SV.split(".")[4])
			left_pos_index=left_genome_pos_range.index(left_genome_pos)
			right_pos_index=right_genome_pos_range.index(right_genome_pos)
			return [curated_SV,SNPs]

def parse_invSV_reads(candidate_invSV_reads_fn,SV,posA_genome_seq,posA_genome_revcmpseq,posB_genome_seq,posB_genome_revcmpseq):
	chrA,posA,strandA,chrB,posB,strandB,insertion=SV.split(".")
	SV_candidate_reads={}
	with open (candidate_invSV_reads_fn,"r") as candidate_invSV_fh:
		for line in candidate_invSV_fh:
			lst=line.strip().split("\t")
			read_id=lst[1]+"_"+lst[2]
			cigar=lst[6]
			read_seq=lst[10]
			read_qual=lst[11]
			if read_seq not in SV_candidate_reads:
				SV_candidate_reads[read_seq]={
					"read_id":[read_id],
					"cigar":[cigar],
					"read_qual":[read_qual]
				}
			else:
				SV_candidate_reads[read_seq]["read_id"].append(read_id)
				SV_candidate_reads[read_seq]["cigar"].append(cigar)
				SV_candidate_reads[read_seq]["read_qual"].append(read_qual)

	curated_SVs={}
	for read_seq in SV_candidate_reads:
		if strandA=="+" and strandB=="-":
			if (read_seq[0:20] in posA_genome_seq) and (read_seq[-20:] in posB_genome_revcmpseq):
				curated_SV,SNPs=locate_invSV(read_seq,chrA,posA_genome_seq,strandA,posA_genome_pos_range,chrB,posB_genome_revcmpseq,strandB,posB_genome_revpos_range)

		if strandA=="-" and strandB=="+":
			if (read_seq[0:20] in posA_genome_revcmpseq) and (read_seq[-20:] in posB_genome_seq):
				curated_SV,SNPs=locate_invSV(read_seq,chrA,posA_genome_revcmpseq,strandA,posA_genome_revpos_range,chrB,posB_genome_seq,strandB,posB_genome_pos_range)
		if curated_SV=="noSV":
			continue
		elif curated_SV in curated_SVs:
			curated_SVs[curated_SV]["read_id"]+=SV_candidate_reads[read_seq]["read_id"]
			curated_SVs[curated_SV]["SNPs"]+=SNPs
		else:
			curated_SVs[curated_SV]={
				"read_id":SV_candidate_reads[read_seq]["read_id"],
				"SNPs":SNPs
			}
	return curated_SVs


def get_SV_padding(SV):
	chrA,posA,strandA,chrB,posB,strandB,insertion=SV.split(".")
	if "chr" not in chrA:
		chrA="chr"+chrA
	if "chr" not in chrB:
		chrA="chr"+chrB
	if strandA=="+":
		posA_padding_command=f"{SAMTOOLS} faidx "+bit_fn+" "+chrA+":"+str(int(posA)-SV_padding+1)+"-"+str(posA)
		# posA_padding_command="samtools faidx "+bit_fn+" "+chrA+":"+str(int(posA)-SV_padding+1)+"-"+str(posA)
		posA_padding_seq=("").join(str(subprocess.check_output(posA_padding_command,shell=True),"utf-8").strip().split("\n")[1:]).upper()
	if strandA=="-":
		posA_padding_command=f"{SAMTOOLS} faidx "+bit_fn+" "+chrA+":"+str(posA)+"-"+str(int(posA)+SV_padding-1)
		# posA_padding_command="samtools faidx "+bit_fn+" "+chrA+":"+str(posA)+"-"+str(int(posA)+SV_padding-1)
		posA_padding_seq=("").join(str(subprocess.check_output(posA_padding_command,shell=True),"utf-8").strip().split("\n")[1:]).upper()
		posA_padding_seq=get_revcmp(posA_padding_seq)
	
	if strandB=="+":
		posB_padding_command=f"{SAMTOOLS} faidx "+bit_fn+" "+chrB+":"+str(posB)+"-"+str(int(posB)+SV_padding-1)
		# posB_padding_command="samtools faidx "+bit_fn+" "+chrB+":"+str(posB)+"-"+str(int(posB)+SV_padding-1)
		posB_padding_seq=("").join(str(subprocess.check_output(posB_padding_command,shell=True),"utf-8").strip().split("\n")[1:]).upper()
	if strandB=="-":
		posB_padding_command=f"{SAMTOOLS} faidx "+bit_fn+" "+chrB+":"+str(int(posB)-SV_padding+1)+"-"+str(posB)
		# posB_padding_command="samtools faidx "+bit_fn+" "+chrB+":"+str(int(posB)-SV_padding+1)+"-"+str(posB)
		posB_padding_seq=("").join(str(subprocess.check_output(posB_padding_command,shell=True),"utf-8").strip().split("\n")[1:]).upper()
		posB_padding_seq=get_revcmp(posB_padding_seq)
	SV_with_padding="%s.%s.%s.%s.%s.%s.%s.%s.%s"%(chrA,posA,strandA,chrB,posB,strandB,posA_padding_seq,insertion,posB_padding_seq)
	return SV_with_padding
		
def get_depth(SV,bam,chr_in_bam):
	chrA,posA,strandA,chrB,posB,strandB,insertion=SV.split(".")
	if chr_in_bam=="yes":
		depth_posA_command=f"{SAMTOOLS} depth -r "+chrA+":"+posA+"-"+posA+" "+bam+"|cut -f3"
		depth_posB_command=f"{SAMTOOLS} depth -r "+chrB+":"+posB+"-"+posB+" "+bam+"|cut -f3"
  		# depth_posA_command="samtools depth -r "+chrA+":"+posA+"-"+posA+" "+bam+"|cut -f3"
		# depth_posB_command="samtools depth -r "+chrB+":"+posB+"-"+posB+" "+bam+"|cut -f3"
	else:
		depth_posA_command=f"{SAMTOOLS} depth -r "+chrA.replace("chr","")+":"+posA+"-"+posA+" "+bam+"|cut -f3"
		depth_posB_command=f"{SAMTOOLS} depth -r "+chrB.replace("chr","")+":"+posB+"-"+posB+" "+bam+"|cut -f3"
  		# depth_posA_command="samtools depth -r "+chrA.replace("chr","")+":"+posA+"-"+posA+" "+bam+"|cut -f3"
		# depth_posB_command="samtools depth -r "+chrB.replace("chr","")+":"+posB+"-"+posB+" "+bam+"|cut -f3"
	try:
		depth_posA=str(subprocess.check_output(depth_posA_command,shell=True),"utf-8").strip()
	except:
		depth_posA="0"
	try:
		depth_posB=str(subprocess.check_output(depth_posB_command,shell=True),"utf-8").strip()
	except:
		depth_posB="0"
	#added line 1418-1421 Dec5_2023
	if depth_posA=="":
		depth_posA="0"
	if depth_posB=="":
		depth_posB="0"
	return (depth_posA,depth_posB)


out_SV_fn=SV_fn+"_curated_SV.txt"
out_SV_fh=open(out_SV_fn,"w")
out_SV_fh.write("sample\tinitial_SV\tcurated_SV\tgenome_version\tcurated_SV_reads\tdepth_posA\tdepth_posB\tcurated_SV_with_padding\n")

out_SV_reads_fn=SV_fn+"_curated_SV_reads.txt"
out_SV_reads_fh=open(out_SV_reads_fn,"w")
out_SV_reads_fh.write("curated_SV_reads\tread_id\tbam\n")
#with open ("SV_test.txt","r") as SV_fh:
with open (SV_fn,"r") as SV_fh:
	for line in SV_fh:
		lst=line.strip().split("\t")
		initial_SV=(".").join(lst[0:6])
		initial_chrA=lst[0]
		initial_chrB=lst[3]
		bam=lst[6]
		bit_fn,genome_version,chr_in_bam=get_bit_fn(bam)
		sample=os.path.basename(bam).replace(".bam", "")
		files_to_delete=[]
		if not ((initial_chrA in centromere_regions[genome_version]) and (initial_chrB in centromere_regions[genome_version])):
			out_SV="%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n"%(sample,initial_SV,"ND","ND","ND","ND","ND","ND")
			out_SV_fh.write(out_SV)
			continue
		if initial_SV.count(".")==5:
			initial_SV=initial_SV+"."
		SV=adjust_SV(initial_SV)
		adjust_chrA,adjust_posA,adjust_strandA,adjust_chrB,adjust_posB,adjust_strandB,adjust_insertion=SV.split(".")
		if not ((adjust_chrA in centromere_regions[genome_version])and (adjust_chrB in centromere_regions[genome_version])):
			out_SV="%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n"%(sample,initial_SV,"ND","ND","ND","ND","ND","ND")
			out_SV_fh.write(out_SV)
			continue
		if (centromere_regions[genome_version][adjust_chrA][0]<int(adjust_posA)<centromere_regions[genome_version][adjust_chrA][1]) or (centromere_regions[genome_version][adjust_chrB][0]<int(adjust_posB)<centromere_regions[genome_version][adjust_chrB][1]):
			out_SV="%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n"%(sample,initial_SV,"ND","ND","ND","ND","ND","ND")
			out_SV_fh.write(out_SV)
			continue

		breakpoint_distance=abs(int(adjust_posB)-int(adjust_posA))

		if (adjust_chrA!=adjust_chrB) or (breakpoint_distance)>300:
			posA_ref_fn,posB_ref_fn,genomeA_pos,genomeB_pos=make_ref(SV,bam,bit_fn)
			files_to_delete+=[posA_ref_fn,posB_ref_fn]
		else:
			SV_ref_fn,genome_pos=make_focalSV_ref(SV,bam,bit_fn)	
			files_to_delete+=[SV_ref_fn]
		sam_fn=extract_sam(SV,bam,chr_in_bam)
		sam_fn_check_command="wc -l "+sam_fn
		sam_fn_check_result=str(subprocess.check_output(sam_fn_check_command,shell=True),"utf-8")
		if int(sam_fn_check_result.split(" ")[0])>sam_reads_limit:
			out_SV="%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n"%(sample,initial_SV,"ND","ND","ND","ND","ND","ND")
			continue
		seq_data,seqid_data,candidate_reads_fn,candidate_fa_fn=extract_fasta(sam_fn)
		files_to_delete+=[candidate_reads_fn,candidate_fa_fn,sam_fn]
		candidate_fa_fn_check_command="wc -l "+candidate_fa_fn
		candidate_fa_fn_check_result=str(subprocess.check_output(candidate_fa_fn_check_command,shell=True),"utf-8")
		if candidate_fa_fn_check_result.split(" ")[0]=="0":
			depth_posA,depth_posB=get_depth(initial_SV,bam,chr_in_bam)
			curated_SV="ND"
			curated_SV_with_padding="ND"
			curated_SV_reads_num=0
			out_SV="%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n"%(sample,initial_SV,curated_SV,genome_version,str(curated_SV_reads_num),depth_posA,depth_posB,curated_SV_with_padding)
			out_SV_fh.write(out_SV)
			for fn_to_delete in set(files_to_delete):
				os.remove(fn_to_delete)
			continue

		if (adjust_chrA!=adjust_chrB) or (breakpoint_distance)>300:
			posA_blast_fn,posB_blast_fn=do_blast(posA_ref_fn,posB_ref_fn,candidate_fa_fn)
			posA_blast_data=parse_blast(posA_blast_fn,seqid_data,SV)
			posB_blast_data=parse_blast(posB_blast_fn,seqid_data,SV)
			files_to_delete+=[posA_blast_fn,posB_blast_fn]
			curated_SVs={}
			for seqid in seqid_data:
				if posA_blast_data[seqid]["hsps"]=="ND" or posB_blast_data[seqid]["hsps"]=="ND":
					continue
				curated_SV,curated_reads=locate_SV(SV,seqid,seqid_data,seq_data,posA_blast_data,posB_blast_data)
				if "ND" in curated_SV:
					continue
				if curated_SV not in curated_SVs:
					curated_SVs[curated_SV]=curated_reads
				else:
					curated_SVs[curated_SV]+=curated_reads
			if len(curated_SVs)==0:
				depth_posA,depth_posB=get_depth(initial_SV,bam,chr_in_bam)
				curated_SV="ND"
				curated_SV_with_padding="ND"
				curated_SV_reads_num=0
				out_SV="%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n"%(sample,initial_SV,curated_SV,genome_version,str(curated_SV_reads_num),depth_posA,depth_posB,curated_SV_with_padding)
				out_SV_fh.write(out_SV)
			else:
				curated_SVs_list=[]		
				for curated_SV in curated_SVs:
					curated_SVs_list.append((curated_SV,len(set(curated_SVs[curated_SV])),set(curated_SVs[curated_SV])))
				sorted_curated_SVs_list=sorted(curated_SVs_list,key=lambda x:x[1],reverse=True)
				if len(sorted_curated_SVs_list)>3:
					sorted_curated_SVs_list=sorted_curated_SVs_list[0:3]
				for curated_SV_info in sorted_curated_SVs_list:
					curated_SV=curated_SV_info[0]
					curated_SV_reads_num=curated_SV_info[1]
					curated_SV_reads=curated_SV_info[2]
					curated_SV_with_padding=get_SV_padding(curated_SV)
					depth_posA,depth_posB=get_depth(curated_SV,bam,chr_in_bam)
					out_SV="%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n"%(sample,initial_SV,curated_SV,genome_version,str(curated_SV_reads_num),depth_posA,depth_posB,curated_SV_with_padding)
					out_SV_fh.write(out_SV)
					for read in curated_SV_reads:
						out_SV_read="%s\t%s\t%s\n"%(curated_SV,read,bam)
						out_SV_reads_fh.write(out_SV_read)
		else:
			if adjust_strandA==adjust_strandB:
				SV_blast_fn=do_focalSV_blast(SV_ref_fn,candidate_fa_fn)
				files_to_delete+=[SV_blast_fn]
				curated_SVs=parse_focalSV_blast(SV_blast_fn,seqid_data,seq_data,SV,genome_pos)
				if len(curated_SVs)==0:
					depth_posA,depth_posB=get_depth(initial_SV,bam,chr_in_bam)
					curated_SV="ND"
					curated_SV_with_padding="ND"
					curated_SV_reads_num=0
					out_SV="%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n"%(sample,initial_SV,curated_SV,genome_version,str(curated_SV_reads_num),depth_posA,depth_posB,curated_SV_with_padding)
					out_SV_fh.write(out_SV)
				else:
					curated_SVs_list=[]		
					for curated_SV in curated_SVs:
						if "ND" in curated_SV:
							continue
						else:
							curated_SVs_list.append((curated_SV,len(set(curated_SVs[curated_SV])),set(curated_SVs[curated_SV])))
					sorted_curated_SVs_list=sorted(curated_SVs_list,key=lambda x:x[1],reverse=True)
					if len(sorted_curated_SVs_list)==0:
						depth_posA,depth_posB=get_depth(initial_SV,bam,chr_in_bam)
						curated_SV="ND"
						curated_SV_with_padding="ND"
						curated_SV_reads_num=0
						out_SV="%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n"%(sample,initial_SV,curated_SV,genome_version,str(curated_SV_reads_num),depth_posA,depth_posB,curated_SV_with_padding)
						out_SV_fh.write(out_SV)
					else:
						if len(sorted_curated_SVs_list)>3:
							sorted_curated_SVs_list=sorted_curated_SVs_list[0:3]
						for curated_SV_info in sorted_curated_SVs_list:
							curated_SV=curated_SV_info[0]
							curated_SV_reads_num=curated_SV_info[1]
							curated_SV_reads=curated_SV_info[2]
							curated_SV_with_padding=get_SV_padding(curated_SV)
							depth_posA,depth_posB=get_depth(curated_SV,bam,chr_in_bam)
							out_SV="%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n"%(sample,initial_SV,curated_SV,genome_version,str(curated_SV_reads_num),depth_posA,depth_posB,curated_SV_with_padding)
							out_SV_fh.write(out_SV)
							for read in curated_SV_reads:
								out_SV_read="%s\t%s\t%s\n"%(curated_SV,read,bam)
								out_SV_reads_fh.write(out_SV_read)
			else:
				posA_genome_seq,posA_genome_pos_range,posA_genome_revcmpseq,posA_genome_revpos_range=make_reference_seq_4inv(adjust_chrA,adjust_posA)
				posB_genome_seq,posB_genome_pos_range,posB_genome_revcmpseq,posB_genome_revpos_range=make_reference_seq_4inv(adjust_chrB,adjust_posB)
				candidate_invSV_reads_fn=collect_SV_reads_4inv(candidate_reads_fn,SV,posA_genome_seq,posA_genome_revcmpseq,posB_genome_seq,posB_genome_revcmpseq)
				curated_SVs=parse_invSV_reads(candidate_invSV_reads_fn,SV,posA_genome_seq,posA_genome_revcmpseq,posB_genome_seq,posB_genome_revcmpseq)
				files_to_delete+=[candidate_invSV_reads_fn]
		
				if len(curated_SVs)==0:
					depth_posA,depth_posB=get_depth(initial_SV,bam,chr_in_bam)
					curated_SV="ND"
					curated_SV_with_padding="ND"
					curated_SV_reads_num=0
					out_SV="%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n"%(sample,initial_SV,curated_SV,genome_version,str(curated_SV_reads_num),depth_posA,depth_posB,curated_SV_with_padding)
					out_SV_fh.write(out_SV)
				else:
					curated_SVs_list=[]		
					for curated_SV in curated_SVs:
						if "ND" in curated_SV:
							continue
						curated_SVs_list.append((curated_SV,len(set(curated_SVs[curated_SV]["read_id"])),set(curated_SVs[curated_SV]["read_id"])))
					sorted_curated_SVs_list=sorted(curated_SVs_list,key=lambda x:x[1],reverse=True)
					if len(sorted_curated_SVs_list)==0:
						depth_posA,depth_posB=get_depth(initial_SV,bam,chr_in_bam)
						curated_SV="ND"
						curated_SV_with_padding="ND"
						curated_SV_reads_num=0
						out_SV="%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n"%(sample,initial_SV,curated_SV,genome_version,str(curated_SV_reads_num),depth_posA,depth_posB,curated_SV_with_padding)
						out_SV_fh.write(out_SV)
						
					else:
						if len(sorted_curated_SVs_list)>3:
							sorted_curated_SVs_list=sorted_curated_SVs_list[0:3]
						for curated_SV_info in sorted_curated_SVs_list:
							curated_SV=curated_SV_info[0]
							curated_SV_reads_num=curated_SV_info[1]
							curated_SV_reads=curated_SV_info[2]
							curated_SV_with_padding=get_SV_padding(curated_SV)
							depth_posA,depth_posB=get_depth(curated_SV,bam,chr_in_bam)
							out_SV="%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n"%(sample,initial_SV,curated_SV,genome_version,str(curated_SV_reads_num),depth_posA,depth_posB,curated_SV_with_padding)
							out_SV_fh.write(out_SV)
							for read in curated_SV_reads:
								out_SV_read="%s\t%s\t%s\n"%(curated_SV,read,bam)
								out_SV_reads_fh.write(out_SV_read)
		for fn_to_delete in set(files_to_delete):
			os.remove(fn_to_delete)

out_SV_fh.close()
out_SV_reads_fh.close()


