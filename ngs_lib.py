"""all utility funcitons for the lab"""

import os
import re
import time 
import subprocess 
import random,string
import sys
from operator import itemgetter

is_debugging_phase = True

#PATH = "/research_jude/rgs01_jude/groups/xmagrp/projects/SensitiveDetection/xmagrp/code/mRNA/ZL/reference/"
PATH="/research/groups/xmagrp/projects/SensitiveDetection/xmagrp/code/dipseek/"

SAMBAMBA      = PATH + 'apps/sambamba'
SAMTOOLS      = PATH + 'apps/samtools'
BWA           = PATH + 'apps/bwa'
BLAT          = PATH + "apps/blat"
sequencerr    = PATH + 'apps/sequencerr_2.09 -mincov=1 -minmapq=15 -qcut=30 -qCutHard=20 -fcut=0.05 '
hg19_fa       = PATH + 'refs/hg19.gz'
hg38_fa       = PATH + 'refs/hg38.gz'
bwaindex_hg19 = PATH + 'refs/BWA-hg19/GRCh37-lite.fa'
bwaindex_hg38 = PATH + 'refs/BWA-hg38/GRCh38_no_alt_chr.fa'
bit_fa_hg38   = PATH + "refs/hg38.2bit"
bit_fa_hg19   = PATH + "refs/hg19.2bit"
gnomad_file   = PATH + 'refs/gnomAD_v3.1.2_hg38_genome/slim_gnomad.genomes.whole.hg38.hg19.concordant.min5.tsv'
exac_file     = PATH + 'refs/gnomAD_v3.1.2_hg38_genome/slim_gnomad.exomes.v4.1.sites.hg38.hg19.final.txt'
pecan_file    = PATH + 'refs/whitelist/lift/pecan_snv.lst_slim_dedup'
annotator     = PATH + 'script/annotation/IndelAnnotation_withNeoSpPred.v2.develop.py' # mnv anno update 03292024

def aa01_is_nucleotide(let: str) -> bool:
    """ determine if a letter is a nucleotide YF:XM:11202023"""
    """ YF updated to consider lower case nuc 11/22/23; tested """
    """ XM: this is great; however, we may need to make sure all nucs are uppercase before processing"""
    nuc_list = ['A', 'C', 'G', 'T', 'a', 'c', 'g', 't']
    result = False
    if(let in nuc_list): result = True
    return (result)

def aa02_marker_classification(marker: str) -> str:
    """ classify marker as snv, mnv, indel and sv (to do: valid nuc check); tested """
    """ TODO need to make sure all letters are valid for indels and for SV nts"""
    result = 'unknown'
    m = marker.split('.')
    if len(m) == 4:
        if('-' in m[2] or '-' in m[3]): result = 'indel'
        if len(m[2])==0 and len(m[3]) > 0: result = 'indel'
        if len(m[3])==0 and len(m[2]) > 0: result = 'indel'
        if len(m[2]) > 1 or len(m[3]) > 1:
            if len(m[2]) != len(m[3]): result = 'indel'
            else:
                if(m[2] != m[3]): result = 'mnv'
        if len(m[2]) == 1 and len(m[3]) == 1:
            if (aa01_is_nucleotide(m[2]) and aa01_is_nucleotide(m[3])): result = 'snv'
    if len(m) == 6 or len(m) == 7: result = 'sv'
    return (result)

def aa03_make_error_bin(lvl: int) -> list:
    """ to make bins for error rate levels YF:XM:11202023; tested """
    bin_edges = [0]
    nbin = 1
    if lvl>0:
        for i in range(lvl,0,-1):
            if(i>1):
                the_one = 1 / pow(10,i-1)
                the_two = 2 * the_one
                the_five= 5 * the_one
                bin_edges += [the_one,the_two,the_five]
            else:
                bin_edges += [1,2,5]
            nbin += 3
    bin_edges += [10, 30, 70, 100]
    nbin += 4
    return ([nbin, bin_edges])

def aa04_determine_error_bin_of_af(nM: int, nT: int, bin_edges: list) -> int:
    """to determine the bin of error ladder for an AF YF:XM:11202023"""
    """ tested; shall we remove the nbin input? As nbin = len(bin_edges)"""
    bin_number = -1
    nbin = len(bin_edges)
    if (nM==0):  bin_number = 0
    elif(nM>=nT): bin_number = nbin-1
    else:
        af = 100*nM/(nT+0.0001)
        for i in range(0,nbin-1,1):
            if(bin_edges[i]<af and af<=bin_edges[i+1]):
                bin_number = i+1
                break
    return (bin_number)

def aa05_does_bam_file_have_chr_tag(bam: str) -> bool:
    """ parse bam file to get chromosome format and size; tested """
    ret = False
    rlt = subprocess.check_output(["samtools view -H {} | grep SN | grep LN | head -1".format(bam)], shell=True).decode('utf-8').strip().split('\t')[1]
    chrom = rlt.split(':')[1]
    if chrom.startswith('chr'): ret = True
    return (ret)

def aa06_determine_genome_build_of_bam(bam: str) -> str:
    """ parse bam header to determine genome build, hg19/hg38; tested in hg38/hg19 """
    ret = 'unknown'
    rlt = subprocess.check_output(["samtools view -H {} | grep SN | grep LN | head -1".format(bam)], shell=True).decode('utf-8').strip().split('\t')[2]
    chrom_length = rlt.split(':')[1]
    if(chrom_length == '248956422'): ret = 'hg38'
    elif(chrom_length == '249250621'): ret = 'hg19'
    elif(chrom_length == '248387328'): ret = 'T2T'
    elif(chrom_length == '247249719'): ret = 'hg18'
    return(ret)

def aa07_check_jobs_status_one_round(job_dict:dict) -> list:
    """ to check job status, one round; tested 'is_debugging_phase' is not defined"""
    """ TODO to account for scenarios when memory exception happens"""
    """ TODO need to take care of SV/Indel Curation/Genotyping intermediate file issue"""
    success_counter = 0
    failure_counter = 0
    finish_counter  = 0
    good_job_dict   = {}
    for job in job_dict:
        oo_job = 'oo.' + job + '.oo'
        if os.path.exists(oo_job):
            ifh = open(oo_job,'r')
            #if is_debugging_phase: print(oo_job)
            success_of_this_job = 'unknown'
            finish_of_this_job  = 'unknown'
            for line in ifh:
                if 'Successfully completed' in line: 
                    success_of_this_job = 'yes'
                    finish_of_this_job  = 'yes'
                if 'Exited with' in line:
                    success_of_this_job = 'no'
                    finish_of_this_job  = 'yes'
                if 'job killed' in line:
                    success_of_this_job = 'no'
                    finish_of_this_job  = 'no'
            ifh.close()
            if success_of_this_job == 'yes':
                success_counter += 1
                good_job_dict[job] = 1
            if success_of_this_job == 'no':
                failure_counter += 1
            if finish_of_this_job == 'yes':
                finish_counter += 1
    return([failure_counter, success_counter, finish_counter, good_job_dict])

def aa08_job_manager_succeeded(job_dict:dict, only_check_finish: bool, initial_sleep_sec: int, wait_sec:int) -> list:
    """to monitor job status; they are submitted outside this function """
    #####os.system(f'sh {jobs_fn}')
    overall_sucess = True
    good_job_dict = {}    
    time.sleep(initial_sleep_sec)
    keep_checking = True
    check_round = 0
    n_jobs = len(job_dict)
    while keep_checking:
        time.sleep(wait_sec)
        check_round += 1
        failure_counter, success_counter, finish_counter, good_jobs = aa07_check_jobs_status_one_round(job_dict)
        if is_debugging_phase:
            print(f'job checking round: {check_round}; finish: {finish_counter}/{n_jobs}') 
        for job in good_jobs: good_job_dict[job] = 1
        if success_counter == n_jobs:
            keep_checking = False
            overall_sucess = True
            break
        if only_check_finish:
            if finish_counter == n_jobs:
                keep_checking = False
                overall_sucess = True
                break
        else:
            if failure_counter>0:
                keep_checking = False
                overall_sucess = False
                break
    time.sleep(wait_sec)
    return([overall_sucess, good_job_dict])

def aa09_name_map_parser(namemapFile: str) -> list:
    """ process the 5-column namemap file: """
    """ IntakeID|bam|TruID|UseTag|TruTag """
    Intake_Sample_Duplicate = 0
    Perfect_Sample_Duplicate = 0
    sample_tag_dict = {}
    intake2perfectname_dict = {}
    perfectname2intake_dict = {}
    ifh = open(namemapFile,'r')
    #if is_debugging_phase: print(namemapFile)
    for line in ifh:
        if '.bam' in line:
            intake_sample, bampath, perfect_pid, perfect_sample_tag = line.strip().split('\t')[:4]
            if perfect_sample_tag not in sample_tag_dict:
                sample_tag_dict[perfect_sample_tag] = 1
            if intake_sample not in intake2perfectname_dict:
                intake2perfectname_dict[intake_sample] = {}
                intake2perfectname_dict[intake_sample]['bam_path'] = bampath
                intake2perfectname_dict[intake_sample]['perfect_pid'] = perfect_pid
                intake2perfectname_dict[intake_sample]['perfect_sample_tag'] = perfect_sample_tag
            else:
                Intake_Sample_Duplicate += 1
            if perfect_pid not in perfectname2intake_dict:
                perfectname2intake_dict[perfect_pid] = {}
            if perfect_sample_tag not in perfectname2intake_dict[perfect_pid]:
                perfectname2intake_dict[perfect_pid][perfect_sample_tag]={}
                perfectname2intake_dict[perfect_pid][perfect_sample_tag]['intake_sample'] = intake_sample
                perfectname2intake_dict[perfect_pid][perfect_sample_tag]['bam_path'] =  bampath
            else:
                Perfect_Sample_Duplicate += 1
    ifh.close()
    return ([Intake_Sample_Duplicate, Perfect_Sample_Duplicate, sample_tag_dict, intake2perfectname_dict, perfectname2intake_dict])

def aa10_table_joining(to_table: str, to_idx_cols: list, from_table: str, from_idx_cols: list, wanted_from_table_cols: str, delim_for_indexing: str, output_file: str) -> list:
    """ join two tables based on specified columns and write the output to a file."""
    """ now accounted for the possibility where from_table has duplicates """
    """ if to_table has too many rows (>10M), need to split it before using this function. """
    """ wanted_from_table_cols input is ':' delimited; use -1 if all cols wanted """
    success = True
    message = ''
    to_table_dict = {}
    wanted_from_table_cols = [int(item) for item in wanted_from_table_cols.split(':') if item]
    if -1 in wanted_from_table_cols:
        ifn = open(from_table, 'r')
        wanted_from_table_cols = list(range(len(ifn.readline().strip().split('\t'))))
        ifn.close()
    for col_idx in from_idx_cols:
        if col_idx in wanted_from_table_cols: wanted_from_table_cols.remove(col_idx)
    min_needed_column_number_of_from_table = max(from_idx_cols + wanted_from_table_cols)
    min_needed_column_number_of_to_table = max(to_idx_cols)
    if(min_needed_column_number_of_to_table == None or min_needed_column_number_of_from_table==None):
        success = False
        message = 'index columns not defined well'
    if success:
        filler_data = '\t'.join(['-1' for i in range(len(wanted_from_table_cols))])
        ifh = open(to_table, 'r')
        for line in ifh:
            rec = line.strip().split('\t')
            if len(rec) >= min_needed_column_number_of_to_table:
                to_idx_key = delim_for_indexing.join([rec[idx] for idx in to_idx_cols])
                if to_idx_key not in to_table_dict:
                    to_table_dict[to_idx_key] = {}
                if line.strip() not in to_table_dict[to_idx_key]:
                    to_table_dict[to_idx_key][line.strip()] = ''
            else: 
                success = False
                message = 'to_table has fewer columns than needed'
                break
        ifh.close()
        nn = 0
        ifh = open(from_table, 'r')
        for line in ifh:
            rec = line.strip().split('\t')
            nn += 1
            if len(rec)>=min_needed_column_number_of_from_table:
                from_idx_key = delim_for_indexing.join([rec[idx] for idx in from_idx_cols])
                if from_idx_key in to_table_dict:
                    wanted_data_in_from_table = '\t'.join([rec[idx] for idx in wanted_from_table_cols])
                    for to_table_line in to_table_dict[from_idx_key]:
                        if to_table_dict[from_idx_key][to_table_line] == '':
                            to_table_dict[from_idx_key][to_table_line] = wanted_data_in_from_table
                        else:
                            to_table_dict[from_idx_key][to_table_line] += ';' + wanted_data_in_from_table   #in case from_table has duplicates
            else:
                success = False
                message = 'from_table has fewer columns than needed'
                break
        ifh.close()
        ofh = open(output_file, 'w')
        for idx_key in to_table_dict:
            for to_table_line in to_table_dict[idx_key]:
                new_data = to_table_dict[idx_key][to_table_line]
                if(new_data == ''): new_data = filler_data
                ofh.write(to_table_line + '\t' + new_data + '\n')
        ofh.close()
    return([success, message])

def aa11_obtain_error_prfile_for_all_marker(ifile: str, ofile: str, lvl: int, minDepth: int)-> None:
    """parse the count file in a block fasion"""
    """marker file has 4 columns: SampleID|Marker|nM|nT"""
    n_bin, bin_edge = aa03_make_error_bin(lvl)
    bin_counts = [0 for i in range(n_bin)]
    n_good_sample = 0
    current_marker = "Nothing"
    for i in range(0,n_bin): bin_counts[i] = 0
    ofh = open(ofile, 'w')
    ofh.write('Curated.Marker\tnGoodSample\tP0')
    for i in range(0,n_bin-1):
        ofh.write('\tP_'+str(bin_edge[i])+'_'+str(bin_edge[i+1]))
    ofh.write('\n')
    ifh = open(ifile, 'r')
    for line in ifh:
        info = line.strip().split('\t')
        mkr = ''
        if (len(info)==4):
            smp, mkr, nM, nT  = info
        elif (len(info)==5):
            mkr, smp, nRa, nRb, nM = info
            nT = int((int(nRa) + int(nRb)) / 2) + int(nM)
        if mkr != '':
            if(mkr != current_marker):
                if(current_marker != "Nothing"):
                    ofh.write(current_marker + '\t' + str(n_good_sample))
                    for i in range(0,n_bin): ofh.write('\t'+str(bin_counts[i]))
                    ofh.write('\n')
                for i in range(0,n_bin): bin_counts[i] = 0
                n_good_sample = 0
            if (int(nT)>=minDepth):
                the_bin_number = aa04_determine_error_bin_of_af(int(nM), int(nT), bin_edge)
                if the_bin_number>=0: bin_counts[the_bin_number] += 1
                n_good_sample += 1
            current_marker = mkr

    ifh.close()
    ofh.write(current_marker + '\t' + str(n_good_sample))
    for i in range(0,n_bin): ofh.write('\t'+str(bin_counts[i]))
    ofh.write('\n')
    ofh.close()

def aa12_sorted_sample_tag(tag_list: list, germline_tags: list = ['G', 'H']):
    """Sort the tag order alphabetically and move specified tags behind, also sorting numerically within each tag group."""
    def tag_sort_key(item: str):
        match = re.match(r"([a-zA-Z]+)([0-9]+)", item)
        if match: letters, number = match.groups()
        else: letters, number = item, 0
        idx = any(item.startswith(tag) for tag in germline_tags)
        return (int(idx), letters, int(number))
    return sorted(tag_list, key=tag_sort_key)

def aa13_marker_is_good_curated_indel(marker: str) -> str:
    """ classify marker as snv, mnv, indel and sv (to do: valid nuc check); tested """
    """ TODO need to make sure all letters are valid for indels and for SV nts"""
    result = False
    m = marker.split('.')
    if len(m) == 4:
        if len(m[2])>=2 and len(m[3])>=2 and abs(len(m[2])-len(m[3]))>0:
            if (not 'N' in m[2]+m[3]) and (not 'n' in m[2]+m[3]): result = True
    return (result)

def aa14_split_big_file_to_small_files(bigfile_name: str, batchSize: int, split_prefix: str) -> dict:
    """ to split Marker or Sample|Marker into smaller batch files and return file list """    
    os.system(f'split -l {batchSize} {bigfile_name} {split_prefix}')
    tmp_split_collection = bigfile_name + '.split'
    os.system(f'ls {split_prefix}* > {tmp_split_collection}')
    split_filename_dict = {}
    ifh = open(tmp_split_collection, 'r')
    for line in ifh:
        rcd = line.strip()
        split_filename_dict[rcd]=1
    ifh.close()
    os.system(f'rm {tmp_split_collection}')
    return (split_filename_dict)

def aa15_merge_small_files_to_bigfile(split_fn_dict: dict, big_out_fn: str, has_header: bool = False):
    ofh = open(big_out_fn, 'w')
    file_counter = 1
    for ifn in split_fn_dict:
        ifh = open(ifn, 'r')
        if has_header:
            header = ifh.readline().strip()
            if file_counter == 1: ofh.write(header + '\n')
        for line in ifh:
           rcd = line.strip()
           if len(rcd)>1: 
               ofh.write(rcd + '\n')
        ifh.close()
        file_counter += 1
    ofh.close()

def aa16_extract_one_region_from_one_bam(bam_path, chrom, start, end, nLine_buffer_cut, min_read_len, o_fastq_handle):
    """ using samtools. when too many reads, will write to tmp file first """
    if int(start)<=int(end):
        sstart = str(start)
        eenddd = str(end)
        sample_name = bam_path.split('/')[-1].split('.')[0]
        if int(start)<1: sstart = '1'
        cmd = f'{SAMTOOLS} view ' + bam_path + ' ' + chrom + ':' + sstart + '-' + eenddd + ' | wc -l'
        proc = subprocess.run(cmd, shell=True, capture_output=True, text=True, check=True)
        nLine = proc.stdout.splitlines()[0]
        if int(nLine)<nLine_buffer_cut:
            cmd = f'{SAMTOOLS} view ' + bam_path + ' ' + chrom + ':' + sstart + '-' + eenddd
            proc = subprocess.run(cmd, shell=True, capture_output=True, text=True, check=True)
            for line in proc.stdout.splitlines():
                rcd = line.strip().split('\t')
                read_name,flag,chrom,pos,mapq,cigar,same_seg,next_start,distance,read_seq,read_phred=rcd[:11]
                if(len(read_seq)==len(read_phred) and len(read_seq)>=int(min_read_len)):
                    read1or2 = 'read1' # determine_read1or2(flag)
                    new_name='@' + sample_name + ':' + read1or2 + ':' + read_name
                    o_fastq_handle.write(new_name+'\n'+read_seq+'\n+\n'+read_phred+'\n')
        else:
            ran_string=''.join([random.choice(string.ascii_letters) for x in range(8)])
            temp_sam_fn = '.'.join([sample_name,chrom,sstart,eenddd,ran_string])
            cmd = f'{SAMTOOLS} view ' + bam_path + ' ' + chrom + ':' + sstart + '-' + eenddd + ' > ' + temp_sam_fn
            subprocess.run(cmd, shell=True, capture_output=True, text=True, check=True)
            ifh = open(temp_sam_fn,'r')
            for line in ifh:
                rcd = line.strip().split('\t')
                read_name,flag,chrom,pos,mapq,cigar,same_seg,next_start,distance,read_seq,read_phred=rcd[:11]
                if(len(read_seq)==len(read_phred) and len(read_seq)>=int(min_read_len)):
                    read1or2 = 'read1' # determine_read1or2(flag)
                    new_name='@' + sample_name + ':' + read1or2 + ':' + read_name
                    o_fastq_handle.write(new_name+'\n'+read_seq+'\n+\n'+read_phred+'\n')
            ifh.close()
            os.system(f'rm {temp_sam_fn}')

def aa17_reverse_complement_nucleotides(dnaseq: str) -> str:
    """ to obtain the reverse compliment of nucleotide sequence"""
    cpm = {'A':'T','T':'A','C':'G','G':'C','N':'N','':'','-':'-'}
    rev = dnaseq[::-1]
    revcmp = ''
    for i in rev: revcmp += cpm[i]
    return (revcmp)

def aa18_extract_seq_segment_from_genome(genome_fa_file, chrom, start, end):
    """to extract sequence from genome using locations"""
    cmd = f'{SAMTOOLS} faidx '+genome_fa_file+' '+chrom+':'+str(start)+'-'+str(end)
    ret = ''.join(str(subprocess.check_output(cmd,shell=True),"utf-8").strip().split("\n")[1:]).upper()
    return (ret)

def aa19_snvindelmnv_to_contig(marker,len_flank,genome_fa_file) -> dict:
    """ make ref and alt contigs using snv/indel/mnv """
    chrom,pos,ref_allele,mut_allele=marker.split('.')
    chrom='chr'+chrom.replace('chr','')
    startA =int(pos)-len_flank
    endA   =int(pos)-1
    startB =int(pos)+len(ref_allele)
    endB   =int(pos)+len(ref_allele)+len_flank-1        
    left_ref_seq    =aa18_extract_seq_segment_from_genome(genome_fa_file, chrom, startA, endA)
    right_ref_seq   =aa18_extract_seq_segment_from_genome(genome_fa_file, chrom, startB, endB)
    ret={'left_ref_seq':left_ref_seq,'ref_allele':ref_allele,'mut_allele':mut_allele,'right_ref_seq':right_ref_seq,'left_chrom':chrom,'left_start':startA,'left_end':endA,'right_chrom':chrom,'right_start':startB,'right_end':endB}
    return(ret)

def aa20_sv_itd_to_contig(marker, len_flank, genome_fa_file) -> dict:
    """ make ref and alt contigs using sv/itd """
    chrA,posA,strandA,chrB,posB,strandB,insertion=marker.split('.')
    chrA='chr'+chrA.replace('chr','')
    chrB='chr'+chrB.replace('chr','')
    if strandA=='+':
        startA=int(posA)-len_flank+1
        endA=int(posA)
    else:
        startA=int(posA)
        endA=int(posA)+len_flank-1
    if strandB=="+":
        startB=int(posB)
        endB=int(posB)+len_flank-1
    else:
        startB=int(posB)-len_flank+1
        endB=int(posB)
    left_ref_seq    =aa18_extract_seq_segment_from_genome(genome_fa_file, chrA, startA, endA)
    right_ref_seq   =aa18_extract_seq_segment_from_genome(genome_fa_file, chrB, startB, endB)
    left_genome_pos=[startA,endA]
    right_genome_pos=[startB,endB]
    if strandA=="-":
        left_ref_seq =aa17_reverse_complement_nucleotides(left_ref_seq)
        left_genome_pos=[endA,startA]
    if strandB=="-":
        right_ref_seq =aa17_reverse_complement_nucleotides(right_ref_seq)
        right_genome_pos=[endB,startB]
    ret={'left_ref_seq':left_ref_seq,'ref_allele':'','mut_allele':insertion,'right_ref_seq':right_ref_seq,'left_chrom':chrA,'left_start':left_genome_pos[0],'left_end':left_genome_pos[1],'right_chrom':chrB,'right_start':right_genome_pos[0],'right_end':right_genome_pos[1]}
    return(ret)

def aa21_marker_to_contig_omni(marker: str, len_flank: int, genome_fa_file: str) -> list:
    """ to convert marker to contig representation for SNV/Indel/MNV/SV """
    marker_type = aa02_marker_classification(marker)
    success_status = False
    result = 'unknown'
    if marker_type in ['snv', 'indel', 'mnv']:
        result = aa19_snvindelmnv_to_contig(marker, len_flank, genome_fa_file)
        success_status = True
    if marker_type in ['sv']:
        result = aa20_sv_itd_to_contig(marker, len_flank, genome_fa_file)
        success_status = True
    return([success_status, result])

def aa22_is_raw_call_file_valid(rawcall_fn: str):
    validity = True
    message  = ''
    nGood = 0
    if(os.path.exists(rawcall_fn)):
        ifh = open(rawcall_fn,'r')
        for line in ifh:
            rcd = line.strip().split('\t')
            if(len(rcd)>=2):
                intake_sample, ori_marker = rcd[:2]
                if (ori_marker.startswith('chr')): nGood += 1                
            else:
                validity = False                
                ifh.close()
                message = 'some rows have less then 2 columns'
                break
        ifh.close()
    else: 
        validity = False
        message = 'file does not exist'
    if nGood == 0: message = 'no valid records'
    return([validity, message])

def aa23_is_curated_call_file_valid(rawcall_fn: str):
    validity = True
    message = ''
    nGood = 0
    if(os.path.exists(rawcall_fn)):
        ifh = open(rawcall_fn,'r')
        for line in ifh:
            rcd = line.strip().split('\t')
            if(len(rcd)>=3):
                intake_sample, ori_marker, curated_marker = rcd[:3]
                if (ori_marker.startswith('chr')): nGood += 1
            else:
                validity = False
                ifh.close()
                message = 'some rows have less then 2 columns'
                break
        ifh.close()
    else: 
        validity = False
        message = 'file does not exist'
    if nGood == 0: message = 'no valid records'
    return([validity,message])


def aa24_self_genotyping(marker_file: str, sample_column: int, is_sample_perfect_pid: bool, snv4column: int, big_count_file: str, platform: str, namemapFile: str, meta_out_file:str)->bool:
    """ using marker file (perfect_sample|marker) or (intake_sample|marker), convert big count file to horizontal format"""
    """ samples not show up in namemap are dropped """
    success = True
    intake_dup, perfect_dup, tag_dict, intake_dict, perfect_dict = aa09_name_map_parser(namemapFile)
    tag_list = list(tag_dict.keys())
    tag_list = aa12_sorted_sample_tag(tag_list) #
    meta_dict = {}
    ifh = open(marker_file, 'r')
    marker_query_dict = {} #
    for line in ifh:
        if 'Marker' not in line.strip(): 
            samplename = line.strip().split('\t')[sample_column]
            marker     = line.strip().split('\t')[snv4column]
            if snv4column == 2: #
                orig_marker = line.strip().split('\t')[1] #
            else: orig_marker = marker #
            if samplename not in marker_query_dict: marker_query_dict[samplename] = {} #
            if marker not in marker_query_dict[samplename]: #
                marker_query_dict[samplename][marker] = orig_marker #
            elif snv4column == 2: #
                old_dist = aa29_marker_distance(marker, marker_query_dict[samplename][marker]) #
                new_dist = aa29_marker_distance(marker, orig_marker) #
                if new_dist < old_dist: marker_query_dict[samplename][marker] = orig_marker #
            if is_sample_perfect_pid:
                perfect_pid = samplename
                if perfect_pid not in meta_dict:
                    meta_dict[perfect_pid] = {}
                if marker not in meta_dict[perfect_pid]:
                    meta_dict[perfect_pid][marker] = {}
            else:
                if samplename in intake_dict:
                    perfect_pid = intake_dict[samplename]['perfect_pid']
                    if perfect_pid not in meta_dict:
                        meta_dict[perfect_pid] = {}
                    if marker not in meta_dict[perfect_pid]:
                        meta_dict[perfect_pid][marker] = {}
    ifh.close()

    ifh = open(big_count_file, 'r')
    for line in ifh:
        info = line.strip().split('\t')
        if len(info) == 4 or len(info)==5: #need to be well formated
            if len(info) == 4: #SNV/Indel
                intake_sample, marker, nM, nT = info
            elif len(info) == 5: #SV
                marker, intake_sample, nRa, nRb, nM = info
                nT = int((int(nRa) + int(nRb)) / 2) + int(nM)
            if intake_sample in intake_dict:
                perfect_pid = intake_dict[intake_sample]['perfect_pid']
                perfect_tag = intake_dict[intake_sample]['perfect_sample_tag']
                if perfect_pid in meta_dict:
                    if marker in meta_dict[perfect_pid]:
                        meta_dict[perfect_pid][marker][perfect_tag] = str(nM) + '\t' + str(nT)
    ifh.close()

    ofh = open(meta_out_file,'w')
    ofh.write('PID\tOriginal.Marker\tCurated.Marker')
    for tag in tag_list:
        ofh.write('\tnM_'+tag+'_'+platform+'\tnT_'+tag+'_'+platform)
    ofh.write('\n')

    filler_text_no_sample = '-1\t-1'
    filler_text_has_sample = '0\t0'
    for pid in meta_dict:
        for marker in meta_dict[pid]:
            orig_marker = marker_query_dict[pid][marker] #
            ofh.write('\t'.join([pid, orig_marker, marker])) #
            for tag in tag_list:
                output_str = '\t' + filler_text_no_sample
                if pid in perfect_dict:
                    if tag in perfect_dict[pid]:
                        output_str = '\t'  + filler_text_has_sample
                if tag in meta_dict[pid][marker]: 
                    output_str = '\t'  + meta_dict[pid][marker][tag]
                ofh.write(output_str)
            ofh.write('\n')
    ofh.close()
    return(success)

def aa25_make_horizontal_file_from_big_count_files(
    metafile:str,
    marker_type: str,
    annotationfile: str,
    selfgeno_all2all_fn: str, 
    controlgeno_all2all_fn: str, 
    nameMapFile: str, 
    outputFile: str, 
    platform: str,
    hg_build: str
    ) -> None:
    """need to check file existance"""
    self_horizontal_file = outputFile + '.self_hori'
    control_horizontal_file = outputFile + '.ctrl_hori'
    genotyping_result_file = outputFile + '.geno_hori'
    anno_horizontal_file = outputFile + '.anno_hori'
    gnomad_horizontal_file = outputFile + '.gnomad_hori'
    exac_horizontal_file = outputFile + '.exac_hori'
     
    print('self genotyping')
    aa24_self_genotyping(metafile, 1-1, True, 3-1, selfgeno_all2all_fn, platform, nameMapFile, self_horizontal_file)
    print('error profiling')
    if controlgeno_all2all_fn.upper() != 'NULL':
        aa11_obtain_error_prfile_for_all_marker(controlgeno_all2all_fn, control_horizontal_file, 1, 20)
        print('genotyping table joining')
        aa10_table_joining(self_horizontal_file, [3-1], control_horizontal_file, [1-1], '-1', '\t', genotyping_result_file)
    else:
        os.system(f'cp {self_horizontal_file} {genotyping_result_file}')
    
    print('adding annotation')
    if marker_type.upper() in ['INDEL', 'MNV', 'SNV']:
        aa10_table_joining(genotyping_result_file, [3-1], annotationfile, [2-1], '2:3:4:5:6', ';', anno_horizontal_file)
    elif marker_type.upper() == 'SV':
        aa10_table_joining(genotyping_result_file, [3-1], annotationfile, [1-1], ':'.join(list(map(str,range(3,27)))), ';', anno_horizontal_file)
    print('checking gnomad')
    gnomad_cols = [1] if hg_build == 'hg19' else [0]
    aa10_table_joining(anno_horizontal_file, [3-1], gnomad_file, gnomad_cols, '2:3', ';', gnomad_horizontal_file)
    print('checking ExAC')
    exac_cols = [1] if hg_build == 'hg19' else [0]
    aa10_table_joining(gnomad_horizontal_file, [3-1], exac_file, exac_cols, '2:3', ';', exac_horizontal_file)
    print('checking PeCan whitelist')
    pecan_cols = [1] if hg_build == 'hg19' else [0]
    aa10_table_joining(exac_horizontal_file, [3-1], pecan_file, pecan_cols, '2:3', ';', outputFile)
    
    os.system(f'rm {self_horizontal_file} {control_horizontal_file} {genotyping_result_file} {anno_horizontal_file} {gnomad_horizontal_file} {exac_horizontal_file}')

def aa26_run_snvindelmnv_annotation(marker_file: str, snv4column: int, hg_build: str, annotation_ofile: str, projectID: str) -> None:
    """run snv/indel/mnv annotation using Zhikai code"""
    before_sort = annotation_ofile + '_' +  projectID + '.before_sort'
    after_sort =  annotation_ofile + '_' +  projectID + '.after_sort'
    marker_dict = {}
    ifh = open(marker_file,'r')
    for line in ifh:
        if 'Marker' not in line:
            marker = line.strip().split('\t')[snv4column]
            marker_dict[marker] = 1
    ifh.close()
    ofh = open(before_sort,'w')
    for marker in marker_dict:
        ofh.write(f'{marker}\n')
    ofh.close()
    os.system(f'cat {before_sort} | sort -t. -k1,1 -k2,2n > {after_sort}; rm {before_sort}')
    ifh = open(after_sort,'r')
    ofh = open(before_sort,'w')
    for line in ifh:
        ofh.write('\t'.join(['SJID',line.strip()])+'\n')
    ifh.close()
    ofh.close()
    os.system(f'mv {before_sort} {after_sort}')
    print('doing annotation')
    os.system(f'python {annotator} {hg_build} {after_sort} > {annotation_ofile}')
    os.system(f'rm {after_sort}')
    print('finished annotation')

def aa27_get_revcmp(strs: str) -> str:
    """ reverse complement the nucleotide sequence """
    cpm = {"A":"T","T":"A","C":"G","G":"C","N":"N","":"","-":"-"}
    rev = strs[::-1]
    revcmp = ""
    for i in rev: revcmp += cpm[i]
    return (revcmp)

def aa28_adjust_SV(SV: str) -> str:
    """ correct the format of sv """
    if SV.count(".")==5: SV=SV+"."
    cmp_strand={"+":"-","-":"+"}
    SV_lst=SV.split(".")
    chrA, posA, strandA, chrB, posB, strandB, insertion = SV_lst[:7]
    if "chr" not in chrA: chrA="chr"+chrA
    if chrA=="chrX": chrA="chr23"
    if chrA=="chrY": chrA="chr24"
    if "chr" not in chrB: chrB="chr"+chrB
    if chrB=="chrX": chrB="chr23"
    if chrB=="chrY": chrB="chr24"

    adjusted_chrA, adjusted_posA, adjusted_strandA = [chrA, posA, strandA]
    adjusted_chrB, adjusted_posB, adjusted_strandB = [chrB, posB, strandB]
    adjusted_insertion = insertion

    if (int(chrA.replace("chr",""))>int(chrB.replace("chr",""))) or (int(chrA.replace("chr",""))==int(chrB.replace("chr","")) and int(posA)>int(posB)):
        adjusted_chrA, adjusted_posA, adjusted_strandA = [chrB, posB, cmp_strand[strandB]]
        adjusted_chrB, adjusted_posB, adjusted_strandB = [chrA, posA, cmp_strand[strandA]]
        adjusted_insertion = aa27_get_revcmp(insertion)

    if adjusted_chrA=="chr23": adjusted_chrA="chrX"
    if adjusted_chrA=="chr24": adjusted_chrA="chrY"
    if adjusted_chrB=="chr23": adjusted_chrB="chrX"
    if adjusted_chrB=="chr24": adjusted_chrB="chrY"

    adjusted_SV="%s.%s.%s.%s.%s.%s.%s"%(adjusted_chrA,adjusted_posA,adjusted_strandA,adjusted_chrB,adjusted_posB,adjusted_strandB,adjusted_insertion)
    return (adjusted_SV)

def aa29_marker_distance(marker1: str, marker2: str) -> int:
    """ calculate the distance between 2 snv4 markers """
    distance = -1
    marker1_type = aa02_marker_classification(marker1)
    marker2_type = aa02_marker_classification(marker2)
    if (marker1_type in ['indel', 'mnv', 'snv']) and (marker2_type in ['indel', 'mnv', 'snv']):
        marker1_info = marker1.split('.')
        marker2_info = marker2.split('.')
        if marker1_info[0] == marker2_info[0]:
            distance = abs(int(marker1_info[1]) - int(marker2_info[1]))
    elif (marker1_type == 'sv') and (marker2_type == 'sv'):
        marker1_adj = aa28_adjust_SV(marker1)
        marker2_adj = aa28_adjust_SV(marker2)
        marker1_info = marker1_adj.split('.')
        marker2_info = marker2_adj.split('.')
        chrom1a, pos1a, strand1a, chrom1b, pos1b, strand1b = marker1_info[:6]
        chrom2a, pos2a, strand2a, chrom2b, pos2b, strand2b = marker2_info[:6]
        if chrom1a == chrom2a and chrom1b == chrom2b:
            distance = abs(int(pos1a) - int(pos2a)) + abs(int(pos1b) - int(pos2b))
    return (distance)

def ba01_bam_slicer(input_bam_file: str, bed_file: str, sliced_bam: str, remove_singletons : bool = False) -> None: 
    """ slice the original bam file based on bed coordinates """
    os.system(f'{SAMBAMBA} view -h -L {bed_file} -f bam -o {sliced_bam} {input_bam_file}')
    if remove_singletons:
        os.system(f"{SAMTOOLS} view {sliced_bam} | cut -f1 | sort | uniq -c | awk '$1 == 1 {{print $2}}' > {sliced_bam}_unpaired_read_names.txt")
        os.system(f'{SAMTOOLS} view -h {sliced_bam} | fgrep -v -f {sliced_bam}_unpaired_read_names.txt > {sliced_bam}_tmp')
        os.system(f'{SAMTOOLS} sort -n -o {sliced_bam}_tmp_sorted {sliced_bam}_tmp')
        os.system(f'mv {sliced_bam}_tmp_sorted {sliced_bam}')
        os.system(f'rm {sliced_bam}_unpaired_read_names.txt {sliced_bam}_tmp')
        
def ba02_bamsam2fastq(sliced_bamsam: str, run_prefix: str, mode: str = 'single-end') -> list: 
    """ use samtools to convert sorted bam|sam file to single-end fastq file """
    output_filenames = []
    if mode == 'single-end':
        R0 = f'{run_prefix}.R0'
        os.system(f'{SAMTOOLS} fastq {sliced_bamsam} > {R0}.fastq')
        os.system(f'gzip {R0}.fastq')
        output_filenames = [f'{R0}.fastq.gz']
    elif mode == 'paired-end':
        R1, R2 = [f'{run_prefix}.R1', f'{run_prefix}.R2']
        os.system(f'{SAMTOOLS} fastq -1 {R1}.fastq -2 {R2}.fastq -s {run_prefix}_singleton.fastq -0 /dev/null -n {sliced_bamsam}')
        os.system(f'gzip {R1}.fastq')
        os.system(f'gzip {R2}.fastq')
        output_filenames = [f'{R1}.fastq.gz', f'{R2}.fastq.gz']
    return (output_filenames)

def ba03_slice_bam_2_fastq(input_bam: str, bed_file: str, run_prefix: str, mode: str = 'paired-end') -> list:
    """ extract reads from bam file and convert to fastq """
    sliced_bam = f'{run_prefix}_sliced.bam'
    ba01_bam_slicer(input_bam, bed_file, sliced_bam, remove_singletons = True)
    fastqfile_list = ba02_bamsam2fastq(sliced_bam, run_prefix, mode = mode)
    cmd = f'rm {sliced_bam} {sliced_bam}.bai'
    os.system(cmd)
    return (fastqfile_list)

def ba04_map_fastq_unpaired(fastq_R0: str, bwaindex: str, oprefix: str, sort_n_index: bool = True) -> None:
    " use bwa to align unpaired fastq file and return sorted & indexed bam file "
    cmd = f'{BWA} aln -t 4 {bwaindex} {fastq_R0} > {oprefix}_0.sai'
    os.system(cmd)
    cmd = f'{BWA} samse {bwaindex} {oprefix}_0.sai {fastq_R0} > {oprefix}.sam; rm {oprefix}_0.sai'
    os.system(cmd)
    if sort_n_index:  
        cmd = f'{SAMTOOLS} view -b {oprefix}.sam -o {oprefix}.bam; rm {oprefix}.sam'
        os.system(cmd)
        cmd = f'{SAMTOOLS} sort {oprefix}.bam -o {oprefix}_sorted.bam'
        os.system(cmd)
        cmd = f'mv {oprefix}_sorted.bam {oprefix}.bam'
        os.system(cmd)
        cmd = f'{SAMTOOLS} index {oprefix}.bam'
        os.system(cmd)

def ba05_map_fastq_paired(fastq_R1: str, fastq_R2: str, bwaindex: str, oprefix: str, sort_n_index: bool = True) -> None:
    " use bwa to align paired-end fastq file and return sorted & indexed bam file "
    cmd = f'{BWA} aln -t 4 {bwaindex} {fastq_R1} > {oprefix}_1.sai'
    os.system(cmd)
    cmd = f'{BWA} aln -t 4 {bwaindex} {fastq_R2} > {oprefix}_2.sai'
    os.system(cmd)
    cmd = f'{BWA} sampe {bwaindex} {oprefix}_1.sai {oprefix}_2.sai {fastq_R1} {fastq_R2} > {oprefix}.sam'
    os.system(cmd)
    cmd = f"rm {oprefix}_1.sai {oprefix}_2.sai {fastq_R1.replace('.R1.fastq.gz', '_singleton.fastq')}"
    os.system(cmd)
    if sort_n_index:
        cmd = f'{SAMTOOLS} view -b {oprefix}.sam -o {oprefix}.bam; rm {oprefix}.sam'
        os.system(cmd)
        cmd = f'{SAMTOOLS} sort {oprefix}.bam -o {oprefix}_sorted.bam'
        os.system(cmd)
        cmd = f'mv {oprefix}_sorted.bam {oprefix}.bam'
        os.system(cmd)
        cmd = f'{SAMTOOLS} index {oprefix}.bam'
        os.system(cmd)

def ba06_run_sequencerr(bam: str, ofn: str) -> None:
    """ run sequencerr genotyping with an input bam file """
    cmd = f'{sequencerr} {bam} {ofn}'
    os.system(cmd)

def ba07_if_sv_relevant(sv1: str, sv2: str, distance_cutoff: int = 20) -> bool:
    """ check if two sv markers are relevant (same) """
    is_relevant = False
    strand_dict = {'+': '-', '-': '+'}
    chrom1a, pos1a, strand1a, chrom1b, pos1b, strand1b = sv1.split('.')[:6]
    chrom2a, pos2a, strand2a, chrom2b, pos2b, strand2b = sv2.split('.')[:6]
    if chrom1a == chrom2a and chrom1b == chrom2b and abs(int(pos1a)-int(pos2a))<distance_cutoff and abs(int(pos1b)-int(pos2b))<distance_cutoff and strand1a==strand2a and strand1b==strand2b:
        is_relevant = True
    elif chrom1a == chrom2b and chrom1b == chrom2a and abs(int(pos1a)-int(pos2b))<distance_cutoff and abs(int(pos1b)-int(pos2a))<distance_cutoff and strand1a==strand_dict[strand2b] and strand1b==strand_dict[strand2a]:
        is_relevant = True
    return (is_relevant)

def ba08_bwa_index_reference(fasta_file: str) -> list:
    """ index the reference file for bwa sequence alignment """
    cmd = f'{BWA} index {fasta_file}'
    os.system(cmd)
    return ([f'{fasta_file}.amb', f'{fasta_file}.ann', f'{fasta_file}.bwt', f'{fasta_file}.fai', f'{fasta_file}.pac', f'{fasta_file}.sa'])

def ba09_cigar_match_only(cigar: str) -> bool:
    """ check if the cigar string only contains match (M)"""
    pattern = r"[^\dM]"
    return re.search(pattern, cigar) is None

def ba10_seq_to_singleend_fastq(seq_dict: dict, output_fastq: str) -> None:
    """ create a pseudo fastq file for sequences """
    ofn = open(output_fastq, 'w')
    for name in seq_dict:
        header = '@' + name
        seq = seq_dict[name]
        qual = 'E' * len(seq)
        ofn.write('\n'.join([header, seq, '+', qual]) + '\n')
    ofn.close()

def ba11_seq_to_fasta(seq_dict: dict, output_fasta: str) -> None:
    """ create a fasta file for sequences """
    ofn = open(output_fasta, 'w')
    for name in seq_dict:
        header = '>' + name
        seq = seq_dict[name]
        ofn.write('\n'.join([header, seq]) + '\n')
    ofn.close()
