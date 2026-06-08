import os
import sys
import glob
import random
import time
import argparse
import subprocess
import pandas as pd
from ngs_lib import *


CODE_DIR = '/research/groups/xmagrp/projects/SensitiveDetection/xmagrp/code/dipseek/'

indel_curator  = CODE_DIR + 'script/curation/indel_curation_with_MNV.py'
sv_curator     = CODE_DIR + 'script/curation/SV_curation_v7.py'
annotator      = CODE_DIR + 'script/annotation/IndelAnnotation_withNeoSpPred.v3.develop.py' # mnv annotation updated 03292024
sv_annotator   = CODE_DIR + 'script/annotation/svCnvAnnotation.py'
svindel_typer  = CODE_DIR + 'script/genotyper/SNVInDelMNVSVGenotyping_2modes_DNA_v2.py'
chromsize_hg19 = CODE_DIR + 'script/ai_check/chr_sizes_hg19.txt'
chromsize_hg38 = CODE_DIR + 'script/ai_check/chr_sizes_hg20.txt'
ai_plot_bam    = CODE_DIR + 'script/ai_check/ai_plot_bam.r'
check_depth    = CODE_DIR + 'script/ai_check/check_depth.r'
dup_check      = CODE_DIR + 'script/ai_check/identify_genetic_dup.r'

sequencerr = CODE_DIR + 'apps/sequencerr_2.09 -mincov=1 -minmapq=15 -qcut=30 -qCutHard=20 -fcut=0.05 '

gnomad_file             = CODE_DIR + 'refs/gnomAD_v3.1.2_hg38_genome/slim_gnomad.genomes.whole.hg38.hg19.concordant.min5.tsv'
gnomad_hg19_marker_meta = CODE_DIR + 'refs/gnomAD_v3.1.2_hg38_genome/slim_gnomad.genomes.hg38.hg19.concordant.min_0.3_PK_hg19_marker_meta'
gnomad_hg38_marker_meta = CODE_DIR + 'refs/gnomAD_v3.1.2_hg38_genome/slim_gnomad.genomes.hg38.hg19.concordant.min_0.3_PK_hg38_marker_meta'

black_list_file         = CODE_DIR + 'refs/genomic_black_list.bed'
blacklist_hg38_sv       = CODE_DIR + 'refs/marker_black_list/blacklist_hg38_sv.txt'
blacklist_hg38_indelmnv = CODE_DIR + 'refs/marker_black_list/blacklist_hg38_indelmnv.txt'
blacklist_hg19_sv       = CODE_DIR + 'refs/marker_black_list/blacklist_hg19_sv.txt'
blacklist_hg19_indelmnv = CODE_DIR + 'refs/marker_black_list/blacklist_hg19_indelmnv.txt'

# def run_sv_annotation(metafile: str, annotation_ofile: str) -> None:
#     """ create pseudo sv annotation output file """
#     curated_marker_list = []
#     ifn = open(metafile, 'r')
#     ofn = open(annotation_ofile, 'w')
#     ofn.write('\t'.join(['Sample', 'Curated.Marker', 'GeneName', 'RefSeq', 'Region', 'Class', 'AAChange']) + '\n')
    
#     for row in ifn:
#         if 'Marker' not in row:
#             curated_marker = row.strip().split('\t')[2]
#             if curated_marker not in curated_marker_list: curated_marker_list.append(curated_marker)
#     ifn.close()
    
#     for curated_marker in curated_marker_list:
#         new_row = '\t'.join(['SJID', curated_marker] + ['-1'] * 5)
#         ofn.write(new_row + '\n')
#     ofn.close()

def prepare_jobs_4_annotation_for_sv(metafile: str, batchSize: int, hg_build: str, sv_annotator: str, queue: str, jobs_ofn: str):
    """ prepare jobs for sv annotation """
    print('preparing jobs for sv annotation')
    curated_marker_dict = {}
    ifh = open(metafile, 'r')
    for line in ifh:
        if 'Curated.Marker' not in line:
            perfect_pid, ori_marker, curated_marker = line.strip().split('\t')[:3]
            if curated_marker not in curated_marker_dict: curated_marker_dict[curated_marker] = 1
    ifh.close()
    big_annotation_file = metafile + '_4svannotation'
    ofh = open(big_annotation_file,'w')
    for curated_marker in curated_marker_dict: ofh.write(curated_marker+'\n')
    ofh.close()
    
    tmp_annotation_bigfile = big_annotation_file + '.tmp'
    os.system(f'shuf {big_annotation_file} > {tmp_annotation_bigfile}')
    os.system(f'mv {tmp_annotation_bigfile} {big_annotation_file}')
    
    split_prefix = big_annotation_file + '_split_'
    split_fn_dict = aa14_split_big_file_to_small_files(big_annotation_file, batchSize, split_prefix)
    os.system(f'rm {big_annotation_file}')
    
    job_list_dict = {}
    ofh = open(jobs_ofn,'w')
    run_root_dir = os.getcwd()
    temp_header = "temp_svanno_header.txt"
    os.system(f'echo "Curated.Marker" > {temp_header}')
    for smallbatch in split_fn_dict:
        subdir_name = 'batch_' + smallbatch
        os.system(f'mkdir {subdir_name}; cat {temp_header} {smallbatch} > ./{subdir_name}/{smallbatch};')
        bsub_op = f'bsub -q {queue} -M10000 -eo eo.{smallbatch}.eo -oo oo.{smallbatch}.oo '
        dir_op = f'cd {run_root_dir}/{subdir_name}'
        anno_op = f'python {sv_annotator} -varFile {smallbatch} -varType SV -varColumns Curated.Marker -genome {hg_build} -outFile {smallbatch}_annotated'
        cmd = f'{bsub_op} "{dir_op}; {anno_op}"'
        ofh.write(cmd + '\n')
        job_list_dict[smallbatch] = subdir_name
    ofh.close()
    os.remove(temp_header)
    return (job_list_dict)
    
def run_sv_annotation(metafile: str, annotation_ofile: str, batchSize: int, hg_build: str, queue: str, projectID: str, initial_sleep_sec: int, wait_sec: int) -> None:
    """ conduct sv annotation by split job """
    print('running annotation for sv')
    sv_jobs_fn = projectID + '_anno_SV.jobs'
    sv_job_dict = prepare_jobs_4_annotation_for_sv(metafile, batchSize, hg_build, sv_annotator, queue, sv_jobs_fn)
    os.system(f'sh {sv_jobs_fn}')
    success = aa08_job_manager_succeeded(sv_job_dict, True, initial_sleep_sec, wait_sec)
    if success:
        split_fn_dict = {}
        for smallbatch in sv_job_dict:
            annotated_sv_small_fn = sv_job_dict[smallbatch] + f'/{smallbatch}_annotated'
            split_fn_dict[annotated_sv_small_fn] = 1               
        aa15_merge_small_files_to_bigfile(split_fn_dict, annotation_ofile, has_header = True)
        for smallbatch in sv_job_dict:
            subdir = sv_job_dict[smallbatch]
            os.system(f'rm {subdir}/*; rmdir {subdir}; rm eo.{smallbatch}.eo oo.{smallbatch}.oo {smallbatch}')
        os.system(f'rm {sv_jobs_fn}')

def read_black_list(file_path: str, hg_build: str) -> dict:
    """ read a 4-column bed file and save to a dict """
    result = []
    ifn = open(file_path, 'r')
    for row in ifn:
        chrom, start, end, version = row.strip().split('\t')
        if version == hg_build: 
            result.append([chrom, int(start), int(end)])
    return (result)

def make_svindel_curation_big_file_with_filtering(raw_internal_sv_indel_big_file: str, marker_type: str, namemapFile: str, output_file: str, hg_build:str) -> None:
    """ prepare input for the curation pipeline """
    """ to accomadate the formatting difference between Indel and SV genotyping """
    """ indels with length of ref/alt >40 are dropped (due to Kohei's correction)"""
    """ to retire later """
    print('preparing input file from internal raw calls')
    Intake_Sample_Duplicate, Perfect_Sample_Duplicate, sample_tag_list, intake2perfectname_dict, perfectname2intake_dict = aa09_name_map_parser(namemapFile)

    black_list = read_black_list(black_list_file,hg_build)
    
    """ read the raw indel calls, filter long alleles, filter count 0"""
    marker_type_of_interest = ''
    if marker_type.upper() in ['INDEL', 'MNV']: marker_type_of_interest = ['indel', 'mnv']
    if marker_type.upper() == 'SV':    marker_type_of_interest = ['sv']
    marker_sample_dict = {}
    ifh = open(raw_internal_sv_indel_big_file, 'r')
    for line in ifh:
        rcd = line.strip().split('\t')
        intake_sample, orig_marker = rcd[:2]
        if (aa02_marker_classification(orig_marker) in marker_type_of_interest) and (intake_sample in intake2perfectname_dict):
            bam_path = intake2perfectname_dict[intake_sample]['bam_path']
            toUse = False
            if marker_type.upper() in ['INDEL', 'MNV']:
                intake_sample, orig_marker, nM, nT = rcd[:4]
                chrom, pos, ref, alt = orig_marker.split('.')
                if(not(len(ref)>=40 and len(alt)>=40) and len(ref)<70 and len(alt)<70) and int(nM)>0:
                    toUse = True
                for interval in black_list:
                    chro, start, end = interval
                    if chrom == chro and int(pos)>=start and int(pos)<=end: toUse = False
            if marker_type.upper() == 'SV': 
                intake_sample,orig_marker, nMA, nMB = rcd[:4]
                chromA,posA,ortA,chromB,posB,ortB = orig_marker.split('.')[:6]
                if(int(nMA)+int(nMB)>0): toUse = True
                for interval in black_list:
                    chrom, start, end = interval
                    if chromA == chrom and int(posA)>=start and int(posA)<=end: toUse = False
                    if chromB == chrom and int(posB)>=start and int(posB)<=end: toUse = False
            if toUse:
                if intake_sample not in marker_sample_dict: marker_sample_dict[intake_sample] = {}
                marker_sample_dict[intake_sample][orig_marker] = bam_path
    ifh.close()
    
    """ put the uniq sample-marker pairs to big file for curation"""
    ofh = open(output_file, 'w')
    for intake_sample in marker_sample_dict:
        for orig_marker in marker_sample_dict[intake_sample]:
            bam_path = marker_sample_dict[intake_sample][orig_marker]
            if marker_type.upper() in ['INDEL', 'MNV']:
                ofh.write('\t'.join([orig_marker, bam_path]) + '\n')
            if marker_type.upper() == 'SV':
                m = orig_marker.split('.')
                ofh.write('\t'.join([m[0],m[1],m[2],m[3],m[4],m[5],bam_path])+'\n')
    ofh.close()

def filter_low_counts_long_alleles_in_curated_svindel_file(curated_snvindel_marker_big_file: str, hg_build: str):
    """ to filter low counts and long alleles from our curation, or the curated indel is 40nt away from original"""    
    print('filtering low counts and long alleles from curated markers')
    os.system(f'cp {curated_snvindel_marker_big_file} before_filter.{curated_snvindel_marker_big_file}')
    tmp_ofn = curated_snvindel_marker_big_file + '.filtered'

    black_list = read_black_list(black_list_file,hg_build)

    ifh = open(curated_snvindel_marker_big_file,'r')
    ofh = open(tmp_ofn,'w')
    for line in ifh:
        rcd = line.strip().split('\t')
        intake_sample, orig_marker, curated_marker = rcd[:3]
        toUse = False
        if(aa02_marker_classification(curated_marker)=='sv'):
            intake_sample, orig_marker, curated_marker,hg_build,nM,nA,nB = rcd[:7]
            if nB == '': nB = '0'
            if (int(nM)>1 and int(nM)>0.01*(int(nA)+int(nB))): toUse = True
            chromA,posA,ortA,chromB,posB,ortB = curated_marker.split('.')[:6]
            if int(nM)>0: toUse = True
            for interval in black_list:
                chro, start, end = interval
                if chro == chromA and int(posA)>=start and int(posA)<=end: toUse = False
                if chro == chromB and int(posB)>=start and int(posB)<=end: toUse = False
            
        if(aa02_marker_classification(curated_marker) in ['indel','mnv']):
            intake_sample, orig_marker, curated_marker,hg_build,nM,nT = rcd[:6]
            chrom, pos, ref, alt = curated_marker.split('.')
            chrom_ori, pos_ori, ref_ori, alt_ori = orig_marker.split('.')            
            if(not(len(ref)>=40 and len(alt)>=40) and len(ref)<70 and len(alt)<70) and int(nM)>1: toUse = True
            if(len(ref)<2 or len(alt)<2): toUse = False
            if(abs(int(pos)-int(pos_ori))>=40): toUse = False
            for interval in black_list:
                chro, start, end = interval
                if chro == chrom and int(pos)>=start and int(pos)<=end: toUse = False
        if toUse:
            ofh.write(line.strip()+'\n')
    ifh.close()
    ofh.close()
    os.system(f'mv {tmp_ofn} {curated_snvindel_marker_big_file}')

def prepare_jobs_4_svindel_curation(raw_internal_sv_indel_big_file: str, marker_type: str, curator: str, batchSize: int, projectID: str, namemapFile: str, queue: str, jobs_ofn: str, hg_build:str): 
    """ to perform SV Indel curation, separately; """ 
    """ now has folder structure that will retire """
    """ 1. make curation input, filter long allele and low counts """
    print('preparing jobs for split sv/indel curation')
    curation_input_big_file = projectID + '_' + marker_type + '_4curate'
    make_svindel_curation_big_file_with_filtering(raw_internal_sv_indel_big_file, marker_type, namemapFile, curation_input_big_file, hg_build)
    
    """ 2. split to smaller batches"""
    split_prefix = projectID + '_' + marker_type + '_4curate_split_'
    split_fn_dict = aa14_split_big_file_to_small_files(curation_input_big_file, batchSize, split_prefix)
    #os.system(f'rm {curation_input_big_file}')
    
    """ 3. make subdirs, copy files, prepare jobs"""
    job_list_dict = {}
    ofh = open(jobs_ofn,'w')
    run_root_dir = os.getcwd()
    for smallbatch in split_fn_dict:
        subdir_name = 'batch_' + smallbatch
        os.system(f'mkdir {subdir_name}; cp {smallbatch} ./{subdir_name}/')
        bsub_op = f'bsub -q {queue} -M12000 -eo eo.{smallbatch}.eo -oo oo.{smallbatch}.oo '
        dir_op = f'cd {run_root_dir}/{subdir_name}'
        curation_op = f'python {curator} {smallbatch}'
        cmd = f'{bsub_op} "{dir_op}; {curation_op}"'
        ofh.write(cmd + '\n')
        job_list_dict[smallbatch] = subdir_name
    ofh.close()
    return (job_list_dict)

def convert_curation_file_to_meta_file(curated_svindel_bigfile: str, meta_svindel_big_ofile: str, namemapFile: str): 
    """curation file: intake_sample|ori|curate to perfectID|ori|curated"""
    print('coverting curated data to 3-column meta file; diff orig may result in same curated marker')
    Intake_Sample_Duplicate, Perfect_Sample_Duplicate, sample_tag_list, intake2perfectname_dict, perfectname2intake_dict = aa09_name_map_parser(namemapFile)
    meta_dict = {}
    ifh = open(curated_svindel_bigfile,'r')
    for line in ifh:
        intake_sample, ori_marker, curated_marker = line.strip().split('\t')[:3]
        if intake_sample in intake2perfectname_dict:
            perfect_pid = intake2perfectname_dict[intake_sample]['perfect_pid']
            if perfect_pid not in meta_dict: meta_dict[perfect_pid] = {}
            if curated_marker not in meta_dict[perfect_pid]: meta_dict[perfect_pid][curated_marker] = {}
            if ori_marker not in meta_dict[perfect_pid][curated_marker]: meta_dict[perfect_pid][curated_marker][ori_marker] = 1
    ifh.close()
    ofh = open(meta_svindel_big_ofile,'w')
    ofh.write('PID\tOriginal.Marker\tCurated.Marker\n')
    for perfect_pid in meta_dict:
        for curated_marker in meta_dict[perfect_pid]:
            for ori_marker in meta_dict[perfect_pid][curated_marker]:
                ofh.write('\t'.join([perfect_pid,ori_marker,curated_marker])+'\n')
    ofh.close()

def run_sv_curation(raw_internal_sv_bigfile: str, sv_curator: str, batchSize: int, projectID: str, namemapFile: str, sv_meta_ofile: str, queue: str, initial_sleep_sec: int, wait_sec: int, hg_build: str): 
    """ 1. perform SV curation """
    """ 2. filter low counts and long indel alleles"""
    print('running sv/indel curation')
    curated_sv_bigfile    = sv_meta_ofile    + '_big_curatedfile'
    sv_jobs_fn = projectID + '.sv_curation.jobs'
    
    sv_job_dict= prepare_jobs_4_svindel_curation(raw_internal_sv_bigfile, 'sv', sv_curator,batchSize, projectID, namemapFile, queue, sv_jobs_fn, hg_build)
    os.system(f'sh {sv_jobs_fn}')
    print('checking sv curation jobs')
    success = aa08_job_manager_succeeded(sv_job_dict, True, initial_sleep_sec, wait_sec)
    if success:
        split_fn_dict = {}
        for smallbatch in sv_job_dict:
            curated_sv_small_fn = sv_job_dict[smallbatch] + '/' + smallbatch + '_curated_SV.txt'
            split_fn_dict[curated_sv_small_fn] = 1
        aa15_merge_small_files_to_bigfile(split_fn_dict, curated_sv_bigfile)
        for smallbatch in sv_job_dict:
            subdir = sv_job_dict[smallbatch]
            os.system(f'rm {subdir}/*; rmdir {subdir}; rm eo.{smallbatch}.eo oo.{smallbatch}.oo {smallbatch}')
        os.system(f'rm {sv_jobs_fn}')
        filter_low_counts_long_alleles_in_curated_svindel_file(curated_sv_bigfile, hg_build)
        convert_curation_file_to_meta_file(curated_sv_bigfile, sv_meta_ofile, namemapFile)
        remove_blacklist_from_meta_file(sv_meta_ofile, 'sv', hg_build)  

def run_indel_curation(raw_internal_indel_bigfile: str, indel_curator: str, batchSize: int, projectID: str, namemapFile: str, indel_meta_ofile: str, queue: str,initial_sleep_sec: int, wait_sec: int, hg_build:str): 
    """ 1. perform Indel curation """
    """ 2. filter low counts and long indel alleles"""
    print('running sv/indel curation')
    curated_indel_bigfile = indel_meta_ofile + '_big_curatedfile'
    indel_jobs_fn = projectID + '.indel_curation.jobs'
    
    indel_job_dict = prepare_jobs_4_svindel_curation(raw_internal_indel_bigfile, 'indel', indel_curator,batchSize, projectID, namemapFile, queue, indel_jobs_fn,hg_build)
    os.system(f'sh {indel_jobs_fn}')
    print('checking indel curation jobs')
    success = aa08_job_manager_succeeded(indel_job_dict, True, initial_sleep_sec, wait_sec)
    if success:
        split_fn_dict = {}
        for smallbatch in indel_job_dict:
            curated_indel_small_fn = indel_job_dict[smallbatch] + '/' + smallbatch + '_curated_indel'
            split_fn_dict[curated_indel_small_fn] = 1               
        aa15_merge_small_files_to_bigfile(split_fn_dict, curated_indel_bigfile)
        for smallbatch in indel_job_dict:
            subdir = indel_job_dict[smallbatch]
            os.system(f'rm {subdir}/*; rmdir {subdir}; rm eo.{smallbatch}.eo oo.{smallbatch}.oo {smallbatch}')
        os.system(f'rm {indel_jobs_fn}')
        filter_low_counts_long_alleles_in_curated_svindel_file(curated_indel_bigfile, hg_build)
        convert_curation_file_to_meta_file(curated_indel_bigfile, indel_meta_ofile, namemapFile)
        remove_blacklist_from_meta_file(indel_meta_ofile, 'indel', hg_build)           

def remove_blacklist_from_meta_file(meta_file: str, marker_type: str, hg_build: str) -> None:
    """ remove black list markers from meta file """
    tmp_meta = f"{meta_file}_rmblacklist"
    black_list_file = ''
    
    if marker_type.upper() == 'SV':
        if hg_build == 'hg38': black_list_file = blacklist_hg38_sv
        elif hg_build == 'hg19': black_list_file = blacklist_hg19_sv
    elif marker_type.upper() in ['INDEL', 'MNV']:
        if hg_build == 'hg38': black_list_file = blacklist_hg38_indelmnv
        elif hg_build == 'hg19': black_list_file = blacklist_hg19_indelmnv
    
    if black_list_file != '':
        cmd = f'fgrep -vw -f {black_list_file} {meta_file} > {tmp_meta}'
        os.system(cmd)
        os.system(f'mv {tmp_meta} {meta_file}')

def run_sv_indel_curation_together(raw_internal_indel_bigfile: str, raw_internal_sv_bigfile: str, indel_curator: str, sv_curator: str, batchSize: int, projectID: str, namemapFile: str, sv_meta_ofile: str, indel_meta_ofile: str, queue: str, initial_sleep_sec: int, wait_sec: int, hg_build:str): 
    """ 1. perform SV and Indel curation """
    """ 2. filter low counts and long indel alleles"""
    print('running sv/indel curation')
    curated_indel_bigfile = indel_meta_ofile + '_big_curatedfile'
    curated_sv_bigfile    = sv_meta_ofile    + '_big_curatedfile'
    indel_jobs_fn = projectID + '.indel_curation.jobs'
    sv_jobs_fn = projectID + '.sv_curation.jobs'
    indel_job_dict = {}
    sv_job_dict = {}
    
    if(not(raw_internal_indel_bigfile == '')):
        indel_job_dict= prepare_jobs_4_svindel_curation(raw_internal_indel_bigfile, 'indel', indel_curator,batchSize, projectID, namemapFile, queue, indel_jobs_fn,hg_build)
        os.system(f'sh {indel_jobs_fn}')
    if(not(raw_internal_sv_bigfile == '')):
        sv_job_dict= prepare_jobs_4_svindel_curation(raw_internal_sv_bigfile, 'sv', sv_curator,batchSize, projectID, namemapFile, queue, sv_jobs_fn,hg_build)
        os.system(f'sh {sv_jobs_fn}')
    if(not(raw_internal_indel_bigfile == '')):
        print('checking indel curation jobs')
        success = aa08_job_manager_succeeded(indel_job_dict, True, initial_sleep_sec, wait_sec)
        if success:
            split_fn_dict = {}
            for smallbatch in indel_job_dict:
                curated_indel_small_fn = indel_job_dict[smallbatch] + '/' + smallbatch + '_curated_indel'
                split_fn_dict[curated_indel_small_fn] = 1               
            aa15_merge_small_files_to_bigfile(split_fn_dict, curated_indel_bigfile)
            for smallbatch in indel_job_dict:
                subdir = indel_job_dict[smallbatch]
                os.system(f'rm {subdir}/*; rmdir {subdir}; rm eo.{smallbatch}.eo oo.{smallbatch}.oo {smallbatch}')
            os.system(f'rm {indel_jobs_fn}')
            filter_low_counts_long_alleles_in_curated_svindel_file(curated_indel_bigfile,hg_build)
            convert_curation_file_to_meta_file(curated_indel_bigfile, indel_meta_ofile, namemapFile)            
            
    if(not(raw_internal_sv_bigfile == '')):
        print('checking sv curation jobs')
        success = aa08_job_manager_succeeded(sv_job_dict, True, initial_sleep_sec, wait_sec)
        if success:
            split_fn_dict = {}
            for smallbatch in sv_job_dict:
                curated_sv_small_fn = sv_job_dict[smallbatch] + '/' + smallbatch + '_curated_SV.txt'
                split_fn_dict[curated_sv_small_fn] = 1
            aa15_merge_small_files_to_bigfile(split_fn_dict, curated_sv_bigfile)
            for smallbatch in sv_job_dict:
                subdir = sv_job_dict[smallbatch]
                os.system(f'rm {subdir}/*; rmdir {subdir}; rm eo.{smallbatch}.eo oo.{smallbatch}.oo {smallbatch}')
            os.system(f'rm {sv_jobs_fn}')
            filter_low_counts_long_alleles_in_curated_svindel_file(curated_sv_bigfile,hg_build)
            convert_curation_file_to_meta_file(curated_sv_bigfile, sv_meta_ofile, namemapFile)

def prepare_jobs_4_self_genotyping_for_svindel(meta_svindel_bigfile: str, marker_type: str, projectID: str, batchSize: int, hg_build: str, genotyper: str, namemapFile: str, queue: str, jobs_ofn: str, use_perfect_sample: bool = False):
    """ 1. obtain patient/marker unique data """
    print('preparing jobs for sv/indel self-genotyping')
    Intake_Sample_Duplicate, Perfect_Sample_Duplicate, sample_tag_list, intake2perfectname_dict, perfectname2intake_dict = aa09_name_map_parser(namemapFile)
    curated_marker_dict = {}
    ifh = open(meta_svindel_bigfile, 'r')
    for line in ifh:
        perfect_pid, ori_marker, curated_marker = line.strip().split('\t')[:3]
        if perfect_pid not in curated_marker_dict:
            curated_marker_dict[perfect_pid] = {}
        if curated_marker not in curated_marker_dict[perfect_pid]:
            curated_marker_dict[perfect_pid][curated_marker] = 1
    ifh.close()
    big_genotyping_file = meta_svindel_bigfile.split('/')[-1] + '_4selfgeno_bigfile'
    ofh = open(big_genotyping_file,'w')
    for perfect_pid in curated_marker_dict:
        for curated_marker in curated_marker_dict[perfect_pid]:
            if perfect_pid in perfectname2intake_dict:
                for perfect_tag in perfectname2intake_dict[perfect_pid]:
                    bam_path = perfectname2intake_dict[perfect_pid][perfect_tag]['bam_path']
                    intake_sample = perfectname2intake_dict[perfect_pid][perfect_tag]['intake_sample']
                    if use_perfect_sample:
                        ofh.write('\t'.join([bam_path, curated_marker, perfect_pid+'_'+perfect_tag])+'\n')
                    else: ofh.write('\t'.join([bam_path, curated_marker, intake_sample])+'\n')
    ofh.close()

    """ 2. split to smaller batches"""
    split_prefix = projectID + '_' + marker_type + '_4selfgeno_split_'
    split_fn_dict = aa14_split_big_file_to_small_files(big_genotyping_file, batchSize, split_prefix)
    os.system(f'rm {big_genotyping_file}')
    
    """ 3. make subdirs, copy files, prepare jobs"""
    job_list_dict = {}
    ofh = open(jobs_ofn,'w')
    run_root_dir = os.getcwd()
    for smallbatch in split_fn_dict:
        subdir_name = 'batch_' + smallbatch
        os.system(f'mkdir {subdir_name}; cp {smallbatch} ./{subdir_name}/')
        bsub_op = f'bsub -q {queue} -M12000 -eo eo.{smallbatch}.eo -oo oo.{smallbatch}.oo '
        dir_op = f'cd {run_root_dir}/{subdir_name}'
        geno_op = f'python {genotyper} -action one-to-one -build {hg_build} -variantfile {smallbatch}'
        cmd = f'{bsub_op} "{dir_op}; {geno_op}"'
        ofh.write(cmd + '\n')
        job_list_dict[smallbatch] = subdir_name
    ofh.close()
    return (job_list_dict)

def run_self_genotyping_for_svindel_together(meta_sv_bigfile: str, meta_indel_bigfile: str, batchSize: int, hg_build: str, genotyper: str, namemapFile: str, queue: str, projectID: str, sv_genotyped_bigfile: str, indel_genotyped_bigfile: str, initial_sleep_sec: int, wait_sec: int):
    """ perform self genotyping by split job"""
    print('running self genotyping for sv/indel')
    sv_jobs_fn = projectID + '_selfgeno_SV.jobs'
    indel_jobs_fn = projectID + '_selfgeno_Indel.jobs'
    sv_job_dict = {}
    indel_job_dict = {}
    if not(meta_sv_bigfile==''):
        sv_job_dict = prepare_jobs_4_self_genotyping_for_svindel(meta_sv_bigfile, 'sv', projectID, batchSize, hg_build, genotyper, namemapFile, queue,sv_jobs_fn)
        os.system(f'sh {sv_jobs_fn}')
    if not(meta_indel_bigfile==''):
        indel_job_dict = prepare_jobs_4_self_genotyping_for_svindel(meta_indel_bigfile, 'indel', projectID, batchSize, hg_build, genotyper, namemapFile, queue,indel_jobs_fn)
        os.system(f'sh {indel_jobs_fn}')
    if not(meta_sv_bigfile==''):
        success = aa08_job_manager_succeeded(sv_job_dict, True, initial_sleep_sec, wait_sec)
        if success:
            split_fn_dict = {}
            for smallbatch in sv_job_dict:
                genotyped_sv_small_fn = sv_job_dict[smallbatch] + '/counts_' + smallbatch + '.onetoone.txt'
                split_fn_dict[genotyped_sv_small_fn] = 1
            print('merging genotyped SVs')
            aa15_merge_small_files_to_bigfile(split_fn_dict, sv_genotyped_bigfile)
            for smallbatch in sv_job_dict:
                subdir = sv_job_dict[smallbatch]
                os.system(f'rm {subdir}/*; rmdir {subdir}; rm eo.{smallbatch}.eo oo.{smallbatch}.oo {smallbatch}')
            os.system(f'rm {sv_jobs_fn}')
    if not(meta_indel_bigfile==''):
        success = aa08_job_manager_succeeded(indel_job_dict, True, initial_sleep_sec, wait_sec)
        if success:
            split_fn_dict = {}
            for smallbatch in indel_job_dict:
                genotyped_indel_small_fn = indel_job_dict[smallbatch] + '/counts_' + smallbatch + '.onetoone.txt'
                split_fn_dict[genotyped_indel_small_fn] = 1               
            print('merging genotyped Indels')
            aa15_merge_small_files_to_bigfile(split_fn_dict, indel_genotyped_bigfile)
            for smallbatch in indel_job_dict:
                subdir = indel_job_dict[smallbatch]
                os.system(f'rm {subdir}/*; rmdir {subdir}; rm eo.{smallbatch}.eo oo.{smallbatch}.oo {smallbatch}')
            os.system(f'rm {indel_jobs_fn}')

def run_self_genotyping_for_sv(meta_sv_bigfile: str, batchSize: int, hg_build: str, genotyper: str, namemapFile: str, queue: str, projectID: str, sv_genotyped_bigfile: str, initial_sleep_sec: int, wait_sec: int):
    """ perform self genotyping by split job"""
    print('running self genotyping for sv')
    sv_jobs_fn = projectID + '_selfgeno_SV.jobs'
    sv_job_dict = prepare_jobs_4_self_genotyping_for_svindel(meta_sv_bigfile, 'sv', projectID, batchSize, hg_build, genotyper, namemapFile, queue,sv_jobs_fn)
    os.system(f'sh {sv_jobs_fn}')
    success = aa08_job_manager_succeeded(sv_job_dict, True, initial_sleep_sec, wait_sec)
    if success:
        split_fn_dict = {}
        for smallbatch in sv_job_dict:
            genotyped_sv_small_fn = sv_job_dict[smallbatch] + '/counts_' + smallbatch + '.onetoone.txt'
            split_fn_dict[genotyped_sv_small_fn] = 1
        print('merging genotyped SVs')
        aa15_merge_small_files_to_bigfile(split_fn_dict, sv_genotyped_bigfile)
        for smallbatch in sv_job_dict:
            subdir = sv_job_dict[smallbatch]
            os.system(f'rm {subdir}/*; rmdir {subdir}; rm eo.{smallbatch}.eo oo.{smallbatch}.oo {smallbatch}')
        os.system(f'rm {sv_jobs_fn}')

def run_self_genotyping_for_indel(meta_indel_bigfile: str, batchSize: int, hg_build: str, genotyper: str, namemapFile: str, queue: str, projectID: str, indel_genotyped_bigfile: str, initial_sleep_sec: int, wait_sec: int):
    """ perform self genotyping by split job"""
    print('running self genotyping for indel')
    indel_jobs_fn = projectID + '_selfgeno_Indel.jobs'
    indel_job_dict = prepare_jobs_4_self_genotyping_for_svindel(meta_indel_bigfile, 'indel', projectID, batchSize, hg_build, genotyper, namemapFile, queue,indel_jobs_fn)
    os.system(f'sh {indel_jobs_fn}')
    success = aa08_job_manager_succeeded(indel_job_dict, True, initial_sleep_sec, wait_sec)
    if success:
        split_fn_dict = {}
        for smallbatch in indel_job_dict:
            genotyped_indel_small_fn = indel_job_dict[smallbatch] + '/counts_' + smallbatch + '.onetoone.txt'
            split_fn_dict[genotyped_indel_small_fn] = 1               
        print('merging genotyped Indels')
        aa15_merge_small_files_to_bigfile(split_fn_dict, indel_genotyped_bigfile)
        for smallbatch in indel_job_dict:
            subdir = indel_job_dict[smallbatch]
            os.system(f'rm {subdir}/*; rmdir {subdir}; rm eo.{smallbatch}.eo oo.{smallbatch}.oo {smallbatch}')
        os.system(f'rm {indel_jobs_fn}')

def prepare_jobs_4_all2all_genotyping_for_svindel(meta_svindel_bigfile: str, marker_type: str, batchSize: int, hg_build: str, genotyper: str, projectID: str, bamlist_filename: str, queue: str, jobs_ofn: str, big_count_fn: str):
    """ 1. obtain unique marker """
    print('preparing jobs for all2all genotyping for sv/indel')
    curated_marker_dict = {}
    ifh = open(meta_svindel_bigfile, 'r')
    for line in ifh:
        if 'Curated.Marker' not in line:
            perfect_pid, ori_marker, curated_marker = line.strip().split('\t')[:3]
            if curated_marker not in curated_marker_dict: curated_marker_dict[curated_marker] = 1
    ifh.close()
    big_genotyping_file = big_count_fn + '_' + marker_type + '_4all2allgeno'
    ofh = open(big_genotyping_file,'w')
    for curated_marker in curated_marker_dict: ofh.write(curated_marker+'\n')
    ofh.close()
    
    """ 1.5 sort the marker to achieve randomness for Zhikai code"""
    tmp_meta_svindel_bigfile = big_genotyping_file + '.tmp'
    os.system(f'shuf {big_genotyping_file} > {tmp_meta_svindel_bigfile}')
    os.system(f'mv {tmp_meta_svindel_bigfile} {big_genotyping_file}')

    """ 2. split to smaller batches"""
    split_prefix = big_genotyping_file + '_split_'
    split_fn_dict = aa14_split_big_file_to_small_files(big_genotyping_file, batchSize, split_prefix)
    os.system(f'rm {big_genotyping_file}')
    
    """ 3. make subdirs, copy files, prepare jobs"""
    job_list_dict = {}
    ofh = open(jobs_ofn,'w')
    run_root_dir = os.getcwd()
    for smallbatch in split_fn_dict:
        subdir_name = 'batch_' + smallbatch
        os.system(f'mkdir {subdir_name}; cp {smallbatch} ./{subdir_name}/; cp {bamlist_filename} ./{subdir_name}/;')
        bsub_op = f'bsub -q {queue} -M12000 -eo eo.{smallbatch}.eo -oo oo.{smallbatch}.oo '
        dir_op = f'cd {run_root_dir}/{subdir_name}'
        geno_op = f'python {genotyper} -action all-to-all -build {hg_build} -bamfile {bamlist_filename} -variantlst {smallbatch}'
        cmd = f'{bsub_op} "{dir_op}; {geno_op}"'
        ofh.write(cmd + '\n')
        job_list_dict[smallbatch] = subdir_name
    ofh.close()
    return (job_list_dict)

def run_all2all_genotyping_for_sv_indel_together(meta_sv_bigfile: str, meta_indel_bigfile: str, batchSize: int, hg_build: str, genotyper: str, bamlist_filename: str, queue: str, projectID: str, sv_genotyped_bigfile: str, indel_genotyped_bigfile: str, initial_sleep_sec: int, wait_sec: int):
    """ perform self genotyping by split job"""
    print('running all2all genotyping for sv/indel')
    sv_jobs_fn = projectID + '_all2allgeno_SV.jobs'
    indel_jobs_fn = projectID + '_all2allgeno_Indel.jobs'
    sv_job_dict = {}
    indel_job_dict = {}
    if not(meta_sv_bigfile==''):
        sv_job_dict = prepare_jobs_4_all2all_genotyping_for_svindel(meta_sv_bigfile, 'sv', batchSize, hg_build, genotyper, projectID, bamlist_filename, queue,sv_jobs_fn,sv_genotyped_bigfile)
        os.system(f'sh {sv_jobs_fn}')
    if not(meta_indel_bigfile==''):
        indel_job_dict = prepare_jobs_4_all2all_genotyping_for_svindel(meta_indel_bigfile, 'indel', batchSize, hg_build, genotyper, projectID, bamlist_filename, queue,indel_jobs_fn,indel_genotyped_bigfile)
        os.system(f'sh {indel_jobs_fn}')
    if not(meta_sv_bigfile==''):
        success = aa08_job_manager_succeeded(sv_job_dict, True, initial_sleep_sec, wait_sec)
        if success:
            split_fn_dict = {}
            for smallbatch in sv_job_dict:
                genotyped_sv_small_fn = sv_job_dict[smallbatch] + '/counts_' + smallbatch + '_' + bamlist_filename +'.txt'
                split_fn_dict[genotyped_sv_small_fn] = 1               
            aa15_merge_small_files_to_bigfile(split_fn_dict, sv_genotyped_bigfile)
            for smallbatch in sv_job_dict:
                subdir = sv_job_dict[smallbatch]
                os.system(f'rm {subdir}/*; rmdir {subdir}; rm eo.{smallbatch}.eo oo.{smallbatch}.oo {smallbatch}')
            os.system(f'rm {sv_jobs_fn}')
    if not(meta_indel_bigfile==''):
        success = aa08_job_manager_succeeded(indel_job_dict, True, initial_sleep_sec, wait_sec)
        if success:
            split_fn_dict = {}
            for smallbatch in indel_job_dict:
                genotyped_indel_small_fn = indel_job_dict[smallbatch] + '/counts_' + smallbatch + '_' + bamlist_filename +'.txt'
                split_fn_dict[genotyped_indel_small_fn] = 1               
            aa15_merge_small_files_to_bigfile(split_fn_dict, indel_genotyped_bigfile)
            for smallbatch in indel_job_dict:
                subdir = indel_job_dict[smallbatch]
                os.system(f'rm {subdir}/*; rmdir {subdir}; rm eo.{smallbatch}.eo oo.{smallbatch}.oo {smallbatch}')
            os.system(f'rm {indel_jobs_fn}')

def run_all2all_genotyping_for_sv(meta_sv_bigfile: str, batchSize: int, hg_build: str, genotyper: str, bamlist_filename: str, queue: str, projectID: str, sv_genotyped_bigfile: str, initial_sleep_sec: int, wait_sec: int):
    """ perform self genotyping by split job"""
    print('running all2all genotyping for sv/indel')
    sv_jobs_fn = projectID + '_all2allgeno_SV.jobs'
    sv_job_dict = prepare_jobs_4_all2all_genotyping_for_svindel(meta_sv_bigfile, 'sv', batchSize, hg_build, genotyper, projectID, bamlist_filename, queue,sv_jobs_fn,sv_genotyped_bigfile)
    os.system(f'sh {sv_jobs_fn}')
    success = aa08_job_manager_succeeded(sv_job_dict, True, initial_sleep_sec, wait_sec)
    if success:
        split_fn_dict = {}
        for smallbatch in sv_job_dict:
            genotyped_sv_small_fn = sv_job_dict[smallbatch] + '/counts_' + smallbatch + '_' + bamlist_filename +'.txt'
            split_fn_dict[genotyped_sv_small_fn] = 1               
        aa15_merge_small_files_to_bigfile(split_fn_dict, sv_genotyped_bigfile)
        for smallbatch in sv_job_dict:
            subdir = sv_job_dict[smallbatch]
            os.system(f'rm {subdir}/*; rmdir {subdir}; rm eo.{smallbatch}.eo oo.{smallbatch}.oo {smallbatch}')
        os.system(f'rm {sv_jobs_fn}')

def run_all2all_genotyping_for_indel(meta_indel_bigfile: str, batchSize: int, hg_build: str, genotyper: str, bamlist_filename: str, queue: str, projectID: str, indel_genotyped_bigfile: str, initial_sleep_sec: int, wait_sec: int):
    """ perform self genotyping by split job"""
    print('running all2all genotyping for sv/indel')
    indel_jobs_fn = projectID + '_all2allgeno_Indel.jobs'
    indel_job_dict = prepare_jobs_4_all2all_genotyping_for_svindel(meta_indel_bigfile, 'indel', batchSize, hg_build, genotyper, projectID, bamlist_filename, queue,indel_jobs_fn,indel_genotyped_bigfile)
    os.system(f'sh {indel_jobs_fn}')
    success = aa08_job_manager_succeeded(indel_job_dict, True, initial_sleep_sec, wait_sec)
    if success:
        split_fn_dict = {}
        for smallbatch in indel_job_dict:
            genotyped_indel_small_fn = indel_job_dict[smallbatch] + '/counts_' + smallbatch + '_' + bamlist_filename +'.txt'
            split_fn_dict[genotyped_indel_small_fn] = 1               
        aa15_merge_small_files_to_bigfile(split_fn_dict, indel_genotyped_bigfile)
        for smallbatch in indel_job_dict:
            subdir = indel_job_dict[smallbatch]
            os.system(f'rm {subdir}/*; rmdir {subdir}; rm eo.{smallbatch}.eo oo.{smallbatch}.oo {smallbatch}')
        os.system(f'rm {indel_jobs_fn}')

def make_blat_contig_from_marker(snv4_marker: str, flanking_distance: int, input_fasta: str):
    """ build ref/a and mut/b contigs for markers in snv4 format """
    marker_type = aa02_marker_classification(snv4_marker)
    ref_a, ref_b, coord_a, coord_b = ['', '', '', '']
    success_status, contig_dict = aa21_marker_to_contig_omni(snv4_marker, flanking_distance, input_fasta)
    if marker_type.upper() in ['SV']:
        ref_a, ref_b = [contig_dict['left_ref_seq'], contig_dict['right_ref_seq']]
        coord_a = f"{contig_dict['left_chrom']}.{str(contig_dict['left_start'])}.{str(contig_dict['left_end'])}.a"
        coord_b = f"{contig_dict['right_chrom']}.{str(contig_dict['right_start'])}.{str(contig_dict['right_end'])}.b"
    elif marker_type.upper() in ['SNV','INDEL','MNV']:
        ref_a = contig_dict['left_ref_seq'] + contig_dict['ref_allele'] + contig_dict['right_ref_seq'] # ref
        ref_b = contig_dict['left_ref_seq'] + contig_dict['mut_allele'] + contig_dict['right_ref_seq'] # mut
        coord_a = f"{contig_dict['left_chrom']}.{str(contig_dict['left_start'])}.{str(contig_dict['right_end'])}.a" # ref 
        coord_b = f"{contig_dict['right_chrom']}.{str(contig_dict['left_start'])}.{str(contig_dict['right_end'])}.b" # mut
    return (marker_type, ref_a, ref_b, coord_a, coord_b)

def run_blat(input_fasta: str, bit_fa_fn: str, output_file: str) -> None:
    """ run blat for fasta file """
    blat_command = f"{BLAT} {bit_fa_fn} {input_fasta} -out=blast8 {output_file}"
    os.system(blat_command)

def parse_blat_blast8(blat_result: str, output_file: str) -> None:
    """ identify top-2 hits from blat (blast8 format) result """
    blat_dict = {}
    
    ifn = open(blat_result, 'r')
    for row in ifn:
        info = row.strip().split('\t')
        refname, chrom, identity, alignedlen = info[:4]
        if '_' not in chrom: # remove alt/random chr scaffolds 
            query_start, query_end, genome_start, genome_end, e_val, score = info[6:12]
            identity, score = [float(identity), float(score)]
            alignedlen, genome_start, genome_end = [int(alignedlen), int(genome_start), int(genome_end)]
            marker, ref_info = refname.split('_')[:2]
            marker_type = aa02_marker_classification(marker)
            ref_chrom, ref_start, ref_end, ref_type = ref_info.split('.')
            ref_start, ref_end = [int(ref_start), int(ref_end)]
            
            if marker not in blat_dict: 
                blat_dict[marker] = {'a': {'self': {}, 'others': {}}, 'b': {'self': {}, 'others': {}}}
            classification = ''
            if marker_type.upper() in ['SV']:
                if chrom == ref_chrom and ref_start == genome_start and ref_end == genome_end: classification = 'self'
                else: classification = 'others'
            elif marker_type.upper() in ['SNV','INDEL','MNV']:
                hit_center_pos = 0.5*(int(genome_start)+int(genome_end))
                ref_hit_center_pos = 0.5*(int(ref_start)+int(ref_end))
                if chrom == ref_chrom and abs(hit_center_pos-ref_hit_center_pos) < 3: classification = 'self'
                elif chrom == ref_chrom and (genome_start == ref_start or genome_end == ref_end): classification = 'self'
                else: classification = 'others'
            if classification != '':
                if score not in blat_dict[marker][ref_type][classification]: blat_dict[marker][ref_type][classification][score] = []
                blat_dict[marker][ref_type][classification][score].append(identity)
    ifn.close()
    
    ofn = open(output_file, 'w')
    ofn.write('\t'.join(['Curated.Marker', 'Ref/A_score1', 'Ref/A_idty1', 'Ref/A_score2', 'Ref/A_idty2', 
                         'Mut/B_score1', 'Mut/B_idty1', 'Mut/B_score2', 'Mut/B_idty2']) + '\n')
    for marker in blat_dict:
        a_results = ['-1', '-1', '-1', '-1']
        a_scores = list(blat_dict[marker]['a']['others'].keys())
        a_scores.sort(reverse=True)
        idx = 0
        for score in a_scores:
            for identity in blat_dict[marker]['a']['others'][score]:
                if idx < 2:
                    a_results[idx*2] = str(score)
                    a_results[idx*2+1] = str(identity)
                    idx += 1
                else: break
            if idx >= 2: break
        b_results = ['-1', '-1', '-1', '-1']
        b_scores = list(blat_dict[marker]['b']['others'].keys())
        b_scores.sort(reverse=True)
        idx = 0
        for score in b_scores:
            for identity in blat_dict[marker]['b']['others'][score]:
                if idx < 2:
                    b_results[idx*2] = str(score)
                    b_results[idx*2+1] = str(identity)
                    idx += 1
                else: break
            if idx >= 2: break
        ofn.write('\t'.join([marker] + a_results + b_results) + '\n')
    ofn.close()

def run_blat_for_meta(meta_file: str, flanking_distance: int, hg_build: str, output_file: str) -> None:
    """ run blat for markers in meta file and identify top-2 non-self hits in genome """
    blat_fasta = f'{meta_file}_blat.fasta'
    blat_output = f'{meta_file}_blat.fasta.blast8'
    
    hg_fasta = bwaindex_hg38 if hg_build == 'hg38' else bwaindex_hg19
    bit_fa_fn = bit_fa_hg38 if hg_build == 'hg38' else bit_fa_hg19
    
    ifn = open(meta_file, 'r')
    curatedmarker_list = [row.strip().split('\t')[2] for row in ifn if 'Marker' not in row]
    ifn.close()
    curatedmarker_list = list(set(curatedmarker_list))
    
    ofn = open(blat_fasta, 'w')
    for marker in curatedmarker_list:
        marker_type, ref_a, ref_b, coord_a, coord_b = make_blat_contig_from_marker(marker, flanking_distance, hg_fasta)
        ofn.write(f'>{marker}_{coord_a}' + '\n' + ref_a + '\n')
        ofn.write(f'>{marker}_{coord_b}' + '\n' + ref_b + '\n')
    ofn.close()
    
    run_blat(blat_fasta, bit_fa_fn, blat_output)
    
    parse_blat_blast8(blat_output, output_file)

    fn_to_delete=[blat_fasta, blat_output]
    for fn in fn_to_delete:
        cmd = f"rm {fn}"
        os.system(cmd)

def run_sv_indel_omni_together(raw_sv_filename: str, raw_indel_filename: str, namemapFile: str, control_bam_list: str, batchSize: int, all2all_batchSize: int, projectID: str, hg_build: str, indel_curator: str, sv_curator: str, svindel_typer: str, queue: str, only_check_finish: bool, initial_sleep_sec: int, wait_sec: int):
    """ omni procedure for sv/indel"""
    sv_meta_ofile = projectID    + '_SV_meta'
    indel_meta_ofile = projectID + '_Indel_meta'
    sv_selfgeno_ofile = projectID    + '_SV.self_geno.bigcount'
    indel_selfgeno_ofile = projectID + '_Indel.self_geno.bigcount'
    control_sv_count_bigfile = projectID    + '_SV.control_geno.bigcount'
    control_indel_count_bigfile = projectID + '_Indel.control_geno.bigcount'
    incohort_sv_count_bigfile = projectID    + '_SV.incohort_geno.bigcount'
    incohort_indel_count_bigfile = projectID + '_Indel.incohort_geno.bigcount'
    
    run_sv_indel_curation_together(raw_indel_filename, raw_sv_filename, indel_curator, sv_curator, batchSize, projectID, 
        namemapFile,sv_meta_ofile, indel_meta_ofile, queue,initial_sleep_sec, wait_sec,hg_build)

    run_self_genotyping_for_svindel_together(sv_meta_ofile, indel_meta_ofile, batchSize, hg_build, svindel_typer, 
        namemapFile, queue, projectID, sv_selfgeno_ofile,indel_selfgeno_ofile, initial_sleep_sec, wait_sec)

    run_all2all_genotyping_for_sv_indel_together(sv_meta_ofile, indel_meta_ofile, all2all_batchSize, hg_build, 
        svindel_typer, control_bam_list, queue, projectID,control_sv_count_bigfile,control_indel_count_bigfile, initial_sleep_sec, wait_sec)
    
    Intake_Sample_Duplicate, Perfect_Sample_Duplicate, sample_tag_list, intake2perfectname_dict, perfectname2intake_dict = aa09_name_map_parser(namemapFile)
    tmp_incohort_bamlist = projectID + '_incohort_bam4svindel.lst'
    ofh = open(tmp_incohort_bamlist,'w')
    for intake_sample in intake2perfectname_dict:
        bam_path = intake2perfectname_dict[intake_sample]['bam_path']
        ofh.write(bam_path+'\n')
    ofh.close()
    
    run_all2all_genotyping_for_sv_indel_together(sv_meta_ofile, indel_meta_ofile, all2all_batchSize, hg_build, 
        svindel_typer, tmp_incohort_bamlist, queue, projectID,incohort_sv_count_bigfile,incohort_indel_count_bigfile, initial_sleep_sec, wait_sec)
    os.system(f'rm {tmp_incohort_bamlist}')

def run_sv_omni(raw_sv_filename: str, namemapFile: str, platform: str, control_bam_list: str, batchSize: int, all2all_batchSize: int, projectID: str, hg_build: str, sv_curator: str, svindel_typer: str, if_incohort: str, queue: str, initial_sleep_sec: int, wait_sec: int):
    """ omni procedure for sv"""
    sv_meta_ofile = projectID    + '_SV_meta'
    sv_selfgeno_ofile = projectID    + '_SV.self_geno.bigcount'
    control_sv_count_bigfile = projectID    + '_SV.control_geno.bigcount'
    incohort_sv_count_bigfile = projectID    + '_SV.incohort_geno.bigcount'
    annotation_ofile = projectID    + '_SV_meta.annotated'
    merged_outputFile = projectID + platform + '.SV_for_review.txt'
    
    run_sv_curation(raw_sv_filename, sv_curator, batchSize, projectID, namemapFile,sv_meta_ofile, queue, initial_sleep_sec, wait_sec, hg_build)
    run_sv_annotation(sv_meta_ofile, annotation_ofile, batchSize, hg_build, queue, projectID, initial_sleep_sec, wait_sec)
    run_self_genotyping_for_sv(sv_meta_ofile, batchSize, hg_build, svindel_typer, namemapFile, queue, projectID, sv_selfgeno_ofile, initial_sleep_sec, wait_sec)
    run_all2all_genotyping_for_sv(sv_meta_ofile, all2all_batchSize, hg_build, svindel_typer, control_bam_list, queue, projectID, control_sv_count_bigfile, initial_sleep_sec, wait_sec)

    aa25_make_horizontal_file_from_big_count_files(sv_meta_ofile, 'sv', annotation_ofile, sv_selfgeno_ofile, control_sv_count_bigfile, namemapFile, merged_outputFile, platform, hg_build)

    if if_incohort.upper() in ['TRUE', 'YES', 'INCOHORT', 'IN_COHORT']:
        Intake_Sample_Duplicate, Perfect_Sample_Duplicate, sample_tag_list, intake2perfectname_dict, perfectname2intake_dict = aa09_name_map_parser(namemapFile)
        tmp_incohort_bamlist = projectID + '_incohort_bam4_sv.lst'
        ofh = open(tmp_incohort_bamlist,'w')
        for intake_sample in intake2perfectname_dict:
            bam_path = intake2perfectname_dict[intake_sample]['bam_path']
            ofh.write(bam_path+'\n')
        ofh.close()

        run_all2all_genotyping_for_sv(sv_meta_ofile, all2all_batchSize, hg_build, svindel_typer, tmp_incohort_bamlist, queue, projectID, incohort_sv_count_bigfile, initial_sleep_sec, wait_sec)
        os.system(f'rm {tmp_incohort_bamlist}')

def run_indel_omni(raw_indel_filename: str, namemapFile: str, platform: str, control_bam_list: str, batchSize: int, all2all_batchSize: int, projectID: str, hg_build: str, indel_curator: str, svindel_typer: str, if_incohort: str, queue: str, initial_sleep_sec: int, wait_sec: int):
    """ omni procedure for indel"""
    indel_meta_ofile = projectID    + '_Indel_meta'
    indel_selfgeno_ofile = projectID    + '_Indel.self_geno.bigcount'
    control_indel_count_bigfile = projectID    + '_Indel.control_geno.bigcount'
    incohort_indel_count_bigfile = projectID   + '_Indel.incohort_geno.bigcount'
    annotation_ofile = projectID    + '_Indel_meta.annotated'
    merged_outputFile = projectID + platform + '.Indel_for_review.txt'
    
    run_indel_curation(raw_indel_filename, indel_curator, batchSize, projectID, namemapFile, indel_meta_ofile, queue, initial_sleep_sec, wait_sec, hg_build)
    aa26_run_snvindelmnv_annotation(indel_meta_ofile, 3-1, hg_build, annotation_ofile, projectID)
    run_self_genotyping_for_indel(indel_meta_ofile, batchSize, hg_build, svindel_typer, namemapFile, queue, projectID, indel_selfgeno_ofile, initial_sleep_sec, wait_sec)
    run_all2all_genotyping_for_indel(indel_meta_ofile, all2all_batchSize, hg_build, svindel_typer, control_bam_list, queue, projectID, control_indel_count_bigfile, initial_sleep_sec, wait_sec)

    aa25_make_horizontal_file_from_big_count_files(indel_meta_ofile, 'indel', annotation_ofile, indel_selfgeno_ofile, control_indel_count_bigfile, namemapFile, merged_outputFile, platform, hg_build)

    if if_incohort.upper() in ['TRUE', 'YES', 'INCOHORT', 'IN_COHORT']:
        Intake_Sample_Duplicate, Perfect_Sample_Duplicate, sample_tag_list, intake2perfectname_dict, perfectname2intake_dict = aa09_name_map_parser(namemapFile)
        tmp_incohort_bamlist = projectID + '_incohort_bam4_indel.lst'
        ofh = open(tmp_incohort_bamlist,'w')
        for intake_sample in intake2perfectname_dict:
            bam_path = intake2perfectname_dict[intake_sample]['bam_path']
            ofh.write(bam_path+'\n')
        ofh.close()
        
        run_all2all_genotyping_for_indel(indel_meta_ofile, all2all_batchSize, hg_build, svindel_typer, tmp_incohort_bamlist, queue, projectID, incohort_indel_count_bigfile, initial_sleep_sec, wait_sec)
        os.system(f'rm {tmp_incohort_bamlist}')

def run_sv_or_indel_omni(raw_filename: str, marker_type: str, namemapFile: str, platform: str, control_bam_list: str, batchSize: int, all2all_batchSize: int, projectID: str, hg_build: str, if_incohort: str, queue: str, initial_sleep_sec: int, wait_sec: int):
    if(marker_type.upper() == 'SV'): 
        run_sv_omni(raw_filename, namemapFile, platform, control_bam_list, batchSize, all2all_batchSize, projectID, hg_build, sv_curator, svindel_typer, if_incohort, queue, initial_sleep_sec, wait_sec)
    elif(marker_type.upper() in ['INDEL', 'MNV']):
        run_indel_omni(raw_filename, namemapFile, platform, control_bam_list, batchSize, all2all_batchSize,projectID, hg_build, indel_curator, svindel_typer, if_incohort, queue, initial_sleep_sec, wait_sec)

def run_sv_or_indel_curation(raw_internal_bigfile: str, marker_type: str, batchSize: int, projectID: str, namemapFile: str, curated_meta_ofile: str, queue: str,initial_sleep_sec: int, wait_sec: int, hg_build: str):
    message = marker_type + ' curation success!'
    if(aa22_is_raw_call_file_valid(raw_internal_bigfile)[0]):
        if(marker_type.upper() in ['INDEL', 'MNV']):
            run_indel_curation(raw_internal_bigfile, indel_curator, batchSize, projectID, namemapFile, curated_meta_ofile, queue, initial_sleep_sec, wait_sec, hg_build) 
        elif(marker_type.upper() == 'SV'):
            run_sv_curation(raw_internal_bigfile, sv_curator, batchSize, projectID, namemapFile, curated_meta_ofile, queue, initial_sleep_sec, wait_sec, hg_build) 
        else: message = 'marker_type is should be one of sv|indel, while I got: '+marker_type
    else: message = 'check existance of file '+raw_internal_bigfile
    print(message +'\n')

def run_sv_or_indel_self_genotyping(meta_bigfile: str, marker_type: str, batchSize: int, projectID: str, namemapFile: str, selfgenotyping_ofile: str, queue: str, initial_sleep_sec: int, wait_sec: int, hg_build: str):
    message = marker_type + ' self genotyping success!'
    if(aa23_is_curated_call_file_valid(meta_bigfile)[0]):
        if marker_type.upper() == 'SV':
            run_self_genotyping_for_sv(meta_bigfile, batchSize, hg_build, svindel_typer, namemapFile, queue, projectID, selfgenotyping_ofile, initial_sleep_sec, wait_sec)
        elif(marker_type.upper() in ['INDEL', 'MNV']):
            run_self_genotyping_for_indel(meta_bigfile, batchSize, hg_build, svindel_typer, namemapFile, queue, projectID, selfgenotyping_ofile, initial_sleep_sec, wait_sec)
        else: message = 'marker_type is should be one of sv|indel, while I got: '+marker_type
    else: message = 'check existance of file '+meta_bigfile
    print(message +'\n')

def run_sv_or_indel_control_genotyping(meta_bigfile: str, marker_type: str, batchSize: int, projectID: str, control_bamlist_filename: str, control_genotyped_bigofile: str, queue: str, initial_sleep_sec: int, wait_sec: int, hg_build: str):
    message = marker_type + ' control genotyping success!'
    if(aa23_is_curated_call_file_valid(meta_bigfile)[0]):
        if marker_type.upper() == 'SV':
            run_all2all_genotyping_for_sv(meta_bigfile, batchSize, hg_build, svindel_typer, control_bamlist_filename, queue, projectID, control_genotyped_bigofile, initial_sleep_sec, wait_sec)
        elif(marker_type.upper() in ['INDEL', 'MNV']):
            run_all2all_genotyping_for_indel(meta_bigfile, batchSize, hg_build, svindel_typer, control_bamlist_filename, queue, projectID, control_genotyped_bigofile, initial_sleep_sec, wait_sec)
        else: message = 'marker_type is should be one of sv|indel, while I got: '+marker_type
    else: message = 'check existance of file '+meta_bigfile
    print(message +'\n')

def run_sv_or_indel_cohort_genotyping(meta_bigfile: str, marker_type: str, batchSize: int, projectID: str, namemapFile: str, cohort_genotyped_bigofile: str, queue: str, initial_sleep_sec: int, wait_sec: int, hg_build: str):
    message = marker_type + ' control genotyping success!'
    if(aa23_is_curated_call_file_valid(meta_bigfile)[0]):
        Intake_Sample_Duplicate, Perfect_Sample_Duplicate, sample_tag_list, intake2perfectname_dict, perfectname2intake_dict = aa09_name_map_parser(namemapFile)
        tmp_incohort_bamlist = projectID + '_incohort_bam4_' + marker_type + '.lst'
        ofh = open(tmp_incohort_bamlist,'w')
        for intake_sample in intake2perfectname_dict:
            bam_path = intake2perfectname_dict[intake_sample]['bam_path']
            ofh.write(bam_path+'\n')
        ofh.close()
        
        if marker_type.upper() == 'SV':
            run_all2all_genotyping_for_sv(meta_bigfile, batchSize, hg_build, svindel_typer, tmp_incohort_bamlist, queue, projectID, cohort_genotyped_bigofile, initial_sleep_sec, wait_sec)
        elif(marker_type.upper() in ['INDEL', 'MNV']):
            run_all2all_genotyping_for_indel(meta_bigfile, batchSize, hg_build, svindel_typer, tmp_incohort_bamlist, queue, projectID, cohort_genotyped_bigofile, initial_sleep_sec, wait_sec)
        else: message = 'marker_type is should be one of sv|indel, while I got: '+marker_type
        os.system(f'rm {tmp_incohort_bamlist}')
    else: message = 'check existance of file '+meta_bigfile
    print(message +'\n')

def make_job_sequencerr_one_bam(bam: str, snv_fgrepfile: str, count_ofn: str) -> str:
    """ make sequencerr jobs per each bam, run by chromosomes, extract markers using fgrep to avoid large file """
    cmd = '" '    
    rom_list = [str(i) for i in list(range(22,0,-1))] + ['X', 'Y']
    has_chr = aa05_does_bam_file_have_chr_tag(bam)
    if os.path.exists(count_ofn): os.system(f'rm {count_ofn}')
    for rom in rom_list:
        chrom = rom
        if has_chr:
            chrom = 'chr' + rom
            cmd += f'{sequencerr} -STDOUT={chrom}:1-500000000 {bam} | fgrep -w -f {snv_fgrepfile} >> {count_ofn};    '
        else:
            cmd += f'{sequencerr} -STDOUT={chrom}:1-500000000 -prefix=chr {bam} | fgrep -w -f {snv_fgrepfile} >> {count_ofn};    '
    cmd += ' "'
    return (cmd)

def make_jobs_sequencerr_many_bam(bam_dict: dict, snv_fgrepfile: str, jobs_ofn: str, projectID: str, queue: str) -> list:
    """ make job file for a list of bams, each job will count chromosome-by-chromosome"""
    job_name_dict = {}
    count_file_list = []
    intake_sample_list = []
    ofh = open(jobs_ofn,'w')
    for bam in bam_dict:
        sample_name = bam.split('/')[-1].split('.')[0]
        job_name = projectID + '.snv.' + sample_name
        count_ofn = job_name + '.count'
        cmd = make_job_sequencerr_one_bam(bam, snv_fgrepfile, count_ofn)
        job_cmd = f'bsub -q {queue} -M8000 -eo eo.{job_name}.eo -oo oo.{job_name}.oo {cmd}'
        ofh.write(job_cmd + '\n')
        job_name_dict[job_name] = 1
        count_file_list.append(count_ofn)
        intake_sample_list.append(sample_name)
    ofh.close()
    return ([job_name_dict, count_file_list, intake_sample_list])

def make_snv_fgrep_file_using_snvlistfile(snvfile: str, snv4column: int, snvfgrep_ofile: str)->bool:
    """to convert snv4 (0-based column index) to snvfgrep file for sequencErr with fgrep"""
    n_good_marker = 0
    fgrep_loci = {}
    ifh = open(snvfile, 'r')
    for line in ifh:
        rcd = line.strip().split('\t')[snv4column]
        if aa02_marker_classification(rcd) == 'snv': 
            chrom = rcd.split('.')[0]
            if not chrom.startswith('chr'): chrom = 'chr'+ chrom
            pos = rcd.split('.')[1]
            fgrep_loci[chrom + '\t' + str(pos)] = 1
            n_good_marker += 1
    ifh.close()
    if n_good_marker>0:
        ofh = open(snvfgrep_ofile, 'w')
        for rcd in fgrep_loci:
            ofh.write(rcd +'\n')
        ofh.close()
    return (n_good_marker)

def parse_and_merge_snv_from_sequencerr(snvfile: str, snv4column: int, count_file_list:list, intake_sample_list:list, ofile: str):
    """ to parse and merge snv calls from sequencerr output"""
    success, message, n_good_record, n_bad_record, n_good_file, n_bad_file = [True,'',0,0,0,0]

    snv_dict = {}
    ifh = open(snvfile, 'r')
    for line in ifh:
        snv4 = line.strip().split('\t')[snv4column]
        if(len(snv4.split('.'))==4):
            chrom, posit, ref, alt = snv4.split('.')
            if chrom not in snv_dict: snv_dict[chrom] = {}
            if posit not in snv_dict[chrom]: snv_dict[chrom][posit] = {}
            if alt not in snv_dict[chrom][posit]: snv_dict[chrom][posit][alt] = ref
    ifh.close()
    print ('snv read done')

    n_files = len(count_file_list)
    ofh = open(ofile, 'w')
    for count_file,intake_sample in zip(count_file_list,intake_sample_list):
        if os.path.exists(count_file):
            n_good_file += 1
            print(f'processing {n_good_file}/{n_files}')
            ifh = open(count_file, 'r')
            for line in ifh:
                rcd = line.strip().split('\t')
                if len(rcd)==7:
                    n_good_record += 1
                    chrom, pos, nA, nC, nG, nT, nN = rcd
                    if chrom in snv_dict:
                        if pos in snv_dict[chrom]:
                            for alt in snv_dict[chrom][pos]:
                                mkr = '.'.join([chrom,pos,snv_dict[chrom][pos][alt],alt])
                                nM = '0'
                                if   alt == 'A': nM = nA
                                elif alt == 'C': nM = nC
                                elif alt == 'G': nM = nG
                                elif alt == 'T': nM = nT
                                ofh.write('\t'.join([intake_sample,mkr,str(nM),str(nN)])+'\n')
                else: 
                    n_bad_record += 1
                    message = f'record did not have 7 columns in {count_file}: {rcd}'
            ifh.close()
        else:
            n_bad_file += 1
            message = f'count file {count_file} did not exist'
    ofh.close()
    return([success, message, n_good_record, n_bad_record, n_good_file, n_bad_file])

def bamlist_2_bamdict(bamlist_filename: str) -> dict:
    """ process both namemap and bamlist file """
    ifn = open(bamlist_filename, 'r')
    cont_list = [row.strip() for row in ifn]
    ifn.close()
    bam_dict = {}
    if len(cont_list[0].split('\t')) > 1:
        Intake_Sample_Duplicate, Perfect_Sample_Duplicate, sample_tag_list, intake2perfectname_dict, perfectname2intake_dict = aa09_name_map_parser(bamlist_filename)
        for intake_sample in intake2perfectname_dict:
            bampath = intake2perfectname_dict[intake_sample]['bam_path']
            bam_dict[bampath] = intake_sample
    elif len(cont_list[0].split('\t')) == 1:
        for row in cont_list:
            bampath = row.strip()
            if '.bam' in bampath:
                sample = bampath.split('/')[-1].replace('.bam', '')
                bam_dict[bampath] = sample
    return (bam_dict)

def sort_largefile(largefile: str, sort_cols: list, output_name: str, batch_lines: int = 30000000) -> None:
    """ sort large file by splitting it into smaller sort-able parts and combining them """
    os.system(f'time split -l {batch_lines} {largefile} {largefile}_part_')
    sort_idx = ' '.join([f'-k{col},{col}' for col in sort_cols])
    os.system(f'for file in {largefile}_part_*; do time sort {sort_idx} "$file" > "sorted_$file"; done')
    os.system(f'time sort {sort_idx} -m sorted_{largefile}_part_* > {output_name}')
    os.system(f'rm {largefile}_part_* sorted_{largefile}_part_*')

def run_all2all_genotyping_for_snv(meta_bigfile: str, batchSize: int, hg_build: str, sequencerr: str, control_bamlist_filename: str, queue: str, projectID: str, control_genotyped_bigofile: str, initial_sleep_sec: int, wait_sec: int) -> None:
    """ perform all2all genotyping """
    print('running all2all genotyping for snv')
    snv4column = 2
    all_success = True
    only_check_finish = True
    errors = ''
    tmp_big_snv_output_file = projectID + '.tmp.' + control_genotyped_bigofile
    snvfgrep_ofile = projectID + '.snvfgrep'
    snvfgrep_success = make_snv_fgrep_file_using_snvlistfile(meta_bigfile, snv4column, snvfgrep_ofile)
    if snvfgrep_success:
        jobs_ofn = projectID + '.snv_by_sequencer.jobs'
        bam_dict = bamlist_2_bamdict(control_bamlist_filename)
        job_name_dict, count_files, intake_samples = make_jobs_sequencerr_many_bam(bam_dict, snvfgrep_ofile, jobs_ofn, projectID, queue)
        os.system(f'sh {jobs_ofn}')
        all_jobs_success, success_job_dict = aa08_job_manager_succeeded(job_name_dict, only_check_finish, initial_sleep_sec, wait_sec)
        if all_jobs_success:
            os.system(f'rm {jobs_ofn}')
            parse_success, err_message, n_good_record, n_bad_record, n_good_file, n_bad_file = parse_and_merge_snv_from_sequencerr(meta_bigfile, snv4column, count_files, intake_samples, tmp_big_snv_output_file)
            if parse_success:
                my_post_sort = tmp_big_snv_output_file + '.sorted'
                sort_largefile(tmp_big_snv_output_file, [2], my_post_sort)
                os.system(f'mv {my_post_sort} {tmp_big_snv_output_file}')
                # os.system(f'sort -k2,2 {tmp_big_snv_output_file} > {my_post_sort};  mv {my_post_sort} {tmp_big_snv_output_file}')
                for job in job_name_dict: os.system(f'rm eo.{job}.eo oo.{job}.oo')
                for fn in count_files: os.system(f'rm {fn}')
            else: 
                all_success = False
                errors += err_message
        else: 
            all_success = False
            errors += 'job management failure'
    else: 
        all_success = False
        errors += 'snvfgrep failure'
    os.system(f'mv {tmp_big_snv_output_file} {control_genotyped_bigofile}')
    os.system(f'rm {snvfgrep_ofile}')
    return([all_success, errors])

def run_snv_curation(raw_internal_bigfile: str, meta_snv_big_ofile: str, nameMapFile: str) -> None:
    """ current snv will not do curation, but will correct the sample namemap info """
    Intake_Sample_Duplicate, Perfect_Sample_Duplicate, sample_tag_list, intake2perfectname_dict, perfectname2intake_dict = aa09_name_map_parser(nameMapFile)
    
    meta_dict = {}
    ifh = open(raw_internal_bigfile, 'r')
    for line in ifh:
        intake_sample, ori_marker = line.strip().split('\t')[:2]
        if (aa02_marker_classification(ori_marker) == 'snv') and (intake_sample in intake2perfectname_dict):
            perfect_pid = intake2perfectname_dict[intake_sample]['perfect_pid']
            if perfect_pid not in meta_dict: meta_dict[perfect_pid] = {}
            if ori_marker not in meta_dict[perfect_pid]: meta_dict[perfect_pid][ori_marker] = 1
    ifh.close()
    ofh = open(meta_snv_big_ofile, 'w')
    ofh.write('PID\tOriginal.Marker\tCurated.Marker\n')
    for perfect_pid in meta_dict:
        for ori_marker in meta_dict[perfect_pid]:
                ofh.write('\t'.join([perfect_pid, ori_marker, ori_marker])+'\n')
    ofh.close()     
    
def run_snv_control_genotyping(meta_bigfile: str, batchSize: int, projectID: str, control_bamlist_filename: str, control_genotyped_bigofile: str, queue: str, initial_sleep_sec: int, wait_sec: int, hg_build: str) -> None:
    """ perform all2all snv genotyping """
    message = 'snv control genotyping success!'
    if(aa23_is_curated_call_file_valid(meta_bigfile)[0]):
        run_all2all_genotyping_for_snv(meta_bigfile, batchSize, hg_build, sequencerr, control_bamlist_filename, queue, projectID, control_genotyped_bigofile, initial_sleep_sec, wait_sec)
    else: message = 'check existance of file ' + meta_bigfile
    print(message +'\n')

def run_bam2ai_count(batchSize: int, projectID: str, bamlist_filename: str, genotyped_bigofile: str, queue: str, initial_sleep_sec: int, wait_sec: int, hg_build: str) -> None:
    """ perform bam2ai snv genotyping """
    message = 'bam2ai snv genotyping success!'
    if hg_build == 'hg19': meta_bigfile=gnomad_hg19_marker_meta
    elif hg_build == 'hg38': meta_bigfile=gnomad_hg38_marker_meta
    else: meta_bigfile=''
    if(aa23_is_curated_call_file_valid(meta_bigfile)[0]):
        run_all2all_genotyping_for_snv(meta_bigfile, batchSize, hg_build, sequencerr, bamlist_filename, queue, projectID, genotyped_bigofile, initial_sleep_sec, wait_sec)
    else: message = 'check existance of gnomad marker file ' + meta_bigfile
    print(message +'\n')

def plot_geno2ai(geno_all2all_fn: str, projectID: str, output_dir: str, queue: str, initial_sleep_sec: int, wait_sec: int, hg_build: str) -> list:
    """ parse the all-to-all genotyping result and plot ai check """
    message = 'ai check plotting success!'
    if hg_build == 'hg19': chromsize_file=chromsize_hg19
    elif hg_build == 'hg38': chromsize_file=chromsize_hg38
    else: chromsize_file=''

    if not os.path.exists(output_dir): os.makedirs(output_dir)
    
    if chromsize_file != '':
        aicheck_file_list = make_aicheck_files(geno_all2all_fn, projectID, output_dir)
        manage_aiplot_jobs(aicheck_file_list, chromsize_file, projectID, queue, initial_sleep_sec, wait_sec)
    else: message = 'invalid genome'
    print(message +'\n')
    
    return (aicheck_file_list)

def ai_dupcheck(aichecklist_filename: str, projectID: str, queue: str, initial_sleep_sec: int, wait_sec: int, batchSize: int = 50) -> str:
    """ """
    all_pair_file = f'{projectID}_all_pairs.txt'
    depth_stat_file = f'{projectID}_depth_stats.txt'
    dupcheck_output = f'{projectID}_rlt_dup_check.txt'
    
    cmd = f'Rscript {check_depth} {aichecklist_filename} {projectID}'
    os.system(cmd)
    
    split_prefix =  f'{projectID}_all_pairs_split_'
    split_fn_dict = aa14_split_big_file_to_small_files(all_pair_file, batchSize, split_prefix)
    dupcheck_file_list = list(split_fn_dict.keys())
    manage_dupcheck_jobs(dupcheck_file_list, projectID, dupcheck_output, queue, initial_sleep_sec, wait_sec)
    for dupcheck_file in dupcheck_file_list:
        print (f'removing {dupcheck_file}')
        os.remove(dupcheck_file)    
    return (dupcheck_output)

def manage_dupcheck_jobs(dupcheck_file_list: list, projectID: str, dupcheck_output: str, queue: str, initial_sleep_sec: int, wait_sec: int) -> None:
    """ manage ai check plotting progress """
    print('submit dup check jobs')
    dupcheck_jobs_fn = projectID + '_dupcheck.jobs'
    dupcheck_job_dict = prepare_jobs_4_dupcheck(dupcheck_file_list, queue, dupcheck_jobs_fn)
    os.system(f'sh {dupcheck_jobs_fn}')
    success = aa08_job_manager_succeeded(dupcheck_job_dict, True, initial_sleep_sec, wait_sec)
    if success:
        print('dup check jobs finished')
        os.system(f'rm {dupcheck_jobs_fn}')
        dup_check_file_list = list(dupcheck_job_dict.values())
        os.system(f"head -1 {dup_check_file_list[0]} > header")
        os.system(f"cat {' '.join(dup_check_file_list)} | grep -v SampleA > {dupcheck_output}_tmp")
        os.system(f"cat header {dupcheck_output}_tmp > {dupcheck_output}")
        os.remove(f"{dupcheck_output}_tmp")
        for job in dupcheck_job_dict:
            print (f'removing eo.{job}.eo oo.{job}.oo {dupcheck_job_dict[job]}')
            os.remove(f'eo.{job}.eo')
            os.remove(f'oo.{job}.oo')
            os.remove(dupcheck_job_dict[job])
        os.system(f'rm {dupcheck_jobs_fn}')

def prepare_jobs_4_dupcheck(dupcheck_file_list: list, queue: str, dupcheck_jobs_fn: str):
    job_list_dict = {}
    ofn = open(dupcheck_jobs_fn, 'w')
    for dupcheck_file in dupcheck_file_list:
        dupcheck_output = dupcheck_file+'_dup_check.txt'
        jobname = dupcheck_file.split('/')[-1].split('.')[0]
        ofn.write(f'bsub -q {queue} -M8000 -eo eo.{jobname}.eo -oo oo.{jobname}.oo Rscript {dup_check} {dupcheck_file} {dupcheck_output}' + '\n')
        job_list_dict[jobname] = dupcheck_output
    ofn.close()
    return (job_list_dict)

def run_bam2ai_omni(batchSize: int, projectID: str, bamlist_filename: str, output_dir: str, queue: str, initial_sleep_sec: int, wait_sec: int, hg_build: str) -> None:
    """ omni procedure for bam2ai """
    genotyped_bigofile = f'{projectID}_gnomad_{hg_build}.incohort_geno.bigcount'
    aichecklist_filename = f'{projectID}_ai_check_files.lst'

    run_bam2ai_count(batchSize, projectID, bamlist_filename, genotyped_bigofile, queue, initial_sleep_sec, wait_sec, hg_build)
    aicheck_file_list = plot_geno2ai(genotyped_bigofile, projectID, output_dir, queue, initial_sleep_sec, wait_sec, hg_build)
    
    ofn = open(aichecklist_filename, 'w')
    ofn.write('\n'.join(aicheck_file_list) + '\n')
    ofn.close()
    ai_dupcheck(aichecklist_filename, projectID, queue, initial_sleep_sec, wait_sec)

def prepare_jobs_4_aicheck(aicheck_file_list: list, chromsize_file: str, queue: str, aicheck_jobs_fn: str):
    job_list_dict = {}
    ofn = open(aicheck_jobs_fn, 'w')
    for aicheck_file in aicheck_file_list:
        aicheck_figure = aicheck_file.replace('.txt', '.jpg')
        jobname = aicheck_file.split('/')[-1].split('.')[0]
        ofn.write(f'bsub -q {queue} -M8000 -eo eo.{jobname}.eo -oo oo.{jobname}.oo Rscript {ai_plot_bam} {chromsize_file} {aicheck_file} {aicheck_figure}' + '\n')
        job_list_dict[jobname] = 1
    ofn.close()
    return (job_list_dict)

def manage_aiplot_jobs(aicheck_file_list: list, chromsize_file: str, projectID: str, queue: str, initial_sleep_sec: int, wait_sec: int) -> None:
    """ manage ai check plotting progress """
    print('submit ai check plotting jobs')
    aicheck_jobs_fn = projectID + '_aicheck.jobs'
    aicheck_job_dict = prepare_jobs_4_aicheck(aicheck_file_list, chromsize_file, queue, aicheck_jobs_fn)
    os.system(f'sh {aicheck_jobs_fn}')
    success = aa08_job_manager_succeeded(aicheck_job_dict, True, initial_sleep_sec, wait_sec)
    if success:
        print('ai check plotting jobs finished')
        for job in aicheck_job_dict:
            os.remove(f'eo.{job}.eo')
            os.remove(f'oo.{job}.oo')
        os.system(f'rm {aicheck_jobs_fn}')

def make_aicheck_files(geno_all2all_fn: str, projectID: str, output_dir: str) -> list:
    """ create ai check input files for each sample """
    if not os.path.exists(output_dir): os.makedirs(output_dir)
    aicheck_files = []
    cmd = f"""awk '{{split($2, markerParts, "."); outputFile = "{output_dir}/ai_check_{projectID}_"$1".txt"; if (!headersWritten[outputFile]++) {{print "Chr\\tPos\\tRef\\tAlt\\tnM\\tnT" > outputFile;}}  print markerParts[1]"\\t"markerParts[2]"\\t"markerParts[3]"\\t"markerParts[4]"\\t"$3"\\t"$4 > outputFile;}}' {geno_all2all_fn}"""
    os.system(cmd)
    aicheck_files = glob.glob(f"{output_dir}/ai_check_{projectID}_*.txt")
    return (aicheck_files)

def run_snv_omni(raw_snv_filename: str, namemapFile: str, platform: str, control_bam_list: str, batchSize: int, all2all_batchSize: int, projectID: str, hg_build: str, queue: str, initial_sleep_sec: int, wait_sec: int) -> None:
    """ omni procedure for snv """
    snv_meta_ofile = projectID + '_SNV_meta'
    annotation_ofile = projectID + '_SNV_meta.annotated'
    # snv_selfgeno_ofile = projectID + '_SNV.self_geno.bigcount'
    control_snv_count_bigfile = projectID + '_SNV.control_geno.bigcount'
    incohort_snv_count_bigfile = projectID + '_SNV.incohort_geno.bigcount'
    outputFile = projectID + platform + '.SNV_for_review.txt'
    
    run_snv_curation(raw_snv_filename, snv_meta_ofile, namemapFile)
    aa26_run_snvindelmnv_annotation(snv_meta_ofile, 3-1, hg_build, annotation_ofile, projectID)
    run_snv_control_genotyping(snv_meta_ofile, batchSize, projectID, namemapFile, incohort_snv_count_bigfile, queue, initial_sleep_sec, wait_sec, hg_build)
    run_snv_control_genotyping(snv_meta_ofile, batchSize, projectID, control_bam_list, control_snv_count_bigfile, queue, initial_sleep_sec, wait_sec, hg_build)
    aa25_make_horizontal_file_from_big_count_files(snv_meta_ofile, 'snv', annotation_ofile, incohort_snv_count_bigfile, control_snv_count_bigfile, namemapFile, outputFile, platform, hg_build)

def meta_union(table1: str, table2: str, idx_cols: list, wanted_col: int, delim_for_indexing: str, output_table: str) -> list:
    """ take the union of 2 meta-tables with the same column numbers and orders """
    success = True
    message = ''
    table_dict = {}
    
    ifh = open(table1, 'r')
    for line in ifh:
        rec = line.strip().split('\t')
        if len(rec) >= max(idx_cols):
            table1_idx_key = delim_for_indexing.join([rec[idx] for idx in idx_cols])
            if table1_idx_key not in table_dict: table_dict[table1_idx_key] = {}
            wanted_info = rec[wanted_col]
            if wanted_info not in table_dict[table1_idx_key]: table_dict[table1_idx_key][wanted_info] = ''
        else:
            success = False
            message = 'table1 has fewer columns than needed'
            break
    ifh.close()
    
    ifh = open(table2, 'r')
    for line in ifh:
        rec = line.strip().split('\t')
        if len(rec) >= max(idx_cols):
            table2_idx_key = delim_for_indexing.join([rec[idx] for idx in idx_cols])
            if table2_idx_key not in table_dict: table_dict[table2_idx_key] = {}
            wanted_info = rec[wanted_col]
            if wanted_info not in table_dict[table2_idx_key]: table_dict[table2_idx_key][wanted_info] = ''
        else:
            success = False
            message = 'table2 has fewer columns than needed'
            break
    ifh.close()

    ofn = open(output_table, 'w')
    for idx_key in table_dict:
        wanted_info = ';'.join(list(table_dict[idx_key].keys()))
        final_cols = idx_key.split(delim_for_indexing)
        final_cols.insert(wanted_col, wanted_info)
        ofn.write('\t'.join(final_cols) + '\n')
    ofn.close()

    return ([success, message])

def identify_genotyping_cols(hori_file: str, platform: str) -> list:
    """ ifentify the genotyping result columns in horizontal table """
    geno_col_list = []
    ifn = open(hori_file, 'r')
    header_line = ifn.readline()
    col_list = header_line.split('\t')
    for idx in range(len(col_list)):
        col = col_list[idx]
        if ('nT_' in col or 'nM_' in col) and (f'_{platform}' in col): geno_col_list.append(str(idx))
    ifn.close()
    return (geno_col_list)

def run_crossplatform_genotyping(
    meta_bigfile: str, 
    marker_type: str, 
    batchSize: int, 
    projectID: str, 
    meta_platformfile: str, 
    controlgeno_all2all_fn: str,
    outputFile: str, 
    queue: str, 
    initial_sleep_sec: int, 
    wait_sec: int, 
    hg_build: str
    ) -> None:
    """ conduct self-genotyping, annotation, and combine genotyping result across all platforms """
    annotation_fn = f'{meta_bigfile}_annotated'
    if marker_type.upper() in ['SNV', 'INDEL', 'MNV']:
        snv4column = 2
        aa26_run_snvindelmnv_annotation(meta_bigfile, snv4column, hg_build, annotation_fn, projectID)
    elif marker_type.upper() == 'SV':
        run_sv_annotation(meta_bigfile, annotation_fn, batchSize, hg_build, queue, projectID, initial_sleep_sec, wait_sec)

    hori_file_dict = {}
    ifn = open(meta_platformfile, 'r')
    for row in ifn:
        platform, nameMapFile = row.strip().split('\t')
        platform_projectID = f'{platform}_{projectID}_{marker_type}'
        selfgenotyping_ofile = f'{platform_projectID}.self_geno.bigcount'
        hori_ofile = f'review_{platform_projectID}.txt'
        if marker_type.upper() in ['INDEL', 'SV', 'MNV']:
            run_sv_or_indel_self_genotyping(meta_bigfile, marker_type, batchSize, platform_projectID, nameMapFile, selfgenotyping_ofile, queue, initial_sleep_sec, wait_sec, hg_build)
        elif marker_type.upper() == 'SNV':
            run_snv_control_genotyping(meta_bigfile, batchSize, platform_projectID, nameMapFile, selfgenotyping_ofile, queue, initial_sleep_sec, wait_sec, hg_build)
        if marker_type.upper() in ['INDEL', 'SNV', 'MNV']:
            aa25_make_horizontal_file_from_big_count_files(meta_bigfile, 'snv', annotation_fn, selfgenotyping_ofile, controlgeno_all2all_fn, nameMapFile, hori_ofile, platform, hg_build)
            hori_file_dict[platform] = hori_ofile
        elif marker_type.upper() == 'SV':
            aa25_make_horizontal_file_from_big_count_files(meta_bigfile, 'sv', annotation_fn, selfgenotyping_ofile, controlgeno_all2all_fn, nameMapFile, hori_ofile, platform, hg_build)
            hori_file_dict[platform] = hori_ofile
    ifn.close()
    
    platform_list = list(hori_file_dict.keys())
    platform_0 = platform_list[0]
    to_table = hori_file_dict[platform_0]
    tmp_file_list = []
    for platform in platform_list[1:]:
        platform_hori_file = hori_file_dict[platform]
        idx_cols = [0, 1, 2]
        from_table = platform_hori_file
        wanted_from_table_cols = ':'.join(identify_genotyping_cols(platform_hori_file, platform))
        delim_for_indexing = ';'
        output_file = f'tmp_{platform_hori_file}'
        aa10_table_joining(to_table, idx_cols, from_table, idx_cols, wanted_from_table_cols, delim_for_indexing, output_file)
        to_table = output_file
        tmp_file_list.append(output_file)
    
    df = pd.read_csv(to_table, sep='\t')
    first_three_cols = df.columns[:3]
    geno_cols = [col for col in df.columns if (col.startswith('nM_') or col.startswith('nT_')) and col not in first_three_cols]
    other_cols = [col for col in df.columns if col not in first_three_cols and col not in geno_cols]
    new_column_order = list(first_three_cols) + geno_cols + other_cols
    df = df[new_column_order]
    df.to_csv(outputFile, sep='\t', index=False)
    os.system(f"rm {' '.join(tmp_file_list)}")

def horiz_to_vertical(ifn: str, good_marker_ifn: str, ofn: str, min_nM: int = 3, min_nT: int = 3) -> None:
    """ convert horizontal genotyping table to vertical format; good_marker: PID.marker"""
    gfh = open(good_marker_ifn, 'r')
    good_mkr = {}
    for line in gfh:
        good_mkr[line.strip()] = 1
    gfh.close()
    print('good marker read done')
    
    ifh = open(ifn, 'r')
    Platforms = {}
    thedata = {}
    anno_idx_dict = {}
    marker_anno_dict = {}
    for line in ifh:
        rcd = line.strip().split('\t')
        if 'Marker' in line:
            for i in range(3, len(rcd)):
                if rcd[i].startswith('nM_'):
                    detail = rcd[i].split('_')
                    if detail[2] not in Platforms: Platforms[detail[2]] = {}
                    if detail[1] not in Platforms[detail[2]]: Platforms[detail[2]][detail[1]] = {}
                    Platforms[detail[2]][detail[1]]['nM'] = i
                elif rcd[i].startswith('nT_'):
                    detail = rcd[i].split('_')
                    if detail[2] not in Platforms: Platforms[detail[2]] = {}
                    if detail[1] not in Platforms[detail[2]]: Platforms[detail[2]][detail[1]] = {}
                    Platforms[detail[2]][detail[1]]['nT'] = i
                elif rcd[i] == 'GeneName': anno_idx_dict['GeneName'] = i
                elif rcd[i] == 'RefSeq': anno_idx_dict['RefSeq'] = i
                elif rcd[i] == 'Region': anno_idx_dict['Region'] = i
                elif rcd[i] == 'Class': anno_idx_dict['Class'] = i
                elif rcd[i] == 'AAChange': anno_idx_dict['AAChange'] = i
        else:
            PID, orig, cura = rcd[:3]
            idx = PID + '.' + cura
            if idx in good_mkr:
                if PID not in thedata: thedata[PID] = {}
                if cura not in thedata[PID]: thedata[PID][cura] = {}
                thedata[PID][cura]['Original.Marker'] = orig
                for plat in Platforms:
                    for tag in Platforms[plat]:
                        nM_col = Platforms[plat][tag]['nM']
                        nT_col = Platforms[plat][tag]['nT']
                        if(int(rcd[nT_col])>=min_nT and int(rcd[nM_col])>=min_nM):
                            if plat not in  thedata[PID][cura]: thedata[PID][cura][plat] = {}
                            if tag not in thedata[PID][cura][plat]: thedata[PID][cura][plat][tag] = {}
                            thedata[PID][cura][plat][tag]['nM'] = rcd[nM_col]
                            thedata[PID][cura][plat][tag]['nT'] = rcd[nT_col]
                if cura not in marker_anno_dict: marker_anno_dict[cura] = {}
                marker_anno_dict[cura]['GeneName'] = rcd[anno_idx_dict['GeneName']] if 'GeneName' in anno_idx_dict else '-1'
                marker_anno_dict[cura]['RefSeq'] = rcd[anno_idx_dict['RefSeq']] if 'RefSeq' in anno_idx_dict else '-1'
                marker_anno_dict[cura]['Region'] = rcd[anno_idx_dict['Region']] if 'Region' in anno_idx_dict else '-1'
                marker_anno_dict[cura]['Class'] = rcd[anno_idx_dict['Class']] if 'Class' in anno_idx_dict else '-1'
                marker_anno_dict[cura]['AAChange'] = rcd[anno_idx_dict['AAChange']] if 'AAChange' in anno_idx_dict else '-1'
    ifh.close()
    print('data read done')

    ofh = open(ofn, 'w')
    Ostring = '\t'.join(['Sample','Original.Marker','Curated.Marker','Platform','nM','nT','Chrom','Pos','Ref','Alt','Origin','Valid','AAChange','GeneName','mRNA','Class','Driver'])
    ofh.write(Ostring + '\n')
    for PID in thedata:
        for cura in thedata[PID]:
            for plat in thedata[PID][cura]:
                if (not plat == 'Original.Marker'):
                    for tag in thedata[PID][cura][plat]:
                        if 'H' not in tag and 'G' not in tag: # 03192024
                            Sample = PID + '_' + tag
                            Original_Marker = thedata[PID][cura]['Original.Marker']
                            Marker = cura
                            Chrom, Pos, Ref, Alt = cura.split('.')[:4]
                            nM = thedata[PID][cura][plat][tag]['nM']
                            nT = thedata[PID][cura][plat][tag]['nT']
                            Origin, Valid, Driver  = ['Somatic', 'Putative', 'ToAdd']
                            AAChange = marker_anno_dict[cura]['AAChange']
                            GeneName = marker_anno_dict[cura]['GeneName']
                            mRNA = marker_anno_dict[cura]['RefSeq']
                            Class = marker_anno_dict[cura]['Class']
                            Ostring = '\t'.join([Sample,Original_Marker,Marker,plat,nM,nT,Chrom,Pos,Ref,Alt,Origin,Valid,AAChange,GeneName,mRNA,Class,Driver])
                            ofh.write(Ostring + '\n')
    ofh.close()
