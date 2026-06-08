import re
import os
import sys
import argparse
import pandas as pd


def qc_file_existence(file_path: str) -> bool:
    """ Check if file exists and is empty by reading its size """
    file_is_good = False
    message = ''
    if os.path.isfile(file_path):
        if os.path.getsize(file_path) != 0: file_is_good = True
        else: message = f"size is 0 for file {file_path}"
    else: message = f"non-exist for file {file_path}"
    return ([file_is_good, message])

def read_bed_file(file_path: str) -> dict:
    """ read a 4-column bed file and save to a dict """
    result = {}
    ifn = open(file_path, 'r')
    for row in ifn:
        chrom, start, end, version = row.strip().split('\t')
        if version not in result: result[version] = []
        result[version].append([chrom, int(start), int(end)])
    return (result)

def transform_marker(marker: str) -> str:
    """ replace empty string in marker with '-' """
    if not isinstance(marker, str): marker = str(marker)
    transformed = re.sub(r'\.\.', '.-.', marker)
    if transformed.endswith('.'):
        transformed = transformed + '-'
    return transformed

def filter_internal_indel(file_path: str, allele_len: int)->None:
    tmpofn = file_path + '.filtered'
    ifh = open(file_path, 'r')
    ofh = open(tmpofn,    'w')
    for line in ifh:
        toUse = True
        if ('Marker' not in line.strip()) and ('Sample' not in line.strip()): 
            marker = line.strip().split('\t')[1]
            chrom, pos, ref, alt = marker.split('.')
            if (len(ref)>=allele_len and len(alt)>=allele_len): toUse = False
            nM = int(line.strip().split('\t')[2])
            nT = int(line.strip().split('\t')[3])
            if (nM<1): toUse = False
            if (nM<0.01*nT): toUse = False
        if toUse:
            ofh.write(line.strip() + '\n')
    ifh.close()
    ofh.close()
    os.system(f'mv {tmpofn} {file_path} ')

def filter_internal_snv(file_path: str)->None:
    tmpofn = file_path + '.filtered'
    ifh = open(file_path, 'r')
    ofh = open(tmpofn,    'w')
    for line in ifh:
        toUse = True
        if ('Marker' not in line.strip()) and ('Sample' not in line.strip()): 
            marker = line.strip().split('\t')[1]
            nM = int(line.strip().split('\t')[2])
            nT = int(line.strip().split('\t')[3])
            if(nM<2): toUse = False
            if(nM<0.01*nT): toUse = False
        if toUse:
            ofh.write(line.strip() + '\n')
    ifh.close()
    ofh.close()
    os.system(f'mv {tmpofn} {file_path} ')


def filter_internal_sv(file_path: str, black_list_dict: dict, hg_build: str = 'hg38') -> None:
    """ filter sv markers """
    tmpofn = file_path + '.filtered'
    ifh = open(file_path, 'r')
    ofh = open(tmpofn,    'w')
    interval_list = black_list_dict[hg_build] if hg_build in black_list_dict else []
    
    for line in ifh:
        toUse = True
        if ('Marker' not in line.strip()) and ('Sample' not in line.strip()): 
            intake_sample, marker, nA, nB, nAT, nBT = line.strip().split('\t')[:6]
            if((int(nA)+int(nB))<1): toUse = False
            if((int(nA)+int(nB))<0.01*(int(nAT)+int(nBT))): toUse = False
            chromA,posA,ortA,chromB,posB,ortB = marker.split('.')[:6]
            for interval in interval_list:
                chrom, start, end = interval
                if chromA == chrom and int(posA)>=start and int(posA)<=end: toUse = False
                if chromB == chrom and int(posB)>=start and int(posB)<=end: toUse = False
        if toUse:
            ofh.write(line.strip() + '\n')
    ifh.close()
    ofh.close()
    os.system(f'mv {tmpofn} {file_path} ')

def snv_summary(input_file: str, sample: str, chrom_filter: bool) -> pd.DataFrame:
    """ summarize a single snv mutation file """
    df_out = pd.DataFrame()
    if os.path.exists(input_file) and os.path.getsize(input_file) > 0:
        df = pd.read_csv(input_file, sep = "\t", low_memory=False)
        columns = df.columns.tolist()
        is_tier23 = False
        if 'ReferenceAllele' not in columns: # old tier 2 and 3 format
            is_tier23 = True
            column_names = ['GeneName', 'SJQuality', 'Sample', 'Chr', 'WU_HG19_Pos', '#Mutant_In_Tumor', '#Total_In_Tumor',
                            '#Mutant_In_Normal', '#Total_In_Normal', 'ReferenceAllele', 'MutantAllele', 'Flanking']
            df = pd.read_csv(input_file, sep = "\t", low_memory=False, header = None, names=column_names)
            
        columns = df.columns.tolist()
        is_new_format = False if '#Mutant_In_Tumor' in columns else True
        df_out = pd.DataFrame()
        
        df['Chr'] = df['Chr'].astype(str)
        mask = ~df['Chr'].str.startswith('chr')
        df.loc[mask, 'Chr'] = 'chr' + df.loc[mask, 'Chr']
        if 'WU_HG18_Pos' in columns: df['WU_HG19_Pos'] = df['WU_HG18_Pos']
        if is_new_format:
            if chrom_filter:
                chrom_list = ['chr'+str(i) for i in list(range(22,0,-1))] + ['chrX', 'chrY']
                df = df[df['Chr'].isin(chrom_list)]
            df_out['Marker'] = df['Chr'] + '.' + df['WU_HG19_Pos'].astype(str) + '.' + df['ReferenceAllele'] + '.' + df['MutantAllele']
            df_out['nM_D'] = df['alternative_tumor_count'].fillna(0).astype(int)
            df_out['nT_D'] = df['alternative_tumor_count'].fillna(0).astype(int) + df['reference_tumor_count'].fillna(0).astype(int)
            df_out['nM_G'] = df['alternative_normal_count'].fillna(0).astype(int)
            df_out['nT_G'] = df['alternative_normal_count'].fillna(0).astype(int) + df['reference_normal_count'].fillna(0).astype(int)
            df_out['Annotation'] = df['GeneName'].astype(str) + '|' + df['mRNA_acc'].astype(str) + '|' + df['Class'].astype(str) + '|' + df['AAChange'].astype(str)
            df_out['QualTag'] = df['QUALTAG'].astype(str) + '|' + df['TRIMTAG'].astype(str) + '|' + df['BLATTAG'].astype(str)
            df_out['Somatic_GSBClass'] = df['Somatic_GSBClass'] # added 8/9/2024
            df_out['GSBClass'] = df['GSBClass']
            df_out['Tier'] = df['TIER']
        else:
            if chrom_filter:
                chrom_list = ['chr'+str(i) for i in list(range(22,0,-1))] + ['chrX', 'chrY']
                df = df[df['Chr'].isin(chrom_list)]
            df_out['Marker'] = df['Chr'] + '.' + df['WU_HG19_Pos'].astype(str) + '.' + df['ReferenceAllele'] + '.' + df['MutantAllele']
            df_out['nM_D'] = df['#Mutant_In_Tumor'].fillna(0).astype(int)
            df_out['nT_D'] = df['#Total_In_Tumor'].fillna(0).astype(int)
            df_out['nM_G'] = df['#Mutant_In_Normal'].fillna(0).astype(int)
            df_out['nT_G'] = df['#Total_In_Normal'].fillna(0).astype(int)
            df_out['Somatic_GSBClass'] = 'OldFormatNoData' # added 8/9/2024
            df_out['GSBClass'] = 'OldFormatNoData'
            if is_tier23:
                df_out['Annotation'] = df['GeneName'].astype(str) + '|' + 'ND' + '|' + 'ND'
                df_out['QualTag'] = df['SJQuality'].astype(str) + '|' + 'ND' + '|' + 'ND'
                df_out['Tier'] = 'SNV_T23'
            else:
                df_out['Annotation'] = df['GeneName'].astype(str) + '|' + df['mRNA_acc'].astype(str) + '|' + df['Class'].astype(str) + '|' + df['AAChange'].astype(str)
                df_out['QualTag'] = df['SJQuality'].astype(str) + '|' + df['State'].astype(str) + '|' + 'ND'
                df_out['Tier'] = 'SNV_T1'
        df_out['Sample'] = sample
        columns = df_out.columns.tolist()
        columns.remove('Sample')
        new_columns = ['Sample'] + columns
        df_out = df_out[new_columns]
    return (df_out)

def indel_summary(input_file: str, sample: str, chrom_filter: bool) -> pd.DataFrame:
    """ summarize a single indel mutation file """
    df_out = pd.DataFrame()
    if os.path.exists(input_file) and os.path.getsize(input_file) > 0:
        df = pd.read_csv(input_file, sep = "\t", low_memory=False)
        
        columns = df.columns.tolist()
        is_new_format = False if '#Mutant_In_Tumor' in columns else True
        df_out = pd.DataFrame()
        
        if 'WU_HG18_Pos' in columns: df['WU_HG19_Pos'] = df['WU_HG18_Pos']
        if is_new_format:
            df['CHROM'] = df['CHROM'].astype(str)
            mask = ~df['CHROM'].str.startswith('chr')
            df.loc[mask, 'CHROM'] = 'chr' + df.loc[mask, 'CHROM']
            if chrom_filter:
                chrom_list = ['chr'+str(i) for i in list(range(22,0,-1))] + ['chrX', 'chrY']
                df = df[df['CHROM'].isin(chrom_list)]
            df_out['Marker'] = df['CHROM'] + '.' + df['POS'].astype(str) + '.' + df['REF'] + '.' + df['ALT']
            df_out['nM_D'] = df['ALT_F_CNT'].fillna(0).astype(int) + df['ALT_R_CNT'].fillna(0).astype(int)
            df_out['nT_D'] = df['REF_F_CNT'].fillna(0).astype(int) + df['REF_R_CNT'].fillna(0).astype(int) + df['ALT_F_CNT'].fillna(0).astype(int) + df['ALT_R_CNT'].fillna(0).astype(int)
            df_out['nM_G'] = df['GL_ALT_F_CNT'].fillna(0).astype(int) + df['GL_ALT_R_CNT'].fillna(0).astype(int)
            df_out['nT_G'] = df['GL_REF_F_CNT'].fillna(0).astype(int) + df['GL_REF_R_CNT'].fillna(0).astype(int) + df['GL_ALT_F_CNT'].fillna(0).astype(int) + df['GL_ALT_R_CNT'].fillna(0).astype(int)
            df_out['Annotation'] = df['GeneName'].astype(str) + '|' + df['mRNA_acc'].astype(str) + '|' + df['Class'].astype(str) + '|' + df['AAChange'].astype(str)
            df_out['QualTag'] =  df['QUALTAG'].astype(str)
            df_out['Somatic_GSBClass'] = df['Somatic_GSBClass'] # added 8/9/2024
            df_out['GSBClass'] = df['GSBClass']
            df_out['Tier'] = 'Indel'
            df_out['ORIGINAL_INDEL'] = df['ORIGINAL_INDEL']
            df_out['ORIGINAL_INDEL'] = df_out['ORIGINAL_INDEL'].apply(transform_marker)
        else:
            df['Chr'] = df['Chr'].astype(str)
            mask = ~df['Chr'].str.startswith('chr')
            df.loc[mask, 'Chr'] = 'chr' + df.loc[mask, 'Chr']
            if chrom_filter:
                chrom_list = ['chr'+str(i) for i in list(range(22,0,-1))] + ['chrX', 'chrY']
                df = df[df['Chr'].isin(chrom_list)]
            df_out['Marker'] = df['Chr'] + '.' + df['WU_HG19_Pos'].astype(str) + '.' + df['ReferenceAllele'] + '.' + df['MutantAllele']
            df_out['nM_D'] = df['#Mutant_In_Tumor'].fillna(0).astype(int)
            df_out['nT_D'] = df['#Total_In_Tumor'].fillna(0).astype(int)
            df_out['nM_G'] = df['#Mutant_In_Normal'].fillna(0).astype(int)
            df_out['nT_G'] = df['#Total_In_Normal'].fillna(0).astype(int)
            df_out['Annotation'] = df['GeneName'].astype(str) + '|' + df['mRNA_acc'].astype(str) + '|' + df['Class'].astype(str) + '|' + df['AAChange'].astype(str)
            df_out['QualTag'] =  df['SJQuality'].astype(str) + '|' + df['State'].astype(str)
            df_out['Somatic_GSBClass'] = 'OldFormatNoData' # added 8/9/2024
            df_out['GSBClass'] = 'OldFormatNoData'
            df_out['Tier'] = 'Indel'
            df_out['ORIGINAL_INDEL'] = df['Chr'] + '.' + df['WU_HG19_Pos'].astype(str) + '.' + df['ReferenceAllele'] + '.' + df['MutantAllele']
            df_out['ORIGINAL_INDEL'] = df_out['ORIGINAL_INDEL'].apply(transform_marker)
        df_out['Sample'] = sample
        columns = df_out.columns.tolist()
        columns.remove('Sample')
        new_columns = ['Sample'] + columns
        df_out = df_out[new_columns]
    return (df_out)

def crest_post_processing(crest_post_file: str) -> pd.DataFrame:
    """ extract sv marker information from crest-post file """
    df = pd.read_csv(crest_post_file, sep = "\t", low_memory=False)
    df_out = pd.DataFrame()
    df_out['Marker'] = df['ChrA'] + '.' + df['PosA'].astype(str) + '.' + df['OrientA'] + '.' + df['ChrB'] + '.' + df['PosB'].astype(str) + '.' + df['OrientB'] + '.'
    df_out['Fusion_Gene'] = df['Fusion Gene']
    df_out['rating'] = df['rating']
    df_out['Tsim_Result'] = df['Tsim_Result']
    df_out['strpnpolymer'] = df['strpnpolymer']
    df_out['GeneA'] = df['GeneA']
    df_out['GeneB'] = df['GeneB']
    return (df_out)

def sv_summary(input_file: str, sample: str, chrom_filter: bool, crest_post_file: str = None) -> pd.DataFrame:
    """ summarize a single sv mutation file """
    df_out = pd.DataFrame()
    if os.path.exists(input_file) and os.path.getsize(input_file) > 0:
        df = pd.read_csv(input_file, sep = "\t", low_memory=False)
        
        df['#chrA'] = df['#chrA'].astype(str)
        mask = ~df['#chrA'].str.startswith('chr')
        df.loc[mask, '#chrA'] = 'chr' + df.loc[mask, '#chrA']
        df['chrB'] = df['chrB'].astype(str)
        mask = ~df['chrB'].str.startswith('chr')
        df.loc[mask, 'chrB'] = 'chr' + df.loc[mask, 'chrB']
        if chrom_filter:
            chrom_list = ['chr'+str(i) for i in list(range(22,0,-1))] + ['chrX', 'chrY']
            df = df[df['#chrA'].isin(chrom_list) & df['chrB'].isin(chrom_list)]
        df_out['Marker'] = df['#chrA'] + '.' + df['posA'].astype(str) + '.' + df['ortA'] + '.' + df['chrB'] + '.' + df['posB'].astype(str) + '.' + df['ortB'] + '.'
        df_out['nM_A'] = df['SClipA'].fillna(0).astype(int)
        df_out['nM_B'] = df['SClipB'].fillna(0).astype(int)
        df_out['nT_A'] = df['CoverA'].fillna(0).astype(int)
        df_out['nT_B'] = df['CoverB'].fillna(0).astype(int)
        df_out['Type'] = df['Type']
        df_out['Annotation'] = df['GeneA'].astype(str) + '|' + df['GeneB'].astype(str)
        df_out['QualTag'] =  'ND'
        df_out['GSBClass'] = 'ND'
        df_out['Tier'] = 'SV'
        df_out['Sample'] = sample
        columns = df_out.columns.tolist()
        columns.remove('Sample')
        new_columns = ['Sample'] + columns
        df_out = df_out[new_columns]
        
        if crest_post_file is not None:
            df_crest_post = crest_post_processing(crest_post_file)
            merged_df = pd.merge(df_out, df_crest_post, on="Marker")
            df_out = merged_df
        
    return (df_out)

if __name__ == '__main__':
    
    parser = argparse.ArgumentParser(description="Prepare input data for dipseek from internal pipeline output",
        add_help=False)

    parser.add_argument('-h', '--help', action='help', default=argparse.SUPPRESS,
        help='works for snv | indel | sv ')

    parser.add_argument('-s', '--samples', type=str, 
        help='a 2-column (new format) or 1-column(old format) input file with a list of samples and the path to their mut record')

    parser.add_argument('-m', '--mut_type', type=str, 
        help='the type of mutation, select from snv | indel | sv, case insensitive')
    
    parser.add_argument('-x', '--exclude', type=str,
        help='a list of samples to exclude')
    
    parser.add_argument('-f', '--chrom_filter', type=bool, default=True,
        help='if only output chrom 1-22 + X + Y')
    
    parser.add_argument('-b', '--black_list', type=str,
        help='a 4-column bed file showing the genomic regions to remove')
    
    parser.add_argument('-r', '--ref', type=str, default='',
        help='version of the reference genome: hg19 or hg38')
    
    parser.add_argument('-p', '--project_id', type=str, default='Proj',
        help='project_id used to label the output files')
    
    parser.add_argument('--crest_post_list', type=str, default=None,
        help='for sv only, a list of crest-post event_fusion.txt files to get extra information')
    
    args = parser.parse_args()
    
    samples = args.samples
    mut_type = args.mut_type
    exclude = args.exclude
    chrom_filter = args.chrom_filter
    black_list = args.black_list
    ref = args.ref
    project_id = args.project_id
    crest_post_list = args.crest_post_list
    
    step = 'step0_0.txt' if exclude is None else 'step0_1.txt'
    exclude_list = []
    if exclude is None:
        step = 'step0_0.txt'
    else:
        step = 'step0_1.txt'
        ifn = open(exclude, 'r')
        for row in ifn: exclude_list.append(row.strip())
        ifn.close()
    output_file = project_id + '_' + mut_type.upper() + '_' + step
    
    if black_list is not None:
        black_list_dict = read_bed_file(black_list)
    else: black_list_dict = {}
    
    crest_post_dict = {}
    if crest_post_list is not None:
        ifn = open(crest_post_list, 'r')
        for row in ifn:
            path = row.strip()
            sample = '_'.join(path.split('/')[-1].split('_')[:2])
            crest_post_dict[sample] = path
        ifn.close()
    
    big_df = pd.DataFrame()
    ifn = open(samples, 'r')
    for row in ifn:
        cols = row.strip().split('\t')
        if len(cols) == 2: # new format
            sample, path = cols
            sample = '_'.join(sample.split('_')[:2])
        elif len(cols) == 1: # old format
            path = cols[0]
            sample = '_'.join(path.split('/')[-1].split('_')[:2])
        else: sys.exit(f'Exiting: invalid input file')
        
        if sample not in exclude_list:
            if qc_file_existence(path):
                if mut_type.upper() == 'SNV':
                    df = snv_summary(input_file=path, sample=sample, chrom_filter=chrom_filter)
                elif mut_type.upper() == 'INDEL':
                    df = indel_summary(input_file=path, sample=sample, chrom_filter=chrom_filter)
                elif mut_type.upper() == 'SV':
                    crest_post_file = crest_post_dict[sample] if sample in crest_post_dict else None
                    df = sv_summary(input_file=path, sample=sample, chrom_filter=chrom_filter, crest_post_file=crest_post_file)
                else: df = pd.DataFrame()
                big_df = pd.concat([big_df, df], ignore_index=True)
    
    ifn.close()
    
    big_df.to_csv(output_file, sep = '\t', index=False)
    if mut_type.upper() == 'INDEL':
        filter_internal_indel(output_file, 40)
    if mut_type.upper() == 'SV':
        filter_internal_sv(output_file, black_list_dict, ref)
    if mut_type.upper() == 'SNV':
        filter_internal_snv(output_file)
