from dipseek_lib import *

batchSize = 500
all2all_batchSize = 200
queue = 'standard' 
initial_sleep_sec = 300
wait_sec = 50

def run_analysis_via_sys_arg(args):
    spaces = '      '
    help_raw_internal    = spaces + 'raw_internal_bigfile:   2-column: intake_sample|original_marker'
    help_namemap         = spaces + 'nameMapFile:            5-column: intake_sample|bam_path|perfect_pid|perfect_sample_tag|true_sample_tag'
    help_prefix          = spaces + 'run_prefix:             a short tag, used to label all the intermediate files to avoid clashing between different runs'
    help_curate_meta     = spaces + 'curated_meta:           the curated marker file. It is organized into 3 columns: perfect_pid|original_marker|curated_marker'
    help_ctrl_bam        = spaces + 'control_bamlist_fn:     1-column: bam_path'
    help_ctrl_count_fn   = spaces + 'ofn_control_counts:     the final output file. It is organized into 4 columns: sample_ID|curated_marker|nM|nT'
    help_cohort_count_fn = spaces + 'ofn_cohort_counts:      the final output file. It is organized into 4 columns: sample_ID|curated_marker|nM|nT'
    help_self_count_fn   = spaces + 'ofn_self_counts:        the final output file. It is organized into 4 columns: sample_ID|curated_marker|nM|nT'
    help_horiz_fn        = spaces + 'ofn_horizontal_data:    the genotyping file in horizontal format'
    help_annotat_fn      = spaces + 'annotation_result:      the annotated marker file. It is organized into 8 columns: perfect_pid|original_marker|curated_marker|5_anno_cols'
    help_platform_meta   = spaces + 'meta_platformfile:      2-column: platform|path_to_namemap'
    help_good_marker     = spaces + 'good_marker_ifn:        1-column: PID.marker'
    help_aicheck_fn      = spaces + 'aichecklist_filename:   1-column: ai check files'
    help_flankingdist    = spaces + 'flanking distance (bp): integer; used to build contig'
    help_blat_fn         = spaces + 'blat_ofile:             the final output file. It is organized into columns: curated_marker and blat score|identify from top-2 '

    if len(args)==1:
        print('usage: python ' + args[0] + ' curation|genotyping|annotation|omni|vertical|bam2ai ' + '\n')
        
    if len(args)==2:
        if (args[1]=='curation'):
            print('usage: python ' + args[0] + ' curation <sv|indel|snv> <raw_internal_bigfile> <nameMapFile> <run_prefix> <outputfile_curated_meta> <hg19|hg38>\n')
            print('\n'.join([help_raw_internal, help_namemap, help_prefix, help_curate_meta]) + '\n')
        if (args[1]=='genotyping'): 
            print('usage: python ' + args[0] + ' genotyping self|cohort|control|merge|crossplatform\n')
        if (args[1]=='omni'): 
            print('usage: python ' + args[0] + ' omni <sv|indel|snv> <raw_internal_bigfile> <nameMapFile> <control_bam_list> <run_prefix> <hg19|hg38> <platform:WGS|WES|VAL|RNA>\n')
            print('\n'.join([help_raw_internal, help_namemap, help_ctrl_bam, help_prefix]) + '\n')
        if (args[1]=='annotation'): 
            print('usage: python ' + args[0] + ' annotation <indel|snv> <meta_bigfile> <run_prefix> <outputfile_annotated_meta> <hg19|hg38>\n')
            print('\n'.join([help_curate_meta, help_namemap, help_prefix, help_annotat_fn]) + '\n')
        if (args[1]=='vertical'): 
            print('usage: python ' + args[0] + ' vertical <horizontal_file> <good_marker_ifn> <outputFile> <min_nM> <min_nT>\n')
            print('\n'.join([help_horiz_fn, help_good_marker]) + '\n')
        if (args[1]=='bam2ai'): 
            print('usage: python ' + args[0] + ' bam2ai count|plot|dupcheck|omni\n')

           
    if len(args)>=3:
        if (args[1]=='curation'):
            if  len(args)==8:
                marker_type = args[2]
                raw_internal_bigfile = args[3]
                nameMapFile = args[4]
                projectID = args[5]
                curated_meta_ofile = args[6]
                hg_build = args[7]
                if marker_type.upper() in ['INDEL', 'SV', 'MNV']:
                    run_sv_or_indel_curation(raw_internal_bigfile, marker_type, batchSize, projectID, nameMapFile, curated_meta_ofile, queue, initial_sleep_sec, wait_sec, hg_build)
                elif marker_type.upper() == 'SNV':
                    run_snv_curation(raw_internal_bigfile, curated_meta_ofile, nameMapFile)
            else:
                print('usage: python ' + args[0] + ' curation <sv|indel|snv> <raw_internal_bigfile> <nameMapFile> <run_prefix> <outputfile_curated_meta> <hg19|hg38>\n')
                print('\n'.join([help_raw_internal, help_namemap, help_prefix, help_curate_meta]) + '\n')

        if (args[1]=='annotation'):
            if  len(args)==7:
                marker_type = args[2]
                metafile = args[3]
                projectID = args[4]
                annotation_ofile = args[5]
                hg_build = args[6]
                if marker_type.upper() in ['SNV', 'INDEL', 'MNV']:
                    snv4column = 2
                    aa26_run_snvindelmnv_annotation(metafile, snv4column, hg_build, annotation_ofile, projectID)
                elif marker_type.upper() == 'SV':
                    run_sv_annotation(metafile, annotation_ofile, batchSize, hg_build, queue, projectID, initial_sleep_sec, wait_sec)
            else:
                print('usage: python ' + args[0] + ' annotation <indel|snv> <meta_bigfile> <run_prefix> <outputfile_annotated_meta> <hg19|hg38>\n')
                print('\n'.join([help_curate_meta, help_namemap, help_prefix, help_annotat_fn]) + '\n')
        
        if (args[1]=='genotyping'):
            if (args[2]=='crossplatform'): 
                if (len(args)==10):
                    marker_type = args[3]
                    meta_bigfile = args[4]
                    meta_platformfile = args[5]
                    controlgeno_all2all_fn = args[6] # 'NULL' if skip error profiling step
                    projectID = args[7]
                    hg_build = args[8]
                    outputFile = args[9]
                    run_crossplatform_genotyping(meta_bigfile, marker_type, batchSize, projectID, meta_platformfile, controlgeno_all2all_fn, outputFile, queue, initial_sleep_sec, wait_sec, hg_build)
                else:
                    print('usage: python ' + args[0] + ' genotyping crossplatform <sv|indel|snv> <meta_bigfile> <meta_platformfile> <ctrlgeno_bigfile> <run_prefix> <hg19|hg38> <outputfile>\n')
                    print('\n'.join([help_curate_meta, help_platform_meta, help_ctrl_count_fn, help_prefix]) + '\n')
            
            elif (args[2]=='merge'):
                if (len(args)==12):
                    marker_type = args[3] # 05/13/2024
                    metafile = args[4]
                    annotation_fn = args[5]
                    selfgeno_all2all_fn = args[6]
                    controlgeno_all2all_fn = args[7]
                    nameMapFile = args[8]
                    outputFile = args[9]
                    platform = args[10]
                    hg_build = args[11]
                    aa25_make_horizontal_file_from_big_count_files(metafile, marker_type, annotation_fn, selfgeno_all2all_fn, controlgeno_all2all_fn, nameMapFile, outputFile, platform, hg_build)
                else:
                    print('usage: python ' + args[0] + ' genotyping merge <sv|indel|snv> <meta_bigfile> <annotation_result> <selfgeno_bigfile> <ctrlgeno_bigfile> <nameMapFile> <horiz_file> <WES|WGS|Val> <hg19|hg38>\n')
                    print('\n'.join([help_curate_meta, help_annotat_fn, help_cohort_count_fn, help_ctrl_count_fn, help_namemap, help_horiz_fn]) + '\n')

            elif (args[2]=='control'):
                if (len(args)==9):
                    marker_type = args[3]
                    meta_bigfile = args[4]
                    control_bamlist_filename = args[5]
                    projectID = args[6]
                    ofn_control_counts = args[7]
                    hg_build = args[8]
                    if marker_type.upper() in ['INDEL', 'MNV']:
                        run_all2all_genotyping_for_indel(meta_bigfile, all2all_batchSize, hg_build, svindel_typer, control_bamlist_filename, queue, projectID, ofn_control_counts, initial_sleep_sec, wait_sec)
                    elif marker_type.upper() == 'SV':
                        run_all2all_genotyping_for_sv(meta_bigfile, all2all_batchSize, hg_build, svindel_typer, control_bamlist_filename, queue, projectID, ofn_control_counts, initial_sleep_sec, wait_sec)
                    elif marker_type.upper() == 'SNV':
                        run_snv_control_genotyping(meta_bigfile, all2all_batchSize, projectID, control_bamlist_filename, ofn_control_counts, queue, initial_sleep_sec, wait_sec, hg_build)
                else:
                    print('usage: python ' + args[0] + ' genotyping control <sv|indel|snv> <meta_bigfile> <control_bamlist_fn> <run_prefix> <ofn_control_counts> <hg19|hg38>\n')
                    print('\n'.join([help_curate_meta, help_ctrl_bam, help_prefix, help_ctrl_count_fn]) + '\n')
                    
            elif (args[2]=='cohort'):
                if (len(args)==9):
                    marker_type = args[3]
                    meta_bigfile = args[4]
                    nameMapFile = args[5]
                    projectID = args[6]
                    ofn_cohort_counts = args[7]
                    hg_build = args[8]
                    if marker_type.upper() in ['INDEL', 'SV', 'MNV']:
                        run_sv_or_indel_cohort_genotyping(meta_bigfile, marker_type, all2all_batchSize, projectID, nameMapFile, ofn_cohort_counts, queue, initial_sleep_sec, wait_sec, hg_build)
                    elif marker_type.upper() == 'SNV':
                        run_snv_control_genotyping(meta_bigfile, all2all_batchSize, projectID, nameMapFile, ofn_cohort_counts, queue, initial_sleep_sec, wait_sec, hg_build)
                else:
                    print('usage: python ' + args[0] + ' genotyping cohort <sv|indel|snv> <meta_bigfile> <nameMapFile> <run_prefix> <ofn_cohort_counts> <hg19|hg38>\n')
                    print('\n'.join([help_curate_meta, help_namemap, help_prefix, help_cohort_count_fn]) + '\n')
                    
            elif (args[2]=='self'): 
                if (len(args)==9):
                    marker_type = args[3]
                    meta_bigfile = args[4]
                    nameMapFile = args[5]
                    projectID = args[6]
                    selfgenotyping_ofile = args[7]
                    hg_build = args[8]
                    if marker_type.upper() in ['INDEL', 'SV', 'MNV']:
                        run_sv_or_indel_self_genotyping(meta_bigfile,marker_type, batchSize, projectID, nameMapFile, selfgenotyping_ofile, queue,initial_sleep_sec, wait_sec,hg_build)
                    elif marker_type.upper() == 'SNV':
                        run_snv_control_genotyping(meta_bigfile, batchSize, projectID, nameMapFile, selfgenotyping_ofile, queue, initial_sleep_sec, wait_sec, hg_build)
                else:
                    print('usage: python ' + args[0] + ' genotyping self <sv|indel|snv> <meta_bigfile> <nameMapFile> <run_prefix> <self_genotyping_counts> <hg19|hg38>\n')
                    print('\n'.join([help_curate_meta, help_namemap, help_prefix, help_self_count_fn]) + '\n')
        
        if (args[1]=='blat'):
            if  len(args)==6:
                metafile = args[2]
                flanking_distance = int(args[3])
                blat_ofile = args[4]
                hg_build = args[5]
                run_blat_for_meta(metafile, flanking_distance, hg_build, blat_ofile)
            else:
                print('usage: python ' + args[0] + ' blat <meta_bigfile> <flanking distance (bp)> <outputfile_blat> <hg19|hg38>\n')
                print('\n'.join([help_curate_meta, help_flankingdist, help_blat_fn]) + '\n')
                    
        if (args[1]=='omni'):
            if (len(args)==9) or (len(args)==10):
                marker_type = args[2]
                raw_filename= args[3]
                namemapFile = args[4]
                control_bam_list = args[5]
                projectID = args[6]
                hg_build = args[7]
                platform = args[8]
                if_incohort = sys.argv[9] if len(args)==10 else 'TRUE'
                if marker_type.upper() in ['INDEL', 'SV', 'MNV']:
                    run_sv_or_indel_omni(raw_filename, marker_type, namemapFile, platform, control_bam_list, batchSize, all2all_batchSize, projectID, hg_build, if_incohort, queue, initial_sleep_sec, wait_sec)
                elif marker_type.upper() == 'SNV':
                    run_snv_omni(raw_filename, namemapFile, platform, control_bam_list, batchSize, all2all_batchSize, projectID, hg_build, queue, initial_sleep_sec, wait_sec)
            else:
                print('usage: python ' + args[0] + ' omni <sv|indel|snv> <raw_internal_bigfile> <nameMapFile> <control_bam_list> <run_prefix> <hg19|hg38> <platform:WGS|WES|VAL|RNA> <if_incohort:TRUE|FALSE>\n')
                print('\n'.join([help_raw_internal, help_namemap, help_ctrl_bam, help_prefix]) + '\n')

        if (args[1]=='vertical'):
            if len(args)==7:
                horizontal_file = args[2]
                good_marker_ifn = args[3]
                outputFile = args[4]
                min_nM = args[5]
                min_nT = args[6]
                horiz_to_vertical(horizontal_file, good_marker_ifn, outputFile, int(min_nM), int(min_nT))
            else:
                print('usage: python ' + args[0] + ' vertical <horizontal_file> <good_marker_ifn> <outputFile> <min_nM> <min_nT>\n')
                print('\n'.join([help_horiz_fn, help_good_marker]) + '\n')

        if (args[1]=='bam2ai'):
            if (args[2]=='count'): 
                if (len(args)==7):
                    bamlist_filename = args[3]
                    projectID = args[4]
                    ofn_counts = args[5]
                    hg_build = args[6]
                    run_bam2ai_count(all2all_batchSize, projectID, bamlist_filename, ofn_counts, queue, initial_sleep_sec, wait_sec, hg_build)
                else:
                    print('usage: python ' + args[0] + ' bam2ai count <bamlist_fn> <run_prefix> <ofn_counts> <hg19|hg38>\n')
                    print('\n'.join([help_ctrl_bam, help_prefix, help_ctrl_count_fn]) + '\n')
            
            elif (args[2]=='plot'): 
                if (len(args)==7):
                    geno_all2all_fn = args[3]
                    projectID = args[4]
                    hg_build = args[5]
                    output_dir = args[6]
                    plot_geno2ai(geno_all2all_fn, projectID, output_dir, queue, initial_sleep_sec, wait_sec, hg_build)
                else:
                    print('usage: python ' + args[0] + ' bam2ai plot <all2all_geno_fn> <run_prefix> <hg19|hg38> <output_dir>\n')
                    print('\n'.join([help_ctrl_count_fn, help_prefix]) + '\n')
            
            elif (args[2]=='dupcheck'): 
                if (len(args)==5):
                    aichecklist_filename = args[3]
                    projectID = args[4]
                    ai_dupcheck(aichecklist_filename, projectID, queue, initial_sleep_sec, wait_sec)
                else:
                    print('usage: python ' + args[0] + ' bam2ai ai_dupcheck <aichecklist_filename> <run_prefix>\n')
                    print('\n'.join([help_aicheck_fn, help_prefix]) + '\n')
            
            elif (args[2]=='omni'): 
                if (len(args)==7):
                    bamlist_filename = args[3]
                    projectID = args[4]
                    hg_build = args[5]
                    output_dir = args[6]
                    run_bam2ai_omni(all2all_batchSize, projectID, bamlist_filename, output_dir, queue, initial_sleep_sec, wait_sec, hg_build)
                else:
                    print('usage: python ' + args[0] + ' bam2ai omni <bamlist_fn> <run_prefix> <hg19|hg38> <output_dir>\n')
                    print('\n'.join([help_ctrl_bam, help_prefix]) + '\n')
                      
run_analysis_via_sys_arg(sys.argv)
