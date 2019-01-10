#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os, subprocess, argparse, re, time

parser = argparse.ArgumentParser(description='Enter file preperation stage number & window directory')
parser.add_argument('stage', metavar='Stage Number', choices=range(1,8), type=int, help='1: create bed files\n2: create fasta files\n3: run trimal\n4: create trees\n5:mark binned branches\n6:run baseml')
parser.add_argument('window', metavar='Target Window Directory', type=str, help='50kb; ASHCE')
args = parser.parse_args()


# SPECIFYING STAGE OF FILE PREPARATION (1-bed, 2-fasta, 3-trimal, 4-trees, 5-binned, 6-baseml)

if args.stage == int(1):
    preparation_stage = '1-bed'
if args.stage == int(2):
    preparation_stage = '2-fasta'
if args.stage == int(3):
    preparation_stage = '3-trimal'
if args.stage == int(4):
    preparation_stage = '4-trees'
if args.stage == int(5):
    preparation_stage = '5-binned'
if args.stage == int(6):
    preparation_stage = '6-baseml'


# SPECIFYING SIZE OF TARGET i.e. floating window (50kb) or individual (ASHCE)

if args.window == str('50kb'):
    window_size = '50kb_targets'
if args.window == str('ASHCE'):
    window_size = 'individual_ASHCE_regions'


# FUNCTION: QSUB FEED LOOP FOR CONTINUOUSLY SUBMITTING JOBS ON SHARC HPCC (FINAL STEP)

def qsub_queue(cmd, job_limit):
    check_qstat_cmd = 'qstat | grep -c "chr"'
    submitted_num = subprocess.Popen(check_qstat_cmd, shell=True, stdout=subprocess.PIPE).communicate()[0].decode("utf-8").split('\n')[0]
    submitted_num = int(submitted_num)
    if submitted_num >= job_limit:
        time.sleep(60)
        qsub_queue(cmd, job_limit)
    else:
        subprocess.call(cmd, shell=True)
    return


# SPECIFYING DIRECTORIES FOR (i) ROOT, (ii) CURRENT TARGETS WINDOW, (iii) CHROMOSOME ALIGNMENT FILES, (iv) TARGETS INPUT DATA, (v) TEMPORARY BATCH SCRIPTS & (vi) OUTPUT

root = '/fastdata/bo4mhe'
window_dir = '{}/{}'.format(root, window_size)
chr_files_dir = '{}/chr_files'.format(root)
region_files_dir = '{}/target_regions/ASHCE'.format(root)
script_dir = '{}/batch_scripts'.format(root)
output_dir = '{}/hpc_output'.format(window_dir)


# CREATING TOP OUT/ERROR & OUTPUT DIRECTORIES

if not os.path.isdir('{}/out_&_error'.format(window_dir)):
    os.mkdir('{}/out_&_error'.format(window_dir))
if not os.path.isdir(output_dir):
    os.mkdir(output_dir)
if not os.path.isdir('{}/out_&_error'.format(output_dir)):
    os.mkdir('{}/out_&_error'.format(output_dir))


# LISTING ALL/SPECIFIC CHROMOSOME DIRECTORIES

chr_dir_list = [chr_dir for chr_dir in os.listdir(window_dir) if chr_dir.startswith('chr')]
#chr_dir_list = [chr_dir for chr_dir in os.listdir(window_dir) if re.search('chr[345]-', chr_dir)]
#chr_dir_list = ['chr{}'.format(i) for i in range(20, 29, 1)]

def chr_tuple(chr_dir):
    chr_dir = chr_dir.strip('chr').partition('-part')
    return chr_dir[0], chr_dir[2]
sorted_chr_dir_list = sorted(chr_dir_list, key=lambda x: (int(chr_tuple(x)[0]) if chr_tuple(x)[0].isdigit() else float('inf'), int(chr_tuple(x)[1]) if chr_tuple(x)[1].isdigit() else 0))


for chr_dir in sorted_chr_dir_list:
    current_chr_dir = '{}/{}'.format(window_dir, chr_dir)
    os.chdir(current_chr_dir)

    # CREATING INDIVIDUAL CHROMOSOME DIRECTORIES FOR (i) TEMPORARY BATCH SCRIPT & (ii) OUT/ERROR

    batch_dir_list = ['batch_scripts', 'batch_scripts/out_&_error']
    for batch_dir in batch_dir_list:
        if not os.path.isdir('{}/{}'.format(current_chr_dir, batch_dir)):
            os.mkdir('{}/{}'.format(current_chr_dir, batch_dir))


    # LISTING INDIVIDUAL CHROMOSOME TARGET DIRECTORIES

    tar_dir_list = [tar_dir for tar_dir in os.listdir(current_chr_dir) if tar_dir.startswith('target_')]
    sorted_tar_dir_list = sorted(tar_dir_list, key=lambda tar: tar.split('_')[1])


    # CREATING INDIVIDUAL CHROMOSOME CREATE_BED_FROM_TARGET SCRIPT (STAGE 1)

    if args.stage == int(1):
        create_bed_from_target_script = open('{}/create_bed_from_target_{}.py'.format(script_dir, chr_dir), 'w')
        create_bed_from_target_script.write('import os, sys\n\n')
        create_bed_from_target_script.write('bedtools="bedtools"\n')
        create_bed_from_target_script.write('bed_maf="{}/Gallus_gallus3_{}.all_birds.JetzHackett_Genbank.clean.wga.bed.gz"\n'.format(chr_files_dir, chr_dir))
        create_bed_from_target_script.write('target="target.txt"\n')
        create_bed_from_target_script.write('name=target.split(".")[0]\n\n')
        
        if args.window == str('ASHCE'):
            create_bed_from_target_script.write('os.system(bedtools+" intersect -a "+bed_maf+ " -b "+target+" | gzip -c > ./bed/target_region.bed.gz")\n\n')
        
        if 'kb' in args.window:
            create_bed_from_target_script.write('os.system(bedtools+" intersect -a "+bed_maf+ " -b "+target+" > tmp_test.bed")\n\n')
            create_bed_from_target_script.write('bed_maf="tmp_test.bed"\n')
            create_bed_from_target_script.write('target="{}/ASHCE.bed"\n'.format(region_files_dir))
            create_bed_from_target_script.write('os.system(bedtools+" intersect -a "+bed_maf+ " -b "+target+" | gzip -c > ./bed/"+name+"_region.bed.gz")\n\n')
        create_bed_from_target_script.close()


    # LISTING TARGET SUBLISTS TO SUBMIT AS BATCH JOBS

    if args.stage < int(6):
        tar_sublist_list = [sorted_tar_dir_list[i:i+50] for i in range(0, len(sorted_tar_dir_list), 50)]
    if args.stage == int(6):
        tar_sublist_list = [sorted_tar_dir_list[i:i+10] for i in range(0, len(sorted_tar_dir_list), 10)]

    counter = 0
    for tar_sublist in tar_sublist_list:
        tar_sublist_num = tar_sublist_list.index(tar_sublist)+1


        # CREATING INDIVIDUAL TARGET SUBLIST BASH SCRIPT

        pipe_sublist_qsub = open('{}/batch_scripts/pipe_sublist_{}_{}.sh'.format(current_chr_dir, preparation_stage, tar_sublist_num), 'w')
        pipe_sublist_qsub.write('#!/bin/bash\n')
        pipe_sublist_qsub.write('#$ -l h_rt=08:00:00\n')
        pipe_sublist_qsub.write('#$ -l rmem=2G\n')
        pipe_sublist_qsub.write('#$ -N {}_{}\n'.format(chr_dir, tar_sublist_num))
        pipe_sublist_qsub.write('#$-o {}/batch_scripts/out_&_error/pipe_sublist_{}_{}.out\n'.format(current_chr_dir, preparation_stage, tar_sublist_num))
        pipe_sublist_qsub.write('#$-e {}/batch_scripts/out_&_error/pipe_sublist_{}_{}.error\n'.format(current_chr_dir, preparation_stage, tar_sublist_num))
        pipe_sublist_qsub.write('#$-V\n')
        pipe_sublist_qsub.write('cd {}/batch_scripts\n'.format(current_chr_dir))
        pipe_sublist_qsub.write('python {}/batch_scripts/pipe_sublist_{}_{}.py'.format(current_chr_dir, preparation_stage, tar_sublist_num))
        pipe_sublist_qsub.close()


        # CREATING INDIVIDUAL TARGET SUBLIST PYTHON SCRIPT

        pipe_sublist_script = open('{}/batch_scripts/pipe_sublist_{}_{}.py'.format(current_chr_dir, preparation_stage, tar_sublist_num), 'w')
        pipe_sublist_script.write('import os, sys, shutil, subprocess\n\n')
        pipe_sublist_script.write('window_dir = "{}"\n'.format(current_chr_dir))
        pipe_sublist_script.write('script_dir = "{}"\n'.format(script_dir))
        pipe_sublist_script.write('script_name = "{}"\n'.format(chr_dir))
        pipe_sublist_script.write('tar_dir_sublist = {}\n'.format(tar_sublist))
        pipe_sublist_script.write('sublist_len = len(tar_dir_sublist)\n')
        pipe_sublist_script.write('sorted_tar_dir_sublist = sorted(tar_dir_sublist)\n\n')
        pipe_sublist_script.write('counter = 0\n')
        pipe_sublist_script.write('for tar_dir in sorted_tar_dir_sublist:\n')
        pipe_sublist_script.write('\tcurrent_tar_dir = "{}/{}".format(window_dir, tar_dir)\n')
        pipe_sublist_script.write('\tos.chdir(current_tar_dir)\n')
        pipe_sublist_script.write('\tif counter == 0:\n')
        pipe_sublist_script.write('\t\tstart_time = subprocess.check_output("date", shell=True).decode("utf-8").rstrip("\\n")\n')
        pipe_sublist_script.write('\t\tprint("{} commenced | {} of {} targets processed in total | {}".format(tar_dir, str(counter), str(sublist_len), start_time))\n')
        pipe_sublist_script.write('\t\tsys.stdout.flush()\n')
        pipe_sublist_script.write('\tcounter +=1\n')
        
        # N.B. STEPS 1 TO 5 NEED TO BE RUN IN AN ENVIRONMENT WITH BOTH PYTHON 2 AND RELEVANT MODULES!

        # STAGE 1 CALLING TONI'S SCRIPT FOR CREATING BED FILES FROM TARGET COORDINATES
        if args.stage == int(1):
            pipe_sublist_script.write('\tsubprocess.call("python {}/create_bed_from_target_{}.py".format(script_dir, script_name), shell=True)\n')
       
        # STAGE 2 CALLING TONI'S SCRIPT FOR CREATING FASTA FILES FROM BED FILES
        if args.stage == int(2):
            pipe_sublist_script.write('\tsubprocess.call("python {}/create_fasta_files_from_bed.py".format(script_dir), shell=True)\n')
       
        # STAGE 3 CALLING TONI'S SCRIPT FOR RUNNING TRIMAL
        if args.stage == int(3):
            pipe_sublist_script.write('\tsubprocess.call("python {}/run_trimal.py".format(script_dir), shell=True)\n')
        
        # CHANGING TO SUBDIRECTORY FOR SUBSEQUENT STAGES
        if args.stage == int(4) or args.stage == int(5):
            pipe_sublist_script.write('\tos.chdir("{}/4tree".format(current_tar_dir))\n')

        # STAGE 4 CALLING TONI'S SCRIPT FOR CREATING TREES
        if args.stage == int(4):
            pipe_sublist_script.write('\tsubprocess.call("python {}/4tree/create_trees.py".format(script_dir), shell=True)\n')
       
        # STAGE 5 CALLING TONI'S SCRIPT FOR MARKING BINNED BRANCHES
        if args.stage == int(5):
            pipe_sublist_script.write('\tsubprocess.call("python {}/4tree/mark_binned_branch.py".format(script_dir), shell=True)\n')

        # STAGE 6 RUNNING BASEML
        if args.stage == int(6):
            pipe_sublist_script.write('\tphylip_file = "{}/paml_input/target_region.fa".format(current_tar_dir)\n')
            pipe_sublist_script.write('\tbinned_tree = "{}/4tree/binned_tree/target_region.fa.tre".format(current_tar_dir)\n')
            pipe_sublist_script.write('\tunbinned_tree = "{}/4tree/final_trees/target_region.fa.tre".format(current_tar_dir)\n')
        
            #BINNED TREE
            pipe_sublist_script.write('\tif os.path.isfile(phylip_file) and os.path.isfile(binned_tree):\n')
            pipe_sublist_script.write('\t\tshutil.copy2(phylip_file, "{}/4baseml_binned/".format(current_tar_dir))\n')
            pipe_sublist_script.write('\t\tshutil.copy2(binned_tree, "{}/4baseml_binned/".format(current_tar_dir))\n')
            pipe_sublist_script.write('\t\tos.chdir("{}/4baseml_binned".format(current_tar_dir))\n')
            pipe_sublist_script.write('\t\tsubprocess.call("baseml", shell=True)\n')

            #UNBINNED TREE
            pipe_sublist_script.write('\tif os.path.isfile(phylip_file) and os.path.isfile(unbinned_tree):\n')
            pipe_sublist_script.write('\t\tshutil.copy2(phylip_file, "{}/4baseml_unbinned/".format(current_tar_dir))\n')
            pipe_sublist_script.write('\t\tshutil.copy2(unbinned_tree, "{}/4baseml_unbinned/".format(current_tar_dir))\n')
            pipe_sublist_script.write('\t\tos.chdir("{}/4baseml_unbinned".format(current_tar_dir))\n')
            pipe_sublist_script.write('\t\tsubprocess.call("baseml", shell=True)\n')

            #FREE TREE
            pipe_sublist_script.write('\tif os.path.isfile(phylip_file) and os.path.isfile(unbinned_tree):\n')
            pipe_sublist_script.write('\t\tshutil.copy2(phylip_file, "{}/4baseml_free/".format(current_tar_dir))\n')
            pipe_sublist_script.write('\t\tshutil.copy2(unbinned_tree, "{}/4baseml_free/".format(current_tar_dir))\n')
            pipe_sublist_script.write('\t\tos.chdir("{}/4baseml_free".format(current_tar_dir))\n')
            pipe_sublist_script.write('\t\tsubprocess.call("baseml", shell=True)\n')

        # WRITING PROGRESS SUMMARY TO OUT
        pipe_sublist_script.write('\tend_time = subprocess.check_output("date", shell=True).decode("utf-8").rstrip("\\n")\n')
        pipe_sublist_script.write('\tprint("{} processed | {} of {} targets processed in total | {}".format(tar_dir, str(counter), str(sublist_len), end_time))\n')
        pipe_sublist_script.write('\tsys.stdout.flush()\n')
        pipe_sublist_script.write('print("All targets completed | {}".format(end_time))\n')
        pipe_sublist_script.close()

        # SUBMITING INDIVIDUAL TARGET SUBLIST BASH SCRIPT TO HPC
        qsub_cmd = 'qsub {}/batch_scripts/pipe_sublist_{}_{}.sh'.format(current_chr_dir, preparation_stage, tar_sublist_num)
        qsub_queue(qsub_cmd, 200)
        #subprocess.call('qsub {}/batch_scripts/pipe_sublist_{}_{}.sh'.format(current_chr_dir, preparation_stage, tar_sublist_num), shell=True)

