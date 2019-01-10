#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os, subprocess, argparse, re

parser = argparse.ArgumentParser(description='Enter file preperation stage number & window directory')
parser.add_argument('window', metavar='Target Window Directory', type=str, help='10kb; 25kb; 50kb; 100kb; ASHCE; repeats')
parser.add_argument('chr', metavar='Chromosome Directory', type=str, help='chr#')
args = parser.parse_args()


# SPECIFYING SIZE OF TARGET i.e. floating window (50kb) or individual (ASHCE)

if args.window == str('50kb'):
    window_size = '50kb_targets'
if args.window == str('ASHCE'):
    window_size = 'individual_regions_ASHCE'


# SPECIFYING DIRECTORIES FOR (i) ROOT & (ii) OUTPUT

root = '/fastdata/bo4mhe'
window_dir = '{}/{}'.format(root, window_size)
output_dir = '{}/hpc_output'.format(window_size)


# FUNCTION: CALCULATING REGION SIZE

def calc_region_size(region_size_list):
    total_size = int(0)
    for line in region_size_list:
        split_line = line.split('\t')
        try:
            start = int(split_line[0])
            end = int(split_line[1])
            size = end - start
            total_size = total_size + size
        except ValueError:
            total_size = 'EMPTY'
    return total_size


# LISTING ALL CHROMOSOME DIRECTORIES

chr_dir_list = [chr_dir for chr_dir in os.listdir(window_dir) if chr_dir.startswith('{}'.format(args.chr))]

def chr_tuple(chr_dir):
    chr_dir = chr_dir.strip('chr').partition('-part')
    return chr_dir[0], chr_dir[2]
sorted_chr_dir_list = sorted(chr_dir_list, key=lambda x: (int(chr_tuple(x)[0]) if chr_tuple(x)[0].isdigit() else float('inf'), int(chr_tuple(x)[1]) if chr_tuple(x)[1].isdigit() else 0))

for chr_dir in sorted_chr_dir_list:


    # CREATING INDIVIDUAL CHROMOSOME OUTPUT FILE

    pipeline_output_file = open('{}/{}_pipeline_output_{}.txt'.format(output_dir, args.window, chr_dir), 'w')
    heading_list = ['target', 'start', 'end', 'target_length', 'phylip_size', 'free_np', 'free_lnL', 'unbinned_np', 'unbinned_lnL', 'binned_np', 'binned_lnL']
    [heading_list.append('#{}'.format(str(i))) for i in range(1,17+1)]
    heading_list = '\t'.join(heading_list)
    print(heading_list, file=pipeline_output_file)


    # LISTING INDIVIDUAL CHROMOSOME TARGET DIRECTORIES

    current_chr_dir = '{}/{}'.format(window_dir, chr_dir)
    os.chdir(current_chr_dir)    
    tar_dir_list = [tar_dir for tar_dir in os.listdir(current_chr_dir) if tar_dir.startswith('target_')]
    sorted_tar_dir_list = sorted(tar_dir_list, key=lambda tar: tar.split('_')[1])


    # EXAMINING INDIVIDUAL TARGET OUTPUT FILES

    for tar_dir in sorted_tar_dir_list:
        
        current_tar_dir = '{}/{}'.format(current_chr_dir, tar_dir)
        os.chdir(current_tar_dir)
        output_line = [tar_dir]


        # FINDING TARGET WINDOW SIZE

        find_tar_size_cmd = 'cat ./target.txt | cut -f 2,3'
        tar_size_list = subprocess.Popen(find_tar_size_cmd, shell=True, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()[0].decode('utf-8').rstrip('\n').split('\t')
        tar_coordinates_list = tar_size_list[0:2]
        output_line.extend(tar_coordinates_list)


        # FINDING ASHCE & REPEAT REGION SIZE

        bed_file_list = ['./bed/target_region.bed.gz']
        for bed_file in bed_file_list:
            if not os.path.isfile(bed_file):
                region_size = 'ABSENT'
            else:
                find_region_size_cmd = 'zcat {} | cut -f 2,3'.format(bed_file)
                region_size_list = subprocess.Popen(find_region_size_cmd, shell=True, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()[0].decode('utf-8').rstrip('\n').split('\n')    
                region_size = str(calc_region_size(region_size_list))
            output_line.append(region_size)


        # FINDING PHYLIP SIZE

        if not os.path.isfile('./paml_input/target_ASHCE.fa'):
            phylip_size = 'ABSENT'
        else:
            find_phylip_size_cmd = 'head -n 1 ./paml_input/target_ASHCE.fa | cut -d " " -f 2'
            phylip_size = subprocess.Popen(find_phylip_size_cmd, shell=True, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()[0].decode('utf-8').rstrip('\n')
        output_line.append(phylip_size)
        

        # FINDING BASEML RESULTS
        
        tree_file_list = ['./4baseml_free/est', './4baseml_unbinned/est', './4baseml_binned/est']
        for tree_file in tree_file_list:
            if not os.path.isfile(tree_file):
                np = lnL = 'ABSENT'
                np_lnL_list = [np, lnL]
            else:
                find_lnL_cmd = 'grep lnL {} | tr -d " " | tr -d "(" | tr -d ")"'.format(tree_file)
                lnL_output = subprocess.Popen(find_lnL_cmd, shell=True, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()[0].decode('utf-8').rstrip('\n')
                lnL_output_list = re.split(':|\+', str(lnL_output))
                np_lnL_list = lnL_output_list[2:4]
            output_line.extend(np_lnL_list)


        # FINDING BRANCH RATES FOR BINNED TREE

        if os.path.isfile('./4baseml_binned/est'):
            find_binned_rates_cmd = 'grep branches ./4baseml_binned/est | tr -s " "'
            binned_rates_list = subprocess.Popen(find_binned_rates_cmd, shell=True, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()[0].decode('utf-8').rstrip('\n').split(' ')
            binned_rates = [i for i in binned_rates_list[4:]]
            output_line.extend(binned_rates)


        # WRITING RESULTS TO OUTPUT FILE

        print(*output_line, sep='\t', file=pipeline_output_file, flush=True)
    
    pipeline_output_file.close()
