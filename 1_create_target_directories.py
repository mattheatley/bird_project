#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os, shutil, argparse, re, math


# SPECIFYING ASHCE OR REPEAT REGIONS (INPUT WHEN EXECUTING SCRIPT)

parser = argparse.ArgumentParser(description='Enter window directory')
parser.add_argument('region', metavar='Target Regions Directory', type=str, help='ASHCE; repeats')
args = parser.parse_args()

if args.region == str('ASHCE'):
    region_type = 'ASHCE'
if args.region == str('repeats'):
    region_type = 'repeats'


# SPECIFYING DIRECTORIES FOR (i) ROOT, (ii) TEMPLATE TARGET & (iii) TARGETS INPUT DATA

root = '/fastdata/bo4mhe'
template_dir = '{}/template'.format(root)
regions_input_dir = '{}/target_regions/{}'.format(root, region_type)


# SPECIFYING & CREATING TOP OUTPUT DIRECTORY

window_dir = '{}/individual_{}_regions'.format(root, region_type)
regions_output_dir = '{}/regions'.format(window_dir)

if not os.path.isdir(window_dir):
    os.mkdir(window_dir)
if not os.path.isdir(regions_output_dir):
    os.mkdir(regions_output_dir)


# IMPORTING CHROMOSOME DATA

region_file_list = [region_file for region_file in os.listdir(regions_input_dir) if region_file.startswith('{}_chr'.format(region_type))]
region_file_list = sorted(region_file_list)

for region_file in region_file_list:

    # CREATING INDIVIDUAL OUTPUT DIRECTORY FOR CURRENT CHROMOSOME/PART

    chr_tuple = region_file.strip('{}_chr'.format(region_type)).rstrip('.bed').partition('-part_')
    # e.g. region_file; ASHCE_chr1-part_1.bed
    chr_num = chr_tuple[0]
    chr_part = chr_tuple[2]

    if chr_part.isdigit():
    # i.e. 1 - 5
        chr_part = '-part{}'.format(chr_part)        
    if chr_part.isalpha():
    # i.e. 'all'
        chr_part = ''

    chr_dir = '{}/chr{}{}'.format(window_dir, str(chr_num), str(chr_part))
    if not os.path.isdir(chr_dir):
        os.mkdir(chr_dir)


    # IMPORTING LIST OF ASHCE/REPEAT TARGETS FOR CURRENT CHROMOSOME/PART

    region_file = open('{}/{}'.format(regions_input_dir, region_file), 'r')
    

    # CREATING REFERENCE LIST OF TARGETS FOR CURRENT CHROMOSOME/PART
    
    chr_targets_file = open('{}/chr{}{}_targets.txt'.format(regions_output_dir, chr_num, chr_part), 'w')
    
    counter=0
    for input_line in region_file:
        counter +=1

        column = input_line.rstrip('\n').split('\t')
        # e.g. chr1	511743	511786	chr1.1-1000000.477	110
        target_num = '%.5i' %(int(counter))
        start = int(column[1])
        end = int(column[2])
        target_dir = '{}/target_{}'.format(chr_dir, target_num)


        # CREATING INDIVIDUAL OUTPUT DIRECTORY FOR CURRENT ASCHE/REPEAT TARGET FROM TEMPLATE

        if not os.path.isdir(target_dir):
            shutil.copytree('{}/'.format(template_dir), '{}/'.format(target_dir))
        target_coords_file = open('{}/target.txt'.format(target_dir), 'w')
        print('chr{}'.format(chr_num), str(start), str(end), sep='\t', file=target_coords_file)


        # PRINTING CURRENT TARGET INFORMATION TO REFERENCE LIST

        print('chr{}'.format(chr_num), 'target_{}'.format(target_num), str(start), str(end), sep='\t', file=chr_targets_file)

    target_coords_file.close()
    chr_targets_file.close()
