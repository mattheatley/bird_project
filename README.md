# bird_project
Exploring the influence of non-coding regions on evolutionary phenotypic diversity in birds
See Yusuf et al, 2020: https://genome.cshlp.org/content/early/2020/03/19/gr.255752.119

Chromosome alignments files produced by henryjuho are either complete or split into several parts depending on size. 
The ASCHE region coordinates (265,984 in total) for corresponding chromosome files are found in target_regions/ASHCE.
Just specify root directory in the script and enter the relevant arguments listed below when promted.


1_create_target_directories.py will create a directory for each target ASHCE region from the template.


2_batch_pipeline.py will process individual targets in 6 different stages. Just specify which the stage number & also the nature of the target i.e. ‘ASHCE’. 
N.B. overlapping 50kb windows i.e. '50kb' that contain multiple individual targets were attempted initially but subsequently abandoned.

Stages:

1 - create corresponding target bed files

2 - convert bed file to fasta file

3 - run trimal on fasta file

4 - create corresponding trees using script by tonig-evo

5 - bin trees using script by tonig-evo

6 - run baseml on data for (i) free model, (ii) unbinned model & (iii) binned model.


3_check_output.py will write the following analysis results for individual targets to an output file for individual chromosomes:

- target_number
- target_start_coordinates
- target_end_coordinates
- ASHCE_region_length
- phylip_file_species_number
- free_model_parameter_number
- free_model_lnL
- unbinned_model_parameter_number
- unbinned_model_lnL
- binned_model_parameter_number
- binned_model_lnL
- binned_tree_rates (#1 to #17)

 
