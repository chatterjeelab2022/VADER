# -*- coding: utf-8 -*-
"""
@author: Rachel Kelemen
@deputy_wizard: Joshua Brown
Finalized April 2019

Main code to process .fastq files from Illumina sequencing of libraries and selections,
filter by quality, and generate a .csv output file containing counts, enrichments,
and other stats for each sequence in a library.

Main steps
---
1. Read in parameter file
2. Identify the .fastq files for each biosample
3. Characterize run quality
4. Quality filter
5. Mismatch filter
6. Generate a list of all possible randomized sequence variants
7. Count sequences
8. Calculate enrichments
9. Generate plots showing the results of the selection
"""

import time
import os
import csv
import datetime as dt
from read_parameters import read_parameter_file
from identify_files import find_fastq_files, find_text_files
from characterize_run_quality import characterize_quality
from Q_score_filter import Q_score_filter
from mismatch_filter import mismatch_filter
from generate_possible_sequence_dict import generate_randomized_sequence_dict, generate_twist_sequence_dict
from count_randomized_sequences import count_randomized_sequences, counts_from_csv
from calculate_enrichments import calculate_enrichments
from plot_results import plot_results


# SETUP

# Copy the parameter file name and enter it here as a raw string (r'string')
parameter_file_name = r'/Users/racheltroyan/BaseSpace/Pol_III_LtR_2020_1_30/Illumina_processing_parameters_LtR_7.csv'

# Which parts of the program do you want to run?

# find_fastq_files identifies the files associated with each biosample
# Always set to True
run_find_fastq_files = True

# Run characterization
# characterize_quality outputs plots describing the overall quality of each fastq file
# limit describes how many sequence records to consider when characterizing quality
# Recommended value is 10000 - the 1st 1000 or so records often have lower quality, but
# using more than 10000 records takes significantly more time
run_characterize_quality = True
limit = 1000

# Filtering
# Q_score_filter and mismatch_filter read an input file and filter by data quality
# If run_Q_score_filter is set to False but run_mismatch_filter is set to True, the
# program looks for previously Q score filtered text files and uses those as the input
# for mismatch filtering
# If both are set to False but a subsequent step is set to True, the program looks for
# previously mismatch filtered text files and uses those as the input for future steps
run_Q_score_filter = True
run_mismatch_filter = True

# Sequence counting
# count_sequences generates a dictionary of all possible sequences and their counts
# in each biosample data set
# from_csv indicates whether the counting should be done from mismatch filtered text files
# (False) or from csv files from Jon (True)
run_count_sequences = True
from_csv = False

# Enrichment counting and plotting
# These steps require run_count_sequences to be set to True
# calculate_enrichments calculates the results for each sequence, sorts by enrichment, and
# writes the results to a .csv file.
# plot_results generates various plots describing the results of the selection(s)
run_calculate_enrichments = True
run_plot_results = True


# Track start time
start_time = time.time()


# 1. READ IN PARAMETER FILE
# read_parameter_file returns p : {'parameter_name': parameter_value}
print('/n---/nReading parameters/n---')
p = read_parameter_file(parameter_file_name)

# Make sure all required folders exist
# These will be used by subsequent functions
required_folders = ['Plots', 'Results']
for folder in required_folders:
    full_path = p['working_folder_name'] + '//' + folder
    if os.path.exists(full_path) == False:
        os.makedirs(full_path)

# Generate a list of biosamples, based on the runs listed in the parameter file
print('/n---/nIdentifying biosamples/n---')
biosamples = [p['lib_run']]
for run in p['sel_1_runs']:
    biosamples.append(run)
if p['sel_2_runs']:
    for run in p['sel_2_runs']:
        biosamples.append(run)
print('> ' + biosamples[0] + ' (library reference sample)')
for biosample in biosamples[1:]:
    print('> ' + biosample)


# 2. IDENTIFY THE .FASTQ FILES FOR EACH BIOSAMPLE
# Returns fastq_files: [[biosample_1_file_1, biosample_1_file_2...],
# [biosample_2_file_1, biosample_2_file_2...]...]
# This accomodates cases where one biosample is associated with multiple .fastq files.
if run_find_fastq_files:
    print('/n---/nLooking for .fastq files in:/n' + p['working_folder_name'])
    fastq_files = find_fastq_files(p['working_folder_name'], biosamples)


# 3. CHARACTERIZE RUN QUALITY
# Generates plots in the folder Plots \ Characterize quality plots
# Defaults to a limit of 10,000 sequences to consider
if run_characterize_quality:
    for fastq_file_list in fastq_files:
        for fastq_file in fastq_file_list:
            characterize_quality(fastq_file, p['working_folder_name'], limit,
                                 p['expected_sequence'], p['exclude_at_start'],
                                 p['exclude_at_end'],
                                 randomized_bases=p['randomized_bases'])


# 4. QUALITY FILTER
if run_Q_score_filter:
    print('/n---/nQuality filtering/n---')
    print('/nQ score filtering parameters for the entire sequence')
    print('/nNo more than {} bases below Q{}'.format(p['QF'], p['Q1']))
    print('and no bases below Q{}'.format(p['Q0']))

    if p['randomized_bases']:
        print('Q score filtering parameters for randomized bases')
        print('No randomized bases below Q{}'.format(p['Q2']))

    # Set up a list of quality-filtered files
    Q_score_filtered_files = []

    # Set up a .csv file to hold the Q score filtering stats
    stats_file_name = str(dt.date.today()) + ' Q score filtering stats.csv'
    stats_file_address = p['working_folder_name'] + '//Results//' + stats_file_name
    with open(stats_file_address, 'w', newline='') as stats_file:
        # Clear the whole file
        stats_file.truncate(0)
        header_writer = csv.writer(stats_file)
        headers = ['Biosample', 'Pct pass', 'Total reads', 'Passing reads',
                   'Short sequences', 'Q0 fail', 'Q1 fail', 'Q2 fail']
        header_writer.writerow(headers)

    # Call Q_score_filter for each biosample and its associated .fastq file(s)
    # Returns the name of the Q score filtered file and writes stats to the .csv file
    for biosample, to_filter in zip(biosamples, fastq_files):
        print('/nQ score filtering:/n' + biosample)
        filtered_file_name = Q_score_filter(p['working_folder_name'], to_filter, biosample,
                                            p['Q1'], p['QF'], p['Q0'], p['Q2'],
                                            p['expected_sequence'], p['randomized_bases'],
                                            p['exclude_at_start'], p['exclude_at_end'],
                                            stats_file_address)
        Q_score_filtered_files.append(filtered_file_name)

else:
    # If you're not running this part of the program, presumably you've already generated
    # filtered files and want to use those.
    # This calls find_text_files, which returns files from your working folder which end
    # with Q score filtered.txt
    if run_mismatch_filter:
        print('/n---/nUsing previously Q score filtered files')
        Q_score_filtered_files = find_text_files(p['working_folder_name'], biosamples,
                                                 'Q score filtered.txt')


# 6. MISMATCH FILTER
if run_mismatch_filter:
    print('/n---/nMismatch filtering/n---')
    print('Filters set to:')
    for m in range(len(p['allowed_mismatch_counts'])):
        print('> Constant region {}: allow {} mismatches'.format(
                m+1, p['allowed_mismatch_counts'][m]))

    # Set up a list of mismatch-filtered files
    mismatch_filtered_files = []

    # Set up a .csv file to hold the Q score filtering stats
    stats_file_name = str(dt.date.today()) + ' mismatch filtering stats.csv'
    stats_file_address = p['working_folder_name'] + '//Results//' + stats_file_name
    with open(stats_file_address, 'w', newline='') as stats_file:
        # Clear the whole file
        stats_file.truncate(0)
        # Write the headers
        header_writer = csv.writer(stats_file)
        headers = ['Biosample', 'Pct pass', 'Total reads', 'Passing reads']
        for constant_region_coords in p['constant_region_coords']:
            headers.append(str(constant_region_coords) + ' failed')
        header_writer.writerow(headers)

    # Call mismatch_filter for each biosample and its associated .fastq file(s)
    # Returns the name of the mismatch filtered file and stats describing how many
    # sequences passed.
    # mismatch_filter also automatically generates plots showing base and mismatch
    # frequency at each position and writes the same data to .csv files.
    for biosample, to_filter in zip(biosamples, Q_score_filtered_files):
        print('/nMismatch filtering:/n'+to_filter)
        filtered_file_name = mismatch_filter(p['working_folder_name'], to_filter,
                                             biosample, p['constant_region_coords'],
                                             p['allowed_mismatch_counts'],
                                             p['expected_sequence'], stats_file_address,
                                             tRNA_coords=p['tRNA_coords'],
                                             randomized_bases=p['randomized_bases'])
        mismatch_filtered_files.append(filtered_file_name)

else:
    # If you're not running this part of the program, presumably you've already generated
    # filtered files and want to use those.
    # This calls find_text_files, which returns files from your working folder which end
    # with mismatch_filtered.txt
    print('/n---/nUsing previously filtered files/n---')
    mismatch_filtered_files = find_text_files(p['working_folder_name'], biosamples,
                                             'mismatch filtered.txt')


# 7. COUNT HOW MANY TIMES EACH RANDOMIZED SEQUENCE VARIANT APPEARS IN THE FILTERED DATA
if run_count_sequences:
    # Generate a list of all of the sequences we expect in the library
    print('/n---/nGenerating a list of processed sequences/n---')
    if p['twist_sequence_file']:
        # Read the 'names' of sequences from the Twist sequence CSV file
        # These are the lists of randomized bases
        print('Twist library')
        possible_sequence_dict = generate_twist_sequence_dict(biosamples,
                                                            p['twist_sequence_file'])
    else:
        # Generate a list of all possible sequences given the number of randomized bases
        print('Site saturation library')
        possible_sequence_dict = generate_randomized_sequence_dict(p['randomized_bases'],
                                                                   biosamples)

    # count_randomized_sequences updates updates possible_sequence dict with counts and
    # returns it as sequence_count_dict {Sequence_1:[count_lib, count_sel_1...]...}
    # Also returns read_counts [total_lib_reads, total_sel_1_reads...]
    if from_csv:
        # Identify the .csv files based on biosample names
        csv_file_list = find_text_files(p['working_folder_name'], biosamples, '.csv')

        # Extract counts from csv files
        print('/n---/nReading randomized sequence counts from csv/n---')
        sequence_count_dict, read_counts = \
        counts_from_csv(biosamples, csv_file_list, p['working_folder_name'],
                        possible_sequence_dict)

    else:
        # Count in each mismatch-filtered sequence file
        print('/n---/nCounting randomized sequences/n---')
        sequence_count_dict, read_counts = \
        count_randomized_sequences(biosamples, mismatch_filtered_files,
                                   p['working_folder_name'], possible_sequence_dict,
                                   p['randomized_bases'])


# 8. CALCULATE ENRICHMENTS
# Create sequence_table containing the following information for all sequences:
# COLUMN   HEADER                     CONTENTS
# 0        Randomized sequence        [formatted randomized sequence]
# 1        Full stem sequence         [formatted full stem sequence]
# 2        Paired bases               [number of paired bases]
# 3        Raw counts                 [lib_count, sel_1_count, sel_2_count...]
# 4        Fraction of total          [lib_fraction, sel_1_fraction, sel_2_fraction...]
# 5        Fold enrichment            [sel_1_fold_enrichment, sel_2_fold_enrichment...]
# 6        Enrichment factors         [sel_1_enrichment_factor, sel_2_enrichment_factor...]
# 7        Enrichment ranks           [sel_1_rank, sel_2_rank...]
# 8        Average enrichment factor  [average enrichment factor]
# 9        StDev                      [enrichment factor standard deviation]
# 10       Rank                       [rank based on average enrichment factor]

if run_calculate_enrichments:
    # calculate_enrichments writes sequence_table to a csv file
    # working_folder\Results\lib_name Results.csv
    # Sequences whose counts in the library run were too low to be identified are written
    # to a separate file working_folder\Results\lib_name Low Abundance Sequences.csv
    print('/n---/nCalculating enrichments/n---')
    sequence_table, low_abundance_table = calculate_enrichments(sequence_count_dict,
        read_counts, biosamples, p['sel_1_name'], p['sel_2_name'], p['sel_2_runs'],
        p['output_format'], p['full_seq_format'], p['rand_expected_pairs'],
        p['working_folder_name'], p['lib_name'], p['min_lib_count'])


# 10. GENERATE PLOTS
if run_plot_results:
    print('/n---/nPlotting results')
    plot_results(p['working_folder_name'], sequence_table, low_abundance_table, biosamples,
                 p['randomized_bases'], p['sel_1_name'], p['sel_2_name'], p['sel_2_runs'],
                 p['expected_pairs'])


# How long did the run take?
end_time = time.time()
print ('/n---/nTotal run time: {:.0f} seconds ({:.2f} minutes)/n---'.format(
        end_time-start_time, (end_time-start_time)/60))
