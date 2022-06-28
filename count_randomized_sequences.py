# -*- coding: utf-8 -*-
"""
@author: Rachel Kelemen
@deputy wizard: Joshua Brown
Finalized April 2019

This module contains functions which counts the frequency of each randomized sequence in
the filtered data for each biosample.
count_randomized_sequences is for quality and mismatch filtered text files.
counts_from_csv is for csv files of sequences and counts from Jon.

Returns a dictionary of sequences and counts:
{sequence: [lib_count, sel_1_count, sel_2_count...]}
"""

import numpy as np
import csv

def count_randomized_sequences(biosamples, filtered_file_list, working_folder_name,
                               possible_sequence_dict, randomized_bases):
    """Count the frequency of each randomized sequence in the filtered data for each
    biosample.

    Parameters
    ---
    biosamples : list
        Biosample names
    filtered_file_list : list
        Filtered files to read
    working_folder_name : str
        Working folder
    possible_sequence_dict : dict
        Dictionary of the form ['Randomized sequence': [0, 0...]]
    randomized_bases : list
        Randomized base coordinates

    Returns
    ---
    possible_sequence_dict : dict
        Dictionary of the form ['Randomized sequence': [Lib_count, sel_1_count...]]
    read_count : list
        List of the total identified reads for each biosample.
    stats : list
        Stats summarizing the filtering results
        [[Biosample, processed_count, found, not_found]...]

    Notes
    ---
    Reads filtered files and updates possible_sequence_dict with counts

    Using randomized sequences as dictionary keys is much faster than using a list
    because dictionary keys are hashed, making searching much more efficient.
    """

    # Setup
    # stats holds statistics describing how many read matched expected library sequences
    # for each biosample. [[biosample, processed_count, found, not_found]...]
    stats = []
    # read_counts holds the number of identified reads which matched expected library
    # sequences for each biosample. This will be used later to calculate enrichment
    # factors.
    read_counts = []

    # Iterate through biosamples
    for i in list(np.arange(len(biosamples))):
        # Track how many sequences were found in the expected library sequences
        processed_count = 0
        found = 0
        not_found = 0

        print('/nProcessing', biosamples[i])

        # Open the filtered sequence file and read each sequence.
        filtered_file_address = working_folder_name + '//' + filtered_file_list[i]
        with open(filtered_file_address) as filtered_file:
            for sequence in filtered_file:
                processed_count += 1
                # Extract randomized bases
                extract_random = ''
                for b in randomized_bases:
                    extract_random += sequence[b]
                # Exclude sequences which are not in the expected sequence list from analysis
                if extract_random in possible_sequence_dict:
                    possible_sequence_dict[extract_random][i] += 1
                    found += 1

                else:
                    not_found += 1

        # Counting statistics
        stats.append([biosamples[i], processed_count, found, not_found])
        read_counts.append(found)

        print('{:,} sequences processed'.format(processed_count))
        print('{:,} sequences ({:.2%}) matched possible sequences'.format(found,
              found/processed_count))
        print('{:,} sequences ({:.2%}) did not match possible sequences'.format(not_found,
              not_found/processed_count))


    # Write counting statistics to a .csv file
    headers = ['Biosample', 'Reads processed', 'Matched possible sequences', 'No match']
    stats_file_name = 'Count possible sequence stats.csv'
    stats_file_address = working_folder_name + '//Results///' + stats_file_name
    with open(stats_file_address, 'w', newline='') as stats_file:
        stats_writer = csv.writer(stats_file)
        stats_writer.writerow(headers)
        for biosample_stats in stats:
            stats_writer.writerow(biosample_stats)


    return possible_sequence_dict, read_counts



def counts_from_csv(biosamples, csv_file_list, working_folder_name,
                        possible_sequence_dict):
    """Count the frequency of each randomized sequence in the filtered data for each
    biosample.

    Parameters
    ---
    biosamples : list
        Biosample names
    csv_file_list : list
        CSV files to read
    working_folder_name : str
        Working folder
    possible_sequence_dict : dict
        Dictionary of the form ['Randomized sequence': [0, 0...]]

    Returns
    ---
    possible_sequence_dict : dict
        Dictionary of the form ['Randomized sequence': [Lib_count, sel_1_count...]]
    read_count : list
        List of the total identified reads for each biosample.
    stats : list
        Stats summarizing the filtering results
        [[Biosample, processed_count, found, not_found]...]
    """

    # Setup
    # stats holds statistics describing how many read matched expected library sequences
    # for each biosample. [[biosample, processed_count, found, not_found]...]
    stats = []
    # read_counts holds the number of identified reads which matched expected library
    # sequences for each biosample. This will be used later to calculate enrichment
    # factors.
    read_counts = []

    # Iterate through biosamples:
    # Iterate through biosamples
    for i in list(np.arange(len(biosamples))):
        # Track how many sequences were found in the expected library sequences
        line_count = 0
        processed_count = 0
        found = 0
        not_found = 0

        print('/nProcessing', biosamples[i])

        # Open the csv file and read each sequence.
        csv_file_address = working_folder_name + '//' + csv_file_list[i]
        with open(csv_file_address, 'r', newline='') as csv_file:
            count_reader = csv.reader(csv_file)
            for line in count_reader:
                line_count += 1
                extract_random, count = line[0], int(line[1])
                processed_count += count
                if extract_random in possible_sequence_dict:
                    found += count
                    possible_sequence_dict[extract_random][i] += count
                else:
                    not_found += count


        # Counting statistics
        stats.append([biosamples[i], processed_count, found, not_found])
        read_counts.append(found)

        print('{:,} total reads'.format(processed_count))
        print('{:,} reads ({:.2%}) matched possible sequences'.format(found,
              found/processed_count))
        print('{:,} reads ({:.2%}) did not match possible sequences'.format(not_found,
              not_found/processed_count))


    # Write counting statistics to a .csv file
    headers = ['Biosample', 'Reads processed', 'Matched possible sequences', 'No match']
    stats_file_name = 'Count possible sequence stats.csv'
    stats_file_address = working_folder_name + '//Results//' + stats_file_name
    with open(stats_file_address, 'w', newline='') as stats_file:
        stats_writer = csv.writer(stats_file)
        stats_writer.writerow(headers)
        for biosample_stats in stats:
            stats_writer.writerow(biosample_stats)


    return possible_sequence_dict, read_counts