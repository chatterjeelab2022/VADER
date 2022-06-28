# -*- coding: utf-8 -*-
"""
Created on Sun Dec 23 17:15:07 2018

@author: Rachel Kelemen
@deputy wizard: Joshua Brown
Finalized April 2019

This module contains a function to characterize the quality of a run based on the
Phred quality scores.
"""

import time
import os
import csv
import datetime as dt
import numpy as np
import matplotlib.pyplot as plt
from Bio import SeqIO
from statistics import mean

def characterize_quality(fastq_file_name, working_folder_name, limit, expected_sequence,
                         exclude_at_start, exclude_at_end, randomized_bases):
    """Characterize the quality of a Fastq data set based on Phred quality scores

    Parameters
    ---
    fastq_file_name : str
        File to be characterized
    working_folder_name : str
        Folder containing .fastq files
    limit : int
        How many sequences to characterize
    expected_sequence : str
        The full expected sequences, including excluded positions
    exclude_at_start : int
        Number of bases to exclude from analysis at the beginning of the sequence,
        typically the diversity sequence
    exclude_at_end : int
        Number of bases to exclude from analysis at the end of the sequence, often
        due to decreased quality scores
    randomized_bases : list
        List of randomized base coordinates

    Notes
    ---
    Output plots are placed working_folder_name\Plots\Run Quality

    The number of sequences to be characterized is limited because reading the whole
    file is very slow.
    """

    # Setup
    # Number of sequences read - doesn't include short sequences
    read_count = 0
    # Number of sequences which are shorter than expected
    short_Q_scores = 0

    # Possible Q scores
    Q_score_list = np.arange(41) #Illumina records scores 0 through 41

    # Lists of Q score frequencies [Q0_count, Q1_count...]
    # Lowest scores for each sequence
    min_Q_score_counts = list(np.zeros(len(Q_score_list),dtype=int))
    # Second lowest scores for each sequence
    min_Q_score_counts_2 = list(np.zeros(len(Q_score_list),dtype=int))
    # Average score for each sequence
    avg_Q_score_counts = list(np.zeros(len(Q_score_list),dtype=int))
    # Total occurrences of each Q score
    all_Q_score_counts = list(np.zeros(len(Q_score_list),dtype=int))

    # Average Q score at each position (this includes excluded bases)
    per_base_avg_scores = list(np.zeros(len(expected_sequence)))

    # Lists of Q score frequencies for randomized bases in each sequence
    if randomized_bases:
        # Lowest scores among randomized bases
        min_rand_Q_score_counts = list(np.zeros(len(Q_score_list),dtype=int))
        # Average scores among randomzied bases
        avg_rand_Q_score_counts = list(np.zeros(len(Q_score_list),dtype=int))
        # Total occurrences of each Q score among randomized bases
        rand_Q_score_counts = list(np.zeros(len(Q_score_list),dtype=int))

    start_time = time.time()


#   Iterate over .fastq file using SeqIO.parse()
    print('/n---/n'+fastq_file_name)
    print('Read up to {} sequences'.format(limit))
    print ('Started at',dt.datetime.now().strftime('%I:%M:%S'))
    fastq_file_address = working_folder_name + '//' + fastq_file_name
    for seq_record in SeqIO.parse(fastq_file_address,'fastq'):
        if read_count >= limit:
            break  # exit the loop if the limit has been reached

        # Extract Q scores from the seq_record
        full_Q_scores = seq_record.letter_annotations['phred_quality']
        # Omit excluded bases from the beginning and end of the sequence
        Q_scores = full_Q_scores[exclude_at_start:len(expected_sequence)-exclude_at_end+1]

        # Reject any truncated reads and keep a count of how many were rejected
        if len(full_Q_scores) < len(expected_sequence) - exclude_at_end:
            short_Q_scores += 1

        # Characterize quality for full-length sequences
        else:
            # Average Q scores per base across the whole sequence
            for i in range(len(expected_sequence)):
                per_base_avg_scores[i] += full_Q_scores[i]/limit

            # Randomized base quality
            if randomized_bases:
                # Extract randomized base Q scores
                randomized_base_Q_scores = []
                for coord in randomized_bases:
                    randomized_base_Q_scores.append(full_Q_scores[coord])
                # Minimum scores
                min_rand_Q_score = min(randomized_base_Q_scores)
                if min_rand_Q_score > 40:  # just in case Q scores > 40 are returned
                    min_rand_Q_score = 40
                min_rand_Q_score_counts[min_rand_Q_score] += 1
                # Average scores
                avg_rand_Q_score = int(mean(randomized_base_Q_scores))
                if avg_rand_Q_score > 40: #just in case Q scores > 40 are returned
                    avg_rand_Q_score = 40
                avg_rand_Q_score_counts[avg_rand_Q_score] += 1
                # Total Q score frequency
                for Q_score in randomized_base_Q_scores:
                    rand_Q_score_counts[Q_score] += 1

            # Sort the Q scores in ascending order
            sorted_Q_scores = sorted(Q_scores)
            # Minimum scores
            min_Q_score = sorted_Q_scores[0]
            if min_Q_score > 40:
                min_Q_score = 40
            min_Q_score_counts[min_Q_score] += 1
            # Second lowest scores
            min_Q_score_2 = sorted_Q_scores[1]
            if min_Q_score_2 > 40:
                min_Q_score_2 = 40
            min_Q_score_counts_2[min_Q_score_2] += 1
            # Average scores
            avg_Q_score = int(mean(Q_scores))
            if avg_Q_score > 40:
                avg_Q_score = 40
            avg_Q_score_counts[avg_Q_score] += 1
            # frequency of each Q score over all positions
            for Q_score in Q_scores:
                all_Q_score_counts[Q_score] += 1

            # Increment read count
            read_count += 1

    # Plot results
    # Make plotting results folder
    plot_folder = working_folder_name + '//Plots//Run Quality ' + str(limit)
    if os.path.exists(plot_folder) == False:
        os.makedirs(plot_folder)

    # Plots where x axis is Q scores and y axis is frequency
    # Minimum Q scores, 2nd lowest Q scores, Average Q scores, Q score frequencies,
    # Randomized base minimum Q scores, Randomzied base average Q scores,
    # Randomized base Q score frequencies
    plot_by_Q_scores = [min_Q_score_counts, min_Q_score_counts_2, avg_Q_score_counts, all_Q_score_counts]
    if randomized_bases:
        for a in [min_rand_Q_score_counts, avg_rand_Q_score_counts, rand_Q_score_counts]:
            plot_by_Q_scores.append(a)
    stat_names = ['Lowest Q scores per sequence', '2nd lowest Q scores per sequence',
                  'Average Q scores per sequence', 'Q score frequencies']
    if randomized_bases:
        for a in ['Minimum randomized base Q scores', 'Average randomized base Q scores',
                  'Q score frequences for randomized bases']:
            stat_names.append(a)
    # Plot colors: red for minimum scores, yellow for 2nd lowest, green for average,
    # blue for total frequencies
    colors = ['r', 'y', 'g', 'b', 'r', 'g', 'b']

    # Make plots
    for stat, stat_name, Color in zip(plot_by_Q_scores, stat_names, colors):
        fig, ax = plt.subplots(figsize=(4,3))
        ax.bar(Q_score_list, stat, edgecolor='black', color=Color)
        ax.set_title(fastq_file_name + '\n' + stat_name)
        ax.set_ylabel('Count')
        ax.set_xlabel('Q score')
        ax.minorticks_on()
        plt.tight_layout()
        # Save figures
        plot_file_name = plot_folder + '//' + fastq_file_name[:-6] + ' ' + stat_name + '.png'
        plt.savefig(plot_file_name, dpi=300)
        plt.close()

    # Plot average Q score by position
    # X axis is sequence coordinate and y axis is average Q score
    fig, ax = plt.subplots(figsize=(12, 3))
    ax.bar(np.arange(1, len(expected_sequence)+1), per_base_avg_scores, color='g', edgecolor='w')
    ax.set_title(fastq_file_name + '/n' + 'Average Q scores by position')
    ax.set_ylabel('Q score')
    ax.set_xlabel('Position')
    ax.axvline(exclude_at_start+0.5, color='black', linestyle='--')
    ax.axvline(len(expected_sequence)-exclude_at_end+0.5, color='black', linestyle='--')
    ax.minorticks_on()
    plt.tight_layout()
    # Save figure
    plot_file_name = plot_folder + '//' + fastq_file_name[:-6] + ' Average Q scores by position.png'
    plt.savefig(plot_file_name, dpi=300)
    plt.close()


    end_time=time.time()


    print ('/n{:,} total reads processed in {:.0f} seconds'.format(
            short_Q_scores + read_count, end_time-start_time))
    print('/nResults:\n{:,} sequences had short Q score lists and were rejected'.format(
            short_Q_scores))
    print('{:,} sequences had full-length lists.'.format(read_count))