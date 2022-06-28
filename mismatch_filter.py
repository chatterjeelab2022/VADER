# -*- coding: utf-8 -*-
"""
@author: Rachel Kelemen
@deputy wizard: Joshua Brown
Finalized April 2019

This module contains a function which accomplishes two things:
1. Filters reads based on agreement with the expected constant regions of the sequence
2. Plots the frequency of each base at every position in the sequence, as well as the
frequency of mismatches relative to the original expected sequence
"""

import csv
import os
import time
import datetime as dt
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter

def mismatch_filter(working_folder_name, to_filter, biosample, constant_region_coords,
                    allowed_mismatch_counts, expected_sequence, stats_file_address,
                    tRNA_coords = None, randomized_bases = None):
    """Filters reads based on agreement with expected constant regions and plots base
    frequency at each position

    Parameters
    ---
    working_folder_name : str
        Working folder
    to_filter : str
        Filename to filter e.g. biosample_Q_score_filtered.txt
    biosample : str
        Biosample name
    constant_region_coords : list
        List of onstant region start and end coordinates e.g. [[0, 15], [20, 30]...]
    allowed_mismatch_counts : list
        Maximum allowed mismatches for each constant region e.g. [2, 5...]
    expected_sequence : str
        Full expected sequence, including excluded and randomized regions

    Optional parameters
    ---
    tRNA_coords : list
        First and last bases of the tRNA e.g. [10, 28]
    randomized_bases : list
        List of randomized positions e.g. [15, 17...]

    Returns
    ---
    output_file_name : str
        File to which passing sequences have been written
    stats : list
        Statistics describing number of sequences passing or failing on each constant region

    Notes
    ---
    Passing sequences are written to working_folder_name\Biosample_mismatch_filtered.txt

    For each constant region, failing sequences are written to working_folder_name\
    Biosample_constant_region_coords_failed.txt

    Plots are saved to working_folder_name\Plots\Mismatch filter plots
    """

    # Setup
    # Constant region numbers
    constant_region_numbers = list(np.arange(len(allowed_mismatch_counts), dtype=int))
    # Partial expected sequences
    partial_expected_sequences = []
    for i in constant_region_numbers:
        partial_expected_sequences.append(expected_sequence[constant_region_coords[i][0]:constant_region_coords[i][1]])

    # Count sequences
    read_count = 0
    pass_count = 0
    mismatch_fail_counts = list(np.zeros(len(allowed_mismatch_counts), dtype=int))

    # Holders for sequences to be written to txt files
    passed_sequence_holder = []
    failed_sequence_holder = []

    # Input file
    input_file_address = working_folder_name + '//' + to_filter

    # Output file names
    output_file_name = str(dt.date.today()) + ' ' + biosample + ' mismatch filtered.txt'
    output_file_address = working_folder_name + '//' + output_file_name
    failed_output_file_name = str(dt.date.today()) + ' ' + biosample + ' mismatch failed.txt'
    failed_output_file_address = working_folder_name + '//' + failed_output_file_name

    # Plotting lists
    base_indices = {'A':0, 'C':1, 'G':2, 'T':3, 'N':4}
    # Base counts - will be converted to fractions later
    # [[1A, 2A, 3A...], [1C, 2C, 3C...], [1G, 2G, 3G], [1T, 2T, 3T]]
    base_count_list = []
    for base in ['A', 'C', 'G', 'T', 'N']:
        base_count_list.append(list(np.zeros(len(expected_sequence), dtype=int)))

    # Mismatch counts - will also be converted to fractions later
    mismatch_count_list = list(np.zeros(len(expected_sequence), dtype=int))


    print('Started at', dt.datetime.now().strftime('%I:%M:%S'))
    start_time = time.time()


    # Open input and output files
    with open(input_file_address, 'r', newline='\n') as input_file, \
    open(output_file_address, 'a', newline='\n') as output_file, \
    open(failed_output_file_address, 'a', newline='\n') as failed_output_file:

        # Clear output files from any previous runs
        output_file.truncate(0)
        failed_output_file.truncate(0)

        # Read each line in input file
        for current_sequence in input_file:
            read_count += 1
            mismatch_pass = True

            # Filter
            # Iterate through constant regions
            for coords, max_mismatches, partial_expected_sequence, i in zip(
                    constant_region_coords, allowed_mismatch_counts,
                    partial_expected_sequences, constant_region_numbers):
                mismatches = 0
                # Try a string comparison first = this is true only if the two sequences
                # are completely identical.
                if current_sequence[coords[0]:coords[1]] != partial_expected_sequence:
                    # If the strings are not completely identical, compare each base
                    # individually.
                    for c, e in zip(current_sequence[coords[0]:coords[1]],
                                    partial_expected_sequence):
                        if c != e:
                            mismatches += 1
                            if mismatches > max_mismatches:
                                mismatch_pass = False
                                mismatch_fail_counts[i] += 1
                                failed_sequence_holder.append(current_sequence)
                                # Stop testing individual positions
                                break
                # If failed, break filtering loop instead of iterating through the
                # remaining regions
                if mismatch_pass == False:
                    break

            # If sequence passed on all regions, add to holder and extract plotting information
            # It's more efficient to do the filtering and plotting at the same time
            if mismatch_pass:
                pass_count += 1
                # Add to holder
                passed_sequence_holder.append(current_sequence)
                # Update base counts and mismatch lists
                # Convert bases to list indices
                for i in range(len(expected_sequence)):
                    # Base counts
                    base_index = base_indices[current_sequence[i]]
                    base_count_list[base_index][i] += 1
                    # Mismatch counts
                    if current_sequence[i] != expected_sequence[i]:
                        mismatch_count_list[i] += 1

            # Whenever any holder is full, write the sequences in the holder to the
            # appropriate file and clear the holder.
            # This is done to avoid generating enormous lists which exceed the
            # computer's available memory.
            holder_length = 100000

            if len(passed_sequence_holder) >= holder_length:
                for seq in passed_sequence_holder:
                    output_file.write(seq)
                passed_sequence_holder.clear()
                print('{:,} sequences passed, {:,} failed'.format(pass_count,
                      read_count-pass_count))

            elif len(failed_sequence_holder) >= holder_length:
                for seq in failed_sequence_holder:
                    failed_output_file.write(seq)
                failed_sequence_holder.clear()


        # Write all sequences left in holders at the end of the loop
        # Passing sequence holder
        for seq in passed_sequence_holder:
            output_file.write(seq)
        passed_sequence_holder.clear()
        print('{:,} sequences passed, {:,} failed'.format(pass_count,
              read_count-pass_count))

        for seq in failed_sequence_holder:
            failed_output_file.write(seq)
        failed_sequence_holder.clear()


    end_time = time.time()
    print('\n{:,} reads processed in {:.0f} seconds'.format(read_count, end_time-start_time))
    print('{:,} sequences passed ({:.2%})'.format(pass_count, pass_count/read_count))
    print('Written to ', output_file_name)


    # Write the results to a .csv file
    # Opening the file in 'a' (append) mode means that lines are added to the end of the
    # existing file.

    # The structure of this file (set up in Main) is:
    # ['Biosample', 'Pct pass', 'Total reads', 'Passing reads', 'Region 1 failed'...]

    # Compile filtering statistics into stats
    stats = [biosample, (pass_count/read_count), read_count, pass_count]
    for fail_count in mismatch_fail_counts:
        stats.append(fail_count)
    # Write to file as a new line
    with open(stats_file_address, 'a', newline='') as stats_file:
        stats_writer = csv.writer(stats_file)
        stats_writer.writerow(stats)


    # Plotting

    # Make plotting results folder
    plot_folder = working_folder_name + '//Plots//Mismatch filtering'
    if os.path.exists(plot_folder) == False:
        os.makedirs(plot_folder)

    # Base fractions
    # [[1A, 2A, 3A...], [1C, 2C, 3C...], [1G, 2G, 3G], [1T, 2T, 3T]]
    base_fraction_list = []
    for base in ['A', 'C', 'G', 'T']:
        base_fraction_list.append([])
    # Mismatch fractions
    mismatch_fraction_list = []

    # Convert all count lists to fraction lists
    for N_count_list, i in zip(base_count_list, [0, 1, 2, 3]):
        for count in N_count_list:
            base_fraction_list[i].append(count/pass_count)
    # Mismatch fractions
    for count in mismatch_count_list:
        mismatch_fraction_list.append(count/pass_count)

    # All plots will have tick labels of the format 'Base', or 'coordinate\nBase' every 5.
    x_labels = []
    x_label_colors = []
    base_colors = {'A':'g', 'C':'b', 'G':'0', 'T':'r'}
    for base, index in zip(expected_sequence, np.arange(1, len(expected_sequence)+1)):
        x_label_colors.append(base_colors[base])
        if (index%5 == 0):
            x_labels.append('{}\n{}'.format(base, index))
        else:
            x_labels.append(base)

    # Make mismatch plot
    fig, ax = plt.subplots(figsize=(20, 4))
    ax.set_title(biosample)
    ax.set_ylabel('Mismatch frequency')
    # Shading should be done before plotting the series because otherwise data points
    # within the shaded boxes will also be shaded.
    # Shade tRNA and randomized sequences
    if tRNA_coords:
        # Shade the promoter region green
        ax.axvspan(-0.5, tRNA_coords[0]-0.5, facecolor='xkcd:aqua', alpha=0.15)
        # Shade the tRNA yellow
        ax.axvspan(tRNA_coords[0]-0.5, tRNA_coords[1]+0.5, facecolor='xkcd:yellow',
                   alpha=0.3)
        # Shade the terminator TTTTTT red
        ax.axvspan(tRNA_coords[1]+0.5, tRNA_coords[1]+6.5, facecolor='r', alpha=0.15)
    # Shade randomized bases dull red
    if randomized_bases:
        for coord in randomized_bases:
            ax.axvspan(coord-0.5, coord+0.5, facecolor='xkcd:orange', alpha=0.5)
    # Plot the series
    ax.plot(np.arange(len(mismatch_fraction_list)), mismatch_fraction_list, marker='o',
            markerfacecolor='r', linestyle='', markeredgewidth=0)
    ax.yaxis.set_major_formatter(FuncFormatter('{0:.0%}'.format))
    # Make horizontal gridlines
    plt.grid(axis='y', color='0.5')
    ax.set_xticks(np.arange(1, len(expected_sequence)+1))
    ax.set_xticklabels(x_labels)
    for xtick, color in zip(ax.get_xticklabels(), x_label_colors):
        xtick.set_color(color)
    plt.tight_layout()
    # Save figure
    plot_folder = working_folder_name + '//Plots//Mismatch filtering'
    plot_file_name = plot_folder + '//' + biosample + ' mismatch frequencies.png'
    plt.savefig(plot_file_name, dpi=300)
    plt.close()

    # Make base frequency plot
    fig, ax = plt.subplots(figsize=(20, 4))
    # Plot A, C, G, T
    ax.yaxis.set_major_formatter(FuncFormatter('{0:.0%}'.format))
    ax.set_title(biosample)
    ax.set_ylabel('Base frequency')
    # Shade tRNA and randomized sequences
    if tRNA_coords:
        # Shade the promoter region green
        ax.axvspan(-0.5, tRNA_coords[0]-0.5, facecolor='xkcd:aqua', alpha=0.15)
        # Shade the tRNA yellow
        ax.axvspan(tRNA_coords[0]-0.5, tRNA_coords[1]+0.5, facecolor='xkcd:yellow',
                   alpha=0.3)
        # Shade the terminator TTTTTT red
        ax.axvspan(tRNA_coords[1]+0.5, tRNA_coords[1]+6.5, facecolor='r', alpha=0.15)
    # Shade randomized bases dull red
    if randomized_bases:
        for coord in randomized_bases:
            ax.axvspan(coord-0.5, coord+0.5, facecolor='xkcd:orange', alpha=0.5)
    # Plot series
    for i, color in zip([0, 1, 2, 3], ['g', 'b', '0', 'r']):
        ax.plot(np.arange(len(base_fraction_list[i])), base_fraction_list[i], marker='o',
            markerfacecolor=color, linestyle='', markeredgewidth=0)
    # Make horizontal gridlines
    plt.grid(axis='y', color='0.5')
    ax.set_xticks(np.arange(1, len(expected_sequence)+1))
    ax.set_xticklabels(x_labels)
    for xtick, color in zip(ax.get_xticklabels(), x_label_colors):
        xtick.set_color(color)
    plt.tight_layout()
    # Save figure
    plot_folder = working_folder_name + '//Plots//Mismatch filtering'
    plot_file_name = plot_folder + '//' + biosample + ' base frequencies.png'
    plt.savefig(plot_file_name, dpi=300)
    plt.close()


    # Write base counts to a csv file
    # [['Original base', 'A', 'C', 'A'...], ['A', 0.1, 0.2, 0.95...]]
    csv_file_name = biosample + ' base counts.csv'
    csv_file_address = working_folder_name + '//Results//' + csv_file_name

    with open(csv_file_address, 'w', newline='') as output_file:
        output_writer = csv.writer(output_file)
        row_headers = ['A', 'C', 'G', 'T', 'N']
        # Write the indices and original sequence to the top rows
        row = ['Position']
        for i in np.arange(1, len(expected_sequence)+1):
            row.append(i)
        output_writer.writerow(row)
        row = ['Original base']
        for base in expected_sequence:
            row.append(base)
        output_writer.writerow(row)
        # Write base frequencies
        for header, base_fractions in zip(row_headers, base_fraction_list):
            row = [header]
            for b in base_fractions:
                row.append(b)
            output_writer.writerow(row)


    return output_file_name