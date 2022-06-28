# -*- coding: utf-8 -*-
"""
@author: Rachel Kelemen
@deputy_wizard: Joshua Brown
Finalized April 2019

This stand-alone script allows easy comparison of selection results from an Illumina
sequencing run ("current run") to a previous sequencing run or to a list of previously
tested hits and their activities.

Instructions
---
Enter run names, file names, and folder address below.
You can compare a current run to multiple past data files.
Plots and csv files will be written to the working folder.

Input csv files should have the following structure:
[[full_stem_sequence, value, stdev]...]
"""

import csv
import matplotlib.pyplot as plt

# Enter names and sequences
working_folder_name = r'C:\Users\Rachel Kelemen\Google Drive\Chatterjee lab\AAV\050919NextSeq\2019_05_10 Ac3 processing\Results'
current_data_name = 'Already tested'
current_data_file_name = 'Already tested for comparison.csv'
past_data_names = ['Ac3', 'Ac3del2']
past_data_file_names = ['2019_05_10 Ac3 Results for comparison.csv', '2019_05_10 Ac3 Results for comparison del 2.csv']


def compare(working_folder_name, current_data_name, current_data_file_name,
                     past_data_name, past_data_file_name):
    """
    Compares Illumina results from a selection to a previous sequencing run or a list of
    previously tested hits and activities.

    Parameters
    ---
    working_folder_name : str
        Working folder full address
    current_data_name : str
        Name of the current data set as you want it to appear on plots and output files
    current_data_file_name : str
        csv file containing the current run results
    past_data_name : str
        Name of the past data set as you want it to appear (e.g. Ac2 sequencing or Tested Activity)
    past_data_file_name : str
        csv file containing the past run results

    Returns
    ---
    None

    Notes
    ---
    Plots and results csv files are saved to the working folder.
    """

    # Read in current run results, and output to dictionary current_run
    # Numbers will be read in as strings and must be converted to floats
    # Using dictionaries instead of lists because searching for matches later will be
    # much more efficient.
    current_data_file_address = working_folder_name + '\\' + current_data_file_name
    current_data = {}
    total_current = 0
    # This will have the form {'sequence': [enrichment_factor, stdev]}
    with open(current_data_file_address) as current_data_file:
        current_reader = csv.reader(current_data_file)
        for row in range(0,1): # skip the header row
            next(current_reader)
        for row in current_reader:
            # Column 0 : sequence
            # Column 1 : average enrichment factor
            # Column 2 : standard deviation
            current_data[row[0]] = [float(row[1]), float(row[2])]
            total_current += 1
    print('{:,} sequences read from current data set {}'.format(total_current,
          current_data_name))

    # Do the same for past run
    past_data_file_address = working_folder_name + '\\' + past_data_file_name
    past_data = {}
    total_past = 0
    with open(past_data_file_address) as past_data_file:
        past_reader = csv.reader(past_data_file)
        for row in range(0,1): # skip the header row
            next(past_reader)
        for row in past_reader:
            past_data[row[0]] = [float(row[1]), float(row[2])]
            total_past += 1
    print('{:,} sequences read from past data set {}'.format(total_past, past_data_name))

    # Compare runs
    common_sequences = {}
    # This will have the format {'sequence':
    # [current_enrichment_factor, current_stdev, past_enrichment_factor, past_stdev]}
    # Track how many of the sequences in the current run also occur in the past run
    match_previous = 0
    not_match_previous = 0
    for sequence in current_data:
        if sequence in past_data:
            match_previous += 1
            current_enrichment_factor = current_data[sequence][0]
            current_stdev = current_data[sequence][1]
            past_enrichment_factor = past_data[sequence][0]
            past_stdev = past_data[sequence][1]
            common_sequences[sequence] = [current_enrichment_factor, current_stdev,
                            past_enrichment_factor, past_stdev]
        else:
            not_match_previous += 1
    print('{:,} sequences in the current data set matched the past data set ({:.1%})'.format(
            match_previous, match_previous/total_current))

    # Write results to csv file
    results_file_name = 'Compare ' + current_data_name + ' vs ' + past_data_name + '.csv'
    results_file_address = working_folder_name + '//' + results_file_name
    with open(results_file_address, 'w', newline='') as results_file:
        results_writer = csv.writer(results_file)
        # Write headers
        results_writer.writerow(['Sequence', current_data_name, '', past_data_name])
        # Write subheaders
        results_writer.writerow(['', 'Average', 'Stdev', 'Avg', 'Stdev'])
        for sequence in common_sequences:
            row = [sequence]
            for value in common_sequences[sequence]:
                row.append(value)
            results_writer.writerow(row)

    # Make x and y data lists for plotting
    x_data = [] # Current data
    y_data = [] # Past data
    for sequence in common_sequences:
        x_data.append(common_sequences[sequence][0])
        y_data.append(common_sequences[sequence][2])

    # Plot results
    # Ideally make elipses with error bars as dimensions?
    # For now scatter plot is way simpler.
    fig, ax = plt.subplots(figsize = (5,5))
    ax.scatter(x_data, y_data, c='black', alpha=1)
    ax.set_xlabel(current_data_name)
    ax.set_ylabel(past_data_name)
    plt.tight_layout()
    # Save figure
    plot_file_name = 'Compare ' + current_data_name + ' vs ' + past_data_name + '.png'
    plot_file_address = working_folder_name + '//' + plot_file_name
    plt.savefig(plot_file_address, dpi=300)
    plt.close()

for past_data_name, past_data_file_name in zip(past_data_names, past_data_file_names):
    compare(working_folder_name, current_data_name, current_data_file_name,
                     past_data_name, past_data_file_name)