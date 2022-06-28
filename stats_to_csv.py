# -*- coding: utf-8 -*-
"""
@author: Rachel Kelemen
@deputy wizard: Joshua Brown
Finalized April 2019

A short script to write stats describing the outcome of Q score or mismatch filtering
to a csv file.
"""

import csv

def write_stats(working_folder_name, biosamples,
                  process, stats):
    """Writes stats describing the outcome of Q score or mismatch filtering to a csv file

    Parameters
    ---
    working_folder_name : str
        Working folder
    biosamples : list
        List of biosamples which will be the row headings in the output table
    process : str
        Which process to summarize. Accepts 'Q score filter' or 'Mismatch filter' only
    stats : list
        List of stats to write to the file

    No returns

    Notes
    ---
    Writes results to \Results\Process stats.csv
    """

    output_file_name = process + ' stats.csv'
    output_file_address = working_folder_name + '\\Results\\' + output_file_name

    if process == 'Q score filter':
        headers = ['Biosample', 'Total reads', 'Passing reads', 'Truncated sequences',
                   'Failed on Q1', 'Failed on Q2']

    elif process == 'Mismatch filter':
        headers = ['Biosample', 'Passing reads', 'Mismatch failed']


    # Compile table and write to output file
    with open(output_file_address, 'w', newline = '') as output_file:
        # Clear any previous entries in the output file
        output_file.truncate(0)

        row_writer = csv.writer(output_file)
        # Write the header row
        row_writer.writerow(headers)
        # Write the row for each biosample
        for biosample, stat_row in zip(biosamples, stats):
            row = []
            row.append(biosample)
            for stat in stat_row:
                row.append(stat)
            row_writer.writerow(row)