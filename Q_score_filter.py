# -*- coding: utf-8 -*-
"""
@author: Rachel Kelemen
@deputy wizard: Joshua Brown
Finalized April 2019

This module contains a function to read FASTQ files, filter reads using their Phred
quality scores, and write the passing sequences to a text file.
"""

import time
import csv
import datetime as dt
from Bio import SeqIO

def Q_score_filter(working_folder_name, to_filter, biosample, Q1, QF, Q0, Q2,
                   expected_sequence, randomized_bases, exclude_at_start, exclude_at_end,
                   stats_file_address):
    """Read FASTQ files, filter reads by quality scores, and write passing sequences
    to a text file.

    Parameters
    ---
    working_folder_name : str
        Folder in which FASTQ files can be found
    to_filter : list
        List of FASTQ files for the given biosample
    biosample : string
        Biosample name
    Q1 : int
        Minimum allowed score for the full sequence
    QF : int
        Number of bases below Q1 allowed per sequence
    Q0 : int
        Absolute minimum allowed score for bases allowed by QF
    Q2 : int
        Minimum allowed score for randomized bases
    expected_sequence : str
        Expected sequence for the full read
    randomized_bases : list
        Coordinates of each randomized base
    exclude_at_start : int
        Number of bases to exclude from analysis at the beginning of the sequence,
        typically the diversity sequence
    exclude_at_end : int
        Number of bases to exclude from analysis at the end of the sequence, often
        due to decreased quality scores
    stats_file_address : str
        File address to which filtering stats should be written

    Returns
    ---
    output_file_name : str
        File to which passing sequences have been written
    stats : list
        Statistics describing number of sequences passing or failing on each parameter

    Notes
    ---
    Passing sequences are written to working_folder_name\Biosample_Q_score_filtered.txt

    For each parameter, failing sequences are written to working_folder_name\
    Biosample_parameter_failed.txt
    """

    # Setup

    # Count sequences which fall into each category
    short_seqs = 0
    Q0_fail_count = 0
    Q1_fail_count = 0
    Q2_fail_count = 0
    pass_count = 0
    read_count = 0

    # Holders for sequences to be written to .txt files
    passed_sequence_holder = []
    failed_sequence_holder = []

    start_time = time.time()


    # File names and addresses
    output_file_name = str(dt.date.today()) + ' ' + biosample + ' Q score filtered.txt'
    output_file_address = working_folder_name + '//' + output_file_name
    failed_output_file_name = str(dt.date.today()) + ' ' + biosample + ' Q score failed.txt'
    failed_output_file_address = working_folder_name + '//' + failed_output_file_name

    # Open input and output files
    with open(output_file_address, 'a', newline='\n') as output_file,\
    open(failed_output_file_address, 'a', newline='\n') as failed_output_file:

        # Clear output files from any previous runs
        output_file.truncate(0)
        failed_output_file.truncate(0)

        # Iterate over fastq files using Bio.SeqIO.parse()
        for input_file_name in to_filter:
            input_file = working_folder_name + '//' + input_file_name
            print ('/nProcessing file:/n' + input_file_name)
            print ('Started filtering reads at',dt.datetime.now().strftime('%I:%M:%S'))
            for seq_record in SeqIO.parse(input_file,'fastq'):
                # Track how many sequences records have been read
                read_count += 1

                # Start out by assuming the sequence passes filters
                Q_pass = True

                # Extract the numeric list of Q scores
                Q_scores = seq_record.letter_annotations['phred_quality']

                # Reject any truncated sequences and add them to a holder to be written
                # to a text file
                if len(Q_scores) < len(expected_sequence)-exclude_at_end or\
                len(seq_record.seq) < len(expected_sequence)-exclude_at_end:
                    Q_pass = False
                    short_seqs += 1
                    failed_sequence_holder.append(str(seq_record.seq)+'\n')

                # Filter based on overall sequence quality
                # No more than QF bases below Q1
                # No bases below Q0
                if Q_pass:
                    base_fail_Q1_count = 0
                    for Q in Q_scores[exclude_at_start:len(expected_sequence)-exclude_at_end]:
                        if Q < Q0:
                            Q_pass = False
                            Q0_fail_count += 1
                            failed_sequence_holder.append(str(seq_record.seq)+'\n')
                            break
                        if Q < Q1:
                            base_fail_Q1_count += 1
                            if base_fail_Q1_count > QF:
                                Q_pass = False
                                Q1_fail_count += 1
                                failed_sequence_holder.append(str(seq_record.seq)+'\n')
                                break

                # Filter based on randomized sequence quality
                # No randomized bases below Q2
                if (Q_pass):
                    for rb in randomized_bases:
                        if Q_scores[rb] < Q2:
                            Q_pass = False
                            Q2_fail_count += 1
                            failed_sequence_holder.append(str(seq_record.seq)+'\n')
                            break

                # Write passing sequences to holder
                if Q_pass:
                    passed_sequence_holder.append(str(seq_record.seq)+'\n')
                    pass_count += 1

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



            # At the end of the loop, write the sequences in each holder to the
            # appropriate file
            for seq in passed_sequence_holder:
                output_file.write(seq)
            passed_sequence_holder.clear()
           # print('{:,} sequences passed, {:,} failed ({:.2% pass})'.format(pass_count,
           #      read_count-pass_count, pass_count/read_count))

            for seq in failed_sequence_holder:
                failed_output_file.write(seq)
            failed_sequence_holder.clear()

    end_time = time.time()

    print ('/nQ score filtering results for ' + biosample)
    print ('{:,} reads processed'.format(read_count))
    print('{:,} sequences had short Q_score lists ({:.2%})'.format(short_seqs,
          short_seqs/read_count))
    print('{:,} sequences failed Q0 = {} ({:.2%})'.format(Q0_fail_count, Q0,
          Q0_fail_count/read_count))
    print('{:,} sequences failed Q1 = {} ({:.2%})'.format(Q1_fail_count, Q1,
          Q1_fail_count/read_count))
    print('{:,} sequences failed Q2 = {} ({:.2%})'.format(Q2_fail_count, Q2,
          Q2_fail_count/read_count))
    print('{:,} sequences passed ({:.2%})'.format(pass_count, pass_count/read_count))
    print('/nFiltered sequences written to:/n' + output_file_name)
    print('/nTotal run time = {:.0f} seconds ({:.1f} minutes)'.format(end_time-start_time,
          (end_time-start_time)/60))


    # Write the results to a .csv file
    # Opening the file in 'a' (append) mode means that lines are added to the end of the
    # existing file.

    # The structure of this file (set up in Main) is:
    # ['Biosample', 'Pct pass', 'Total reads', 'Passing reads',
    # 'Short sequences', 'Q0 fail', 'Q1 fail', 'Q2 fail']

    # Compile filtering statistics into stats
    stats = [biosample, (pass_count/read_count), read_count, pass_count, short_seqs,
             Q0_fail_count, Q1_fail_count, Q2_fail_count]
    with open(stats_file_address, 'a', newline='') as stats_file:
        stats_writer = csv.writer(stats_file)
        stats_writer.writerow(stats)

    return output_file_name