# -*- coding: utf-8 -*-
"""
@author: Rachel Kelemen
@deputy wizard: Joshua Brown
Finalized April 2019

This module contains two functions which generate dictionaries where each possible
sequence for a given library is a key.

generate_randomized_sequence_dict is for site-saturation libraries. It generates all
possible sequences, assuming N at all randomized positions.

generate_twist_sequence_dict is for oligo pool libraries. It reads the randomized base
sequences from the sequence names in the Twist oligo pool .csv file.

Both functions return a dictionary of sequences and counts:
{sequence: [0, 0, 0...]}
"""

import numpy as np
import csv

def generate_randomized_sequence_dict(randomized_bases, biosamples):
    """Generate all possible sequences for site-saturation libraries with N at all
    positions.

    Parameters
    ---
    randomized_bases : list
        Randomized positions, used here to indicate how many positions are randomized.
    biosamples : list
        List of biosamples, used here to indicate how many placeholder zeros need to be
        placed in the count list.

    Returns
    ---
    possible_sequence_dict : dict
        Randomized base sequence: count list {sequence: [0, 0, 0...]}
    """

    # Setup
    bases = ['A', 'C', 'G', 'T']

    # Generate a list of possible sequences.
    # We start with a list of sequences of length 1.
    # Each round, we generate new sequences which are one base longer, with each possible
    # base at the new position.
    # These go in holder_2
    # At the end of each round we:
    #   clear holder_1,
    #   write the sequences in holder_2 to holder_1
    #   and then clear holder_2
    holder_1 = []
    holder_2 = []

    # Start with single-base sequences
    for base in bases:
        holder_1.append(base)

    # Loop through and add the next bases one at a time
    for l in range(1,len(randomized_bases)):
        for seq in holder_1:
            for base in bases:
                holder_2.append(seq + base)
        holder_1.clear()
        holder_1 = holder_2.copy()
        holder_2.clear()

    # Sort alphabetically
    holder_1 = sorted(holder_1)

    # Add complete sequences to possible_sequence_dict as keys
    # The values are a list of zeroes, one for each biosample. These will be updated
    # during count_randomized_sequences.
    possible_sequence_dict = {}
    for sequence in holder_1:
        possible_sequence_dict[sequence] = list(np.zeros(len(biosamples),dtype=int))

    print('/nGenerated {:,} possible sequences:'.format(len(holder_1)))
    print('[{}...{}]'.format(holder_1[0], holder_1[-1]))
    return possible_sequence_dict

def generate_twist_sequence_dict(biosamples, twist_sequence_file_address):
    """Generate all possible sequences for Twist libraries, based on the file used to
    create the oligo pool.

    Parameters
    ---
    twist_sequence_file_address : str
        File used to generate the oligo pool. The sequence names are the randomized
        base sequences.
    biosamples : list
        List of biosamples, used here to indicate how many placeholder zeros need to be
        placed in the count list.

    Returns
    ---
    possible_sequence_dict : dict
        Randomized base sequence: count list {sequence: [0, 0, 0...]}
    """

    # Read sequence names from the Twist oligo pool .csv file
    # The sequence names are a list of the randomzied bases for each sequence.
    possible_sequence_dict = {}
    with open (twist_sequence_file_address, 'r', newline='') as sequence_file:
        sequence_file_reader = csv.reader(sequence_file)
        for line in sequence_file_reader:
            possible_sequence_dict[line[0]] = list(np.zeros(len(biosamples), dtype=int))

    print('/nGenerated {:,} possible sequences'.format(len(possible_sequence_dict)))
    return possible_sequence_dict