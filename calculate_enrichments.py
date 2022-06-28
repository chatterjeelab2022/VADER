# -*- coding: utf-8 -*-
"""
@author: Rachel Kelemen
@deputy wizard: Joshua Brown
Finalized April 2019

This module contains a function which calculates the enrichment of each sequence in each
selection and writes the results to a csv file.
"""

import csv
import numpy as np
from statistics import mean, stdev

def calculate_enrichments(sequence_count_dict, read_counts, biosamples, sel_1_name,
                          sel_2_name, sel_2_runs, output_format, full_sequence_format,
                          rand_expected_pairs, working_folder_name, lib_name, min_lib_count):
    """Calculates the enrichment of each sequence in eaach selection relative to the
    library and writes the results to results.csv.

    Parameters
    ---
    sequence_count_dict : dict
    biosample : list
    sel_1_name : str
    sel_2_name : str
    sel_2_runs : list
    output_format : str
    full_sequence_format : str
    rand_expected_pairs : list
    working_folder_name : str
    lib_name : str

    Returns
    ---
    results_table : list
        For each sequence the entry looks like
        [[randomized_sequence], [full_stem_sequence], []]

    Notes
    ---
    results_table structure
    COLUMN   HEADER                     CONTENTS
    0        Randomized sequence        [formatted randomized sequence]
    1        Full stem sequence         [formatted full stem sequence]
    2        Paired bases               [number of paired bases]
    3        Raw counts                 [lib_count, sel_1_count, sel_2_count...]
    4        Fraction of total          [lib_fraction, sel_1_fraction, sel_2_fraction...]
    5        Fold enrichment            [sel_1_fold_enrichment, sel_2_fold_enrichment...]
    6        Enrichment factors         [sel_1_enrichment_factor, sel_2_enrichment_factor...]
    7        Enrichment ranks           [sel_1_rank, sel_2_rank...]
    8        Average enrichment factor  [average enrichment factor]
    9        StDev                      [enrichment factor standard deviation]
    10       Rank                       [rank based on average enrichment factor]

    Results table is written to working_folder\results\lib_name results.csv
    """

    # Setup
    # Set up column indices by name because it's clearer than referring to them by number.
    (Randomized_sequence, Full_stem_sequence, Paired_bases,
    Raw_counts, Fraction_of_total,
    Fold_enrichment, Enrichment_factor, Enrichment_rank,
    Avg_enrichment_factor, Stdev, Avg_rank) = (0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10)

    # Column headers
    header_list = ['Randomized sequence', 'Full stem sequence', 'Paired bases',
                   'Raw counts','Fraction of total',
                   'Fold enrichment', 'Enrichment factors', 'Enrichment ranks',
                   'Average enrichment factor', 'StDev', 'Average rank']

    # For pairing analysis, which pairs are considered valid?
    # Typically GU pairs are "valid" for tRNAs
    # A future goal is to come up with some metric for pairing strength
    # Pi stacking between adjacent bases?
    valid_pairs = {'A':['U'],
                   'C':['G'],
                   'G':['C','U'],
                   'U':['A','G']}

    # Sequences will be added to one of two tables, sequence_table or
    # low_abundance_table.

    # If there are any sequences which are present in selection runs but not present in
    # the library, this will cause divide by 0 errors later in the analysis.
    # Therefore, only certain parts of the analysis - calculating paired bases and % of
    # total reads - can be performed for these sequences.

    # Additionally, sequences which are only observed a very small number of times in the
    # library could lead to inaccurate enrichment factors - if a sequence was seen twice
    # instead of once in the library, for instance, the enrichment factor would be off
    # by a factor of two.

    # sequence_table holds sequences which were found >= min_lib_count times in the
    # library biosample.
    sequence_table = []
    # Enrichment factors will be calculated for these sequences.

    # low_abundance_table holds sequences which were identified < min_lib_count times
    # in the library biosample.
    low_abundance_table = []
    # Enrichment factors will still be calculated for sequences which were identified
    # > 0 times in the library biosample, but these values should be used cautiously.
    # Track how many sequences were not observed at all
    not_observed = 0

    # Lastly, we'll need to find the highest fold enrichment value for each selection
    # before we can calculate the enrichment factors.
    max_fold_enrichments = np.zeros(len(biosamples))


    # Read sequence counts from sequence_count_dict.
    for sequence in sequence_count_dict:
        # Set up table row sequence_data
        sequence_data = []
        for header in header_list:
            sequence_data.append([])

        # 0. Randomized sequence
        # Format sequences for output
        # Replace T's with U's, but keep the original as key_sequence for sequence_count_dict
        key_sequence = sequence
        sequence = sequence.replace('T', 'U')
        if output_format:
            # output_format is represented as 'NNN/NNN' where each N is to be replaced
            # with a base from sequence and all other characters are preserved.
            # Iterate through characters in output_format to generate formatted_sequence.
            formatted_sequence = ''
            # i represents the current position in sequence, starting at the beginning.
            i = 0
            for c in output_format:
                if c == 'N':
                    formatted_sequence += sequence[i]
                    i += 1
                else:
                    formatted_sequence += c
        else:
            formatted_sequence = sequence
        sequence_data[Randomized_sequence].append(formatted_sequence)

        # 1. Full stem sequence
        # Allows visualization of an entire stem, even if only part of the stem was
        # randomized in the lirbary, and facilitates comparison between libraries in the
        # same region.
        # Similar to how formatted_sequence is generated
        if full_sequence_format:
            full_sequence = ''
            # i represents the current position in sequence, starting at the beginning.
            i = 0
            for c in full_sequence_format:
                if c == 'N':
                    full_sequence += sequence[i]
                    i += 1
                else:
                    full_sequence += c
        else:
            full_sequence = formatted_sequence
        sequence_data[Full_stem_sequence].append(full_sequence)

        # 2. Paired bases
        # Count how many bases in the sequence are paired, based on valid_pairs
        # This uses sequence, not formatted_sequence
        # Start with 0 pairs
        paired = 0
        # Iterate through expected pairs using rand_expected_pairs coordinates, which only
        # include the randomized bases (unlike expected_pairs which provides pair
        # coordinates for the full sequence).
        for pair in rand_expected_pairs:
            # Extract coordinates and check for pairs
            a, b = pair[0], pair[1]
            if sequence[a] in valid_pairs[sequence[b]]:
                paired += 1
        sequence_data[Paired_bases].append(paired)

        # 3. Raw counts and
        # 4. Fraction of total
        for raw_count, total in zip(sequence_count_dict[key_sequence], read_counts):
            sequence_data[Raw_counts].append(raw_count)
            sequence_data[Fraction_of_total].append(raw_count/total)

        # Sequences found in library biosample
        # Calculate the change in abundance before and after each selection,
        # relative to the library.
        if sequence_data[Raw_counts][0] > 0:
            # 6. Fold enrichments
            library_abundance = sequence_data[Fraction_of_total][0]
            for selection_abundance, i in zip(sequence_data[Fraction_of_total],
                                              np.arange(len(max_fold_enrichments), dtype=int)):
                fold_enrichment = selection_abundance / library_abundance
                sequence_data[Fold_enrichment].append(fold_enrichment)
                # Update max_fold_enrichments unless the value comes from a low abundance
                # sequence.
                if sequence_data[Raw_counts][0] >= min_lib_count:
                    if fold_enrichment > max_fold_enrichments[i]:
                        max_fold_enrichments[i] = fold_enrichment

            if sequence_data[Raw_counts][0] >= min_lib_count:
                sequence_table.append(sequence_data)
            else:
                low_abundance_table.append(sequence_data)
                not_observed += 1

        # Sequences not found in library biosample
        # Enrichment cannot be calculated because the library abundance is 0.
        # Calculate the value if the library count was 1
        else:
            # What would library_abundance be if the library count was 1?
            library_abundance = 1 / read_counts[0]
            for selection_abundance, i in zip(sequence_data[Fraction_of_total],
                                              np.arange(len(max_fold_enrichments))):
                fold_enrichment = selection_abundance / library_abundance
                sequence_data[Fold_enrichment].append(fold_enrichment)
                # Do not update fold enrichments.
            low_abundance_table.append(sequence_data)


    # 6. Enrichment factors,
    # 8. Average enrichment factors, and
    # 9. StDev
    # Once finished iterating through all sequences and adding to tables, use
    # max_fold_enrichments to determine enrichment factor for each sequence
    # enrichment_factor = fold_enrichment / max_fold_enrichment
    for table in [sequence_table, low_abundance_table]:
        for sequence_data in table:
            for fold_enrichment, max_fold_enrichment in zip(sequence_data[Fold_enrichment],
                                                            max_fold_enrichments):
                enrichment_factor = fold_enrichment / max_fold_enrichment
                sequence_data[Enrichment_factor].append(enrichment_factor)
            # Calculate average enrichments
            # If there are distinctly identified sel_1 and sel_2 runs (e.g. 1x vs 2x TAG)
            # Calculate three means and stdevs: [sel_1_mean, sel_2_mean, all_sels_mean]
            if sel_2_runs:
                sel_1_enrichment_factors = sequence_data[Enrichment_factor][1:-len(sel_2_runs)]
                sel_2_enrichment_factors = sequence_data[Enrichment_factor][-len(sel_2_runs):]
                sequence_data[Avg_enrichment_factor].append(mean(sel_1_enrichment_factors))
                sequence_data[Avg_enrichment_factor].append(mean(sel_2_enrichment_factors))
                if len(sel_1_enrichment_factors) > 1:
                    sequence_data[Stdev].append(stdev(sel_1_enrichment_factors))
                else:
                    sequence_data[Stdev].append(0)
                if len(sel_2_enrichment_factors) > 1:
                    sequence_data[Stdev].append(stdev(sel_2_enrichment_factors))
                else:
                    sequence_data[Stdev].append(0)
            # Calculate the mean enrichment factor and stdev for all selections.
            all_enrichment_factors = sequence_data[Enrichment_factor][1:]
            sequence_data[Avg_enrichment_factor].append(mean(all_enrichment_factors))
            if len(all_enrichment_factors) > 1:
                sequence_data[Stdev].append(stdev(all_enrichment_factors))
            else:
                sequence_data[Stdev].append(0)

    # 7. Rank by enrichment factors and
    # 9. Rank by average enrichment factor(s)
    for table in [sequence_table, low_abundance_table]:
        # Individual selection enrichment factors
        for i in np.arange(len(biosamples)):
            # Sort sequence_table by enrichment factor i, largest values first
            table = sorted(table, key = lambda x: x[Enrichment_factor][i],
                                    reverse=True)
            # Add rank+1 to sequence_table (add 1 to start at 1 not 0)
            for rank in range(len(table)):
                table[rank][Enrichment_rank].append(rank+1)
        # Average enrichment factor(s)
        # How many enrichment factors to expect?
        if sel_2_runs:
            n = 3
        else:
            n = 1
        for i in np.arange(n):
            # Sort sequence_table by enrichment factor i, largest values first
            table = sorted(table, key = lambda x: x[Avg_enrichment_factor][i],
                                    reverse=True)
            # Add rank+1 to sequence_table (add 1 to start at 1 not 0)
            for rank in range(len(table)):
                table[rank][Avg_rank].append(rank+1)

    # Sort tables by average enrichment factor outside of the above loop
    sequence_table = sorted(sequence_table, key = lambda x: x[Avg_enrichment_factor][-1],
                            reverse=True)
    low_abundance_table = sorted(low_abundance_table, key = lambda x: x[Avg_enrichment_factor][-1],
                                 reverse=True)

    # Print results
    print('/nStats for the library run:')
    print('Of {:,} sequences, {:,} were observed >= {} times'.format(
            len(sequence_count_dict), len(sequence_table), min_lib_count))
    print('{:,} sequences were observed < {} times, with {} not observed'.format(
            len(low_abundance_table), min_lib_count, not_observed))

    print('/nMost enriched sequence: ', sequence_table[0][Randomized_sequence])
    print('Least enriched sequence: ', sequence_table[-1][Randomized_sequence])

    print('/nData for the top sequence:')
    for header, value in zip(header_list, sequence_table[0]):
        print(header, value)


    # Write the results to csv files.
    # Tables are already sorted by average enrichment factor.
    # Table structure reminder:
    '''
    COLUMN   HEADER                     CONTENTS
    0        Randomized sequence        [formatted randomized sequence]
    1        Full stem sequence         [formatted full stem sequence]
    2        Paired bases               [number of paired bases]
    3        Raw counts                 [lib_count, sel_1_count, sel_2_count...]
    4        Fraction of total          [lib_fraction, sel_1_fraction, sel_2_fraction...]
    5        Fold enrichment            [sel_1_fold_enrichment, sel_2_fold_enrichment...]
    6        Enrichment factors         [sel_1_enrichment_factor, sel_2_enrichment_factor...]
    7        Enrichment ranks           [sel_1_rank, sel_2_rank...]
    8        Average enrichment factor  [average enrichment factor]
    9        StDev                      [enrichment factor standard deviation]
    10       Rank                       [rank based on average enrichment factor]
    '''

    for table, file_name_end in zip([sequence_table, low_abundance_table], [' Results.csv',
                               ' Low abundance sequences.csv']):
        output_file_address = working_folder_name + '//Results//' + lib_name + file_name_end
        with open(output_file_address, 'w', newline='') as output_file:
            # Clear any previous data
            output_file.truncate(0)

            # Write
            output_writer = csv.writer(output_file)
            # Write headers and subheaders
            header_row = []
            subheader_row = []
            # Randomized sequence, Full stem sequence, and Paired bases each only
            # contain 1 value.
            for header in header_list[Randomized_sequence:Raw_counts]:
                header_row.append(header)
                subheader_row.append('')
            # Raw counts, Fraction of total, Fold enrichment, Enrichment factors, and
            # Enrichment ranks have one value for each biosample.
            for header in header_list[Raw_counts:Avg_enrichment_factor]:
                header_row.append(header)
                for l in range(len(biosamples)-1):
                    header_row.append('')
                for biosample in biosamples:
                    subheader_row.append(biosample)
            # Average enrichment factor and Rank 3 values if there are distinct sel_1 and
            # sel_2 conditions, [Sel_1, Sel_2, and All selections].
            if sel_2_runs:
                for header in header_list[Avg_enrichment_factor:]:
                    header_row.append(header)
                    header_row.append('')
                    header_row.append('')
                for selection_name in [sel_1_name, sel_2_name, 'All selections']:
                    subheader_row.append(selection_name)
            # Otherwise there is just one value for all selections.
            else:
                for header in header_list[Avg_enrichment_factor:]:
                    header_row.append(header)
                    subheader_row.append('All selections')
            output_writer.writerow(header_row)
            output_writer.writerow(subheader_row)

            # Write data
            for sequence_data in table:
                data_row = []
                for column in sequence_data:
                    for value in column:
                        data_row.append(value)
                output_writer.writerow(data_row)

        print('/nData written to')
        print(output_file_address)


    return sequence_table, low_abundance_table