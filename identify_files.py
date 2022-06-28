# -*- coding: utf-8 -*-
"""
@author: Rachel Kelemen
@deputy wizard: Joshua Brown
Finalized April 2019

Short scripts to identify the files associated with each biosample given a processing
stage and extension (e.g. xxx.fastq or xxx_Q_score_filtered.txt)

"""

import os

def find_fastq_files(working_folder_name, biosamples):
    """Finds all fastq files associated with each biosample in a list

    Parameters
    ---
    working_folder_name : str
        full path of the folder in which the files will be found
    biosamples : list
        list of biosamples for which to identify files

    Returns
    ---
    fastq_files : list
        [[biosample_1_file_1, biosample_1_file_2...],[biosample_2_file_1, biosample_2_file_2...]]

    Notes
    ---
    Set up to allow for one biosample having multiple associated .fastq files as is
    true for some NextSeq samples
    """

    # identify all files in working_folder_name and set up a list of .fastq files
    found_file_names = os.listdir(working_folder_name)
    fastq_files = []

    # search for all files which contain the name of a given biosample and end with .fastq
    print('/nIdentified the following files:')
    end_name = '.fastq'
    for biosample in biosamples:
        biosample_files = []
        print(biosample)
        for f in found_file_names:
            if (biosample in f and f[-len(end_name):] == end_name):
                print('>',f)
                biosample_files.append(f)
        fastq_files.append(biosample_files)

    return fastq_files


def find_text_files(working_folder_name, biosamples, end_name):
    """Finds all files with a specified ending associated with each biosample in a list

    Parameters
    ---
    working_folder_name : str
        full path of the folder in which the files will be found
    biosamples : list
        list of biosamples for which to identify files
    end_name : str
        file name ending and extension

    Returns
    ---
    found_files : list
        [biosample_1_file, biosample_2_file...]

    Notes
    ---
    One list entry per biosample
    """

    # identify all files in working_folder_name and set up a list of .fastq files
    found_file_names = os.listdir(working_folder_name)
    found_files = []

    # search for all files which contain the name of a given biosample and have the
    # specified ending
    print('Identified the following files:')
    for biosample in biosamples:
        for f in found_file_names:
            if (biosample in f and f[-len(end_name):] == end_name):
                print('>',biosample,'=',f)
                found_files.append(f)

    return found_files

if __name__ == "__main__":
    working_folder_name = r'C:\Users\Rachel Kelemen\Google Drive\Chatterjee lab\AAV\Illumina\2018_11_18 LJ NextSeq112118\Ac3'
    biosamples = ['Ac3-1-1_S5_L002_R1_001','Ac3-1-1']

    print('/nFind fastq files')
    fastq_files = find_fastq_files(working_folder_name, biosamples)

    print('/nFind text files')
    text_files = find_text_files(working_folder_name, biosamples, 'Q_score_filtered.txt')


