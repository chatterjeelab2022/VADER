# VADER
VADER codes require Bio Python 3.6. 

To run the VADER data processing pipeline:

1. Fill in the file Illumina processing parameters.xlsx with information corresponding to your sequencing sample. 
-In cell 16C, enter the "library name," a key phrase/brief series of characters that shows up in all of the FASTQ files corresponding to that library. 
-In cell 17C, enter the "before selection run," a key phrase/brief series of characters that shows up in all of the FASTQ files containing the input library.
-In row 18, enter the "after selection run(s)," a key phrase/brief series of characters that shows up in all of the FASTQ files containing the input library. If you execute replicate selections, you can enter them in cells 18C, 18D, 18E, and so on.
-Row 19 optionally allows you to assess a second kind of selection simultaneously, but can be left blank.
-In cell 20C, you can optionally give a title to the type of selection entered on line 18.
-In cell 21C, you can optionally give a title to the type of selection entered in row 19, if any.
-In cell 24C, enter the location address of the folder containing the FASTQ files.
-In cell 27C, if the library was generated through a Twist oligo pool, enter the location address of the .csv containing all library sequences. Otherwise, leave blank.
-In cell 29C, enter the sequence that you expect to find in the sequencing FASTQ, beginning after the Illumina sequencing primer, in all capital letters, using T and not U.
-In cell 31C, if the library sample was made with a diversity sequence prefacing the library to enhance the Illumina sequencing run quality, list the length of the diversity sequence.
-In cell 32C, optionally list the number of bases at the end of the read that you would like to exclude from data analysis.
-In cell 34C, list the number of base in cell 29C that corresponds to the first base of the tRNA you are sequencing.
-In cell 35C, list the number of base in cell 29C that corresponds to the last base of the tRNA.
-In rows 37 and 38, enter the bases corresponding to the beginning and end of all constant regions found within the sequencing sample, i.e., all regions that are not randomized parts of your library.
-In row 40, beginning at cell 40C, list the number of base in cell 29C that corresponds to all bases in the tRNA that were randomized in your library. 
-In rows 41 and 42, enter the bases that you expect to be paired with each other, i.e. all bases that fall within a stem region, e.g. the 5' member of a base pair would be entered in 41C and its corresponding 3' member would be entered in 42C.
-In cells 45C, 47C, 49C, and 50C, enter the Q-score filtering parameters of your choice. Q1, the general minimum Q score, denotes the lowest Q score that you permit in a sequencing sample. Q2, the randomized bases minimum Q score, denotes the lowest Q score that you permit in a region of the tRNA that you've randomized in your library. QF denotes the number of Q1 failures allowed per sequence. Q0, the absolute Q floor, denotes how low a Q score associated with a permitted Q1 failure sequence can be while still being quantified for analysis.
-In cell 59C, optionally enter how many counts a library member must have to show up in the data analysis.
-In cell 60C, enter the format in which you want library members to be displayed in analysis. The format must have the same number of Ns as you have randomized bases listed on row 40.
-In cell 61C, optionally enter a second format in which you want library members to be displayed, e.g. if you have only randomized part of a stem but you would like to display the randomized members in the context of the full stem.

2. After filling in the spreadsheet, save a copy of the spreadsheet as a comma delimited .csv file.

3. In Main.py, on line 42, enter the address of the Illumina processing parameteres .csv file.

4. Within Main.py, to assess the quality of the sequencing run by identifying how many sequences in your FASTQ files pass the Q-score filtering criteria defined in the Illumina processing parameters file, set lines 48 and 55 to true, and on line 56, set the limit to 1000.

5. To execute data analysis, set lines 48, 65, 66, 73, 81, and 82 to true, and set line 74 to false. 
