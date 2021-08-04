# Seq_Sort
Sort FASTQ output from Illumina MiSeq into "good", "candidate", and "bad" reads based on the integrity of the adapters and the size of the library. Uses mismatches instead of edit distance.

# Requirements

Developed/tested in Ubuntu 20.04, compiled using gcc. Uses getline() so likely will not work if the operating system does not support it.

# Installation

Download the 'src' folder and move to desired installation directory. Change directory into the 'src' directory and then type 'make'. You should then see a binary executable called 'sort' in the installation directory.

# sort

./sort JOBNAME INPUTFILE OFFSET LIBRARYSIZE HEAD TAIL GOODERRORS MAXERRORS REVERSELIBRARY

* JOBNAME is the output folder name (and prefix for output files).
* INPUTFILE is the filepath to the input file (more information about this below).
* OFFSET is an integer that indicates the number of nucleotides at the start prior to the head adapter.
* HEAD is the sequence of the head adapter.
* TAIL is the sequence of the tail adapter.
* GOODERRORS is the number of cumulative errors allowed in the head and tail adapters to be considered a good read.
* MAXERRORS is the number of cumulative errors allowed in the head and tail adapters to be considered a candidate read.
* REVERSELIBRARY takes either 'true' or 'false' and will reverse the library sequence in the output if set to 'true'.

Sort can also be run without any parameters and it will be the near equivalent (no prefix on the output files) of:
./sort job data.txt 0 15 ACAC CACA 0 1 false

# Additional Instructions

The inputfile is a preprocessed version of the FASTQ files (essentially just the sequences on each line). You can simply write a script that saves every 4n+2 line of the FASTQ file. If necessary, I may upload my preprocessing script (that accounts for complementing and reversing).

# Runtime

In the average case, the runtime should be O(nlogn). The two major time consuming steps are sorting the sequences into hashmaps w/ singly linked lists (average O(nlogn)) and a mergesort (average/worst O(nlogn)). Worst case of hashmap w/ singly linked lists is O(n^2) since no re-hashing method was built in. However, it was incredibly unlikely to be of significant impact as the index uses a modulo value of about ~1 million (thus this uses fairly large spatial complexity) and the sequence to key method creates absolutely unique keys (for library size <= 19) and reasonably unique keys (for library size > 19). Moreover, the MiSeq's upper limit on sequences is ~25 million sequences or so (usually with many duplicates), thus we expect no more than ~25 sequences in any given linked list in the average case.
