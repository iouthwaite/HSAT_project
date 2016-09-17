# HSAT_project
HSATII Subfamily Loci Finder

to use:

run HSAT_hitFinder on your genomes of choice; they should be in FASTA format (.txt files).
Then, create a training dataset of what you think are loci from the hitFinder output - this should be a .txt file as well, formatted as a list of lists (ex: [[a, b, c],[d, e, f]] )
Then, run HSAT_locifinder to group HSAT hits into possible loci of HSATII subfamilies in your genome. See file for command line instructions.
Output of HSAT_locifinder is in the terminal, as well as an excel file.