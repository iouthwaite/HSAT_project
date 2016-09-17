# HSAT_project
HSATII Subfamily Loci Finder

to use:

run HSAT_hitFinder on your genomes of choice; they should be in FASTA format (.txt files).
Then, create a training dataset of what you think are loci from the hitFinder output - this should be a .txt file as well, formatted as a list of lists (ex: [[a, b, c],[d, e, f]] )
Then, run HSAT_locifinder to group HSAT hits into possible loci of HSATII subfamilies in your genome. See file for command line instructions.
Output of HSAT_locifinder is in the terminal, as well as an excel file.

HSATresults.txt is an example results file from hitFinder
HSATLoci.txt is an example curated training dataset for locifinder
Genome Lengths is a useful file with chromosome lengths (in base pairs)
HSATfileforplotproduction_... is an example output file from lociFinder
PlotGenerator is a script to help make plots using plotly. Requires a genome lengths file, as well as a plotly account