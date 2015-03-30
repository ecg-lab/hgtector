# hgtector
Modified HGTector in Python

The following Python packages are required:
Biopython
Pandas
RPY2
numpy

R and blast+ are also required for proper functioning.

The following files/data will are required:
TAXDMP taxonomic information from NCBI
NR BLAST database from NCBI
MultispeciesAutonomousProtein2taxname file from RefSeq

These files can be specified in the Python script as constants.

To use script, place genome fasta files in ./input/ directory, then run script.
You can change parameters by creating a config.txt file in current direcory.

Parameters to adjust:
interactive: 1 (0 to turn off interactive mode)
nThreads: 1 (multithread BLAST)
nHits: 500 (specify number of BLAST hits to keep)
evalue: 1e-5 (specify evalue cutoff)
selfGroup: None (specify TaxIDs of self-group)
closeGroup: None (specify TaxIDs of close group)
modKCO: 1 (0 to turn off conservative cutoffs)
selfLow: 0 (1 to use self-group cutoffs)
