# hgtector
Modified HGTector in Python

The original method can be found here:
Qiyun Zhu, Michael Kosoy, and Katharina Dittmar
HGTector: an automated method facilitating genome-wide
discovery of putative horizontal gene transfers
BMC Genomics 2014, 15: 717
doi: 10.1186/1471-2164-15-717

If you have any questions, please contact:
Timothy.J.Straub.GR@dartmouth.edu
olgazh@dartmouth.edu

The following Python packages are required:
* Biopython
* Pandas
* RPY2
* numpy

R and BLAST+ are also required.

The following files/data are required:
* TAXDMP taxonomic information from NCBI
	ftp://ftp.ncbi.nih.gov/pub/taxonomy
* NR BLAST database from NCBI
	ftp://ftp.ncbi.nlm.nih.gov/blast/db/
* MultispeciesAutonomousProtein2taxname file from RefSeq
	ftp://ftp.ncbi.nlm.nih.gov/refseq/release/release-catalog/

These files can be specified in the Python script as constants.

To use script, place genome fasta file(s) in ./input/ directory, then run script.
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


