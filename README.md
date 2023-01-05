# reptools

A python package which calls V and J gene segments and extracts CDR3
sequences from parallel sequencing data (FASTQ files).

If provided with FASTA files giving V, J and C gene segments, and a CSV file
pinpointing the start/end of the CDR3 region within V/J segments (as defined by
the first base of the codon encoding the conserved cysteine residue at the start
of the CDR3/the first base of  the conserved phenylalanine or tryptophan at the
end of the CDR3), reptools will call genes and extract the CDR3 sequences.

Reptools uses BLASTN to call the C, V and J genes, then Swift-Waterman searches
(performed by SWIPE) to accurately pinpoint the start and
end of the CDR3.

It is also capable of trimming adapters (by calling bbduk), denoising, and
filtering by expected error rate, allowing data to be processed in a single call.
----
I developed the first versions of reptools before good tools for repertoire
calling were widely available, and have continued to use it for my own work
because (a) it facilitates the use of custom databases (in my case, largely
chicken TCR), and (b) in our preliminary testing, it outperformed IMGT/Blast
(significantly) and MiXCR (to a lesser extent) with syntehtic data.  I've yet to
find the time to explore this properly and publish our results, but hope to do so.
----
Reptools requires Python >=3.7, running under Linux (although reptools itself
should run under Windows, the SWIPE aligner is currently only available for linux,
so far as I am aware).

It requires bbduk for adapter trimming and BLAST+ (specifically, blastn and
makeblastdb) and SWIPE for calling genes/extracting CDR3.

These can either be placed in a location specified in the PATH, or their locations
can be given at runtime using the --alignerPaths flag.

They are currently available at:
https://jgi.doe.gov/data-and-tools/software-tools/bbtools/
https://sourceforge.net/projects/bbmap/

https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/

https://github.com/torognes/swipe
----
Reptools is most easily used as a commandline tool, which as of the most recent
version (0.17.0) takes an arbitrarily long set of input files and performs one or 
more of several functions on them: call (call gene segments and extract CDR3),
denoise (collapses rare sequences into very similar more common ones), EEfilter
(filters CDR3 sequences by expected error rate (calculated from the FASTQ qual
string), and VDJtools (outputs a summary in a VDJtools format CSV file).

A call to reptools with default options, performing all functions, might be:

reptools call denoise EEfilter VDJtools -i /path/to/data/*.fastq --adaptors [adaptor
or primer sequences here, seperated by spaces]

to do the same, but omitting adaptor trimming:
reptools call denoise EEfilter VDJtools -i /path/to/data/*.fastq --notrim
----
Full documentation has not yet been developed, but reptools --help details all
available options.
----
Databases are not yet included on the github: get in touch if you would like me to
send mouse or chicken TCR databases, or write up the (very simple) format.
