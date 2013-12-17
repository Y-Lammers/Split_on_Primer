Split\_on\_Primer.py
===================

To increase the efficiency of next generation sequences such as the Illumina HiSeq and the
IonTorrent multiple samples can be pooled into a single run on these machines. To separate
these samples from the mountain of reads produced small barcode sequences (such as MID
labels) are usually attached to the separate samples prior to sequencing. These labels can
be reused or skipped if the samples differ from each other based on the primers used for
sequencing. For example in DNA barcoding projects ITS samples can be mixed with CO1 samples.
The Split\_on\_Primer.py script is desigend to split these next generation reads based on
the primers used for sequencing. Both substitution (hammington) and fuzzy (levenshtein)
string searches are used to accurately determine the origin of a read.

Installation instructions
-------------------------

The script can be cloned from the github repository and be directly ran from the commandline.
The python code should work with both python versions 2.7 and 3.0 or higher.

General usage
-------------

The basic command to run the pipeline is:

`Split_on_Primer.py -f <sequence_file> -p <primer_file>`

This command will split the sequences in the sequences file based on the primers listed
in the primer file. The input sequences can be provided in either FASTA or FASTQ format.
By default the primer file is in .csv format, with column 1 containing the primer name and
column 2 the primer sequences, each primer should be placed on a new line. for example:

	ITS1,CCGTAGGTGAACCTGCGG
	ITS2,CATCGATGAAGAACGCAGC

The script will create an output sequence file for each primer provided in the primer file,
together with one 'unsorted' file. The output file formats are the same as in the input
sequence file format.

Advanced options
----------------

Different separators (such as tabs) can be provided for the primer file with the `-d`,
`--delimiter` argument (for tabs the argument 'tab' or '\t' needs to be provided).

The primer sequences can be trimmed from the reads sequences with the `--trim` argument.

By default the script will only look for exact matches between the reads and the primers.
(hamming distance is 0). Substitutions can be allowed with the `-m` or `--mis` argument.
The number of mismatches allowed is a trade-off between accuracy (a low number of mismatches
allowed) and quantity (a high number of mismatches, but increases the number of false-
positives). The mismatch parameter should depend on the length of the primers used (longer
primers might need higher mismatch settings since there is more room for substitutions)
and the general quality of the sequencing run.

Shifts can occur between the reads and the primers, either due to sequencing errors or due
to post-processing of the reads (trimming of barcode sequences). By default the pipeline
does not look for shifted sequences, so if a shift occurs the number of mismatches shoots
up and the read cannot be assigned to a primer. With the `-s` or `--shift` argument the
maximum number of nucleotides shifted between the primers and the reads can be specified.
A higher shift value will be able to assign more reads to the primers, but might introduce
false-positives.

By default the script will utilize all cores available at the computer. The number of used
cores can be altered with the `--cores` argument. The script memory usage can be influenced
with the `--chunk` argument. The chunk size indicated how many reads it will load into the
memory at any given time. By default the chunk size is set to 25000 * the number of cores,
which, depending on the read format and length, results in 50mb ram usage per core.

Full command information
------------------------

Command line arguments:

	Split\_on\_Primer.py [-h] [-f Sequence file] [-p Primer file]
		[-m Mismatches allowed] [-s Nucleotide shift allowed]
		[-t] [-d CSV delimiter] [-c Number of Cores]
		[--chunk Chunk size]

The following list can be generated with the `-h` flag. All commands can be
used in either short or long form.

	-h, --help
		show this help message and exit
	-f <Sequence file>, --sequence_file <Sequence file>
		The sequence file in either fastq or fasta format.
	-p <Primer file>, --primers <Primer file>
		Separated value file containing the primers. Format =
		primer_name,primer_sequence
	-m <Mismatches allowed>, --mis <Mismatches allowed>
		The maximum number of mismatches allowed between the primer
		and reads (default = 0)
	-s <Nucleotide shift allowed>, --shift <Nucleotide shift allowed>
		The maximum sequence shift allowed between the primer and reads.
		(default = 0)
	-t, --trim
		Trim the primers from the sequences after splitting.
	-d <delimiter>, --delimiter <delimiter>
		CSV delimiter used in the primer file (default = ',')
	-c <Number of Cores>, --cores <Number of Cores>
		The number of CPU cores the script will use (default = max
		number of CPUs available)
	--chunk <Chunk size>
		The maximum number of reads that will be loaded into the memory.
		A higher value will be faster but will take up more RAM space.
		(default = 25.000 * number of CPUs)

