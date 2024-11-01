# Blanche
Blanche is a simple piece of code designed to perform sequence-based cartography for virus sequence data describing within-host patterns of evolution.  It applies a simple method of dimension reduction to genome sequences, and is likely suitable for other applications in which there is a need to create a visual display of differences between genome sequences. 

The code compiles with a make command.  The Makefile likely needs to be edited for your specific system.  In particular the code makes use of the GSL library.  If you have a Mac, GSL can easily be installed using the Homebrew software package (see website...) 

**Inputs:**

Input is taken in the form of a .fasta file, identified using the --input flag.  Sequences in this file should have been pre-aligned using multiple sequence alignment.  There are multiple software packages available for alignment, including Muscle and Clustal Omega.  Manually checking your sequence alignment is recommended.

Sequences may contain ambiguous nucleotide.  In our code these are treated as stochastic objects, given an expected value based upon the composition of non-ambiguous nucleotides at the same position.  For example if the non-ambiguous nucleotides at a position are 3 x A and 5 x T, then an N at this position would be equivalent to 3/8 A and 5/8 T: The distance from a T to this nucleotide would be equal to 5/8.

Running the code:

The code can be run from command line with the commmand

./blanche <flags>

Flags:

--input <file> [Default: Sets.in] Specifies the file from which to read in aligned sequences.

--dim <int> [Default: 2] Specifies the number of dimensions to use in the inference.

--verb <flag> [Default: 0] Flag for more verbose output.

--fix <flag> [Default: 0] Incorporates ambiguous nucleotides into the distance measurement.  If left to zero, only differences in unambiguous A, C, G, and T nucleotides are accounted for in the calculation of Hamming distances.


**Outputs:**

[By default]

Output_points.dat: Inferred points in Euclidean space

Distances.dat: Describes the distance between the Hamming distance matrix and the matrix of Euclidean distances between the inferred Euclidean points.

Subsets.dat: List of points which have identical sequences

[With --verb 1]

Variant_positions.dat: List of sites in the genome at which variants were found

Sequence_distances.dat: Matrix of Hamming distances between sequences


**Example: **

An example can be found in the Example directory.
