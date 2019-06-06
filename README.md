# BioinfoToolbox

Assignments for Algorithms, databases and tools in bioinformatics (NDBI044)

Following assignments from http://bioinformatika.mff.cuni.cz/repository/#/stories/detail?id=bioinformatics_toolbox

1) Processing FASTA files
Upon initialization of the class, a file in FASTA is loaded according to a given filename. The class provides functions for returning the description or sequence of a specific molecule, its sequence length and subsequence.

2) Measuring sequence similarity using Hamming distance
This class contains two functions. GetDist takes two user-given sequences and returns either an error message, if they are of different length, or calls the function Hamming to return their Hamming distance. The sequences can be strings or sequence objects.

3) Sequence alignment using edit distance
The class EditDist takes two sequences as strings or sequence objects in initializationn. It contains the function GetDist, which uses a matrix to calculate the edit distance of two given sequences and return it. The function GetAlignment returns one of the optimal alignments obtained by backtracking. An additional function CalcMin is used whithin GetDist.

4) Processing PDB files
The initialization function of PDBparser loads a .pdb file. The class contains functions, which return different entities whithin the model. See code for details. The function GetWidth calculates the maximum distance between any two atoms in the structure. GetAtmsInRadius and GetResInRadius return a list of atoms or residues, respectively, in a given redius from a ligand (HETATM).

5) Processing multiple sequence alignment
The initialization function of the MSA class loads a multiple sequence alignment from a file in clustal format. It contains functions to returs the number and length of sequences, a specific sequence or a specific column of the alignment. A scoring matrix can be loaded from a file in csv format (LoadMatrix). Then the sum of pairs score can be calculated for a specific column or for the whole alignment.

Determinaion of conservation of a position in a sequence is based on the Jensen-Shannon divergence, steps taken from: 
John A. Capra, Mona Singh; Predicting functionally important residues from sequence conservation, Bioinformatics, Volume 23, Issue 15, 1 August 2007, Pages 1875–1882, https://doi.org/10.1093/bioinformatics/btm270
