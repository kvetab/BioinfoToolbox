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

6) Conservation determination from multiple aligned sequences
Tha class JSD is initialized by an MSAparser (see assignment 5), which already contains a loaded multiple sequence alignment. 
Upon initialization, weights of sequences are calculated (CalcWeights) using methods described in Henikoff and Henikoff (1994). Next, distributions of residues in columns are determined (CalcDistr). For simplicity, a uniform background distribution is assumed. The function JSDivergence returns the conservation score of a given column of the multiple sequence alignment. Determinaion of conservation of a position in a sequence is based on the Jensen-Shannon divergence, steps taken from Capra and Singh (2007). Functions CalcColFreq and RelativeEntropy are used as intermediate steps in CalcWeights and JSDivergence, respectively.

7) Computing structure-related properties
The class StructAnalysis requires a model loaded by a PDBparser and a radius. It calculates the half sphere exposure for each residue. The function GetHSE returns the contact number and HSE beta up, which is supposed to be the best measure of exposure according to Hamelryck (2005). GetCNs then acquires these values for each residue and keeps counts of buried and exposed residues (based on arbitrary thresholds). The distribution can be plotted as a histogram by GetHistogram. PortionPolar returns either the proportion of amino acids in the core, that are polar, or the proportion of amino acids on the surface, that are polar, or the ratio of polar amino acids in the core vs. on the surface. The function GetDiam returns the greatest distance between atoms, as in assignment 4.

References:
Henikoff S, Henikoff J. Position-based sequence weights, J. Mol. Biol., 1994, vol. 243 (pg. 574-578)
John A. Capra, Mona Singh; Predicting functionally important residues from sequence conservation, Bioinformatics, Volume 23, Issue 15, 1 August 2007, Pages 1875–1882, https://doi.org/10.1093/bioinformatics/btm270
Hamelryck, Thomas. "An amino acid has two sides: a new 2D measure provides a different view of solvent exposure." Proteins: Structure, Function, and Bioinformatics 59.1 (2005): 38-48.
