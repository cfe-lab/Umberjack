/*_____________________________________________________________________ */
// Finds the total codon substitutions
function computeScalingFactorB(rateMatrix, baseFreqs)
{
	B = 0;
	for (n1 = 0; n1 < Rows(rateMatrix); n1 = n1+1)
	{
		for (n2 = 0; n2 < Columns(rateMatrix); n2 = n2+1)
		{
			if (n2 != n1)
			{
				B = B + baseFreqs[n1]*baseFreqs[n2]*rateMatrix[n1][n2];
			}
		}
	}
	return B;
}

/*______________________________________________________________________*/
// Writes Branch Names and Branch Lengths to CSV
function writeBranchLenCSV(_tree&, _branchLenCsv)
{
	_flatTreeRep	  = Abs (_tree);  // numeric vector, where the i-th entry stores the post-order traversal index for the parent of the i-th node (also in post-order traversal);
	_branchNames = BranchName(_tree, -1);  //  All branch names in postorder
	_branchLengths = BranchLength(_tree, -1);  // All branch lengths in postorder
	assert(Columns(_branchNames) == Columns(_branchLengths), "Expect same number of branch lengths to branch names");
	
	fprintf (_branchLenCsv, CLEAR_FILE, KEEP_OPEN);  // KEEP_OPEN does not close the file after writing to it
	
	fprintf(_branchLenCsv, "ParentNodeName,ChildNodeName,BranchLen\n");  // Header
	for (iChild=0; iChild < Columns(_branchNames); iChild=iChild+1)
	{
		_iPar = _flatTreeRep[iChild];
		fprintf(_branchLenCsv, _branchNames[_iPar], ",", _branchNames[iChild], ",", _branchLengths[iChild], "\n");
	}
	fprintf (_branchLenCsv, CLOSE_FILE);  // close file
}



/*_____________________________________________________________________ */

NICETY_LEVEL = 3;  // how much cpu to use
VERBOSITY_LEVEL = 6;  // likelihood function verbosity, higher is more verbose

dummy = HYPHY_LIB_DIRECTORY + DIRECTORY_SEPARATOR + "TemplateBatchFiles" + DIRECTORY_SEPARATOR + "TemplateModels" + DIRECTORY_SEPARATOR + "chooseGeneticCode.def";
ExecuteCommands ("#include \""+dummy+"\";");


SetDialogPrompt ("Please specify a codon data file:");
DataSet codon_ds = ReadDataFile (PROMPT_FOR_FILE);
DataSetFilter codon_dsf = CreateFilter (codon_ds,3,"","",GeneticCodeExclusions);
fprintf(stdout, "Codon Data File=", LAST_FILE_PATH, "\n");

/* read in tree from file */
ACCEPT_BRANCH_LENGTHS 	= 1;
ACCEPT_ROOTED_TREES		= 1;
IGNORE_INTERNAL_NODE_LABELS = 0;
SetDialogPrompt ("Please select a file containing a tree with branch lengths: ");
fscanf(PROMPT_FOR_FILE, "String", tree_string);
fprintf(stdout, "Tree File=", LAST_FILE_PATH, "\n");


/* generate codon model (MG94customModel) + GTR */
#include "fit_codon_model.ibf";

Tree	codon_tree = tree_string;	// NB:  Model must be defined before the Tree




/* constrain branch lengths */
// original branch lengths in nucleotide subst/site.  We want branch lengths in syn subst/site, since that is what SLAC uses for scaling expected syn, nonsyn sites.
fprintf(stdout, "\nConstraining branch lengths", "\n");
global scalingB 		= computeScalingFactorB (MG94custom, vectorOfFrequencies);
ReplicateConstraint("this1.?.?:=this2.?.?__/scalingB", codon_tree, codon_tree); 

AUTO_PARALLELIZE_OPTIMIZE = 1;	/* attempt to use MPI */
LikelihoodFunction codon_lf = (codon_dsf, codon_tree);
Optimize (res, codon_lf);

LIKELIHOOD_FUNCTION_OUTPUT = 7;		/* save LF to file */
SetDialogPrompt ("Please specify a file to export likelihood function: ");
fprintf (PROMPT_FOR_FILE, CLEAR_FILE);
fprintf(LAST_FILE_PATH, codon_lf);


// Write out codon tree branch lengths
SetDialogPrompt ("Please specify a csv file to export codon tree branch lengths: ");
fprintf (PROMPT_FOR_FILE, CLEAR_FILE);
outCodonTreeBranchLenCsv = LAST_FILE_PATH;
fprintf(stdout, "\tWriting Branch Length to ", outCodonTreeBranchLenCsv);
writeBranchLenCSV("codon_tree", outCodonTreeBranchLenCsv);




