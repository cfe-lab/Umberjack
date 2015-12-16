

/*
Spits out a tab-separated file containing substitutions for each site at each branch

Columns:
Site, ParentNode, ChildNode,  Nonsynonymous Substitutions, Synonymous Substitutions, Expected nonsynonymous Substitutions, Expected synonymous substitutions

and/or 

Spits out a tab-separate file containing site substituions and site dn-ds.
Columns:
Site, Observed S Changes, Observed NS Changes,  E[S Sites], E[NS Sites], dN, dS, dN-dS, Scaled dN-dS


Will ask for input in this order:
- TSV file to write per-site-branch substitutions to
- TSV file to write per-site dN-dS to
- Genetic Code
- Codon Fasta File
- Tree File
- File to write nucleotide model fit to
- File to write codon model fit to
- CSV File to write codon branch lengths to
- Fasta file to write reconstructed ancestor and leaf sequences


Maintains the same input tree root, input tree node names, input tree branch lengths.



*/

RESOLVE_AMBIGUITIES = 0;  // Do not resolve ambiguous codons when counting substitutions

#include "matrix_site_branch_subst.bf";

SetDialogPrompt ("Please specify Site-Branch Substitution TSV File or -1 to Skip:");
fscanf(stdin, "String", outPerSiteBrSubFile);
if (outPerSiteBrSubFile != "-1")
{
	fprintf(outPerSiteBrSubFile, CLEAR_FILE);
}
fprintf (stdout, "\t[WRITING Site-Branch Substitutions TO ",outPerSiteBrSubFile,"]\n");	


SetDialogPrompt ("Please specify Site dN-dS TSV File or -1 to Skip:");
fscanf(stdin, "String", outPerSiteDnDsFile);
if (outPerSiteDnDsFile != "-1")
{
	fprintf(outPerSiteDnDsFile, CLEAR_FILE);
}
fprintf (stdout, "\t[WRITING Site dN-dS TO ",outPerSiteDnDsFile,"]\n");	

printSubTSV(outPerSiteBrSubFile, outPerSiteDnDsFile)
