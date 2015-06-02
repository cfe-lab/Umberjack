

/*
Spits out a tab-separated file containing substitutions for each site at each branch

Columns:
Site, ParentNode, ChildNode,  Nonsynonymous Substitutions, Synonymous Substitutions, Expected nonsynonymous Substitutions, Expected synonymous substitutions


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


#include "matrix_site_branch_subst.bf";

SetDialogPrompt ("Please specify Site-Branch Substitution TSV File or Empty String to Skip:");
fprintf(PROMPT_FOR_FILE, CLEAR_FILE);
outPerSiteBrSubFile = LAST_FILE_PATH;
fprintf (stdout, "\t[WRITING Site-Branch Substitutions TO ",outPerSiteBrSubFile,"]\n");	


SetDialogPrompt ("Please specify Site dN-dS TSV File or Empty String to Skip:");
fprintf(PROMPT_FOR_FILE, CLEAR_FILE);
outPerSiteDnDsFile = LAST_FILE_PATH;
fprintf (stdout, "\t[WRITING Site dN-dS TO ",outPerSiteDnDsFile,"]\n");	

printSubTSV(outPerSiteBrSubFile, outPerSiteDnDsFile)
