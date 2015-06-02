COUNT_GAPS_IN_FREQUENCIES = 0;  // Do not include gaps (-) in the codon, nucleotide frequencies.  Gaps will not get resolved.  N's still get resolved though.


// Fit nucleotide model, codon model
// Write out tree scaled into units of codon subst/site, write out codon tree branch length csv
ExecuteAFile ( HYPHY_BASE_DIRECTORY + DIRECTORY_SEPARATOR + "fit_codon_model.bf");  // Allow underlying batch files wlil prompt user for input


// Reconstruct ancestors
fprintf (stdout, "\t[RECONSTRUCTING ANCESTORS]\n");	

DataSet dsAnc					= ReconstructAncestors (codon_lf);  //codon_lf defined in fit_codon_model.ibf
DataSet		   dsJoint = Combine(dsAnc,codon_ds);
DataSetFilter filteredDataJoint = CreateFilter (dsJoint,3,"","",GeneticCodeExclusions);  // filteredDataJoint has both leaf sequences and internal node sequences filtered for stop codons

// Output leafs and ancestors to nexus
DATA_FILE_PRINT_FORMAT = 9;  // FASTA SEQUENTIAL
DataSetFilter filteredDataAll = CreateFilter(dsJoint,1);
SetDialogPrompt ("Please specify output nexus file for ancestors and leafs:");
fprintf(PROMPT_FOR_FILE, CLEAR_FILE);
ancestoroutSubstTsvFile = LAST_FILE_PATH;
fprintf (stdout, "\t[WRITING ANCESTORS TO ",ancestoroutSubstTsvFile,"]\n");	
fprintf (ancestoroutSubstTsvFile, filteredDataAll);









