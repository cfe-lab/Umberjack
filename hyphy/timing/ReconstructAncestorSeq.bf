


fscanf  (stdin,  "String", fastaFile);
fscanf  (stdin,  "String", treeFile);
fscanf  (stdin,  "String", outNucModelfit);
fscanf  (stdin,  "String", outCodonModelfit);
fscanf (stdin, "String", outLeafAncFastaFile);
fscanf (stdin, "String", outSubstTsvFile);
fscanf (stdin, "String", outCodonTreeFile);
fscanf (stdin, "String", outCodonTreeBranchLenCsv);

// Fit nucleotide model, codon model
// Write out tree scaled into units of codon subst/site, write out codon tree branch length csv
inputs = {};
inputs[0] = "Universal";   // genetic code  - chooseGeneticCode.def
inputs[1] = fastaFile;  // codon fasta file - fit_codon_model.bf
inputs[2] = treeFile;  // tree - fit_codon_model.bf
inputs[3] = outNucModelfit;  //output nucleotide model fit file - fit_nuc_model.ibf
inputs[4] = outCodonModelfit;  //output codon model fit file - fit_codon_model.ibf
inputs[5] = outCodonTreeFile;  // output codon tree file
inputs[6] = outCodonTreeBranchLenCsv;  //output codon tree branch lengths csv file
ExecuteAFile ( HYPHY_BASE_DIRECTORY + DIRECTORY_SEPARATOR + "fit_codon_model.bf", inputs);


// Reconstruct ancestors
fprintf (stdout, "\t[RECONSTRUCTING ANCESTORS]\n");	
matrix_site_branch_subst_inputs = {};
matrix_site_branch_subst_inputs[0] = outLeafAncFastaFile;
matrix_site_branch_subst_inputs[1] = outSubstTsvFile;
ExecuteAFile ( HYPHY_BASE_DIRECTORY + DIRECTORY_SEPARATOR +  "matrix_site_branch_subst.bf", matrix_site_branch_subst_inputs);










