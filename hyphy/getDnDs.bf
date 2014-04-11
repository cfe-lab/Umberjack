fscanf(stdin, "String",	genCode);
fscanf(stdin, "String",	newRestore);
fscanf(stdin, "String",	codonFasta);
fscanf(stdin, "String",	model);
fscanf(stdin, "String",	treeFile);
fscanf(stdin, "String",	modelFitFile);
fscanf(stdin, "String",	dndsBias);
fscanf(stdin, "String",	ancestorCount);
fscanf(stdin, "String",	branchCorr);
fscanf(stdin, "String",	SLAC);
fscanf(stdin, "String",	treatAmbig);
fscanf(stdin, "String",	testStat);
fscanf(stdin, "String",	sigLevel);
fscanf(stdin, "String",	outputOption);
fscanf(stdin, "String",	outputTSV);
fscanf(stdin, "String",	rateClassEst);
/*fscanf(stdin, "String",	resultsProcess);*/

fprintf(stdout, "00 GeneticCode=" + genCode + "\n");
fprintf(stdout, "01 New/Restore=" + newRestore + "\n");
fprintf(stdout, "02 Codon Fasta=" + codonFasta + "\n");
fprintf(stdout, "03 Model Options=" + newRestore + "\n");
fprintf(stdout, "04 Tree File=" + treeFile + "\n");
fprintf(stdout, "05 Model Fit File=" + modelFitFile + "\n");
fprintf(stdout, "06 DN/DS Bias=" + dndsBias + "\n");
fprintf(stdout, "07 Ancestor Counting Option=" + ancestorCount + "\n");
fprintf(stdout, "08 Branch Correction=" + branchCorr + "\n");
fprintf(stdout, "09 SLAC Option=" + SLAC + "\n");
fprintf(stdout, "10 Treatment of Ambiguities=" + treatAmbig + "\n");
fprintf(stdout, "11 Test Statistic=" + testStat + "\n");
fprintf(stdout, "12 Significance Level=" + sigLevel + "\n");
fprintf(stdout, "13 Output Options=" + outputOption + "\n");
fprintf(stdout, "14 Output DN/DS TSV=" + outputTSV + "\n");
fprintf(stdout, "15 Rate Class Estimator=" + rateClassEst + "\n");
/*fprintf(stdout, "ResultProcessingTool=" + resultsProcess + "\n");*/


stdinRedirect = {};
stdinRedirect["00"] = genCode;
stdinRedirect["01"] = newRestore;
stdinRedirect["02"] = codonFasta;
stdinRedirect["03"] = model;
stdinRedirect["04"] = treeFile;
stdinRedirect["05"] = modelFitFile;
stdinRedirect["06"] = dndsBias;
stdinRedirect["07"] = ancestorCount;
stdinRedirect["08"] = branchCorr;
stdinRedirect["09"] = SLAC;
stdinRedirect["10"] = treatAmbig;
stdinRedirect["11"] = testStat;
stdinRedirect["12"] = sigLevel;
stdinRedirect["13"] = outputOption;
stdinRedirect["14"] = outputTSV;
stdinRedirect["15"] = rateClassEst;
/*stdinRedirect["16"] = resultsProcess;*/


/*fprintf(stdout, HYPHY_LIB_DIRECTORY + DIRECTORY_SEPARATOR + "TemplateBatchFiles" + DIRECTORY_SEPARATOR )*/
ExecuteAFile (HYPHY_LIB_DIRECTORY + DIRECTORY_SEPARATOR + "TemplateBatchFiles" + DIRECTORY_SEPARATOR + "QuickSelectionDetection.bf", stdinRedirect);

