/*
Spits out a tab-separated file containing substitutions for each site at each branch

Columns:
Site, ParentNode, ChildNode,  Nonsynonymous Substitutions, Synonymous Substitutions, Expected nonsynonymous Substitutions, Expected synonymous substitutions

and/or 

Spits out a tab-separate file containing site substituions and site dn-ds.
Columns:
Site, Observed S Changes, Observed NS Changes,  E[S Sites], E[NS Sites], dN, dS, dN-dS, Scaled dN-dS


Tries to duplicate the behaviour QuickSelectonDetection.bf and SGEmulator.bf as much as possible.
Differences: 
- Does not bother finding the null distro for checking if sites are statistically significant selection.
- Does not count branches in which parent or child has completely ambiguous codon NNN.
- Does not bother with CF3x4 where it tries to solve the bias that occurs when estimating 
	equilbrium nucleotide frequencies without taking stop codons into account.
  From http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0011230 
  Correcting the Bias of Empirical Frequency Parameter Estimators in Codon Models, Plos1, 2010
  the inaccuracy doesn't matter when we are calculating branch lengths or dN/dS.



*/


// Columns for Site-Branch Substitutions
__RESULT_SITEBR_COLS_SITE = "Site";
__RESULT_SITEBR_COLS_PARNAME = "ParentNodeName";
__RESULT_SITEBR_COLS_PARORDER = "ParentNodePostOrder";
__RESULT_SITEBR_COLS_CHILDNAME = "ChildNodeName";
__RESULT_SITEBR_COLS_CHILDORDER = "ChildNodePostOrder";
__RESULT_SITEBR_COLS_OS = "ObservedSynSubst";
__RESULT_SITEBR_COLS_ON = "ObservedNonsynSubst";
__RESULT_SITEBR_COLS_ES = "ExpectedSynSubst";
__RESULT_SITEBR_COLS_EN = "ExpectedNonsynSubst";
__RESULT_SITEBR_COLS_BRANCHLEN = "BranchLen";
__RESULT_SITEBR_COLS_PARSEQ = "ParentCodon";
__RESULT_SITEBR_COLS_CHILDSEQ = "ChildCodon";

__RESULT_SITEBR_COLI_SITE = 0;
__RESULT_SITEBR_COLI_PARNAME = 1;
__RESULT_SITEBR_COLI_PARORDER = 2;
__RESULT_SITEBR_COLI_CHILDNAME = 3;
__RESULT_SITEBR_COLI_CHILDORDER = 4;
__RESULT_SITEBR_COLI_OS = 5;
__RESULT_SITEBR_COLI_ON = 6;
__RESULT_SITEBR_COLI_ES = 7;
__RESULT_SITEBR_COLI_EN = 8;
__RESULT_SITEBR_COLI_BRANCHLEN = 9;
__RESULT_SITEBR_COLI_PARSEQ = 10;
__RESULT_SITEBR_COLI_CHILDSEQ = 11;

// Columns for Site Substitutions, dN-dS
__RESULT_SITE_COLS_SITE = "Site";
__RESULT_SITE_COLS_OS = "Observed S Changes";
__RESULT_SITE_COLS_ON = "Observed NS Changes";
__RESULT_SITE_COLS_ES = "E[S Sites]";
__RESULT_SITE_COLS_EN = "E[NS Sites]";
__RESULT_SITE_COLS_DS = "dS";
__RESULT_SITE_COLS_DN = "dN";
__RESULT_SITE_COLS_DN_MINUS_DS = "dN-dS";
__RESULT_SITE_COLS_SCALED_DN_MINUS_DS = "Scaled dN-dS";

__RESULT_SITE_COLI_SITE = 0;
__RESULT_SITE_COLI_OS = 1;
__RESULT_SITE_COLI_ON = 2;
__RESULT_SITE_COLI_ES = 3;
__RESULT_SITE_COLI_EN = 4;
__RESULT_SITE_COLI_DS = 5;
__RESULT_SITE_COLI_DN = 6;
__RESULT_SITE_COLI_DN_MINUS_DS = 7;
__RESULT_SITE_COLI_SCALED_DN_MINUS_DS = 8;


/*___________________________________________________________________________________________________________*/


// Prints Site-Branch substitutions to File
// dataMatrix:  matrix to write
// outSubstTsvFile file to write to
// RETURNS:  total lines written to file
function	PrintTableToFile (dataMatrix, outSubstTsvFile)
{
	totalLinesWritten=0;
	for (i=0; i<Rows(dataMatrix); i = i + 1)
	{
		bufferString = "" +  dataMatrix[i][0];  // first dataMatrix column doesn't need a preceeding delimiter
		for (j = 1; j < Columns (dataMatrix); j = j+1)
		{
		    bufferString = bufferString +  "\t" + dataMatrix[i][j];
		}
		bufferString = bufferString + "\n";	
		fprintf (outSubstTsvFile,bufferString);
		totalLinesWritten = totalLinesWritten + 1;
	}
	
	return totalLinesWritten;
}



// Gets comma delimited list of codons from sparse column vector
// siteInfo:  1 if the codon coded by the element index exists
function commaCodonList(siteInfo, codonCodeToStrMap)
{
	sList = "";
	for (codonCode=0; codonCode<Rows(siteInfo); codonCode=codonCode+1)
	{
		if (siteInfo[codonCode] > 0)
		{
			if (sList != "")
			{
				sList = sList + ",";
			}
			sList = sList + codonCodeToStrMap[codonCode];
		}
	}
	return sList;
}

// Makes a Nucleotide transition Matrix from an Observed Codon Position Frequency Table
// Nucleotide transition rates are relative to AG
// _rNuc_cPos_2fObsFreq:   4 x 3 matrix.
//  _rNuc_cPos_2fObsFreq[Nucleotide 0-3][Position in Codon 0-2] = observed frequency of that nucleotide in that codon position
// Expects that global model parameters AC, AT, CT, ... etc are defined
function makeNucTransMatrixFromCodonPosFreqMatrix(_rNuc_cPos_2fObsFreq)
{
	_pooledFreqs = {4,1};

	for (iNuc=0; iNuc<4; iNuc=iNuc+1)
	{
		_pooledFreqs[iNuc] = (_rNuc_cPos_2fObsFreq[iNuc][0]+_rNuc_cPos_2fObsFreq[iNuc][1]+_rNuc_cPos_2fObsFreq[iNuc][2])/3;  // observedFreq from fit_codon_model.ibf
	}
	
	_matNucTrans = {{1,AC__*_pooledFreqs[1],_pooledFreqs[2],AT__*_pooledFreqs[3]}
				{AC__*_pooledFreqs[0],1,CG__*_pooledFreqs[2],CT__*_pooledFreqs[3]}
				{_pooledFreqs[0],CG__*_pooledFreqs[3],1,GT__*_pooledFreqs[3]}
				{AT__*_pooledFreqs[0],CT__*_pooledFreqs[3],GT__*_pooledFreqs[2],1}};
	
	return _matNucTrans;
}


// Prints out the per-site-branch substitutions and/or per-site dN-dS
function printSubTSV(outPerSiteBrSubFile, outPerSiteDnDsFile)
{
	
	// ReconstructAncestorSeq.bf will fit nucleotide model, codon model, reconstruct ancestral sequences
	// Defines:
	// DataSet dsAnc - ancestral sequences
	// DataSet		   dsJoint - both leaf and ancestral sequences
	// DataSetFiler filteredDataJoint has both leaf sequences and internal node sequences filtered for stop codons
	// codon_lf:  codon model fit
	// observedFreq:  observed nucleotide frequencies per codon position
	//				observedFreq[Nucleotide 0-3][Position in Codon 0-2] = observed frequency of that nucleotide in that codon position
	ExecuteAFile ( HYPHY_BASE_DIRECTORY + DIRECTORY_SEPARATOR +  "ReconstructAncestorSeq.bf");

	

	if (outPerSiteBrSubFile != "-1") 
	{
		//Print column headers for per-site-branch subs
		labelMatrixSiteBr = {1,12};
		labelMatrixSiteBr[__RESULT_SITEBR_COLI_SITE] = __RESULT_SITEBR_COLS_SITE;
		labelMatrixSiteBr[__RESULT_SITEBR_COLI_PARNAME] = __RESULT_SITEBR_COLS_PARNAME;
		labelMatrixSiteBr[__RESULT_SITEBR_COLI_PARORDER] = __RESULT_SITEBR_COLS_PARORDER;
		labelMatrixSiteBr[__RESULT_SITEBR_COLI_CHILDNAME] = __RESULT_SITEBR_COLS_CHILDNAME;
		labelMatrixSiteBr[__RESULT_SITEBR_COLI_CHILDORDER] = __RESULT_SITEBR_COLS_CHILDORDER;
		labelMatrixSiteBr[__RESULT_SITEBR_COLI_OS] = __RESULT_SITEBR_COLS_OS;
		labelMatrixSiteBr[__RESULT_SITEBR_COLI_ON] = __RESULT_SITEBR_COLS_ON;
		labelMatrixSiteBr[__RESULT_SITEBR_COLI_ES] = __RESULT_SITEBR_COLS_ES;
		labelMatrixSiteBr[__RESULT_SITEBR_COLI_EN] = __RESULT_SITEBR_COLS_EN;
		labelMatrixSiteBr[__RESULT_SITEBR_COLI_BRANCHLEN] = __RESULT_SITEBR_COLS_BRANCHLEN;
		labelMatrixSiteBr[__RESULT_SITEBR_COLI_PARSEQ] = __RESULT_SITEBR_COLS_PARSEQ;
		labelMatrixSiteBr[__RESULT_SITEBR_COLI_CHILDSEQ] = __RESULT_SITEBR_COLS_CHILDSEQ;

		fprintf (outPerSiteBrSubFile, CLEAR_FILE, KEEP_OPEN);  // KEEP_OPEN does not close the file after writing to it
		

		// Write header
		PrintTableToFile  (labelMatrixSiteBr, outPerSiteBrSubFile);
	}
	
	if (outPerSiteDnDsFile != "-1") 
	{	
		//Print column headers for per-site dnds
		labelMatrixSite = {1,9};
		labelMatrixSite[__RESULT_SITE_COLI_SITE] = __RESULT_SITE_COLS_SITE;
		labelMatrixSite[__RESULT_SITE_COLI_OS] = __RESULT_SITE_COLS_OS;
		labelMatrixSite[__RESULT_SITE_COLI_ON] = __RESULT_SITE_COLS_ON;
		labelMatrixSite[__RESULT_SITE_COLI_ES] = __RESULT_SITE_COLS_ES;
		labelMatrixSite[__RESULT_SITE_COLI_EN] = __RESULT_SITE_COLS_EN;
		labelMatrixSite[__RESULT_SITE_COLI_DS] = __RESULT_SITE_COLS_DS;
		labelMatrixSite[__RESULT_SITE_COLI_DN] = __RESULT_SITE_COLS_DN;
		labelMatrixSite[__RESULT_SITE_COLI_DN_MINUS_DS] = __RESULT_SITE_COLS_DN_MINUS_DS;
		labelMatrixSite[__RESULT_SITE_COLI_SCALED_DN_MINUS_DS] = __RESULT_SITE_COLS_SCALED_DN_MINUS_DS;

		fprintf (outPerSiteDnDsFile, CLEAR_FILE, KEEP_OPEN);  // KEEP_OPEN does not close the file after writing to it		

		// Write header
		PrintTableToFile  (labelMatrixSite, outPerSiteDnDsFile);
	}
	
	
	/* Call CodonToolsMain.def to calculate the number of synonymous and nonsynonymous substitutions between all possible codon transitions, 
			averaged over all 1-nucleotide-substitution paths from initial to end codon.
	  Also calculates the average number of synonymous and nonsynonymous sites between all possible codon pairs.
	  
	  Will create matrices:
	  _PAIRWISE_S_:    _PAIRWISE_S_[init codon 0-61][end codon 0-61] = ave number of syn sites between codon pair
	  _PAIRWISE_NS_: _PAIRWISE_NS_[init codon 0-61][end codon 0-61] = ave number of nonsyn sites between codon pair
	  _OBSERVED_S_: _OBSERVED_S_[init codon 0-61][end codon 0-61] = ave number of observed syn subst between codon pair
	  _OBSERVED_NS_:  _OBSERVED_NS_[init codon 0-61][end codon 0-61] = ave number of observed nonsyn subst between codon pair
	
	
	  It requires that the following matrices are already defined:
	  _EFV_MATRIX0_ - This is supposed to be the nucleotide transition frequency for 1st codon position
	  _EFV_MATRIX1_	- This is supposed to be the nucleotide transition frequency for 2nd codon position
	  _EFV_MATRIX2_ - This is supposed to be the nucleotide transition frequency for 3rd codon position
	*/
	//We use the same nucleotide transition frequency for all codon positions
	_EFV_MATRIX0_ = makeNucTransMatrixFromCodonPosFreqMatrix(observedFreq);
	_EFV_MATRIX1_ 			= _EFV_MATRIX0_;
	_EFV_MATRIX2_ 			= _EFV_MATRIX0_;
				

	incFileName = HYPHY_LIB_DIRECTORY + DIRECTORY_SEPARATOR + "TemplateBatchFiles" + DIRECTORY_SEPARATOR +  "Distances" + DIRECTORY_SEPARATOR + "CodonToolsMain.def";
	ExecuteCommands  ("#include \""+incFileName+"\";");



	numCdnCodes = Columns(_Genetic_Code);
	branchNames = BranchName (codon_tree,-1);  // all the branch names in post-order traversal.  
	numBranch = Columns (branchNames);	


	/* maps branch post order to sequence name in filteredDataJoint*/
	branchToSeqMap = {numBranch, 1};  // maps post traversal index to sequence index
	for (v=0; v<numBranch; v=v+1)
	{
		for (k=0; k<filteredDataJoint.species && !isFound; k=k+1)  
		{
			GetString (seqName, filteredDataJoint, k);	
			if (branchNames[v] % seqName)  // % means case insensitive string match
			{
				branchToSeqMap[v] = k;
			}
		}
	}

	/* total tree length */
	totalTreeLength = 0;
	
	branchLengths   = BranchLength(codon_tree,-1);  // all branch lengths in post-order traversal

	for (v=Columns(branchLengths)-1; v>=0; v=v-1)
	{
		totalTreeLength = totalTreeLength + branchLengths[v];
	}



	// Multiply matrixTrick by a column vector to extract the index of non-zero values
	matrixTrick  = {1,stateCharCount};
	matrixTrick  = matrixTrick["_MATRIX_ELEMENT_COLUMN_"];
	// Multiply matrixTrickSummer by column vector to sum all the contents of column vector
	matrixTrickSummer = {1,stateCharCount};
	matrixTrickSummer = matrixTrickSummer["1"];  

	/* get codon matrix */
	rSeq_cUniqSite_2iCdn = {filteredDataJoint.species, filteredDataJoint.unique_sites}; // [0-based ancestor index in filteredDataJoint, unique codon site index] = codon code from 0-61
	for (k=0; k<filteredDataJoint.species;k=k+1)
	{
		for (v=0; v<filteredDataJoint.unique_sites;v=v+1)
		{
			GetDataInfo (siteInfo, filteredDataJoint, k, v);
			if ((matrixTrickSummer*siteInfo)[0] > 1)  // >1  possible codon at site v for this sequence
			{
				rSeq_cUniqSite_2iCdn[k][v] = -1;
			}
			else  // only 1 possible codon at site v
			{
				rSeq_cUniqSite_2iCdn[k][v] = (matrixTrick * siteInfo)[0];
			}
		}
	}


	GetDataInfo	   (site2uniqSite, filteredDataJoint); // site2uniqSite[site index] = unique site index for that site


	flatTreeRep	  = Abs (codon_tree);  // numeric vector, where the i-th entry stores the post-order traversal index for the parent of the i-th node (also in post-order traversal);

	// Converts numerical codon code to string, defined in chooseGeneticCode.def
	codonCodeToStrMap = ComputeCodonCodeToStringMap (_Genetic_Code);


	// For each site, traverse each branch.  
	// Count substitutions at each site-branch.
	for (iSite=0; iSite<filteredDataJoint.sites;iSite=iSite+1)
	{
	
		totalSiteOS = 0;
		totalSiteON = 0;
		totalSiteES = 0;
		totalSiteEN = 0;
			
		iUniqSite = site2uniqSite[iSite];
		for (iChildPostOrder=0; iChildPostOrder<numBranch; iChildPostOrder=iChildPostOrder+1) // branches are in Post-Order
		{
			_SITEBR_ES_COUNT = {stateCharCount,stateCharCount};
			_SITEBR_EN_COUNT = {stateCharCount,stateCharCount};
			_SITEBR_OS_COUNT = {stateCharCount,stateCharCount};
			_SITEBR_ON_COUNT = {stateCharCount,stateCharCount};
	
			iParentPostOrder = flatTreeRep[iChildPostOrder];
			if (iParentPostOrder < 0)
			{
				continue;  // The current child node is the root node, there is no parent
			}
		
			// It's possible for the tree to be rooted with a leaf sequence
			iParentSeq = branchToSeqMap[iParentPostOrder]; 
			sParentName = branchNames[iParentPostOrder];
			iParCdn = rSeq_cUniqSite_2iCdn[iParentSeq][iUniqSite];			
		
		
			iChildSeq = branchToSeqMap[iChildPostOrder];
			sChildName = branchNames[iChildPostOrder];
			iChildCdn = rSeq_cUniqSite_2iCdn[iChildSeq][iUniqSite];	
		
			branchFactor = branchLengths[iChildPostOrder]/totalTreeLength;
		
			sParCdn = "";
			sChildCdn = "";
			if (iChildCdn>=0 && iParCdn>=0)  		//no ambiguities
			{			
				_SITEBR_OS_COUNT[iParCdn][iChildCdn] = _SITEBR_OS_COUNT[iParCdn][iChildCdn] + 1;		
				_SITEBR_ON_COUNT[iParCdn][iChildCdn] = _SITEBR_ON_COUNT[iParCdn][iChildCdn] + 1;		
				_SITEBR_ES_COUNT[iParCdn][iChildCdn] = _SITEBR_ES_COUNT[iParCdn][iChildCdn] + branchFactor;		
				_SITEBR_EN_COUNT[iParCdn][iChildCdn] = _SITEBR_EN_COUNT[iParCdn][iChildCdn] + branchFactor;		

				sChildCdn = codonCodeToStrMap[iChildCdn];
				sParCdn = codonCodeToStrMap[iParCdn];
			}	
			else  //ambiguities here
			{

				GetDataInfo    (childAmbInfo, filteredDataJoint, iChildSeq, iUniqSite);
				GetDataInfo    (parAmbInfo, filteredDataJoint, iParentSeq, iUniqSite);		
			
				totalPossChildCodons = (matrixTrickSummer*childAmbInfo)[0];
				totalPossParCodons = (matrixTrickSummer*parAmbInfo)[0];
				totalPossPairs = totalPossParCodons * totalPossChildCodons;
			
				sChildCdn = commaCodonList(childAmbInfo, codonCodeToStrMap);
				sParCdn = commaCodonList(parAmbInfo, codonCodeToStrMap);
			
				// ??? We aren't outputting the site-substitutions for NNN codons, but NNN still affect the codon resolution frequency across the phylo for 2-N or 1-N ambiguous codons.
				// Ignore site-branches where parent or child has completely ambiguous codons.
				// This is because every substitution is possible and will pollute any inference we make at this site.
				if ( totalPossParCodons < stateCharCount && totalPossChildCodons < stateCharCount)  
				{
					// Count codon frequencies for this site only
					siteFilter = ""+(iSite*3)+"-"+(iSite*3+2);
					DataSetFilter filteredDataSite = CreateFilter (dsJoint,3,siteFilter,"",GeneticCodeExclusions);  // Remove any site where either inner node or leaf has stop codon
					// Count all codon frequencies, including stop codons.  
					// Frequencies will total to 1.
					// All the possible codons corresponding to an ambiguous codon count as a fraction of a codon.  EG)  GGN => GGA=0.25, GGC=0.25, GGT=0.25, GGG=0.25
					HarvestFrequencies			  (observedCEFV,filteredDataSite,3,3,0);  
					tempMx = {stateCharCount,1};

					// Remove stop codons from codon frequency table
					hShift = 0;			
					for (k=0; k<numCdnCodes; k=k+1)
					{
						if (IsStop(k, _Genetic_Code))  //IsStop() defined in chooseGeneticCode.def
						{
							hShift = hShift+1;
						}
						else
						{
							tempMx[k-hShift] = observedCEFV[k];
						}
					}	
					observedCEFV = tempMx;		
		
					// childAmbInfo, parAmbInfo only tell us if possible codon exists
					// observedCEFV tells us the frequency of the possible codon
					childAmbCodon_phyloFreqs = childAmbInfo$observedCEFV;  //For each possible codons at this site, we get the frequency of its occurences across the phylo
					parAmbCodon_phyloFreqs = parAmbInfo$observedCEFV;  //For each possible codons at this site, we get the frequency of its occurences across the phylo
				
					// sum frequency of occurences across phylo for possible child codons at this site, including occurences from resolution of ambiguous codons
					totalFreqChildPossCodons = (matrixTrickSummer*childAmbCodon_phyloFreqs)[0];
					// sum frequency of occurences across phylo for possible parent codons at this site, including occurences from resolution of ambiguous codons    
					totalFreqParPossCodons = (matrixTrickSummer*parAmbCodon_phyloFreqs)[0];    
				
				
					// NB:  If COUNT_GAPS_IN_FREQUENCIES = 0, then HarvestFrequencies does not resolve codons with gaps before counting frequencies.
					// BUT even if COUNT_GAPS_IN_FREQUENCIES = 0, GetDataInfo() will resolve codons with gaps!!!
					// Thus it is possible to get possible codons at zero frequency.  Ignore these site-branches.
					if (totalFreqChildPossCodons == 0 || totalFreqParPossCodons == 0)
					{
						continue;
					}
					
					for (iPossParCdn=0; iPossParCdn<stateCharCount; iPossParCdn=iPossParCdn+1)
					{
						if (parAmbInfo[iPossParCdn]>0)
						{		
							parWeightFactor = parAmbCodon_phyloFreqs[iPossParCdn]/totalFreqParPossCodons;  // ambiguous codon count weighted by frequency across phylo
							for (iPossChildCdn=0; iPossChildCdn<stateCharCount; iPossChildCdn=iPossChildCdn+1)
							{
								if (childAmbInfo[iPossChildCdn]>0)
								{		
									childWeightFactor = childAmbCodon_phyloFreqs[iPossChildCdn]/totalFreqChildPossCodons;  // ambiguous codon count weighted by frequency across phylo
									avePairWeightfactor = parWeightFactor *  childWeightFactor;
									_SITEBR_OS_COUNT[iPossParCdn][iPossChildCdn] = _SITEBR_OS_COUNT[iPossParCdn][iPossChildCdn] + avePairWeightfactor;		
									_SITEBR_ON_COUNT[iPossParCdn][iPossChildCdn] = _SITEBR_ON_COUNT[iPossParCdn][iPossChildCdn] + avePairWeightfactor;		
									_SITEBR_ES_COUNT[iPossParCdn][iPossChildCdn] = _SITEBR_ES_COUNT[iPossParCdn][iPossChildCdn] + avePairWeightfactor*branchFactor;		
									_SITEBR_EN_COUNT[iPossParCdn][iPossChildCdn] = _SITEBR_EN_COUNT[iPossParCdn][iPossChildCdn] + avePairWeightfactor*branchFactor;
								}
							} 
						}
					}
				}
			}
		
			// Aggregate each codon-codon transition count for this site-branch
			// Weight each transition by the average substitutions across all 1-substituion paths from initial codon to end codon
			totalSiteBrOS = (matrixTrickSummer*(_OBSERVED_S_$_SITEBR_OS_COUNT)*Transpose(matrixTrickSummer))[0];
			totalSiteBrON = (matrixTrickSummer*(_OBSERVED_NS_$_SITEBR_ON_COUNT)*Transpose(matrixTrickSummer))[0];		
			totalSiteBrES = (matrixTrickSummer*(_PAIRWISE_S_$_SITEBR_ES_COUNT) *Transpose(matrixTrickSummer))[0];
			totalSiteBrEN = (matrixTrickSummer*(_PAIRWISE_NS_$_SITEBR_EN_COUNT) *Transpose(matrixTrickSummer))[0];
			
			// Aggregate each codon-codon transition count for this site across all branches
			totalSiteOS = totalSiteOS + totalSiteBrOS;
			totalSiteON = totalSiteON + totalSiteBrON;
			totalSiteES = totalSiteES + totalSiteBrES;
			totalSiteEN = totalSiteEN + totalSiteBrEN;
		
			if (outPerSiteBrSubFile != "-1")
			{			
				resultMatrixSiteBr = {1, 12};  // Convert to string so that names can be included
				resultMatrixSiteBr[__RESULT_SITEBR_COLI_SITE] = "" + iSite;
				resultMatrixSiteBr[__RESULT_SITEBR_COLI_PARORDER] = "" + iParentPostOrder;
				resultMatrixSiteBr[__RESULT_SITEBR_COLI_PARNAME] = "" + sParentName;
				resultMatrixSiteBr[__RESULT_SITEBR_COLI_PARSEQ] = "" + sParCdn;
				resultMatrixSiteBr[__RESULT_SITEBR_COLI_CHILDORDER] = ""+ iChildPostOrder;
				resultMatrixSiteBr[__RESULT_SITEBR_COLI_CHILDNAME] = ""+ sChildName;
				resultMatrixSiteBr[__RESULT_SITEBR_COLI_CHILDSEQ] = ""+ sChildCdn;
				resultMatrixSiteBr[__RESULT_SITEBR_COLI_OS] = "" + totalSiteBrOS;
				resultMatrixSiteBr[__RESULT_SITEBR_COLI_ON] = "" + totalSiteBrON;
				resultMatrixSiteBr[__RESULT_SITEBR_COLI_ES] = "" + totalSiteBrES;
				resultMatrixSiteBr[__RESULT_SITEBR_COLI_EN] = "" + totalSiteBrEN;
				resultMatrixSiteBr[__RESULT_SITEBR_COLI_BRANCHLEN] = "" + branchLengths[iChildPostOrder];

	
				//Print to CSV: Per-Site-Branch substitutions
				totalLinesWritten = PrintTableToFile  (resultMatrixSiteBr, outPerSiteBrSubFile);
				assert(totalLinesWritten>0, "No site-branch substitutions were written to the file ");
			}
		}
		
		if (outPerSiteDnDsFile != "-1")  //Print to Site dN-dS Results File
		{					
			site_dN = (totalSiteON/totalSiteEN);
			site_dS = (totalSiteOS/totalSiteES);
			site_dN_minus_dS = site_dN - site_dS;
			site_scaled_dN_minus_dS = site_dN_minus_dS/totalTreeLength;
	
			resultMatrixSite = {1, 9};
			resultMatrixSite[__RESULT_SITE_COLI_SITE] = iSite;
			resultMatrixSite[__RESULT_SITE_COLI_OS] = totalSiteOS;
			resultMatrixSite[__RESULT_SITE_COLI_ON] = totalSiteON;
			resultMatrixSite[__RESULT_SITE_COLI_ES] = totalSiteES;
			resultMatrixSite[__RESULT_SITE_COLI_EN] = totalSiteEN;
			resultMatrixSite[__RESULT_SITE_COLI_DN] = site_dN;
			resultMatrixSite[__RESULT_SITE_COLI_DS] = site_dS;
			resultMatrixSite[__RESULT_SITE_COLI_DN_MINUS_DS] = site_dN_minus_dS;
			resultMatrixSite[__RESULT_SITE_COLI_SCALED_DN_MINUS_DS] = site_scaled_dN_minus_dS;

			
			totalLinesWritten = PrintTableToFile  (resultMatrixSite, outPerSiteDnDsFile);
			assert(totalLinesWritten>0, "No site-branch substitutions were written to the file ");
			
		}
	}
	fprintf (stdout, "\nTotal codon tree length (subs/nuc/unit time): ", Format (totalTreeLength,10,5), "\n");
	if (outPerSiteBrSubFile != "-1")
	{
		fprintf (outPerSiteBrSubFile, CLOSE_FILE);  // close file
	}
	if (outPerSiteDnDsFile != "-1")
	{
		fprintf (outPerSiteDnDsFile, CLOSE_FILE);  // close file
	}

}
