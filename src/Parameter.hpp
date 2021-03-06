/********************************************************************************
 *	This file is part of the MOSGWA program code.				*
 *	Copyright ©2011–2015, Bernhard Bodenstorfer.		*
 *										*
 *	This program is free software; you can redistribute it and/or modify	*
 *	it under the terms of the GNU General Public License as published by	*
 *	the Free Software Foundation; either version 3 of the License, or	*
 *	(at your option) any later version.					*
 *										*
 *	This program is distributed in the hope that it will be useful,		*
 *	but WITHOUT ANY WARRANTY; without even the implied warranty of		*
 *	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.			*
 *	See the GNU General Public License for more details.			*
 ********************************************************************************/

#ifndef PARAMETER_HPP
#define PARAMETER_HPP

#include <string>
#include "parser/ConfigParser.hpp"

/** Holds relevant additional parameters, which can be read from a file. */
class Parameter : protected parser::ConfigParser {

public:

//++++++++++++
// Input:

	/** path + name for the plinkfiles .bed .fam .bim */
	std::string in_files_plink;

	/** true = use extra file for Y-values, false use values in .fam file */
	bool y_value_extra_file;

	/** path + name for the HDF5 input file.
	* Use either this or PLink input files, not both.
	*/
	std::string in_file_hdf5;

	/** position of trait in .fam plus .yvm file, starting in .fam with 0, continuing to .yvm */
	int in_values_int;

	/** true = use extra file for covariables */
	bool cov_extra_file;

	/** model to read in a filename with first line length
	* of models
	* and then every line a SNP name but nothing else, check
	* for errors!!
	*/
	std::string models_file;

	/** Cache size limit in number of vectors.
	* It is interpreted as unsigned. Thus, as a side effect, negative means almost limitless.
	*/
	int cache_limit;

//+++++++++++++
// Parameters describing the Data

	/** if Phenotype is affection (case-control) or quantitative.
	* setting is determined by MData::checkYValues() */
	bool affection_status_phenotype;

	/** Helper to set {@link affection_status_phenotype}. */
	int regressionType;

	/** Codes for regressionType. */
	static const int
		regressionType_Linear = 1,
		regressionType_Firth = 3;

	/** 1 = Recessive, 2 = Additive, 3 = Dominant */
	int genetic_model;

	double
		/** which number repressents missing phenotypes */
		missing_phenotype_value,
		/** which number is coding the controls */
		control_value,
		/** which number is coding the cases */
		case_value;

//++++++++++++
// Output:

	/** path + name for the output files .log ... */
	std::string out_file_name;
	std::string singlefile;

	// Logging:

	/** severity threshold for logging */
	int log_level;

	/** true = no output on screen */
	bool silent;

	/** true = all model selection steps are written to logfile */
	bool detailed_selction;

	/** Threshold used in the display of "strongly correlated" SNPs. */
	double correlation_threshold;

//++++++++++++
//ModelSelection:

	/** Constants to interpret {@link #selectionCriterium}. */
	static const int
		selectionCriterium_BIC = 1,
		selectionCriterium_EBIC = 2,
		selectionCriterium_mBIC = 3,
		selectionCriterium_mBIC_firstRound = -3,	// uses different expected_causal_SNPs
		selectionCriterium_mBIC2 = 4;

	/** Which criterium to use. */
	int selectionCriterium;

	/** Constants to interpret {@link #searchStrategy}. */
	static const int
		searchStrategy_greedy = 0,
		searchStrategy_memetic_algorithm = 1;

	/** Which search strategy to use. */
	int searchStrategy;

	/** Parameter in the EBIC criterium. */
	double EBIC_gamma;

	/** Number of Expected SNPs, in the preselection with mBIC */
	int mBIC_expectedCausalSNPs;
	int mBIC_firstRound_expectedCausalSNPs;
	int nSNPKriterium;
	int maximalModelSize;
	/* *the maximal number of SNPs added in the Multi-Forward-Step */
	int ms_MaximalSNPsMultiForwardStep;

	/** the threshold below which SNPs with lower p-Values are considered for adding in the Forward Step */
	double ms_MaximalPValueForwardStep;
	/** a seperate variable for the new multiforward step */
	int ms_forward_step_max;
	bool ms_FastMultipleForwardStep;

	size_t PValueBorder;
	size_t reset;
	int	jump_back;
	int saveguardsteps;
//++++++++++++
// Parameters for Affection/Case-Control

	/** Constants to interpret {@link #singleMarkerTest}. */
	static const int
		singleMarkerTest_CHI_SQUARE = 1,
		singleMarkerTest_COCHRAN_ARMITAGE = 2;

	/** What to use for the initial single-marker tests. */
	int singleMarkerTest;

// Logistic Regression Control Parameters
	int logrC_maxit;
	int logrC_maxhs;
	int logrC_maxstep;
	double logrC_lconv;
	double logrC_gconv;
	double logrC_xconv;

// Memetic algorithm
  int modelsNo;
  int maxNoProgressIter;
  int B;
  double pCross;
  double pMutation;
  int tournamentSize;
  double correlationThreshold;
  int correlationRange;
  int maxPoolSize;
	double regionMinCorrelation;
//	int outNo;
//	std::string modelsFilename;

// TEST
	int replications;

	/** Constructor */
	Parameter ();

	/** Initialise */
	void setParameters ( const int argn, const char* argv[] );
};

/** To be initialised to point to a globally shared parameter object */
extern Parameter * parameter;

#endif
