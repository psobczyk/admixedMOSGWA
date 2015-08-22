/********************************************************************************
 *	This file is part of the MOSGWA program code.				*
 *	Copyright ©2011–2015, Erich Dolejsi, Bernhard Bodenstorfer.		*
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

#include "Parameter.hpp"
#include <cmath>
#include <limits>
#include <cstdlib>
#include <sstream>  //for stringstream
#include <fstream>

#include "logging/Logger.hpp"

using namespace std;
using namespace parser;
using namespace logging;

/** Internally declares the configuration variables. */
Parameter::Parameter () {

	// input file settings
	declare( "input", "plink_files", in_files_plink );
	declare( "input", "hdf5_file", in_file_hdf5 );
	declare( "input", "use_extra_yvalues", y_value_extra_file );
	declare( "input", "use_extra_covariables", cov_extra_file );
	declare(
		"input",
		"missing_phenotype_value",
		missing_phenotype_value = numeric_limits<double>::quiet_NaN()
	);
	declare( "input", "control_value", control_value = 0. );
	declare( "input", "case_value", case_value = 1. );
	declare( "input", "models_file", models_file);
	declare( "input", "cache_limit", cache_limit = 1024 );

	// input data settings
	declare( "data", "trait_index", in_values_int );
	{
		map< const string, int > choice;
		choice[ "recessive" ] = 1;
		choice[ "additive" ] = 2;
		choice[ "dominant" ] = 3;
		declare( "data", "genetic_model", genetic_model, choice );
	}

	// output file settings
	declare( "output", "files", out_file_name );
	declare( "output", "singlefile", singlefile);//for creating the Hlasso genotyp file only once
	declare( "output", "correlation_threshold", correlation_threshold = 0.999 );

	// log settings
	{
		map< const string, int > choice;
		choice[ "DEBUG" ] = 0;
		choice[ "INFO" ] = 1;
		choice[ "WARNING" ] = 2;
		choice[ "ERROR" ] = 3;
		declare( "log", "level", log_level = 1, choice );
	}
	declare( "log", "silent", silent );
	declare( "log", "detailed", detailed_selction );

	// single marker test settings
	{
		map< const string, int > choice;
		choice[ "chi_square" ] = singleMarkerTest_CHI_SQUARE;
		choice[ "cochran_armitage" ] = singleMarkerTest_COCHRAN_ARMITAGE;
		declare(
			"single_marker",
			"test",
			singleMarkerTest = singleMarkerTest_CHI_SQUARE,
			choice
		);
	}

	// model preselection (i.e. first round with mBIC) settings
	declare( "model_preselection", "mBIC_expected_causal_SNPs", mBIC_firstRound_expectedCausalSNPs = 20 );

	// model selection settings
	declare( "model_selection", "mBIC_expected_causal_SNPs", mBIC_expectedCausalSNPs = 1 );
	declare( "model_selection", "EBIC_gamma", EBIC_gamma = ::nan( "default" ) );
	declare( "model_selection", "nSNPKriterium", nSNPKriterium = 0 );
	{
		map< const string, int > choice;
		choice[ "BIC" ] = selectionCriterium_BIC;
		choice[ "EBIC" ] = selectionCriterium_EBIC;
		choice[ "mBIC" ] = selectionCriterium_mBIC;
		choice[ "mBIC2" ] = selectionCriterium_mBIC2;
		declare(
			"model_selection",
			"selection_criterium",
			selectionCriterium = selectionCriterium_mBIC2,
			choice
		);
	}
	declare( "model_selection", "maximalModelSize", maximalModelSize = 3 );
	declare( "model_selection", "multi_forward_step_max", ms_MaximalSNPsMultiForwardStep = 1 );
	declare( "model_selection", "multi_forward_pvalue_max", ms_MaximalPValueForwardStep );
	declare( "model_selection", "forward_step_max", ms_forward_step_max = 1 );
	declare( "model_selection", "fast_multi_forward", ms_FastMultipleForwardStep );
	declare( "model_selection", "PValueBorder", PValueBorder = 350 );
	declare( "model_selection", "reset", reset = 350 );
	declare( "model_selection", "jump_back", jump_back = 350 );
	declare( "model_selection", "saveguardsteps", saveguardsteps = 2 );
	{
		map< const string, int > choice;
		choice[ "greedy" ] = searchStrategy_greedy;
		choice[ "memetic_algorithm" ] = searchStrategy_memetic_algorithm;
		declare(
			"model_selection",
			"search_strategy",
			searchStrategy = searchStrategy_greedy,
			choice
		);
	}

	// log_regression settings
	declare( "log_regression", "max_it", logrC_maxit = 15 );
	declare( "log_regression", "max_hs", logrC_maxhs = 3 );
	declare( "log_regression", "max_step", logrC_maxstep = 3 );
	declare( "log_regression", "l_conv", logrC_lconv = 1e-2 );
	declare( "log_regression", "g_conv", logrC_gconv = 1e-2 );
	declare( "log_regression", "x_conv", logrC_xconv = 1e-2 );

	// memetics algorithm settings
	declare( "memetic_algorithm", "modelsNo", modelsNo = 40 );
	declare( "memetic_algorithm", "maxNoProgressIter", maxNoProgressIter = 1000 );
	declare( "memetic_algorithm", "B", B = 3 );
	declare( "memetic_algorithm", "pCross", pCross = 0.9 );
	declare( "memetic_algorithm", "pMutation", pMutation = 0.05 );
	declare( "memetic_algorithm", "tournamentSize", tournamentSize = 2 );
	declare( "memetic_algorithm", "correlationThreshold", correlationThreshold = 0.5 );
	declare( "memetic_algorithm", "correlationRange", correlationRange = 50 );
	declare( "memetic_algorithm", "causalModelFilename", causalModelFilename = "" );
	declare( "memetic_algorithm", "regionMinCorrelation", regionMinCorrelation = 0.01 );
	declare( "memetic_algorithm", "modelsFilename", modelsFilename );
	declare( "memetic_algorithm", "outNo", outNo );

	// Erich's testcase gearatator
	declare( "TESTING", "replications", replications );	//how many Y vectors should be produced
}

/** Initialise the parameter holder from the files given on the command line. */
void Parameter::setParameters ( const int argn, const char* argv[] ) {

	// Parse command line arguments as config files
	for ( int i = 1; i < argn; ++i ) {
		ifstream configuration;
		configuration.exceptions ( ifstream::badbit );
		try {
			configuration.open( argv[i], ifstream::in );
			// See http://gehrcke.de/2011/06/reading-files-in-c-using-ifstream-dealing-correctly-with-badbit-failbit-eofbit-and-perror/
			if ( ! configuration.is_open() ) {
				logger->error(
					"Cannot open configuration file \"%s\".",
					argv[i]
				);
			} else {
				parse( configuration, argv[i] );
			}
			configuration.close();
		} catch ( ifstream::failure e ) {
			logger->error(
				"Cannot read configuration file \"%s\": %s",
				argv[i],
				e.what()
			);
		}
	}
	logger->info( "y_value_extra_file=%d, position=%d", y_value_extra_file, in_values_int );
	//this creates an unique file for all yvm and of course also in the version without
	//this has to be done before the number is added in the yvm case
	singlefile = out_file_name;	//this is for simulators which need in every run the genotyp file
	if ( y_value_extra_file && 0 < in_values_int ) {
		stringstream number;
		number << in_values_int;
		out_file_name = out_file_name + number.str();
	}
}

Parameter * parameter = NULL;
