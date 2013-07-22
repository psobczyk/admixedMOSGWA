#ifndef PARAMETER_HPP
#define PARAMETER_HPP
#include <vector>  //for SNPs
#include <string>
#include <iostream>
#include <fstream>

#include "parser/ConfigParser.hpp"

using namespace std;

using namespace parser;

/** Holds relevant additional parameters, which can be read from a file. */
class Parameter : protected ConfigParser {

public:

//++++++++++++
// Input:

	/** path + name for the plinkfiles .bed .fam .bim */
	string in_files_plink;

	/** path + name for the plinkfiles .bed (if different) */
	string in_files_plink_bed;

	/** path + name for the plinkfiles .fam (if different) */
	string in_files_plink_fam;

	/** path + name for the plinkfiles .bim (if different) */
	string in_files_plink_bim;

	// Y-value
	/** true = use extra file for Y-values, false use values in .fam file */
	bool y_value_extra_file;

	/** path + name for the HDF5 input file.
	* Use either this or PLink input files, not both.
	*/
	string in_file_hdf5;

	/** path + name for the Y-value matrix file .yvm (if different) */
	string in_files_values_yvm;

	/** position of trait in .yvm file */
	int in_values_int;

	/** use a name to find trait in .yvm file */
	string in_values_name;

	/** name for the trait used for Y-values */
	string y_value_name;

	/** true = use extra file for covariables */
	bool cov_extra_file;

	/** path + name for the .cov file */
	string cov_file_name;

//+++++++++++++
// Parameters describing the Data

	/** if Phenotype is affection (case-control) or quantitative.
	* setting is determined by MData::checkYValues() */
	bool affection_status_phenotype;

	/** the number of covariables */
	int covariables;

	/** 1 = Recessive, 2 = Additive, 3 = Dominant */
	int genetic_model;

	/** which number repressents missing phenotypes */
	double missing_phenotype_code;
        /** which number is coding the controls */
	int control_value; 
	int case_value;

//++++++++++++
// Output:

	/** path + name for the output files .log ... */
	string out_file_name;
        string singlefile;
	/** true = no output on screen */
	bool silent;

	/** true = all model selection steps are written to logfile */
	bool detailed_selction;

//++++++++++++
//Imputation:

	/** the Neighbourhood of SNPs considered for Imputation */
	int imp_Neighbours_No;

	/** the closed matching SNPs considered for Imputation */
	int imp_Best_SNPs_No;

	/** if the Genom-Data is already imputaded */
	bool imp_is_imputated;

//++++++++++++
//ModelSelection:

	/** Number of Expected Causal SNPs, used in Model-Selction Criterias */
	int ms_ExpectedCausalSNPs;
	/** Number of Expected SNPs, in the preselection _*/
	int expected_causal_snps1;
        int maximalModelSize;
	/* *the maximal number of SNPs added in the Multi-Forward-Step */
	int ms_MaximalSNPsMultiForwardStep;

	/** the threshold below which SNPs with lower p-Values are considered for adding in the Forward Step */
	double ms_MaximalPValueForwardStep;
        /** a seperate variable for the new multiforward step */
	int ms_forward_step_max;
	bool ms_FastMultipleForwardStep;
	
        int     PValueBorder;
	int     reset;
	int	jump_back;

//++++++++++++
// Parameters for Affection/Case-Control

	/** use Pearson-Chi-Square for Single Marker Test for affection */
	bool cc_SingleMarkerTest_ChiSq;

	/** use Cochran-Armitage trend test (CATT) for Single Marker Test for affection */
	bool cc_SingleMarkerTest_CATT;

// Logistic Regression Control Parameters
	int logrC_maxit;
	int logrC_maxhs;
	int logrC_maxstep;
	double logrC_lconv;
	double logrC_gconv;
	double logrC_xconv;

// Genetic algorithm
  int modelsNo;
  int maxNoProgressIter;
  double pCross;
  double pMutation;
  int tournamentSize;
  double correlationThreshold;
  int correlationRange;
  
	
// TEST 
	int test;
	bool binary;
	/** upper limit for beta s*/
	double betaup;
	/** lower limit for beta s*/
	double inter;
	double betadown;
	double beta;
	int  nSNP;
	int  replications; 
	vector<string> SNPs;
        vector<string> SNPs1;
        vector<string> SNPs2;
//++++++++++++

	/** Constructor */
	Parameter ();

	/** Initialise */
	void setParameters ( const int argn, const char* argv[] );
};

/** A globally shared parameter object */
extern Parameter parameter;

#endif
