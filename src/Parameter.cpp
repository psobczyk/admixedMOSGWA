#include "Parameter.hpp"
#include <limits>
#include <stdlib.h>
#include <sstream>  //for stringstream
/** Internally declares the configuration variables. */
Parameter::Parameter () {

	// input file settings
	declare( "input", "plink_files", in_files_plink );
	declare( "input", "plink_bed_file", in_files_plink_bed );
	declare( "input", "plink_fam_file", in_files_plink_fam );
	declare( "input", "plink_bim_file", in_files_plink_bim );
	declare( "input", "hdf5_file", in_file_hdf5 );
	declare( "input", "use_extra_yvalues", y_value_extra_file );
	declare( "input", "extra_yvalues_file", in_files_values_yvm );
	declare( "input", "use_extra_covariables", cov_extra_file );
	declare( "input", "extra_covariables_file", cov_file_name );
	declare( "input", "control_value",control_value); //the control are normally 0
	declare( "input", "case_value",case_value); //the case are normally 1

	// input data settings
	declare( "data", "trait_name", y_value_name );
	declare( "data", "trait_position_in_yvm", in_values_int );
	declare( "data", "trait_name_in_yvm", in_values_name );
	declare( "data", "covariables", covariables );
	{
		map< const string, int > choice;
		choice[ "recessive" ] = 1;
		choice[ "additive" ] = 2;
		choice[ "dominant" ] = 3;
		declare( "data", "genetic_model", genetic_model, choice );
	}

	// output file settings
	declare( "output", "files", out_file_name );
	declare( "output", "singlefile",singlefile);//for creating the Hlasso genotyp file only once


	// log settings
	declare( "log", "silent", silent );
	declare( "log", "detailed", detailed_selction );

	// imputation settings
	declare( "imputation", "data_is_imputated",  imp_is_imputated );
	declare( "imputation", "consider_n_neighbours", imp_Neighbours_No );
	declare( "imputation", "consider_n_close_snps",	 imp_Best_SNPs_No );

	// single marker test settings
	declare( "single_marker", "chi_square", cc_SingleMarkerTest_ChiSq );
	declare( "single_marker", "cochran_armitage", cc_SingleMarkerTest_CATT );

	// model selection settings
	declare( "model_selection", "expected_causal_snps", ms_ExpectedCausalSNPs );
	declare( "model_selection", "expected_causal_snps1",expected_causal_snps1 );
	declare( "model_selection", "maximalModelSize",maximalModelSize );
	declare( "model_selection", "multi_forward_step_max", ms_MaximalSNPsMultiForwardStep );
	declare( "model_selection", "multi_forward_pvalue_max", ms_MaximalPValueForwardStep );
        declare( "model_selection", "forward_step_max",ms_forward_step_max);
	declare( "model_selection", "fast_multi_forward",  ms_FastMultipleForwardStep );
        declare( "model_selection","PValueBorder",PValueBorder);
	declare( "model_selection","reset",reset);
        declare( "model_selection","jump_back",jump_back);

	// log_regression settings
	declare( "log_regression", "max_it", logrC_maxit );
	declare( "log_regression", "max_hs", logrC_maxhs );
	declare( "log_regression", "max_step", logrC_maxstep );
	declare( "log_regression", "l_conv", logrC_lconv );
	declare( "log_regression", "g_conv", logrC_gconv );
	declare( "log_regression", "x_conv", logrC_xconv );

// genetics algorithm settings
  declare( "genetic_algorithm", "modelsNo", modelsNo);
  declare( "genetic_algorithm", "maxNoProgressIter", maxNoProgressIter);
  declare( "genetic_algorithm", "pCross", pCross);
  declare( "genetic_algorithm", "pMutation", pMutation);
  declare( "genetic_algorithm", "tournamentSize", tournamentSize);
  declare( "genetic_algorithm", "correlationThreshold", correlationThreshold);
  declare( "genetic_algorithm", "correlationRange", correlationRange);
  
//Erichs testcase gearatator  
declare( "TESTING", "test", test );
declare( "TESTING", "binary", binary);
declare( "TESTING", "betaup", betaup ); //setzt das beta uniform für alle SNP
declare( "TESTING", "inter", inter ); 
declare( "TESTING", "betadown", betadown ); //setzt das beta uniform für alle SNP

declare( "TESTING", "beta", beta ); //setzt das beta uniform für alle SNP 

declare( "TESTING", "nSNP", nSNP ); //setzt die Anzahl der kausalen SNP

declare( "TESTING", "replications", replications ); //how many Y vectors should be produced  
declare ("TESTING", "SNPs", SNPs); //which SNP to select
declare ("TESTING", "SNPs1", SNPs1); //which SNP to select
declare ("TESTING", "SNPs2", SNPs2); //which SNP to select
}

/** Initialise the parameter holder from the files given on the command line. */
void Parameter::setParameters ( const int argn, const char* argv[] ) {
   if (1==argn)
   {cerr<<"You need to provide an input file to use MOSGWA!"<<endl;
           exit(1); //or an better exit code
   }

	missing_phenotype_code = 10000;// numeric_limits<double>::quiet_NaN();

	// Parse command line arguments as config files
	for ( int i = 1; i < argn; ++i ) {
		ifstream configuration;
		configuration.exceptions ( ifstream::badbit );
		try {
			configuration.open( argv[i], ifstream::in );
			// See http://gehrcke.de/2011/06/reading-files-in-c-using-ifstream-dealing-correctly-with-badbit-failbit-eofbit-and-perror/
			if ( ! configuration.is_open() ) {
				cerr << "Cannot open configuration file \"" << argv[i] << "\"." << endl;
				exit( 255 );
			} else {
				parse( configuration, argv[i] );
			}
			configuration.close();
		} catch ( ifstream::failure e ) {
			cerr << "Cannot read configuration file \"" << argv[i] << "\": " << e.what() << endl;
			exit( 255 );
		}
	}
	cerr<<y_value_extra_file<<","<<y_value_name.empty()<<","<<in_values_int<<endl;
//this creates an unique file for all yvm and of course also in the version without 
//this has to be done before the number is added in the yvm case	
singlefile=out_file_name; //this is for simulators which need in every run the genotyp file 
if(y_value_extra_file)
	if(!y_value_name.empty())
             out_file_name=out_file_name+y_value_name;
	else if (0< in_values_int)
	{    stringstream number;
	     number<<in_values_int;
	     out_file_name=out_file_name+number.str();
	     //cerr<<out_file_name<<endl;
	}

	if ( in_file_hdf5.empty() ) {
		//  parameter.out_file_name=parameter.out_file_name+parameter.y_value_name;			
		if ( in_files_plink_bed.empty() ) {
			in_files_plink_bed = in_files_plink+".bed";
		}
		if ( in_files_plink_fam.empty() ) {
			in_files_plink_fam = in_files_plink+".fam";
		}
		if ( in_files_plink_bim.empty() ) {
			in_files_plink_bim = in_files_plink+".bim";
		}
		if ( y_value_extra_file and in_files_values_yvm.empty() ) {
			in_files_values_yvm = in_files_plink+".yvm";
		}
		if ( cov_extra_file and cov_file_name.empty() ) {
			cov_file_name = in_files_plink+".cov";
		}
	}

	//forward step settings
        if(0==ms_MaximalSNPsMultiForwardStep)
                ms_MaximalSNPsMultiForwardStep=1;

        if (0==ms_forward_step_max)
	    	ms_forward_step_max=ms_MaximalSNPsMultiForwardStep; //saveguard against a 0 forward step
//standart settings for logistic regression;

	if ( 0==logrC_maxit)
	    	logrC_maxit=15; 
	if (0==logrC_maxhs)
		logrC_maxhs=3;
        if (0==logrC_maxstep)
		logrC_maxstep=3;
	if (0==logrC_lconv)
		logrC_lconv=10e-3;
	if (0==logrC_gconv) 
		logrC_gconv=10e-3;  
        if (0==logrC_xconv)
		logrC_xconv=10e-3;

        if (0==expected_causal_snps1)
        	expected_causal_snps1=2; //MBIC2 as the standart
	if (0==maximalModelSize)
		maximalModelSize=35;
	if(0==control_value)
	;	//this is ok
        if (0==case_value)
		case_value=1; //this is the normal setting
	if (case_value==control_value)
	{cerr<<"ERROR case und control have the same value";
	 exit(2);}	
	      
}
