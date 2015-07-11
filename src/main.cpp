/********************************************************************************
 *	This file is part of the MOSGWA program code.				*
 *	Copyright ©2011–2013, Erich Dolejsi, Bernhard Bodenstorfer.		*
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

#include "buildinfo.hpp"
#include "Parameter.hpp"
#include "Helpfull.hpp"
#include "Log.hpp"
#include "Exception.hpp"
#include "MData.hpp"
#include "Model.hpp"
#include "GA.hpp"
#include <string>

using namespace std;

/** global variable for logfile */
ofstream LOG;

/** Global parameter holder */
Parameter parameter;

void runGreedy ();
void runGA ( const int outNo = -1, const string &modelsFileName = "" );

string timeStamp ( const time_t t ) {
	int hour, min, sec;
	sec = t;
	hour = sec / 3600;
	sec = sec % 3600;
	min = sec / 60;
	sec = sec % 60;
	return int2strPadWith( hour, 2, '0' )
		+ ":" + int2strPadWith( min, 2, '0' )
		+ ":" + int2strPadWith( sec, 2, '0' );  
}

/** Application entry point */
int main ( const int argc, const char *argv[] ) {

	// print logo on screen
	printStartScreen();

	// set parameters
	parameter.setParameters( argc, argv );

	// init logging
	string logFileName( parameter.out_file_name + ".log" );

	try {
		LOG.open( logFileName.c_str(), fstream::out );
		printLOG( "Start: open log to file \"" + logFileName + "\"" );
		printLOG( "Program version: " + string( buildinfo::version ) );
		printLOG( "Program timestamp: " + string( buildinfo::timestamp ) );
		printLOG( "Program call arguments:" );
		for ( int i = 0; i < argc; ++i ) {
			printLOG( argv[i] );
		}
	} catch ( ofstream::failure e ) {
		cerr << "Could not open logfile \"" + logFileName + "\"" <<endl;
	}

	// checks if logfile can be written.
	LOG.exceptions ( ofstream::eofbit | ofstream::failbit | ofstream::badbit );

	try {
		if ( Parameter::searchStrategy_genetic_algorithm == parameter.searchStrategy ) {
			runGA();
		} else {
			runGreedy();
		}
	} catch ( const Exception e ) {
		printLOG( e.what() );
	}

	try {
		printLOG("end");
		LOG.close();
	} catch ( ofstream::failure e ) {
		cerr << "Could not close logfile \"" + logFileName + "\"" << endl;
	}
}

void runGreedy () {
	// Read in data by generating MData object
	MData data;

	// calculate single marker test (need for modelselection)
	data.calculateIndividualTests();

	// complete model selection process
	//OLD data.selectModel();

	Model model0( data );
	Model firstmodel(data);
	model0.computeRegression();
	data.setLL0M( model0.getMJC() );
	Model *modelin=&model0;
	data.selectModel(
		modelin,
		parameter.PValueBorder,
		parameter.maximalModelSize,
		Parameter::selectionCriterium_mBIC_firstRound
	);
	modelin->printModel( "First round result" );
	data.selectModel(
		modelin,
		5000,
		parameter.maximalModelSize,
		parameter.selectionCriterium
	);
}

/**
 * @brief Runs Genetics Algorithm as a model selection method.
 * @param outNo - it is number of output. It is used to run GA in a loop and save all the results in the different files.
 * @param modelsFileName - it is the file name which contains information about an initial population. 
 *                         If a modelsFileName is empty GA creates a new initial population
 * 
 */
void runGA(int outNo, const string &modelsFileName)
{
  cerr << "runGA" << endl;
  unsigned int modelsNo_ = parameter.modelsNo;
  unsigned int maxNoProgressIter_ = parameter.maxNoProgressIter; 
  double pCross_ = parameter.pCross;                             
  double pMutation_ = parameter.pMutation;                       
  unsigned int tournamentSize_ = parameter.tournamentSize;       
  double correlationThreshold_ = parameter.correlationThreshold; 
  int correlationRange_ = parameter.correlationRange;           
  double regionMinCorrelation_ = parameter.regionMinCorrelation;
  int B_ = parameter.B;
  string old_out_file_name;
  
  cout <<   "modelsNo_ " <<  parameter.modelsNo << endl
       << "maxNoProgressIter_ = " << parameter.maxNoProgressIter << endl
  << "pCross_ = " << parameter.pCross << endl
  << "pMutation_ = " << parameter.pMutation << endl                      
  << "tournamentSize_ = " << parameter.tournamentSize << endl
  << "correlationThreshold_ = " << parameter.correlationThreshold << endl
  << "correlationRange_ = " << parameter.correlationRange << endl
  << "regionMinCorrelation_ = " <<  parameter.regionMinCorrelation << endl
  << "B_ = " << parameter.B << endl
  << "causalModelFilename = \"" << parameter.causalModelFilename << "\"" << endl;
  
  if (outNo >= 0)  
  {
    old_out_file_name = parameter.out_file_name;
  }

  const time_t time_start = time(NULL);
  stringstream ss;
  timespec ts;
  clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &ts);
  time_t now;
  time(&now);
  srand(now);

  map<snp_index_t, int> mapSNPCausal_ind; 
  vector< multiset<long double> > tabCausalPost; 
  vector< multiset<long double> > tabCausalPost_b; 

  GA ga(modelsNo_, maxNoProgressIter_, pCross_, pMutation_, tournamentSize_, B_, modelsFileName, correlationThreshold_, correlationRange_, regionMinCorrelation_);
  ga.run();

  ss.str("");
  ss.clear();
  ss << "GA time: ";
  ga.writePoolToFile(ss);
  
  ga.initCausalPost( mapSNPCausal_ind );
  tabCausalPost.resize( mapSNPCausal_ind.size() );
  tabCausalPost_b.resize( mapSNPCausal_ind.size() );
  
  stringstream ssModels;
  ssModels << "";
  ga.computePosteriorProbability(ssModels, mapSNPCausal_ind, tabCausalPost, tabCausalPost_b);//, minPosterior);
  writePosterior((parameter.out_file_name + "_post_12.txt").c_str(), mapSNPCausal_ind, tabCausalPost, 1);

  const time_t time_end = time(NULL);          //get current calendar time
  ss << (time_end > time_start? sec2time(time_end - time_start) : sec2time(24 * 3600 + time_end - time_start));
  cout << ss.str() << endl;
  ga.writePoolToFile(ss);

  ofstream  Pmi_YsortFile;
  Pmi_YsortFile.exceptions ( ofstream::eofbit | ofstream::failbit | ofstream::badbit ); 
  try
  {
    Pmi_YsortFile.open( ( parameter.out_file_name + "_PjMi_YsortFile" /*+ int2str(parameter.in_values_int)*/ + ".txt" ).c_str(),  fstream::out | fstream::trunc );
    Pmi_YsortFile << ss.str() << endl << ssModels.str() << endl;
    Pmi_YsortFile.flush();
    Pmi_YsortFile.close();
  }
  catch (ofstream::failure e)
  {
    cerr << "Could not write PMi_Y-sorted File" <<endl;
  }
  cout << "Posterior probalibities of SNPs are in the file: " 
       <<  ( parameter.out_file_name + "_PjMi_YsortFile" /*+ int2str(parameter.in_values_int)*/ + ".txt" ).c_str() << endl;
  cout << "A report is in the file: " <<  ( parameter.out_file_name + "_recognized_SNPs" /*+ int2str(parameter.in_values_int)*/ + ".txt" ).c_str() << endl;  
  printLOG(ss.str());   
  if (outNo >= 0)  
  {
    parameter.out_file_name = old_out_file_name;
  }
}
