/********************************************************************************
 *	This file is part of the MOSGWA program code.				*
 *	Copyright ©2011–2013, Artur Gola, Erich Dolejsi, Bernhard Bodenstorfer.	*
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
#include "Helpfull.hpp"
#include "Log.hpp"
#include "MData.hpp"
#include "Model.hpp"
#include "GA.hpp"
#include <iomanip>

/** global variable for logfile */
ofstream LOG;

/** Global parameter holder */
Parameter parameter;


/** Application entry point */
int main ( const int argc, const char *argv[] ) {

	// print logo on screen
	printStartScreen();
	//return 0;
	/*
	const int argcA = 2;
	const char * argvA[argcA] = {"MOSGWA_GA",
                               "HeadRD.conf"
                              };
*/                              
//  for (int i = 0; i < argcA; i++)
//    cout << setw(2) << i + 1 << "] " << argvA[i] << endl;									
	// set parameters
	parameter.setParameters( argc, argv );

	// init logging
	string logFileName( parameter.out_file_name + ".log" );

	try {
		printLOG( "Start: open log to file \"" + logFileName + "\"" );
		LOG.open( logFileName.c_str(), fstream::out );
	} catch ( ofstream::failure e ) {
		cerr << "Could not open logfile \"" + logFileName + "\"" <<endl;
	}

	// checks if logfile can be written.
	LOG.exceptions ( ofstream::eofbit | ofstream::failbit | ofstream::badbit );

	// Read in data by generating MData object
/*  MData data;

// imputate data if nessessary
	if ( ! parameter.imp_is_imputated ) {
		data.imputateMissingValues();
		data.writeBEDfile();
	}
  //data.calculateIndividualTests();
  
	Model model_1(data);
  model_1.computeRegression();
  cout << endl << "Model 1:" << endl;
  cout << "0 snp: " << setprecision(12) << model_1.computeMSC() << endl;
  model_1.addSNPtoModel(1);
  cout << "snp <1>: " << setprecision(12) << model_1.computeMSC() << endl;
  model_1.addSNPtoModel(250);
  cout << "snp <1, 250>: " << setprecision(12) << model_1.computeMSC() << endl;
  cout << endl << "Model 2:" << endl;
  Model model_2(data);
  model_2.computeRegression();
  cout << "0 snp: " << setprecision(12) << model_2.computeMSC() << endl;
  model_2.addSNPtoModel(250);
  cout << "snp <250>: " << setprecision(12) << model_2.computeMSC() << endl;
  model_2.addSNPtoModel(1);
  cout << "snp <250, 1>: " << setprecision(12) << model_2.computeMSC() << endl;
*/	
	
  unsigned int modelsNo_ = parameter.modelsNo;                   // 5
  unsigned int maxNoProgressIter_ = parameter.maxNoProgressIter; //1000;
  double pCross_ = parameter.pCross;                             // 0.65;
  double pMutation_ = parameter.pMutation;                       // 0.05;
  unsigned int tournamentSize_ = parameter.tournamentSize;       //2;
  double correlationThreshold_ = parameter.correlationThreshold; //0.7;
  int correlationRange_ = parameter.correlationRange;            // 500
  /*
  cout << "modelsNo_ = " << parameter.modelsNo << endl
       << "maxNoProgressIter_ = " << parameter.maxNoProgressIter << endl
       << "pCross_ = " << parameter.pCross << endl
       << "pMutation_ = " << parameter.pMutation  << endl
       << "tournamentSize_ = " << parameter.tournamentSize << endl
       << "correlationThreshold_ = " << parameter.correlationThreshold << endl
       << "int correlationRange = " << parameter.correlationRange << endl;;
  exit(1);
  */
  GA ga(modelsNo_, maxNoProgressIter_, pCross_, pMutation_, tournamentSize_, correlationThreshold_, correlationRange_);
//  ga.testRecombination();
  ga.run();
//  ga.tests(); 
  


	// calculate single marker test (need for modelselection)
//	data.calculateIndividualTests();

	// complete model selection process
	//data.selectModel();

	try {
		printLOG("End");
		LOG.close();
	} catch ( ofstream::failure e ) {
		cerr << "Could not close logfile \"" + logFileName + "\"" << endl;
	}
}




/** Application entry point */
/*
int main ( const int argc, const char *argv[] ) {

	// print logo on screen
	printStartScreen();

	// set parameters
	parameter.setParameters( argc, argv );

	// init logging
	string logFileName( parameter.out_file_name + ".log" );

	try {
		printLOG( "Start: open log to file \"" + logFileName + "\"" );
   		LOG.open( logFileName.c_str(), fstream::out );
	} catch ( ofstream::failure e ) {
		cerr << "Could not open logfile \"" + logFileName + "\"" <<endl;
	}

	// checks if logfile can be written.
	LOG.exceptions ( ofstream::eofbit | ofstream::failbit | ofstream::badbit );

	// Read in data by generating MData object
	MData data;

	// imputate data if nessessary
	if ( ! parameter.imp_is_imputated ) {
		data.imputateMissingValues();
		data.writeBEDfile();
	}

	// calculate single marker test (need for modelselection)
	data.calculateIndividualTests();

	// complete model selection process
	data.selectModel();

	try {
		printLOG("End");
		LOG.close();
	} catch ( ofstream::failure e ) {
		cerr << "Could not close logfile \"" + logFileName + "\"" << endl;
	}
}
*/
