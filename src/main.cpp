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
		//modelin->printModel("we are before select model",3);//3 is mBIC
		data.selectModel(
			modelin,
			parameter.PValueBorder,
			parameter.maximalModelSize,
			Parameter::selectionCriterium_mBIC
		);
		modelin->printModel("first result");
		data.selectModel(
			modelin,
			5000
		);
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
