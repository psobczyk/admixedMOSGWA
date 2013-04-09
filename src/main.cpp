#include "Parameter.hpp"
#include "Helpfull.hpp"
#include "Log.hpp"
#include "MData.hpp"
#include "Model.hpp"
#include "GA.hpp"

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
