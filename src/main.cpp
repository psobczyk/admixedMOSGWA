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

#include "buildinfo.hpp"
#include "Parameter.hpp"
#include "Exception.hpp"
#include "logging/FileLogger.hpp"
#include "search/egreedy/ElaboratedGreedy.hpp"
#include "search/genetica/GeneticAlgorithm.hpp"
#include <string>

using namespace std;
using namespace logging;

/** Application entry point */
int main ( const int argc, const char *argv[] ) {

	// set parameters
	parameter.reset( new Parameter );
	parameter->setParameters( argc, argv );

	// init logging
	const string logFileName( parameter->out_file_name + ".log" );
	logger.reset( new FileLogger( logFileName.c_str() ) );
	switch( parameter->log_level ) {
		case 0: logger->setLimit( Logger::DEBUG );
			break;
		case 1: logger->setLimit( Logger::INFO );
			break;
		case 2: logger->setLimit( Logger::WARNING );
			break;
		case 3: logger->setLimit( Logger::ERROR );
			break;
	}

	// log start
	logger->info( " _____ _____ _____ _____ _ _ _ _____ " );
	logger->info( "|     |     |   __|   __| | | |  _  |\tModel Selection" );
	logger->info( "| | | |  |  |__   |  |  | | | |     |\tfor Genome-wide" );
	logger->info( "|_|_|_|_____|_____|_____|_____|__|__|\tAssociations" );
	logger->info( "release %s, %s, called with arguments:", buildinfo::version, buildinfo::timestamp );
	for ( int i = 0; i < argc; ++i ) {
		logger->info( "%s", argv[i] );
	}

	try {
		if ( Parameter::searchStrategy_genetic_algorithm == parameter->searchStrategy ) {
			genetica::GeneticAlgorithm artur;
			artur.run();
		} else {
			egreedy::ElaboratedGreedy erich;
			erich.run();
		}
	} catch ( const Exception e ) {
		logger->error( "%s", e.what() );
	}
	logger->info( "Close logs." );
}
