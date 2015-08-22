/********************************************************************************
 *	This file is part of the MOSGWA program code.				*
 *	Copyright ©2015, Bernhard Bodenstorfer.					*
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

#include "MultiLogger.hpp"
#include <cstddef>	// NULL

using namespace std;

namespace logging {

	void MultiLogger::add ( Logger& logger ) {
		Logger* pointer = &logger;
		if ( NULL == pointer ) {
			error( "Attempt to add a null pointer as logger. Ignored." );
		} else {
			destinations.push_back( pointer );
		}
	}

	void MultiLogger::setThreshold ( const Severity threshold ) {
		Logger::setThreshold( threshold );
		for (
			vector<Logger*>::iterator iterator = destinations.begin();
			iterator != destinations.end();
			++iterator
		) {
			Logger* pointer = *iterator;
			pointer->setThreshold( threshold );
		}
	}

	void MultiLogger::log ( const Severity level, const char * const format, va_list arguments ) {
		if ( getThreshold() <= level ) {
			for (
				vector<Logger*>::iterator iterator = destinations.begin();
				iterator != destinations.end();
				++iterator
			) {
				Logger* pointer = *iterator;
				pointer->log( level, format, arguments );
			}
		}
	}

	void MultiLogger::write ( const char * const text ) {
		for (
			vector<Logger*>::iterator iterator = destinations.begin();
			iterator != destinations.end();
			++iterator
		) {
			Logger* pointer = *iterator;
			pointer->write( text );
		}
	}

}
