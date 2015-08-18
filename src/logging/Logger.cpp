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

#include "Logger.hpp"
#include "../Exception.hpp"	// MAX_MSG_SIZE
#include <cstdio>
#include <ctime>

using namespace std;

namespace logging {

	Logger::Logger () : minimumLevel( INFO ) {
	}

	void Logger::setLimit ( const Severity minimumLevel ) {
		this->minimumLevel = minimumLevel;
	}

	Logger::Severity Logger::getLimit () const {
		return minimumLevel;
	}

	void Logger::log ( const Severity level, const char * const format, va_list arguments ) {
		if ( level >= minimumLevel ) {
			const char* levelText;
			char
				buffer[MAX_MSG_SIZE],
				*cursor;
			const time_t t = time(NULL);
			struct tm time;
			localtime_r( &t, &time );
			switch( level ) {
				case DEBUG: levelText = "DEBUG";
					break;
				case INFO: levelText = "INFO";
					break;
				case WARNING: levelText = "WARNING";
					break;
				case ERROR: levelText = "ERROR";
					break;
				default: levelText = "FATAL";
			}
			const size_t timestampNeeds = snprintf(
				buffer,
				MAX_MSG_SIZE,
				"%04u-%02u-%02u %02u:%02u:%02u\t%s\t",
				1900 + time.tm_year,
				1 + time.tm_mon,
				time.tm_mday,
				time.tm_hour,
				time.tm_min,
				time.tm_sec,
				levelText
			);
			if ( MAX_MSG_SIZE <= timestampNeeds ) {
				write( "logging buffer used up by timestamp" ); // should not occur in near future
			} else {
				vsnprintf(
					buffer + timestampNeeds,
					MAX_MSG_SIZE - timestampNeeds,
					format,
					arguments
				);
				write( buffer );
			}
		}
	}

	void Logger::debug ( const char * const format, ... ) {
		va_list arguments;
		va_start( arguments, format );
		log( DEBUG, format, arguments );
		va_end( arguments );
	}

	void Logger::info ( const char * const format, ... ) {
		va_list arguments;
		va_start( arguments, format );
		log( INFO, format, arguments );
		va_end( arguments );
	}

	void Logger::warning ( const char * const format, ... ) {
		va_list arguments;
		va_start( arguments, format );
		log( WARNING, format, arguments );
		va_end( arguments );
	}

	void Logger::error ( const char * const format, ... ) {
		va_list arguments;
		va_start( arguments, format );
		log( ERROR, format, arguments );
		va_end( arguments );
	}

	Logger::~Logger () {}

}
