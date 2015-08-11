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

#ifndef LOG_FILELOGGER_HPP
#define LOG_FILELOGGER_HPP

#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>

#include "Logger.hpp"

namespace log {

	/** Logs to a file. */
	class FileLogger : public Logger {

		std::ofstream file;

		public:

		/** Opens file for appending logs. */
		FileLogger ( const char * const filename );

		/** Destructor closes the file. */
		virtual ~FileLogger ();

		protected:

		/** Implement this abstract facility to actually log a line. */
		virtual void write ( const char * const text );
	};

}

#endif
