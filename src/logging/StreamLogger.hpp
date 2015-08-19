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

#ifndef LOG_STREAMLOGGER_HPP
#define LOG_STREAMLOGGER_HPP

#include <iostream>
#include <iomanip>
#include <string>

#include "Logger.hpp"

namespace logging {

	/** Logs to a C++ output stream. */
	class StreamLogger : public Logger {

		std::ostream * stream;

		public:

		/** Connect to stream for appending logs. */
		StreamLogger ( std::ostream& stream = std::cerr );

		/** Destructor. */
		virtual ~StreamLogger ();

		protected:

		/** Change the stream to log to. */
		void setStream ( std::ostream& stream = std::cerr );

		/** Implement this abstract facility to actually log a line. */
		virtual void write ( const char * const text );
	};

}

#endif
