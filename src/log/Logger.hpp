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

#ifndef LOG_LOGGER_HPP
#define LOG_LOGGER_HPP

/** A very simple abstract logging interface.
* @author Bernhard Bodenstorfer
*/
namespace log {

	/** Abstract interface to logging. */
	class Logger {

		public:

		/** Message severity encoding. */
		typedef enum { INFO = 0, WARNING = 1, ERROR = 2 } Severity;

		/** Limit from below the severity of messages to be actually logged. */
		void setLimit ( const Severity minimumLevel );

		/** Push a log message, use format similar as in <code>printf</code>. */
		virtual void log ( const Severity level, const char * const format, ... );

		/** Destructor frees any internal resources. */
		virtual ~Logger ();

		protected:

		/** Implement this abstract facility to actually log a line. */
		virtual void write ( const char * const text ) = 0;

		private:

		/** The current minimum log level. */
		Severity minimumLevel;
	};

}

#endif
