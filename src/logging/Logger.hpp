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

#include <cstdarg>

/** A very simple abstract logging interface.
* @author Bernhard Bodenstorfer
*/
namespace logging {

	/** Abstract interface to logging. */
	class Logger {

		friend class MultiLogger;

		public:

		/** Message severity encoding.
		* The numeric values may change, only the ascending order is part of the contract.
		*/
		typedef enum { DEBUG = 0, INFO = 1, WARNING = 2, ERROR = 3 } Severity;

		/** Construct with default severity limit <code>INFO</code> */
		Logger ();

		/** Limit from below the severity of messages to be actually logged.
		* @see getThreshold
		*/
		void setThreshold ( const Severity threshold );

		/** Query the lower severity limit for messages to be actually logged.
		* @see setThreshold
		*/
		Severity getThreshold () const;

		/** Convenience method to {@link log} with severity {@link DEBUG}.
		* Use in analogy to <code>printf</code>.
		* @see log regarding security advice
		*/
		void debug ( const char * const format, ... );

		/** Convenience method to {@link log} with severity {@link INFO}.
		* Use in analogy to <code>printf</code>.
		* @see log regarding security advice
		*/
		void info ( const char * const format, ... );

		/** Convenience method to {@link log} with severity {@link WARNING}.
		* Use in analogy to <code>printf</code>.
		* @see log regarding security advice
		*/
		void warning ( const char * const format, ... );

		/** Convenience method to {@link log} with severity {@link ERROR}.
		* Use in analogy to <code>printf</code>.
		* @see log regarding security advice
		*/
		void error ( const char * const format, ... );

		/** Destructor frees any internal resources. */
		virtual ~Logger ();

		protected:

		/** Push a log message, use format similar as in <code>vprintf</code>.
		* Security advice: never log using code like
		* <code>logger.log( Logger::INFO, variablemessage )</code>
		* unless you have full control over <code>variablemessage</code>,
		* because the second argument will be interpreted as format string.
		* Rather use code like
		* <code>logger.log( Logger::INFO, "%s", variablemessage )</code>
		* for the purpose of logging potentially untrusted strings.
		* TODO<BB>: replace printf style by two stream references (to a dead-end for level<threshold)
		*/
		virtual void log ( const Severity level, const char * const format, va_list arguments );

		/** Implement this abstract facility to actually log a line. */
		virtual void write ( const char * const text ) = 0;

		private:

		/** The current minimum log level. */
		Severity threshold;
	};

	/** To be initialised (in {@link main}) to point to a globally shared logger object. */
	extern Logger* logger;

}

#endif
