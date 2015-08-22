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

#ifndef LOG_MULTILOGGER_HPP
#define LOG_MULTILOGGER_HPP

#include "Logger.hpp"
#include <vector>

/** A very simple abstract logging interface.
* @author Bernhard Bodenstorfer
*/
namespace logging {

	/** Forwards log output to several {@link Logger}s downstream. */
	class MultiLogger : public Logger {

		std::vector<Logger*> destinations;

		public:

		/** Add logger as destination.
		* @param logger to add,
		* it must remain valid, i.e. must not be destructed as long as logs are expected.
		*/
		void add ( Logger& logger );

		/** Set threshold for <code>this</code> and all yet connected destinations.
		* Note that destination loggers may afterwards be set to a different limit.
		* However, it does not make sense to be less restrictive than <code>this</code>,
		* because <code>this</code> will not forward messages below threshold.
		*/
		void setThreshold ( const Severity threshold );

		protected:

		/** Distribute a log message. */
		virtual void log ( const Severity level, const char * const format, va_list arguments );

		/** Distribute a text. This is actually circumvented by {@link log}. */
		virtual void write ( const char * const text );
	};

}

#endif
