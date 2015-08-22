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
#include "../TestSuite.hpp"

namespace test {

	/** Implements {@link Logger} for testing purposes.
	* It compares each log message with the <code>expected</code> one.
	* If <code>NULL</code> is set as expected, an error would be raised if a message were "written".
	* Of course, this depends on the log level.
	*/
	class TestLogger : public logging::Logger, private unitpp::suite {
		const char * expected;
		public:
		TestLogger ();
		void setExpected ( const char * const expected );
		protected:
		virtual void write ( const char * const text );
	};

}
