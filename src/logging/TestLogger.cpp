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

#include "TestLogger.hpp"
#include <cctype>
#include <string>

using namespace std;
using namespace unitpp;
using namespace logging;

namespace test {

	TestLogger::TestLogger () : suite( "TestLogger" ), expected( NULL ) {
	}

	void TestLogger::setExpected ( const char * const expected ) {
		this->expected = expected;
	}

	void TestLogger::write ( const char * const text ) {
		if ( NULL == expected ) {
			assert_fail( "Expected not to forward log" );
		} else {
			const char * cursor = text;
			while( isdigit( *cursor ) ) {
				++cursor;
			}
			assert_true( ">= 4 digits year expected", cursor - text >= 4 );
			assert_eq( "year-month separator", '-', *cursor++ );
			assert_true( "month 10", isdigit( *cursor++ ) );
			assert_true( "month 1", isdigit( *cursor++ ) );
			assert_eq( "month-day separator", '-', *cursor++ );
			assert_true( "day 10", isdigit( *cursor++ ) );
			assert_true( "day 1", isdigit( *cursor++ ) );
			assert_eq( "day time separator", ' ', *cursor++ );
			assert_true( "hour 10", isdigit( *cursor++ ) );
			assert_true( "hour 1", isdigit( *cursor++ ) );
			assert_eq( "hour:minute separator", ':', *cursor++ );
			assert_true( "minute 10", isdigit( *cursor++ ) );
			assert_true( "minute 1", isdigit( *cursor++ ) );
			assert_eq( "minute:second separator", ':', *cursor++ );
			assert_true( "second 10", isdigit( *cursor++ ) );
			assert_true( "second 1", isdigit( *cursor++ ) );
			assert_eq( "timestamp text separator", '\t', *cursor++ );
			assert_eq( "log text", string( expected ), string( cursor ) );
		}
	}

}
