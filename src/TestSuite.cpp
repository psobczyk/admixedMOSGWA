/********************************************************************************
 *	This file is part of the MOSGWA program code.				*
 *	Copyright ©2012–2013, Bernhard Bodenstorfer.				*
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

#include "TestSuite.hpp"
#include <cmath>
#include <cstdio>

using namespace std;
using namespace unitpp;

namespace test {

	TestSuite::TestSuite ( const string &title ) : suite ( title ) {
		suite::main().add( title, this );
	}

	void TestSuite::assert_close(
		const string& msg,
		const double expected,
		const double actual,
		const double tolerance
	) {
		char buffer[1024];
		sprintf(
			buffer,
			" expected close to %f and got %f, difference was %e, tolerance %e",
			expected,
			actual,
			actual-expected,
			tolerance
		);
		string xmsg = msg + buffer;
		assert_true( xmsg, fabs( expected - actual ) < tolerance );
	}

}
