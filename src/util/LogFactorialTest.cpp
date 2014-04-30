/********************************************************************************
 *	This file is part of the MOSGWA program code.				*
 *	Copyright ©2014, Bernhard Bodenstorfer.					*
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

#include "LogFactorial.hpp"
#include "../TestSuite.hpp"
#include <cstdlib>
#include <cstdio>
#include <cmath>

using namespace unitpp;
using namespace util;

namespace test {

	/** Tests the template class {@link util::LogFactorial}. */
	struct LogFactorialTest : public TestSuite {

		LogFactorialTest ();

		void testLogFactorial ();
		void testLogChoose ();

	} * logFactorialTest = new LogFactorialTest();	// automatically freed by unit++

	LogFactorialTest::LogFactorialTest () : TestSuite( "util::LogFactorial Test" ) {
		addTestMethod( "LogFactorialTest::testLogFactorial", this, &LogFactorialTest::testLogFactorial );
		addTestMethod( "LogFactorialTest::testLogChoose", this, &LogFactorialTest::testLogChoose );
	}

	void LogFactorialTest::testLogFactorial () {
		const size_t max = 1000;
		double f[max];
		f[0] = 0;
		for ( size_t i = 1; i < max; ++i ) {
			f[i] = f[i-1] + log( i );
		}
		{
			LogFactorial logFactorial;
			for ( size_t i = 0; i < max; ++i ) {
				char message[ 2 + max ];
				snprintf( message, sizeof( message ), "%u", i );
				assert_close( message, f[i], logFactorial.logFactorial( i ) );
			}
		}
		{
			LogFactorial logFactorial;
			assert_close( "last", f[max-1], logFactorial.logFactorial( max-1 ) );
		}
	}

	/** Test against a Pascal triangle. */
	void LogFactorialTest::testLogChoose () {
		const size_t max = 100;
		double choose[max][max] = { 1.0 };
		LogFactorial logFactorial;
		for ( size_t i = 1; i < max; ++i ) {
			for ( size_t j = 0; j <= i; ++j ) {
				double x = 0.0;
				if ( 0 < j ) {
					x += choose[i-1][j-1];
				}
				if ( j < i ) {
					x += choose[i-1][j];
				}
				choose[i][j] = x;
				{
					char message[ 5 + max ];
					snprintf( message, sizeof( message ), "%u %u", i, j );
					assert_close(
						message,
						log( choose[i][j] ),
						logFactorial.logChoose( i, j ),
						1e-8
					);
				}
			}
		}
	}

}
