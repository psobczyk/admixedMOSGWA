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

#include "Exception.hpp"
#include "TestSuite.hpp"

using namespace std;
using namespace unitpp;

namespace test {

	/** Tests the class {@link Exception}. */
	class ExceptionTest : public TestSuite {

		/** Test {@link Exception::Exception} and {@link Exception::getMessage}. */
		void testMessaging ();

		public:

		ExceptionTest () : TestSuite( "ExceptionTest" ) {
			addTestMethod( "ExceptionTest::testMessaging", this, &ExceptionTest::testMessaging );
		}
	} * exceptionTest = new ExceptionTest();	// automatically freed by unit++

	void ExceptionTest::testMessaging () {
		Exception exception( "Want to report error %d in file \"%s\".", 42, "ExceptionTest" );
		string message = exception.what();
		assert_eq( "Exception message", "Want to report error 42 in file \"ExceptionTest\".", message );
	}
}
