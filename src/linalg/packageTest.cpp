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

#include "package.hpp"
#include "../TestSuite.hpp"

using namespace unitpp;
using namespace linalg;

namespace test {

	/** Tests the package commons. */
	class PackageTest : public TestSuite {

		private:

		/** Test function {@link upperPowerOf2}. */
		void testUpperPowerOf2 () {
			assert_eq( "0", 0, upperPowerOf2( 0 ) );
			assert_eq( "1", 1, upperPowerOf2( 1 ) );
			assert_eq( "2", 2, upperPowerOf2( 2 ) );
			assert_eq( "3", 4, upperPowerOf2( 3 ) );
			assert_eq( "4", 4, upperPowerOf2( 4 ) );
			assert_eq( "5", 8, upperPowerOf2( 5 ) );
			assert_eq( "6", 8, upperPowerOf2( 6 ) );
			assert_eq( "7", 8, upperPowerOf2( 7 ) );
			assert_eq( "8", 8, upperPowerOf2( 8 ) );
			assert_eq( "9", 16, upperPowerOf2( 9 ) );
			assert_eq( "15", 16, upperPowerOf2( 15 ) );
			assert_eq( "1023", 1024, upperPowerOf2( 1023 ) );
			assert_eq( "1024", 1024, upperPowerOf2( 1024 ) );
			assert_eq( "1025", 2048, upperPowerOf2( 1025 ) );
		}

		public:

		PackageTest () : TestSuite( "Test linalg package" ) {
			addTestMethod( "testUpperPowerOf2", this, &PackageTest::testUpperPowerOf2 );
		}
	} * packageTest = new PackageTest();	// automatically freed by unit++

}
