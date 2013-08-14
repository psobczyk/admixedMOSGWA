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

#include "DoubleList.hpp"
#include "Hdf5TestSuite.hpp"
#include <string>

using namespace std;
using namespace unitpp;
using namespace hdf5;

namespace test {

	/** Tests the class {@link hdf5::DoubleList}. */
	struct DoubleListTest : public Hdf5TestSuite {

		/** Construct the test object. */
		DoubleListTest ();

		/** Run the test. */
		void test ();

	} * doubleListTest = new DoubleListTest();	// automatically freed by unit++

	DoubleListTest::DoubleListTest () : Hdf5TestSuite( "DoubleListTest" ) {
		addTestMethod( "DoubleListTest::test", this, &DoubleListTest::test );
	}

	void DoubleListTest::test () {
		const string testFilename = setUp();

		// begin a block in order to avoid tearDown() before object release
		{
			File file( testFilename );

			DoubleList double1dim( file, "/double1dim" );
			assert_eq( "double1dim.countItems", 2, double1dim.countItems() );
			assert_eq( "double1dim.countDimensions", 2, double1dim.countDimensions() );
			double data[3];
			data[2] = 17.0;
			double1dim.readAll( data );
			assert_eq( "double1dim data[0]", 0.0, data[0] );
			assert_eq( "double1dim data[1]", 0.1, data[1] );
			assert_eq( "double1dim data after unchanged", 17.0, data[2] );
		}

		tearDown( testFilename );
	}

}
