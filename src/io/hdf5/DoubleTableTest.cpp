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

#include "DoubleTable.hpp"
#include "Hdf5TestSuite.hpp"
#include <string>

using namespace std;
using namespace unitpp;
using namespace hdf5;

namespace test {

	/** Tests the class {@link hdf5::DoubleTable}. */
	struct DoubleTableTest : public Hdf5TestSuite {

		/** Construct the test object. */
		DoubleTableTest ();

		/** Run the test. */
		void test ();

	} * doubleTableTest = new DoubleTableTest();	// automatically freed by unit++

	DoubleTableTest::DoubleTableTest () : Hdf5TestSuite( "DoubleTableTest" ) {
		addTestMethod( "DoubleTableTest::test", this, &DoubleTableTest::test );
	}

	void DoubleTableTest::test () {
		const string testFilename = setUp();

		// begin a block in order to avoid tearDown() before object release
		{
			File file( testFilename );

			DoubleTable double2dim( file, "/double2dim" );
			assert_eq( "double2dim.countItems", 6, double2dim.countItems() );
			assert_eq( "double2dim.countColumns", 2, double2dim.countColumns() );
			assert_eq( "double2dim.countRows", 3, double2dim.countRows() );
			double data[7];
			data[6] = 17.0;
			double2dim.readAll( data );
			assert_eq( "double2dim data[0]", 0.0, data[0] );
			assert_eq( "double2dim data[1]", 0.1, data[1] );
			assert_eq( "double2dim data[2]", 1.0, data[2] );
			assert_eq( "double2dim data[3]", 1.1, data[3] );
			assert_eq( "double2dim data[4]", 2.0, data[4] );
			assert_eq( "double2dim data[5]", 2.1, data[5] );
		}

		tearDown( testFilename );
	}

}
