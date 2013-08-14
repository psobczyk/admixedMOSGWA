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

#include "StringList.hpp"
#include "Hdf5TestSuite.hpp"
#include <string>

using namespace std;
using namespace unitpp;
using namespace hdf5;

namespace test {

	/** Tests the class {@link hdf5::StringList}. */
	struct StringListTest : public Hdf5TestSuite {

		/** Construct the test object. */
		StringListTest ();

		/** Run the test. */
		void test ();

	} * stringListTest = new StringListTest();	// automatically freed by unit++

	StringListTest::StringListTest () : Hdf5TestSuite( "StringListTest" ) {
		addTestMethod( "StringListTest::test", this, &StringListTest::test );
	}

	void StringListTest::test () {
		const string testFilename = setUp();

		// begin a block in order to avoid tearDown() before object release
		{
			File file( testFilename );

			{
				StringList fixString1dim( file, "/fixString1dim" );
				assert_eq( "fixString1dim.countItems", 3, fixString1dim.countItems() );
				assert_eq( "fixString1dim.countDimensions", 3, fixString1dim.countDimensions() );
				string data[4];
				data[3] = string( "17" );
				fixString1dim.readAll( data );
				assert_eq( "fixString1dim data[0]", string( "," ), data[0] );
				assert_eq( "fixString1dim data[1]", string( ",1" ), data[1] );
				assert_eq( "fixString1dim data[2]", string( ",2" ), data[2] );
				assert_eq( "fixString1dim data after unchanged", string( "17" ), data[3] );
			}

			{
				StringList varString1dim( file, "/varString1dim" );
				assert_eq( "varString1dim.countItems", 2, varString1dim.countItems() );
				assert_eq( "varString1dim.countDimensions", 2, varString1dim.countDimensions() );
				string data[3];
				data[2] = string( "17" );
				varString1dim.readAll( data );
				assert_eq( "varString1dim data[0]", string( "" ), data[0] );
				assert_eq( "varString1dim data[1]", string( "-" ), data[1] );
				assert_eq( "varString1dim data after unchanged", string( "17" ), data[2] );
			}
		}

		tearDown( testFilename );
	}

}
