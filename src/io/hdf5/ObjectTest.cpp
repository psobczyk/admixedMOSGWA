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

#include "Object.hpp"
#include "Hdf5TestSuite.hpp"
#include <string>

using namespace std;
using namespace unitpp;
using namespace hdf5;

namespace test {

	/** Tests the class {@link hdf5::Object}. */
	struct ObjectTest : public Hdf5TestSuite {

		/** Construct the test object. */
		ObjectTest ();

		/** Test the template Object<D>. */
		void test ();

	} * objectTest = new ObjectTest();	// automatically freed by unit++

	ObjectTest::ObjectTest () : Hdf5TestSuite( "ObjectTest" ) {
		addTestMethod( "ObjectTest::test", this, &ObjectTest::test );
	}

	void ObjectTest::test () {
		const string testFilename = setUp();

		// begin a block in order to avoid tearDown() before object release
		{
			File file( testFilename );

			{
				Object<0> double0dim( file, "/double0dim" );
				assert_eq( "double0dim.countItems", 1, double0dim.countItems() );
				double data[2];
				data[1] = 17.0;
				double0dim.readAll( data );
				assert_eq( "double0dim data[0]", 0.0, data[0] );
				assert_eq( "double0dim data after unchanged", 17.0, data[1] );
			}

			{
				Object<1> double1dim( file, "/double1dim" );
				assert_eq( "double1dim.countItems", 2, double1dim.countItems() );
				double data[3];
				data[2] = 17.0;
				double1dim.readAll( data );
				assert_eq( "double1dim data[0]", 0.0, data[0] );
				assert_eq( "double1dim data[1]", 0.1, data[1] );
				assert_eq( "double1dim data after unchanged", 17.0, data[2] );
			}

			{
				Object<2> double2dim( file, "/double2dim" );
				assert_eq( "double2dim.countItems", 6, double2dim.countItems() );
				double data[7];
				data[6] = 17.0;
				double2dim.readAll( data );
				assert_eq( "double2dim data[0]", 0.0, data[0] );
				assert_eq( "double2dim data[1]", 0.1, data[1] );
				assert_eq( "double2dim data[2]", 1.0, data[2] );
				assert_eq( "double2dim data[3]", 1.1, data[3] );
				assert_eq( "double2dim data[4]", 2.0, data[4] );
				assert_eq( "double2dim data[5]", 2.1, data[5] );
				assert_eq( "double2dim data after unchanged", 17.0, data[6] );
			}

			{
				Object<3> double3dim( file, "/double3dim" );
				assert_eq( "double3dim.countItems", 12, double3dim.countItems() );
				double data[13];
				data[12] = 17.0;
				double3dim.readAll( data );
				assert_eq( "double3dim data[0]", 0.0, data[0] );
				assert_eq( "double3dim data[1]", 0.1, data[1] );
				assert_eq( "double3dim data[2]", 1.0, data[2] );
				assert_eq( "double3dim data[3]", 1.1, data[3] );
				assert_eq( "double3dim data[4]", 2.0, data[4] );
				assert_eq( "double3dim data[5]", 2.1, data[5] );
				assert_eq( "double3dim data[6]", -0.0, data[6] );
				assert_eq( "double3dim data[7]", -0.1, data[7] );
				assert_eq( "double3dim data[8]", -1.0, data[8] );
				assert_eq( "double3dim data[9]", -1.1, data[9] );
				assert_eq( "double3dim data[10]", -2.0, data[10] );
				assert_eq( "double3dim data[11]", -2.1, data[11] );
				assert_eq( "double3dim data after unchanged", 17.0, data[12] );
			}

			{
				Object<0> fixString0dim( file, "/fixString0dim" );
				assert_eq( "fixString0dim.countItems", 1, fixString0dim.countItems() );
				string data[2];
				data[1] = string( "17" );
				fixString0dim.readAll( data );
				assert_eq( "fixString0dim data[0]", string( "," ), data[0] );
				assert_eq( "fixString0dim data after unchanged", string( "17" ), data[1] );
			}

			{
				Object<1> fixString1dim( file, "/fixString1dim" );
				assert_eq( "fixString1dim.countItems", 3, fixString1dim.countItems() );
				string data[4];
				data[3] = string( "17" );
				fixString1dim.readAll( data );
				assert_eq( "fixString1dim data[0]", string( "," ), data[0] );
				assert_eq( "fixString1dim data[1]", string( ",1" ), data[1] );
				assert_eq( "fixString1dim data[2]", string( ",2" ), data[2] );
				assert_eq( "fixString1dim data after unchanged", string( "17" ), data[3] );
			}

			{
				Object<2> fixString2dim( file, "/fixString2dim" );
				assert_eq( "fixString2dim.countItems", 6, fixString2dim.countItems() );
				string data[7];
				data[6] = string( "17" );
				fixString2dim.readAll( data );
				assert_eq( "fixString2dim data[0]", string( "," ), data[0] );
				assert_eq( "fixString2dim data[1]", string( ",1" ), data[1] );
				assert_eq( "fixString2dim data[2]", string( ",2" ), data[2] );
				assert_eq( "fixString2dim data[3]", string( "1," ), data[3] );
				assert_eq( "fixString2dim data[4]", string( "1,1" ), data[4] );
				assert_eq( "fixString2dim data[5]", string( "1,2" ), data[5] );
				assert_eq( "fixString2dim data after unchanged", string( "17" ), data[6] );
			}

			{
				Object<0> varString0dim( file, "/varString0dim" );
				assert_eq( "varString0dim.countItems", 1, varString0dim.countItems() );
				string data[2];
				data[1] = string( "17" );
				varString0dim.readAll( data );
				assert_eq( "varString0dim data[0]", string( "" ), data[0] );
				assert_eq( "varString0dim data after unchanged", string( "17" ), data[1] );
			}

			{
				Object<1> varString1dim( file, "/varString1dim" );
				assert_eq( "varString1dim.countItems", 2, varString1dim.countItems() );
				string data[3];
				data[2] = string( "17" );
				varString1dim.readAll( data );
				assert_eq( "varString1dim data[0]", string( "" ), data[0] );
				assert_eq( "varString1dim data[1]", string( "-" ), data[1] );
				assert_eq( "varString1dim data after unchanged", string( "17" ), data[2] );
			}

			{
				Object<2> varString2dim( file, "/varString2dim" );
				assert_eq( "varString2dim.countItems", 6, varString2dim.countItems() );
				string data[7];
				data[6] = string( "17" );
				varString2dim.readAll( data );
				assert_eq( "varString2dim data[0]", string( "" ), data[0] );
				assert_eq( "varString2dim data[1]", string( "-" ), data[1] );
				assert_eq( "varString2dim data[2]", string( "+" ), data[2] );
				assert_eq( "varString2dim data[3]", string( "+-" ), data[3] );
				assert_eq( "varString2dim data[4]", string( "++" ), data[4] );
				assert_eq( "varString2dim data[5]", string( "++-" ), data[5] );
				assert_eq( "varString2dim data after unchanged", string( "17" ), data[6] );
			}
		}

		tearDown( testFilename );
	}

}
