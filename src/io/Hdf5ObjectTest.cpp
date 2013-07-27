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

#include "Hdf5Object.hpp"
#include "../TestSuite.hpp"
#include <unistd.h>	// for unlink(char[]), close()
#include <cstdlib>	// for mkstemp(char[])
#include <string>

using namespace std;
using namespace unitpp;
using namespace io;

namespace test {

	/** Tests the class {@link Hdf5Object}. */
	class Hdf5ObjectTest : public TestSuite {

		/** Prepare test setup.
		* @returns name of the test file.
		*/
		string setUp ();

		/** Remove test setup. */
		void tearDown ( const string& testFilename );

		/** Test the template Hdf5Object<D>. */
		void testGeneric ();

		/** Test the specific Table and List objects. */
		void testSpecific ();

		public:

		/** Construct the test object. */
		Hdf5ObjectTest ();

	} * hdf5ObjectTest = new Hdf5ObjectTest();	// automatically freed by unit++

	Hdf5ObjectTest::Hdf5ObjectTest () : TestSuite( "Hdf5ObjectTest" ) {
		addTestMethod( "Hdf5ObjectTest::testGeneric", this, &Hdf5ObjectTest::testGeneric );
		addTestMethod( "Hdf5ObjectTest::testSpecific", this, &Hdf5ObjectTest::testSpecific );
	}

	string Hdf5ObjectTest::setUp () {
		char buffer[] = "hdf5object-test.hdf5.XXXXXX";
		const int fd = mkstemp( buffer );
		assert_true( "Create test file", 0 <= fd );
		assert_eq( "Close fresh test file", 0, close( fd ) );
		const string testFilename( buffer );

		const hid_t h5FileId = H5Fcreate( testFilename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT );

		const hsize_t doubleDims[3] = {2, 3, 2};
		const double doubleData[2][3][2] = {
			{
				{ 0.0, 0.1 },
				{ 1.0, 1.1 },
				{ 2.0, 2.1 }
			},
			{
				{ -0.0, -0.1 },
				{ -1.0, -1.1 },
				{ -2.0, -2.1 }
			}
		};

		for ( char dims = 0; dims <= 3; ++dims ) {
			const char doubleName[] = { '/', 'd', 'o', 'u', 'b', 'l', 'e', '0'+dims, 'd', 'i', 'm', '\0' };
			const hid_t doubleSpace = H5Screate( H5S_SIMPLE );
			H5Sset_extent_simple( doubleSpace, dims, doubleDims + 3 - dims, NULL );
			const hid_t doubleDataset = H5Dcreate2(
				h5FileId,
				doubleName,
				H5T_NATIVE_DOUBLE,
				doubleSpace,
				H5P_DEFAULT,
				H5P_DEFAULT,
				H5P_DEFAULT
			);
			H5Dwrite( doubleDataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, doubleData );
			H5Dclose( doubleDataset );
			H5Sclose( doubleSpace );
		}

		{
			const hsize_t fixStringDims[2] = {2, 3};
			const char fixStringData[] = ",\000\000,1\000,2\0001,\0001,11,2";

			const hid_t fixStringType = H5Tcopy( H5T_C_S1 );
			H5Tset_strpad( fixStringType, H5T_STR_NULLPAD );
			H5Tset_size( fixStringType, 3 );

			for ( char dims = 0; dims <= 2; ++dims ) {
				const char fixStringName[] = { '/', 'f', 'i', 'x', 'S', 't', 'r', 'i', 'n', 'g', '0'+dims, 'd', 'i', 'm', '\0' };
				const hid_t fixStringSpace = H5Screate( H5S_SIMPLE );
				H5Sset_extent_simple( fixStringSpace, dims, fixStringDims + 2 - dims, NULL );
				const hid_t fixStringDataset = H5Dcreate2(
					h5FileId,
					fixStringName,
					fixStringType,
					fixStringSpace,
					H5P_DEFAULT,
					H5P_DEFAULT,
					H5P_DEFAULT
				);
				H5Dwrite( fixStringDataset, fixStringType, H5S_ALL, H5S_ALL, H5P_DEFAULT, fixStringData );
				H5Dclose( fixStringDataset );
				H5Sclose( fixStringSpace );
			}

			H5Tclose( fixStringType );
		}

		{
			const hid_t varStringType = H5Tcopy( H5T_C_S1 );
			H5Tset_size( varStringType, H5T_VARIABLE );

			const hsize_t varStringDims[2] = {3, 2};
			const char* varStringData[3][2] = {
				{ "", "-" },
				{ "+", "+-" },
				{ "++", "++-" }
			};

			for ( char dims = 0; dims <= 2; ++dims ) {
				const char varStringName[] = { '/', 'v', 'a', 'r', 'S', 't', 'r', 'i', 'n', 'g', '0'+dims, 'd', 'i', 'm', '\0' };
				const hid_t varStringSpace = H5Screate( H5S_SIMPLE );
				H5Sset_extent_simple( varStringSpace, dims, varStringDims + 2 - dims, NULL );
				const hid_t varStringDataset = H5Dcreate2(
					h5FileId,
					varStringName,
					varStringType,
					varStringSpace,
					H5P_DEFAULT,
					H5P_DEFAULT,
					H5P_DEFAULT
				);
				H5Dwrite( varStringDataset, varStringType, H5S_ALL, H5S_ALL, H5P_DEFAULT, varStringData );
				H5Dclose( varStringDataset );
				H5Sclose( varStringSpace );
			}

			H5Tclose( varStringType );
		}

		H5Fclose( h5FileId );

		return testFilename;
	}

	void Hdf5ObjectTest::tearDown ( const string& testFilename ) {
		assert_eq( "Remove test file failed", 0, unlink( testFilename.c_str() ) );
	}

	void Hdf5ObjectTest::testGeneric () {
		const string testFilename = setUp();

		// begin a block in order to avoid tearDown() before object release
		{
			Hdf5FileId fileId( testFilename );

			{
				Hdf5Object<0> double0dim( fileId, "/double0dim" );
				assert_eq( "double0dim.countItems", 1, double0dim.countItems() );
				double data[2];
				data[1] = 17.0;
				double0dim.readAll( data );
				assert_eq( "double0dim data[0]", 0.0, data[0] );
				assert_eq( "double0dim data after unchanged", 17.0, data[1] );
			}

			{
				Hdf5Object<1> double1dim( fileId, "/double1dim" );
				assert_eq( "double1dim.countItems", 2, double1dim.countItems() );
				double data[3];
				data[2] = 17.0;
				double1dim.readAll( data );
				assert_eq( "double1dim data[0]", 0.0, data[0] );
				assert_eq( "double1dim data[1]", 0.1, data[1] );
				assert_eq( "double1dim data after unchanged", 17.0, data[2] );
			}

			{
				Hdf5Object<2> double2dim( fileId, "/double2dim" );
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
				Hdf5Object<3> double3dim( fileId, "/double3dim" );
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
				Hdf5Object<0> fixString0dim( fileId, "/fixString0dim" );
				assert_eq( "fixString0dim.countItems", 1, fixString0dim.countItems() );
				string data[2];
				data[1] = string( "17" );
				fixString0dim.readAll( data );
				assert_eq( "fixString0dim data[0]", string( "," ), data[0] );
				assert_eq( "fixString0dim data after unchanged", string( "17" ), data[1] );
			}

			{
				Hdf5Object<1> fixString1dim( fileId, "/fixString1dim" );
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
				Hdf5Object<2> fixString2dim( fileId, "/fixString2dim" );
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
				Hdf5Object<0> varString0dim( fileId, "/varString0dim" );
				assert_eq( "varString0dim.countItems", 1, varString0dim.countItems() );
				string data[2];
				data[1] = string( "17" );
				varString0dim.readAll( data );
				assert_eq( "varString0dim data[0]", string( "" ), data[0] );
				assert_eq( "varString0dim data after unchanged", string( "17" ), data[1] );
			}

			{
				Hdf5Object<1> varString1dim( fileId, "/varString1dim" );
				assert_eq( "varString1dim.countItems", 2, varString1dim.countItems() );
				string data[3];
				data[2] = string( "17" );
				varString1dim.readAll( data );
				assert_eq( "varString1dim data[0]", string( "" ), data[0] );
				assert_eq( "varString1dim data[1]", string( "-" ), data[1] );
				assert_eq( "varString1dim data after unchanged", string( "17" ), data[2] );
			}

			{
				Hdf5Object<2> varString2dim( fileId, "/varString2dim" );
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

	void Hdf5ObjectTest::testSpecific () {
		const string testFilename = setUp();

		// begin a block in order to avoid tearDown() before object release
		{
			Hdf5FileId fileId( testFilename );

			{
				Hdf5DoubleList double1dim( fileId, "/double1dim" );
				assert_eq( "double1dim.countItems", 2, double1dim.countItems() );
				assert_eq( "double1dim.countDimensions", 2, double1dim.countDimensions() );
				double data[3];
				data[2] = 17.0;
				double1dim.readAll( data );
				assert_eq( "double1dim data[0]", 0.0, data[0] );
				assert_eq( "double1dim data[1]", 0.1, data[1] );
				assert_eq( "double1dim data after unchanged", 17.0, data[2] );
			}

			{
				Hdf5DoubleTable double2dim( fileId, "/double2dim" );
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

			{
				Hdf5StringList fixString1dim( fileId, "/fixString1dim" );
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
				Hdf5StringList varString1dim( fileId, "/varString1dim" );
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
