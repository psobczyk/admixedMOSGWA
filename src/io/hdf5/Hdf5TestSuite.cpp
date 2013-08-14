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

#include "Hdf5TestSuite.hpp"
#include <unistd.h>	// for unlink(char[]), close()
#include <cstdlib>	// for mkstemp(char[])
#include <hdf5.h>

using namespace std;
using namespace unitpp;

namespace test {

	Hdf5TestSuite::Hdf5TestSuite ( const char* name ) : TestSuite( name ) {}

	string Hdf5TestSuite::Hdf5TestSuite::setUp () {
		char buffer[] = "hdf5-object-test.hdf5.XXXXXX";
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

	void Hdf5TestSuite::tearDown ( const string& testFilename ) {
		assert_eq( "Remove test file failed", 0, unlink( testFilename.c_str() ) );
	}

}
