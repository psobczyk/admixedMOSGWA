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

#include "Hdf5Input.hpp"
#include "../TestSuite.hpp"
#include <unistd.h>	// for unlink(char[]), close()
#include <cstdlib>	// for mkstemp(char[])
#include <cstring>
#include <cmath>	// for NaN
#include <vector>

using namespace std;
using namespace io;
using namespace linalg;
using namespace unitpp;

namespace test {

	/** Tests the class {@link Hdf5Input}. */
	class Hdf5InputTest : public TestSuite {

		/** Name template for temporary test data file. */
		static const char * const tmpFilenameTemplate;

		/** Prepare test setup.
		* @returns name of the test file.
		*/
		string setUp ();

		/** Remove test setup. */
		void tearDown ( const string& testFilename );

		/** Test interface methods. */
		void testRead ();

		public:

		/** Construct the test object. */
		Hdf5InputTest ();

	} * hdf5InputTest = new Hdf5InputTest();	// automatically freed by unit++

	const char * const Hdf5InputTest::tmpFilenameTemplate = "hdf5input-test.hdf5.XXXXXX";

	Hdf5InputTest::Hdf5InputTest () : TestSuite( "Hdf5InputTest" ) {
		addTestMethod( "Hdf5InputTest::testRead", this, &Hdf5InputTest::testRead );
	}

	string Hdf5InputTest::setUp () {
		vector<char> tmpFilename( strlen( tmpFilenameTemplate ) + 1 );	// +1 for trailing 0
		strcpy( tmpFilename.data(), tmpFilenameTemplate );
		const int fd = mkstemp( tmpFilename.data() );
		assert_true( "Create test file failed", 0 <= fd );
		assert_eq( "Close fresh test file failed", 0, close( fd ) );
		const string testFilename( tmpFilename.data() );

		const hid_t h5FileId = H5Fcreate( testFilename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT );

		const hid_t fixStringType = H5Tcopy( H5T_C_S1 );
		H5Tset_strpad( fixStringType, H5T_STR_NULLPAD );
		H5Tset_size( fixStringType, 6 );
		const hid_t varStringType = H5Tcopy( fixStringType );
		H5Tset_size( varStringType, H5T_VARIABLE );

		const hsize_t dims[2] = {4, 3};

		{
			const double genomeData[4][3] = {
				{ 0.0, 0.1, 0.2 },
				{ 1.0, 1.1, 1.2 },
				{ 2.0, ::nan("test"), 2.2 },
				{ 3.0, 3.1, 3.2 }
			};
			const hid_t genomeSpace = H5Screate( H5S_SIMPLE );
			H5Sset_extent_simple( genomeSpace, 2, dims, NULL );
			const hid_t genomeDataset = H5Dcreate2(
				h5FileId,
				"/genome_matrix",
				H5T_NATIVE_DOUBLE,
				genomeSpace,
				H5P_DEFAULT,
				H5P_DEFAULT,
				H5P_DEFAULT
			);
			H5Dwrite( genomeDataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, genomeData );
			H5Dclose( genomeDataset );
			H5Sclose( genomeSpace );
		}

		{
			const char * snpsData[4] = {
				"1_2",
				"22_11",
				"0_0",
				"98_97"
			};
			const hid_t snpsSpace = H5Screate( H5S_SIMPLE );
			H5Sset_extent_simple( snpsSpace, 1, dims, NULL );
			const hid_t snpsDataset = H5Dcreate2(
				h5FileId,
				"/single_nucleotide_polymorphisms",
				varStringType,
				snpsSpace,
				H5P_DEFAULT,
				H5P_DEFAULT,
				H5P_DEFAULT
			);
			H5Dwrite( snpsDataset, varStringType, H5S_ALL, H5S_ALL, H5P_DEFAULT, snpsData );
			H5Dclose( snpsDataset );
			H5Sclose( snpsSpace );
		}

		{
			const hid_t individualsSpace = H5Screate( H5S_SIMPLE );
			H5Sset_extent_simple( individualsSpace, 1, dims + 1, NULL );

			const char individualsData[] = "Eric\000\000BenhaoFlo\000\000\000";
			const hid_t individualsDataset = H5Dcreate2(
				h5FileId,
				"/individuals",
				fixStringType,
				individualsSpace,
				H5P_DEFAULT,
				H5P_DEFAULT,
				H5P_DEFAULT
			);
			H5Dwrite( individualsDataset, fixStringType, H5S_ALL, H5S_ALL, H5P_DEFAULT, individualsData );
			H5Dclose( individualsDataset );

			const hid_t phenotypeDataset = H5Dcreate2(
				h5FileId,
				"/phenotypes",
				H5T_NATIVE_FLOAT,
				individualsSpace,
				H5P_DEFAULT,
				H5P_DEFAULT,
				H5P_DEFAULT
			);
			const float phenotypeData[3] = {
				-0.0, ::nan("test"), -0.2
			};
			H5Dwrite( phenotypeDataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, phenotypeData );
			H5Dclose( phenotypeDataset );

			H5Sclose( individualsSpace );
		}

		H5Tclose( varStringType );
		H5Tclose( fixStringType );
		H5Fclose( h5FileId );

		return testFilename;
	}

	void Hdf5InputTest::tearDown ( const string& testFilename ) {
		assert_eq( "Remove test file failed", 0, unlink( testFilename.c_str() ) );
	}

	void Hdf5InputTest::testRead () {
		const string testFilename( setUp() );
		Hdf5Input hdf5Input( testFilename.c_str() );
		AutoVector vector( hdf5Input.countIndividuals() );

		{
			const SNP snp0 = hdf5Input.getSnp( 0 );
			assert_eq( "SNP[0].id", string( "1_2" ), snp0.getSnpId() );
			assert_eq( "SNP[0].chromosome", 1, snp0.getChromosome() );
			assert_eq( "SNP[0].position", 2, snp0.getBasePairPosition() );
		}

		{
			const SNP snp1 = hdf5Input.getSnp( 1 );
			assert_eq( "SNP[1].id", string( "22_11" ), snp1.getSnpId() );
			assert_eq( "SNP[1].chromosome", 22, snp1.getChromosome() );
			assert_eq( "SNP[1].position", 11, snp1.getBasePairPosition() );
		}

		{
			const SNP snp2 = hdf5Input.getSnp( 2 );
			assert_eq( "SNP[2].id", string( "0_0" ), snp2.getSnpId() );
			assert_eq( "SNP[2].chromosome", 0, snp2.getChromosome() );
			assert_eq( "SNP[2].position", 0, snp2.getBasePairPosition() );
		}

		{
			const SNP snp3 = hdf5Input.getSnp( 3 );
			assert_eq( "SNP[3].id", string( "98_97" ), snp3.getSnpId() );
			assert_eq( "SNP[3].chromosome", 98, snp3.getChromosome() );
			assert_eq( "SNP[3].position", 97, snp3.getBasePairPosition() );
		}

		{
			const Individual idv0 = hdf5Input.getIndividual( 0 );
			assert_eq( "Individual[0].id", string( "Eric" ), idv0.getIndividualID() );
			assert_eq( "Individual[0].phenotype", -0.0f, idv0.getPhenotype() );
		}

		{
			const Individual idv1 = hdf5Input.getIndividual( 1 );
			assert_eq( "Individual[1].id", string( "Benhao" ), idv1.getIndividualID() );
			assert_true( "Individual[1].phenotype is NaN", ::isnan( idv1.getPhenotype() ) );
		}

		{
			const Individual idv2 = hdf5Input.getIndividual( 2 );
			assert_eq( "Individual[2].id", string( "Flo" ), idv2.getIndividualID() );
			assert_eq( "Individual[2].phenotype", -0.2f, idv2.getPhenotype() );
		}

		{
			hdf5Input.retrieveGenotypesIntoVector( 0, vector );
			assert_eq( "genotypeMatrixNontransposed[0,0]", 0.0, vector.get( 0 ) );
			assert_eq( "genotypeMatrixNontransposed[1,0]", 0.1, vector.get( 1 ) );
			assert_eq( "genotypeMatrixNontransposed[2,0]", 0.2, vector.get( 2 ) );
		}

		{
			hdf5Input.retrieveGenotypesIntoVector( 1, vector );
			assert_eq( "genotypeMatrixNontransposed[0,1]", 1.0, vector.get( 0 ) );
			assert_eq( "genotypeMatrixNontransposed[1,1]", 1.1, vector.get( 1 ) );
			assert_eq( "genotypeMatrixNontransposed[2,1]", 1.2, vector.get( 2 ) );
		}

		{
			hdf5Input.retrieveGenotypesIntoVector( 2, vector );
			assert_eq( "genotypeMatrixNontransposed[0,2]", 2.0, vector.get( 0 ) );
			assert_true( "genotypeMatrixNontransposed[1,2] is NaN", ::isnan( vector.get( 1 ) ) );
			assert_eq( "genotypeMatrixNontransposed[2,2]", 2.2, vector.get( 2 ) );
		}

		{
			hdf5Input.retrieveGenotypesIntoVector( 3, vector );
			assert_eq( "genotypeMatrixNontransposed[0,3]", 3.0, vector.get( 0 ) );
			assert_eq( "genotypeMatrixNontransposed[1,3]", 3.1, vector.get( 1 ) );
			assert_eq( "genotypeMatrixNontransposed[2,3]", 3.2, vector.get( 2 ) );
		}

		{
			hdf5Input.retrievePhenotypesIntoVector( vector );
			// Mind that in the test data, phenotype has been stored from a float array.
			assert_eq( "phenotypeVector[0]", -0.0f, vector.get( 0 ) );
			assert_true( "phenotypeVector[1] is NaN", ::isnan( vector.get( 1 ) ) );
			assert_eq( "phenotypeVector[2]", -0.2f, vector.get( 2 ) );
		}

		tearDown( testFilename );
	}

}
