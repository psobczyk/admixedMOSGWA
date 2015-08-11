/********************************************************************************
 *	This file is part of the MOSGWA program code.				*
 *	Copyright ©2013, Bernhard Bodenstorfer.					*
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

#include "Hdf5Output.hpp"
#include "Hdf5Constants.hpp"
#include "../TestSuite.hpp"
#include "../Exception.hpp"
#include "../linalg/AutoVector.hpp"
#include "../linalg/AutoMatrix.hpp"
#include <unistd.h>	// for unlink(char[]), close()
#include <cstdlib>	// for mkstemp(char[])
#include <cstring>
#include <cmath>	// for NaN
#include <vector>

using namespace std;
using namespace io;
using namespace io::Hdf5Constants;
using namespace linalg;
using namespace unitpp;

namespace test {

	/** Tests the class {@link io::Hdf5Output}. */
	class Hdf5OutputTest : public TestSuite {

		/** Name template for temporary test data file. */
		static const char * const tmpFilenameTemplate;

		/** Prepare test setup.
		* @returns name of the test file.
		*/
		string setUp ();

		/** Remove test setup. */
		void tearDown ( const string& testFilename );

		/** Test interface methods for writing. */
		void testWrite ();

		public:

		/** Construct the test object. */
		Hdf5OutputTest ();

	} * hdf5OutputTest = new Hdf5OutputTest();	// automatically freed by unit++

	const char * const Hdf5OutputTest::tmpFilenameTemplate = "hdf5output-test.hdf5.XXXXXX";

	Hdf5OutputTest::Hdf5OutputTest () : TestSuite( "Hdf5OutputTest" ) {
		addTestMethod( "Hdf5OutputTest::testWrite", this, &Hdf5OutputTest::testWrite );
	}

	string Hdf5OutputTest::setUp () {
		vector<char> tmpFilename( strlen( tmpFilenameTemplate ) + 1 );	// +1 for trailing 0
		strcpy( tmpFilename.data(), tmpFilenameTemplate );
		const int fd = mkstemp( tmpFilename.data() );
		assert_true( "Create test file failed", 0 <= fd );
		assert_eq( "Close fresh test file failed", 0, close( fd ) );
		const string testFilename( tmpFilename.data() );
		unlink( testFilename.c_str() );
		return testFilename;
	}

/*
		{
			const hid_t genomeSpace = H5Screate( H5S_SIMPLE );
			H5Sset_extent_simple( genomeSpace, 2, dims, NULL );
			const hid_t genomeDataset = H5Dcreate2(
				h5FileId,
				genotypeMatrixPath,
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
		*/

	void Hdf5OutputTest::tearDown ( const string& testFilename ) {
		assert_eq( "Remove test file failed", 0, unlink( testFilename.c_str() ) );
	}

	const Individual individuals[4] = {
		Individual( "Lustig", "Gustl", "Julius", "Calliope", Individual::MISSING ),
		Individual( "Frustlos", "Nana", "Jeff", "Martha", Individual::FEMALE ),
		Individual( "Tucor", "Mira", "Frank", "Mary", Individual::FEMALE ),
		Individual( "Wuspig", "Wunislaus", "Julius", "Mitzi", Individual::MALE )
	};
	const SNP snps [3] = {
		SNP( 1, "1_11", 1.11, 11, 'A', 'T' ),
		SNP( 2, "99_97", ::nan("test"), 0, 'C', 'G' ),
		SNP( 23, "40_0", -1e1, 11, 'G', 'A' )
	};
	const string covariates[2] = {
		"ignorance",
		"irrelevance"
	};

	void Hdf5OutputTest::testWrite () {
		const string testFilename( setUp() );
		Hdf5Output hdf5Output( testFilename.c_str(), 4, 3, 2 );
		hdf5Output.setIndividuals( individuals );
		hdf5Output.setSnps( snps );
		hdf5Output.setCovariates( covariates );
		const double data[16] = {
			0.0, 0.1, 0.2, ::nan("test"),
			1.0, 1.1, 1.2, 1.3,
			2.0, ::nan("test"), 2.2, ::nan("test"),
			nan("test"), 3.1, 3.2, ::nan("test")
		};
		AutoMatrix a( 4, 4 );
		a.fill( data );

		hdf5Output.storePhenotypeVector( a.diagonalVector() );

		hdf5Output.storeGenotypeVector( 0, a.columnVector( 0 ) );
		hdf5Output.storeGenotypeVector( 1, a.columnVector( 1 ) );
		hdf5Output.storeGenotypeVector( 2, a.rowVector( 0 ) );

		hdf5Output.storeCovariateVector( 0, a.columnVector( 3 ) );
		hdf5Output.storeCovariateVector( 1, a.rowVector( 0 ) );

		const hid_t file = H5Fopen( testFilename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT );
		assert_true( "file open", 0 <= file );

		// Test individuals
		{
			const hid_t dataset = H5Dopen2(
				file,
				Hdf5Constants::individualListPath,
				H5P_DEFAULT
			);
			assert_true( "individuals open", 0 <= dataset );

			const hid_t datatype = H5Dget_type( dataset );
			assert_true( "individuals type", 0 <= datatype );
			assert_eq( "individuals class", H5T_STRING, H5Tget_class( datatype ) );
			assert_true( "individuals var string", 0 < H5Tis_variable_str( datatype ) );

			const hid_t dataspace = H5Dget_space( dataset );
			assert_true( "individuals space", 0 <= dataspace );
			assert_eq( "individuals D", 1, H5Sget_simple_extent_ndims( dataspace ) );
			hsize_t dims[1], maxdims[1];
			assert_true( "individuals dims", 0 <= H5Sget_simple_extent_dims( dataspace, dims, maxdims ) );
			assert_eq( "individuals actual dim", 4, dims[0] );
			assert_eq( "individuals maximum dim", 4, maxdims[0] );
			char * ids[3];
			assert_true( "individuals read", 0 <= H5Dread( dataset, datatype, H5S_ALL, H5S_ALL, H5P_DEFAULT, ids ) );
			for ( size_t i = 0; i < 4; ++i ) {
				char sbuf[256];
				snprintf( sbuf, 255, "individuals[%u]", i );
				assert_eq( sbuf, individuals[i].getIndividualID(), ids[i] );
			}

			// Note: This is unit-test code.
			// Production quality code would need to ensure reclaim even after exceptions.
			H5Dvlen_reclaim( datatype, dataspace, H5P_DEFAULT, ids );

			H5Sclose( dataspace );
			H5Tclose( datatype );
			H5Dclose( dataset );
		}

		// Test SNPs
		{
			const hid_t dataset = H5Dopen2(
				file,
				Hdf5Constants::snpListPath,
				H5P_DEFAULT
			);
			assert_true( "snps open", 0 <= dataset );

			const hid_t datatype = H5Dget_type( dataset );
			assert_true( "snps type", 0 <= datatype );
			assert_eq( "snps class", H5T_STRING, H5Tget_class( datatype ) );
			assert_true( "snps var string", 0 < H5Tis_variable_str( datatype ) );

			const hid_t dataspace = H5Dget_space( dataset );
			assert_true( "snps space", 0 <= dataspace );
			assert_eq( "snps D", 1, H5Sget_simple_extent_ndims( dataspace ) );
			hsize_t dims[1], maxdims[1];
			assert_true( "snps dims", 0 <= H5Sget_simple_extent_dims( dataspace, dims, maxdims ) );
			assert_eq( "snps actual dim", 3, dims[0] );
			assert_eq( "snps maximum dim", 3, maxdims[0] );
			char * ids[3];
			assert_true( "snps read", 0 <= H5Dread( dataset, datatype, H5S_ALL, H5S_ALL, H5P_DEFAULT, ids ) );
			for ( size_t i = 0; i < 3; ++i ) {
				char sbuf[256];
				snprintf( sbuf, 255, "snps[%u]", i );
				assert_eq( sbuf, snps[i].getSnpId(), ids[i] );
			}

			// Note: This is unit-test code.
			// Production quality code would need to ensure reclaim even after exceptions.
			H5Dvlen_reclaim( datatype, dataspace, H5P_DEFAULT, ids );

			H5Sclose( dataspace );
			H5Tclose( datatype );
			H5Dclose( dataset );
		}

		// Test covariates
		{
			const hid_t dataset = H5Dopen2(
				file,
				Hdf5Constants::covariateListPath,
				H5P_DEFAULT
			);
			assert_true( "covariates open", 0 <= dataset );

			const hid_t datatype = H5Dget_type( dataset );
			assert_true( "covariates type", 0 <= datatype );
			assert_eq( "covariates class", H5T_STRING, H5Tget_class( datatype ) );
			assert_true( "covariates var string", 0 < H5Tis_variable_str( datatype ) );

			const hid_t dataspace = H5Dget_space( dataset );
			assert_true( "covariates space", 0 <= dataspace );
			assert_eq( "covariates D", 1, H5Sget_simple_extent_ndims( dataspace ) );
			hsize_t dims[1], maxdims[1];
			assert_true( "covariates dims", 0 <= H5Sget_simple_extent_dims( dataspace, dims, maxdims ) );
			assert_eq( "covariates actual dim", 2, dims[0] );
			assert_eq( "covariates maximum dim", 2, maxdims[0] );
			char * ids[3];
			assert_true( "covariates read", 0 <= H5Dread( dataset, datatype, H5S_ALL, H5S_ALL, H5P_DEFAULT, ids ) );
			for ( size_t i = 0; i < 2; ++i ) {
				char sbuf[256];
				snprintf( sbuf, 255, "covariates[%u]", i );
				assert_eq( sbuf, covariates[i], ids[i] );
			}

			// Note: This is unit-test code.
			// Production quality code would need to ensure reclaim even after exceptions.
			H5Dvlen_reclaim( datatype, dataspace, H5P_DEFAULT, ids );

			H5Sclose( dataspace );
			H5Tclose( datatype );
			H5Dclose( dataset );
		}

		// Test phenotype
		{
			const hid_t dataset = H5Dopen2(
				file,
				Hdf5Constants::phenotypeVectorPath,
				H5P_DEFAULT
			);
			assert_true( "phenotype vector open", 0 <= dataset );

			const hid_t datatype = H5Dget_type( dataset );
			assert_true( "phenotype vector type", 0 <= datatype );
			assert_eq( "phenotype vector class", H5T_FLOAT, H5Tget_class( datatype ) );

			const hid_t dataspace = H5Dget_space( dataset );
			assert_true( "phenotype vector space", 0 <= dataspace );
			assert_eq( "phenotype vector D", 1, H5Sget_simple_extent_ndims( dataspace ) );
			hsize_t dims[1], maxdims[1];
			assert_true( "phenotype vector dims", 0 <= H5Sget_simple_extent_dims( dataspace, dims, maxdims ) );
			assert_eq( "phenotype vector actual dim", 4, dims[0] );
			assert_eq( "phenotype vector maximum dim", 4, maxdims[0] );
			double buffer[4];
			assert_true( "phenotype vector read", 0 <= H5Dread( dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, buffer ) );
			for ( size_t i = 0; i < 2; ++i ) {
				char sbuf[256];
				snprintf( sbuf, 255, "phenotype vector[%u]", i );
				const double
					expected = a.get( i, i ),
					actual = buffer[i];
				if ( ::isnan( expected ) ) {
					assert_true( sbuf, ::isnan( actual ) );
				} else {
					assert_eq( sbuf, expected, actual );
				}
			}

			H5Sclose( dataspace );
			H5Tclose( datatype );
			H5Dclose( dataset );
		}

		// Test genotype
		{
			const hid_t dataset = H5Dopen2(
				file,
				Hdf5Constants::genotypeMatrixPath,
				H5P_DEFAULT
			);
			assert_true( "genotype matrix open", 0 <= dataset );

			const hid_t datatype = H5Dget_type( dataset );
			assert_true( "genotype matrix type", 0 <= datatype );
			assert_eq( "genotype matrix class", H5T_FLOAT, H5Tget_class( datatype ) );

			const hid_t dataspace = H5Dget_space( dataset );
			assert_true( "genotype matrix space", 0 <= dataspace );
			assert_eq( "genotype matrix D", 2, H5Sget_simple_extent_ndims( dataspace ) );
			hsize_t dims[2], maxdims[2];
			assert_true( "genotype matrix dims", 0 <= H5Sget_simple_extent_dims( dataspace, dims, maxdims ) );
			assert_eq( "genotype matrix actual dim 0", 3, dims[0] );
			assert_eq( "genotype matrix actual dim 1", 4, dims[1] );
			assert_eq( "genotype matrix maximum dim 0", 3, maxdims[0] );
			assert_eq( "genotype matrix maximum dim 1", 4, maxdims[1] );
			double buffer[3][4];
			assert_true( "genotype matrix read", 0 <= H5Dread( dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, buffer ) );
			// 0 and 1 are column vectors of A
			for ( size_t i = 0; i < 2; ++i ) {
				for ( size_t j = 0; j < 4; ++j ) {
					char sbuf[256];
					snprintf( sbuf, 255, "genotype matrix[%u][%u]", i, j );
					const double
						expected = a.get( j, i ),
						actual = buffer[i][j];
					if ( ::isnan( expected ) ) {
						assert_true( sbuf, ::isnan( actual ) );
					} else {
						assert_eq( sbuf, expected, actual );
					}
				}
			}
			// 2 is row vector of A
			for ( size_t i = 0; i < 4; ++i ) {
				char sbuf[256];
				snprintf( sbuf, 255, "genotype matrix[2][%u]", i );
				const double
					expected = a.get( 0, i ),
					actual = buffer[2][i];
				if ( ::isnan( expected ) ) {
					assert_true( sbuf, ::isnan( actual ) );
				} else {
					assert_eq( sbuf, expected, actual );
				}
			}

			H5Sclose( dataspace );
			H5Tclose( datatype );
			H5Dclose( dataset );
		}

		// Test covariates
		{
			const hid_t dataset = H5Dopen2(
				file,
				Hdf5Constants::covariateMatrixPath,
				H5P_DEFAULT
			);
			assert_true( "covariate matrix open", 0 <= dataset );

			const hid_t datatype = H5Dget_type( dataset );
			assert_true( "covariate matrix type", 0 <= datatype );
			assert_eq( "covariate matrix class", H5T_FLOAT, H5Tget_class( datatype ) );

			const hid_t dataspace = H5Dget_space( dataset );
			assert_true( "covariate matrix space", 0 <= dataspace );
			assert_eq( "covariate matrix D", 2, H5Sget_simple_extent_ndims( dataspace ) );
			hsize_t dims[2], maxdims[2];
			assert_true( "covariate matrix dims", 0 <= H5Sget_simple_extent_dims( dataspace, dims, maxdims ) );
			assert_eq( "covariate matrix actual dim 0", 2, dims[0] );
			assert_eq( "covariate matrix actual dim 1", 4, dims[1] );
			assert_eq( "covariate matrix maximum dim 0", 2, maxdims[0] );
			assert_eq( "covariate matrix maximum dim 1", 4, maxdims[1] );
			double buffer[3][4];
			assert_true( "covariate matrix read", 0 <= H5Dread( dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, buffer ) );
			// 0 is column 3 vector of A
			for ( size_t i = 0; i < 4; ++i ) {
				char sbuf[256];
				snprintf( sbuf, 255, "covariate matrix[0][%u]", i );
				const double
					expected = a.get( i, 3 ),
					actual = buffer[0][i];
				if ( ::isnan( expected ) ) {
					assert_true( sbuf, ::isnan( actual ) );
				} else {
					assert_eq( sbuf, expected, actual );
				}
			}
			// 1 is row 0 vector of A
			for ( size_t i = 0; i < 4; ++i ) {
				char sbuf[256];
				snprintf( sbuf, 255, "covariate matrix[1][%u]", i );
				const double
					expected = a.get( 0, i ),
					actual = buffer[1][i];
				if ( ::isnan( expected ) ) {
					assert_true( sbuf, ::isnan( actual ) );
				} else {
					assert_eq( sbuf, expected, actual );
				}
			}

			H5Sclose( dataspace );
			H5Tclose( datatype );
			H5Dclose( dataset );
		}

		assert_true( "file close", 0 <= H5Fclose( file ) );

		tearDown( testFilename );
	}

}
