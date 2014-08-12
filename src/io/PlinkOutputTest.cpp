/********************************************************************************
 *	This file is part of the MOSGWA program code.				*
 *	Copyright ©2014, Bernhard Bodenstorfer.					*
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

#include "PlinkOutput.hpp"
#include "PlinkConstants.hpp"
#include "../TestSuite.hpp"
#include "../Exception.hpp"
#include "../linalg/AutoVector.hpp"
#include <string>
#include <vector>
#include <bitset>
#include <fstream>
#include <cmath>
#include <cstring>
#include <unistd.h>

using namespace std;
using namespace linalg;
using namespace io;
using namespace io::PlinkConstants;
using namespace unitpp;

namespace test {

	/** Tests the class {@link io::PlinkOutput}. */
	class PlinkOutputTest : public TestSuite {

		/** Name template for temporary test data directory and files. */
		static const char
			* const tmpDirnameTemplate,
			* const filenameTrunc;

		static const Individual individuals[3];
		static const SNP snps[2];
		static const string
			covariates[2],
			traits[2];

		static const double
			genotype[],
			covariate[],
			phenotype[];

		/** Prepare test setup.
		* @param snpMajour determines whether SNP majour mode should be used. Else Individual majour mode.
		* @returns name of the test directory.
		*/
		string setUp ();

		/** Remove test setup. */
		void tearDown (
			const string& testDirname,
			const bool withCovariates,
			const bool withExtraPhenotypes
		 );

		/** Test {@link io::Output} interface methods. */
		void testWrite ();

		/** Test {@link io::OutputCo} interface extension methods. */
		void testCovariates ();

		/** Test {@link io::OutputCo} interface extension methods. */
		void testExtraPhenotypes ();

		public:

		/** Construct the test object. */
		PlinkOutputTest ();

	} * plinkOutputTest = new PlinkOutputTest();	// automatically freed by unit++

	const char
		* const PlinkOutputTest::tmpDirnameTemplate = "plinkoutput-test.XXXXXX",
		* const PlinkOutputTest::filenameTrunc = "data";

	const Individual PlinkOutputTest::individuals[3] = {
		Individual( "Bierwetter", "Bernhard", "Camillo", "Praxedis", Individual::MALE ),
		Individual( "Wetter", "Praxedis", "Gianluca", "Caecilia", Individual::FEMALE ),
		Individual( "Bier", "Camillo", "Gonzago", "Pandora", Individual::MALE )
	};

	const SNP PlinkOutputTest::snps[2] = {
		SNP( 26, string( "SNP26" ), 26.26, 2626, 'A', 'C' ),
		SNP( 17, string( "SNP17" ), 17.17, 1717, 'C', 'T' )
	};

	const string
		PlinkOutputTest::covariates[2] = { "Age", "Sex" },
		PlinkOutputTest::traits[2] = { "Kebbelzahn", "Thirdnosehole" };

	const double
		PlinkOutputTest::genotype[] = {
			1.0,	0.0,
			0.0,	::nan( "missing" ),
			-1.0,	1.0
		},
		PlinkOutputTest::covariate[] = {
			3.1,	::nan( "missing" ),
			::nan( "missing" ),	5e+20,
			3e-20,	::nan( "missing" )
		},
		PlinkOutputTest::phenotype[] = {
			0.0,	0.0,
			1.0,	::nan( "missing" ),
			::nan( "missing" ),	3.5e+15
		};

	const char
		*famTextWithoutPhenotype[] = {
			"Bierwetter Bernhard\tCamillo\tPraxedis\t1\t-",
			"Wetter Praxedis\tGianluca\tCaecilia\t2\t-",
			"Bier Camillo\tGonzago\tPandora\t1\t-"
		},
		*famTextWithPhenotype[] = {
			"Bierwetter Bernhard\tCamillo\tPraxedis\t1\t0",
			"Wetter Praxedis\tGianluca\tCaecilia\t2\t1",
			"Bier Camillo\tGonzago\tPandora\t1\t-"
		},
		*snpText[] = {
			"26 SNP26\t26.26\t2626\tA\tC",
			"17 SNP17\t17.17\t1717\tC\tT"
		},
		*covText[] = {
			"FID IID\tAge\tSex",
			"Bierwetter Bernhard\t3.1\t-",
			"Wetter Praxedis\t-\t5e+20",
			"Bier Camillo\t3e-20\t-"
		},
		*yvmText[] = {
			"FID IID\tThirdnosehole",
			"Bierwetter Bernhard\t0",
			"Wetter Praxedis\t-",
			"Bier Camillo\t3.5e+15"
		},
		genotypeDataClear[] = {
			bedFileMagic[0], bedFileMagic[1], bedFileMagic[2],
			0b00010101, 0b00010101
		},
		genotypeDataSet[] = {
			bedFileMagic[0], bedFileMagic[1], bedFileMagic[2],
			0b00001011, 0b00110110
		};

	PlinkOutputTest::PlinkOutputTest () : TestSuite( "PlinkOutputTest" ) {
		addTestMethod( "PlinkOutputTest::testWrite", this, &PlinkOutputTest::testWrite );
		addTestMethod( "PlinkOutputTest::testCovariates", this, &PlinkOutputTest::testCovariates );
		addTestMethod( "PlinkOutputTest::testExtraPhenotypes", this, &PlinkOutputTest::testExtraPhenotypes );
	}

	string PlinkOutputTest::setUp () {
		return createTmpDir( tmpDirnameTemplate );
	}

	void PlinkOutputTest::tearDown (
		const string& testDirname,
		const bool withCovariates,
		const bool withExtraPhenotypes
	) {
		assert_eq(
			"Remove test SNP file failed",
			0,
			unlink( ( testDirname + "/" + filenameTrunc + snpListExtension ).c_str() )
		);
		assert_eq(
			"Remove test Individual file failed",
			0,
			unlink( ( testDirname + "/" + filenameTrunc + individualListExtension ).c_str() )
		);
		assert_eq(
			"Remove test genome file failed",
			0,
			unlink( ( testDirname + "/" + filenameTrunc + genotypeMatrixExtension ).c_str() )
		);
		if ( withCovariates ) {
			assert_eq(
				"Remove test covariate file failed",
				0,
				unlink( ( testDirname + "/" + filenameTrunc + covariateMatrixExtension ).c_str() )
			);
		}
		if ( withExtraPhenotypes ) {
			assert_eq(
				"Remove test extra phenotype file failed",
				0,
				unlink( ( testDirname + "/" + filenameTrunc + phenotypeMatrixExtension ).c_str() )
			);
		}
		removeTmpDir( testDirname );
	}

	void PlinkOutputTest::testWrite () {
		const string testDirname( setUp() );
		const string testFilenameTrunc( testDirname + "/" + filenameTrunc );
		const char * const tft = testFilenameTrunc.c_str();

		AutoMatrix genotypeMat(
			sizeof( individuals ) / sizeof( Individual ),
			sizeof( snps ) / sizeof( SNP )
		);
		genotypeMat.fill( genotype );

		{
			PlinkOutput plinkOutput(
				tft,
				genotypeMat.countRows(),
				genotypeMat.countColumns(),
				0,
				0
			);
			plinkOutput.setIndividuals( individuals );
			plinkOutput.setSnps( snps );
			plinkOutput.storeGenotypeVector( 0, genotypeMat.columnVector( 0 ) );
			plinkOutput.storeGenotypeVector( 1, genotypeMat.columnVector( 1 ) );
		}

		{
			const string individualFilename( testFilenameTrunc + individualListExtension );
			ifstream individualFile( individualFilename.c_str() );
			assert_true(
				string( "Individual file " ) + individualFilename + " has not been created.",
				individualFile.good()
			);
			assertText(
				"Individual file",
				sizeof( famTextWithoutPhenotype ) / sizeof( famTextWithoutPhenotype[0] ),
				famTextWithoutPhenotype,
				individualFile
			);
			individualFile.close();
		}

		{
			const string snpFilename( testFilenameTrunc + snpListExtension );
			ifstream snpFile( snpFilename.c_str() );
			assert_true(
				string( "SNP file " ) + snpFilename + " has not been created.",
				snpFile.good()
			);
			assertText(
				"SNP file",
				sizeof( snpText ) / sizeof( snpText[0] ),
				snpText,
				snpFile
			);
			snpFile.close();
		}

		{
			const string genotypeFilename( testFilenameTrunc + genotypeMatrixExtension );
			ifstream genotypeFile( genotypeFilename.c_str() );
			assert_true(
				string( "Genotype file " ) + genotypeFilename + " has not been created.",
				genotypeFile.good()
			);
			assertData( "Genotype file", sizeof( genotypeDataSet ), genotypeDataSet, genotypeFile );
			genotypeFile.close();
		}

		tearDown( testDirname, false, false );
	}

	void PlinkOutputTest::testCovariates () {
		const string testDirname( setUp() );
		const string testFilenameTrunc( testDirname + "/" + filenameTrunc );
		const char * const tft = testFilenameTrunc.c_str();

		AutoMatrix covariateMat(
			sizeof( individuals ) / sizeof( Individual ),
			sizeof( covariates ) / sizeof( covariates[0] )
		);
		covariateMat.fill( covariate );

		{
			PlinkOutput plinkOutput(
				tft,
				covariateMat.countRows(),
				sizeof( snps ) / sizeof( SNP ),		// test writing of "missing" if no storeGenotypeVector is called
				covariateMat.countColumns(),
				0
			);
			plinkOutput.setIndividuals( individuals );
			plinkOutput.setSnps( snps );
			plinkOutput.setCovariates( covariates );
			plinkOutput.storeCovariateVector( 0, covariateMat.columnVector( 0 ) );
			plinkOutput.storeCovariateVector( 1, covariateMat.columnVector( 1 ) );
		}

		{
			const string individualFilename( testFilenameTrunc + individualListExtension );
			ifstream individualFile( individualFilename.c_str() );
			assert_true(
				string( "Individual file " ) + individualFilename + " has not been created.",
				individualFile.good()
			);
			assertText(
				"Individual file",
				sizeof( famTextWithoutPhenotype ) / sizeof( famTextWithoutPhenotype[0] ),
				famTextWithoutPhenotype,
				individualFile
			);
			individualFile.close();
		}

		{
			const string genotypeFilename( testFilenameTrunc + genotypeMatrixExtension );
			ifstream genotypeFile( genotypeFilename.c_str() );
			assert_true(
				string( "Genotype file " ) + genotypeFilename + " has not been created.",
				genotypeFile.good()
			);
			assertData( "Genotype file", sizeof( genotypeDataClear ), genotypeDataClear, genotypeFile );
			genotypeFile.close();
		}

		{
			const string covariateFilename( testFilenameTrunc + covariateMatrixExtension );
			ifstream covariateFile( covariateFilename.c_str() );
			assert_true(
				string( "Covariate file " ) + covariateFilename + " has not been created.",
				covariateFile.good()
			);
			assertText(
				"Covariate file",
				sizeof( covText ) / sizeof( covText[0] ),
				covText,
				covariateFile
			);
			covariateFile.close();
		}

		tearDown( testDirname, true, false );
	}

	void PlinkOutputTest::testExtraPhenotypes () {
		const string testDirname( setUp() );
		const string testFilenameTrunc( testDirname + "/" + filenameTrunc );
		const char * const tft = testFilenameTrunc.c_str();

		AutoMatrix phenotypeMat(
			sizeof( individuals ) / sizeof( Individual ),
			sizeof( traits ) / sizeof( traits[0] )
		);
		phenotypeMat.fill( phenotype );

		{
			PlinkOutput plinkOutput(
				tft,
				phenotypeMat.countRows(),
				0,
				0,
				phenotypeMat.countColumns()
			);
			plinkOutput.setIndividuals( individuals );
			plinkOutput.setSnps( snps );
			plinkOutput.setTraits( traits );
			plinkOutput.storePhenotypeVector( 0, phenotypeMat.columnVector( 0 ) );
			plinkOutput.storePhenotypeVector( 1, phenotypeMat.columnVector( 1 ) );
		}

		{
			const string individualFilename( testFilenameTrunc + individualListExtension );
			ifstream individualFile( individualFilename.c_str() );
			assert_true(
				string( "Individual file " ) + individualFilename + " has not been created.",
				individualFile.good()
			);
			assertText(
				"Individual file",
				sizeof( famTextWithPhenotype ) / sizeof( famTextWithPhenotype[0] ),
				famTextWithPhenotype,
				individualFile
			);
			individualFile.close();
		}

		{
			const string genotypeFilename( testFilenameTrunc + genotypeMatrixExtension );
			ifstream genotypeFile( genotypeFilename.c_str() );
			assert_true(
				string( "Genotype file " ) + genotypeFilename + " has not been created.",
				genotypeFile.good()
			);
			assertData( "Genotype file", sizeof( bedFileMagic ), bedFileMagic, genotypeFile );
			genotypeFile.close();
		}

		{
			const string phenotypeFilename( testFilenameTrunc + phenotypeMatrixExtension );
			ifstream phenotypeFile( phenotypeFilename.c_str() );
			assert_true(
				string( "Phenotype file " ) + phenotypeFilename + " has not been created.",
				phenotypeFile.good()
			);
			assertText(
				"Phenotype file",
				sizeof( yvmText ) / sizeof( yvmText[0] ),
				yvmText,
				phenotypeFile
			);
			phenotypeFile.close();
		}

		tearDown( testDirname, false, true );
	}
}
