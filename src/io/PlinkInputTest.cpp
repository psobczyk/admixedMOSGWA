/********************************************************************************
 *	This file is part of the MOSGWA program code.				*
 *	Copyright ©2012–2014, Bernhard Bodenstorfer.				*
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

#include "PlinkInput.hpp"
#include "PlinkConstants.hpp"
#include "../TestSuite.hpp"
#include "../Exception.hpp"
#include "../linalg/AutoVector.hpp"
#include <unistd.h>	// for rmdir(char[])
#include <cstdlib>	// for mkdtemp(char[])
#include <cstring>
#include <cmath>	// for isnan(...)
#include <string>
#include <vector>
#include <bitset>
#include <fstream>

using namespace std;
using namespace linalg;
using namespace io;
using namespace io::PlinkConstants;
using namespace unitpp;

namespace test {

	/** Tests the class {@link io::PlinkInput}. */
	class PlinkInputTest : public TestSuite {

		/** Name template for temporary test data directory and files. */
		static const char
			* const tmpDirnameTemplate,
			* const filenameTrunc;

		/** Prepare test setup.
		* @param snpMajour determines whether SNP majour mode should be used. Else Individual majour mode.
		* @returns name of the test directory.
		*/
		string setUp (
			const bool snpMajour,
			const bool withCovariates,
			const bool withExtraPhenotypes,
			const bool withMoreData = false
		);

		/** Remove test setup. */
		void tearDown (
			const string& testDirname,
			const bool withCovariates,
			const bool withExtraPhenotypes
		 );

		/** Test {@link io::Input} interface methods. */
		void testRead ();

		/** Test {@link io::InputCo} interface extension methods. */
		void testCovariates ();

		/** Test {@link io::InputCo} interface extension methods. */
		void testExtraPhenotypes ();

		/** Test {@link io::InputCo} interface extension methods error reporting for unknown individuals. */
		void testUnknownIndividuals ();

		/** Test {@link io::InputCo} interface extension methods error reporting for duplicate individuals. */
		void testDuplicateIndividuals ();

		/** Test {@link io::InputCo} interface extension methods error reporting for individuals with too short covariate lines. */
		void testIncompleteCovariates ();

		public:

		/** Construct the test object. */
		PlinkInputTest ();

	} * plinkInputTest = new PlinkInputTest();	// automatically freed by unit++

	const char
		* const PlinkInputTest::tmpDirnameTemplate = "plinkinput-test.XXXXXX",
		* const PlinkInputTest::filenameTrunc = "data";

	PlinkInputTest::PlinkInputTest () : TestSuite( "PlinkInputTest" ) {
		addTestMethod( "PlinkInputTest::testRead", this, &PlinkInputTest::testRead );
		addTestMethod( "PlinkInputTest::testCovariates", this, &PlinkInputTest::testCovariates );
		addTestMethod( "PlinkInputTest::testExtraPhenotypes", this, &PlinkInputTest::testExtraPhenotypes );
		addTestMethod( "PlinkInputTest::testUnknownIndividuals", this, &PlinkInputTest::testUnknownIndividuals );
		addTestMethod( "PlinkInputTest::testDuplicateIndividuals", this, &PlinkInputTest::testDuplicateIndividuals );
		addTestMethod( "PlinkInputTest::testIncompleteCovariates", this, &PlinkInputTest::testIncompleteCovariates );
	}

	string PlinkInputTest::setUp (
		const bool snpMajour,
		const bool withCovariates,
		const bool withExtraPhenotypes,
		const bool withMoreData
	) {
		vector<char> tmpDirname( strlen( tmpDirnameTemplate ) + 1 );	// +1 for trailing 0
		strcpy( tmpDirname.data(), tmpDirnameTemplate );
		assert_true( "Create test directory failed", NULL != mkdtemp( tmpDirname.data() ) );
		const string testDirname( tmpDirname.data() );

		{
			const string snpFilename( testDirname + "/" + filenameTrunc + snpListExtension );
			ofstream bim( snpFilename.c_str() );
			bim << "0\tAdenin\t0\t0\tA\tT" << endl;
			bim << "4\tCytosin\t1.7e1\t27\tC\tG" << endl;
			bim << "4\tGuanin\t-7e-1 271717\tG C" << endl;
			bim << "7 Thymin 7.0 1 T A" << endl;	// space instead tab
			if ( withMoreData ) {
				bim << "8\tUracil\t0.001 1\tT A" << endl;
			}
			bim.close();
		}

		{
			const string idvFilename( testDirname + "/" + filenameTrunc + individualListExtension );
			ofstream fam( idvFilename.c_str() );
			fam << "Fohliks\tFlo\tStefan\tLudmilla\t1\t-1.0e-1" << endl;
			fam << "Grün\tGeorg\tKarl\tWaltraud\t1\t2e2" << endl;
			fam << "Skala Esra Reginald Rabia 2 7" << endl;		// space instead tab
			if ( withMoreData ) {
				fam << "Skala\tSophus\tTergonius\tScholastika\t1\t0.0" << endl;
			}
			fam.close();
		}

		{
			const string genFilename( testDirname + "/" + filenameTrunc + genotypeMatrixExtension );
			ofstream bed( genFilename.c_str(), ofstream::binary );
			bed << static_cast<char>( 0x6c );
			bed << static_cast<char>( 0x1b );
			if ( snpMajour ) {
				bed << static_cast<char>( 1 );
				{
					const bitset<8> adeninGenome( string(
						withMoreData ? "00001111" : "001111"
					) );
					bed << static_cast<char>( adeninGenome.to_ulong() );
				}
				{
					const bitset<8> cytosinGenome( string(
						withMoreData ? "01001101" : "001101"
					) );
					bed << static_cast<char>( cytosinGenome.to_ulong() );
				}
				{
					const bitset<8> guaninGenome( string(
						withMoreData ? "10001110" : "001110"
					) );
					bed << static_cast<char>( guaninGenome.to_ulong() );
				}
				{
					const bitset<8> thyminGenome( string(
						withMoreData ? "11001100" : "001100"
					) );
					bed << static_cast<char>( thyminGenome.to_ulong() );
				}
				if ( withMoreData ) {
					const bitset<8> uracilGenome( string( "00011011" ) );
					bed << static_cast<char>( uracilGenome.to_ulong() );
				}
			} else {
				bed << static_cast<char>( 0 );
				{
					const bitset<8> floGenome( string( "00100111" ) );
					bed << static_cast<char>( floGenome.to_ulong() );
					if ( withMoreData ) {
						const bitset<8> floGenome2( string( "11" ) );
						bed << static_cast<char>( floGenome2.to_ulong() );
					}
				}
				{
					const bitset<8> georgGenome( string( "11111111" ) );
					bed << static_cast<char>( georgGenome.to_ulong() );
					if ( withMoreData ) {
						const bitset<8> georgGenome2( string( "10" ) );
						bed << static_cast<char>( georgGenome2.to_ulong() );
					}
				}
				{
					const bitset<8> esraGenome( string( "00000000" ) );
					bed << static_cast<char>( esraGenome.to_ulong() );
					if ( withMoreData ) {
						const bitset<8> esraGenome2( string( "01" ) );
						bed << static_cast<char>( esraGenome2.to_ulong() );
					}
				}
				if ( withMoreData ) {
					const bitset<8> sophusGenome( string( "11100100" ) );
					bed << static_cast<char>( sophusGenome.to_ulong() );
					const bitset<8> sophusGenome2( string( "00" ) );
					bed << static_cast<char>( sophusGenome2.to_ulong() );
				}
			}
			bed.close();
		}

		if ( withCovariates ) {
			const string covFilename( testDirname + "/" + filenameTrunc + covariateMatrixExtension );
			ofstream covStream( covFilename.c_str() );
			covStream << "FID IID\tAge Salary" << endl;	// space+tab mix
			covStream << "Skala Esra 1.5e-1 -1.3" << endl;	// space instead tab
			covStream << "Grün\tGeorg\t10\t2e7" << endl;
			// Flo missing
			covStream.close();
		}

		if ( withExtraPhenotypes ) {
			const string yvmFilename( testDirname + "/" + filenameTrunc + phenotypeMatrixExtension );
			ofstream yvmStream( yvmFilename.c_str() );
			yvmStream << "FID\tIID third_nosehole\tbaldness" << endl;	// space+tab mix
			yvmStream << "Skala Esra - -1.3e-2" << endl;	// space instead tab
			yvmStream << "Fohliks\tFlo - -\t" << endl;	// two missing
			// Georg missing
			yvmStream.close();
		}

		return testDirname;
	}

	void PlinkInputTest::tearDown (
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
		assert_eq( "Remove test directory failed", 0, rmdir( testDirname.c_str() ) );
	}

	void PlinkInputTest::testRead () {
		// Test both arrangements for PLink genotype data.
		for ( int withMoreData = 0; withMoreData < 2; ++withMoreData )
		for ( int snpMajour = 0; snpMajour < 2; ++snpMajour ) {
			const string testDirname( setUp( 0 < snpMajour, false, false, withMoreData ) );
			const string testFilenameTrunc( testDirname + "/" + filenameTrunc );
			const char * const tft = testFilenameTrunc.c_str();
			PlinkInput plinkInput( tft );
			AutoVector vector( plinkInput.countIndividuals() );

			{
				const SNP * snp = plinkInput.getSnps();
				assert_eq(
					"countSnps",
					withMoreData ? 5 : 4,
					plinkInput.countSnps()
				);

				assert_eq( "snp[0].chromosome", 0, snp[0].getChromosome() );
				assert_eq( "snp[0].id", string( "Adenin" ), snp[0].getSnpId() );
				assert_eq( "snp[0].distance", 0.0, snp[0].getGeneticDistance() );
				assert_eq( "snp[0].position", 0, snp[0].getBasePairPosition() );
				assert_eq( "snp[0].allele1", 'A', snp[0].getAllele1() );
				assert_eq( "snp[0].allele2", 'T', snp[0].getAllele2() );

				assert_eq( "snp[1].chromosome", 4, snp[1].getChromosome() );
				assert_eq( "snp[1].id", string( "Cytosin" ), snp[1].getSnpId() );
				assert_eq( "snp[1].distance", 17.0, snp[1].getGeneticDistance() );
				assert_eq( "snp[1].position", 27, snp[1].getBasePairPosition() );
				assert_eq( "snp[1].allele1", 'C', snp[1].getAllele1() );
				assert_eq( "snp[1].allele2", 'G', snp[1].getAllele2() );

				assert_eq( "snp[2].chromosome", 4, snp[2].getChromosome() );
				assert_eq( "snp[2].id", string( "Guanin" ), snp[2].getSnpId() );
				assert_eq( "snp[2].distance", -0.7, snp[2].getGeneticDistance() );
				assert_eq( "snp[2].position", 271717, snp[2].getBasePairPosition() );
				assert_eq( "snp[2].allele1", 'G', snp[2].getAllele1() );
				assert_eq( "snp[2].allele2", 'C', snp[2].getAllele2() );

				assert_eq( "snp[3].chromosome", 7, snp[3].getChromosome() );
				assert_eq( "snp[3].id", string( "Thymin" ), snp[3].getSnpId() );
				assert_eq( "snp[3].distance", 7.0, snp[3].getGeneticDistance() );
				assert_eq( "snp[3].position", 1, snp[3].getBasePairPosition() );
				assert_eq( "snp[3].allele1", 'T', snp[3].getAllele1() );
				assert_eq( "snp[3].allele2", 'A', snp[3].getAllele2() );
			}

			{
				const Individual * idv = plinkInput.getIndividuals();
				assert_eq(
					"countIndividuals",
					withMoreData ? 4 : 3,
					plinkInput.countIndividuals()
				);

				assert_eq( "individual[0].familyId", string( "Fohliks" ), idv[0].getFamilyID() );
				assert_eq( "individual[0].id", string( "Flo" ), idv[0].getIndividualID() );
				assert_eq( "individual[0].paternalId", string( "Stefan" ), idv[0].getPaternalID() );
				assert_eq( "individual[0].maternalId", string( "Ludmilla" ), idv[0].getMaternalID() );
				assert_eq( "individual[0].sex", Individual::MALE, idv[0].getSexCode() );

				assert_eq( "individual[1].familyId", string( "Grün" ), idv[1].getFamilyID() );
				assert_eq( "individual[1].id", string( "Georg" ), idv[1].getIndividualID() );
				assert_eq( "individual[1].paternalId", string( "Karl" ), idv[1].getPaternalID() );
				assert_eq( "individual[1].maternalId", string( "Waltraud" ), idv[1].getMaternalID() );
				assert_eq( "individual[1].sex", Individual::MALE, idv[1].getSexCode() );

				assert_eq( "individual[2].familyId", string( "Skala" ), idv[2].getFamilyID() );
				assert_eq( "individual[2].id", string( "Esra" ), idv[2].getIndividualID() );
				assert_eq( "individual[2].paternalId", string( "Reginald" ), idv[2].getPaternalID() );
				assert_eq( "individual[2].maternalId", string( "Rabia" ), idv[2].getMaternalID() );
				assert_eq( "individual[2].sex", Individual::FEMALE, idv[2].getSexCode() );
			}

			/* Scheme of genotype array is:
					SNP	0	1	2	3	4
				IDV
				0		+	nan	0	-	+
				1		+	+	+	+	0
				2		-	-	-	-	nan
				3		-	nan	0	+	-
			*/
			{
				plinkInput.retrieveGenotypeVector( 0, vector );
				assert_eq( "genotypeMatrixNontransposed[0,0]", 1.0, vector.get( 0 ) );
				assert_eq( "genotypeMatrixNontransposed[1,0]", 1.0, vector.get( 1 ) );
				assert_eq( "genotypeMatrixNontransposed[2,0]", -1.0, vector.get( 2 ) );
				if ( withMoreData ) {
					assert_eq( "genotypeMatrixNontransposed[3,0]", -1.0, vector.get( 3 ) );
				}
			}

			{
				plinkInput.retrieveGenotypeVector( 1, vector );
				assert_true( "genotypeMatrixNontransposed[0,1]", ::isnan( vector.get( 0 ) ) );
				assert_eq( "genotypeMatrixNontransposed[1,1]", 1.0, vector.get( 1 ) );
				assert_eq( "genotypeMatrixNontransposed[2,1]", -1.0, vector.get( 2 ) );
				if ( withMoreData ) {
					assert_true( "genotypeMatrixNontransposed[3,1]", ::isnan( vector.get( 3 ) ) );
				}
			}

			{
				plinkInput.retrieveGenotypeVector( 2, vector );
				assert_eq( "genotypeMatrixNontransposed[0,2]", 0.0, vector.get( 0 ) );
				assert_eq( "genotypeMatrixNontransposed[1,2]", 1.0, vector.get( 1 ) );
				assert_eq( "genotypeMatrixNontransposed[2,2]", -1.0, vector.get( 2 ) );
				if ( withMoreData ) {
					assert_eq( "genotypeMatrixNontransposed[3,2]", 0.0, vector.get( 3 ) );
				}
			}

			{
				plinkInput.retrieveGenotypeVector( 3, vector );
				assert_eq( "genotypeMatrixNontransposed[0,3]", -1.0, vector.get( 0 ) );
				assert_eq( "genotypeMatrixNontransposed[1,3]", 1.0, vector.get( 1 ) );
				assert_eq( "genotypeMatrixNontransposed[2,3]", -1.0, vector.get( 2 ) );
				if ( withMoreData ) {
					assert_eq( "genotypeMatrixNontransposed[3,3]", 1.0, vector.get( 3 ) );
				}
			}

			if ( withMoreData ) {
				plinkInput.retrieveGenotypeVector( 4, vector );
				assert_eq( "genotypeMatrixNontransposed[0,4]", 1.0, vector.get( 0 ) );
				assert_eq( "genotypeMatrixNontransposed[1,4]", 0.0, vector.get( 1 ) );
				assert_true( "genotypeMatrixNontransposed[2,4]", ::isnan( vector.get( 2 ) ) );
				assert_eq( "genotypeMatrixNontransposed[3,4]", -1.0, vector.get( 3 ) );
			}

			{
				plinkInput.retrievePhenotypeVector( 0, vector );
				assert_eq( "phenotypeVector[0]", -0.1, vector.get( 0 ) );
				assert_eq( "phenotypeVector[1]", 200.0, vector.get( 1 ) );
				assert_eq( "phenotypeVector[2]", 7.0, vector.get( 2 ) );
			}

			assert_eq( "No covariates", 0, plinkInput.countCovariates() );

			tearDown( testDirname, false, false );
		}
	}

	void PlinkInputTest::testCovariates () {
		const string testDirname( setUp( false, true, false ) );
		const string testFilenameTrunc( testDirname + "/" + filenameTrunc );
		const char * const tft = testFilenameTrunc.c_str();
		PlinkInput plinkInput( tft );
		assert_eq( "Two covariates", 2, plinkInput.countCovariates() );
		AutoVector vector( plinkInput.countIndividuals() );

		{
			const string * cov = plinkInput.getCovariates();

			assert_eq( "cov[0].name", string( "Age" ), cov[0] );
			assert_eq( "cov[1].name", string( "Salary" ), cov[1] );
		}

		{
			plinkInput.retrieveCovariateVector( 0, vector );
			assert_true( "covariateMatrixNontransposed[0,0] NaN", ::isnan( vector.get( 0 ) ) );
			assert_eq( "covariateMatrixNontransposed[1,0]", 10.0, vector.get( 1 ) );
			assert_eq( "covariateMatrixNontransposed[2,0]", 0.15, vector.get( 2 ) );
		}

		{
			plinkInput.retrieveCovariateVector( 1, vector );
			assert_true( "covariateMatrixNontransposed[0,1]", ::isnan( vector.get( 0 ) ) );
			assert_eq( "covariateMatrixNontransposed[1,1]", 20000000.0, vector.get( 1 ) );
			assert_eq( "covariateMatrixNontransposed[2,1]", -1.3, vector.get( 2 ) );
		}

		tearDown( testDirname, true, false );
	}

	void PlinkInputTest::testExtraPhenotypes () {
		const string testDirname( setUp( false, false, true ) );
		const string testFilenameTrunc( testDirname + "/" + filenameTrunc );
		const char * const tft = testFilenameTrunc.c_str();
		PlinkInput plinkInput( tft );
		assert_eq( "One normal and two extra phenotypes", 3, plinkInput.countTraits() );
		AutoVector vector( plinkInput.countIndividuals() );

		{
			const string * trait = plinkInput.getTraits();

			assert_eq( "trait[0].name from FAM file", string( "" ), trait[0] );
			assert_eq( "trait[1].name", string( "third_nosehole" ), trait[1] );
			assert_eq( "trait[2].name", string( "baldness" ), trait[2] );
		}

		{
			plinkInput.retrievePhenotypeVector( 0, vector );
			assert_eq( "phenotypeVector[0][0]", -0.1, vector.get( 0 ) );
			assert_eq( "phenotypeVector[0][1]", 200.0, vector.get( 1 ) );
			assert_eq( "phenotypeVector[0][2]", 7.0, vector.get( 2 ) );
		}

		{
			plinkInput.retrievePhenotypeVector( 1, vector );
			assert_true( "phenotypeVector[1][0] NaN", ::isnan( vector.get( 0 ) ) );
			assert_true( "phenotypeVector[1][1] NaN", ::isnan( vector.get( 1 ) ) );
			assert_true( "phenotypeVector[1][2] NaN", ::isnan( vector.get( 2 ) ) );
		}

		{
			plinkInput.retrievePhenotypeVector( 2, vector );
			assert_true( "phenotypeVector[2][0] NaN", ::isnan( vector.get( 0 ) ) );
			assert_true( "phenotypeVector[2][1] NaN", ::isnan( vector.get( 1 ) ) );
			assert_eq( "phenotypeVector[2][2]", -0.013, vector.get( 2 ) );
		}

		tearDown( testDirname, false, true );
	}

	void PlinkInputTest::testUnknownIndividuals () {
		const string testDirname( setUp( false, true, true ) );
		const string testFilenameTrunc( testDirname + "/" + filenameTrunc );
		const string covFilename( testDirname + "/" + filenameTrunc + covariateMatrixExtension );
		ofstream covStream( covFilename.c_str(), ios::app );
		covStream << "Grau\tGustav\t1\t0" << endl;
		covStream.close();

		const char * const tft = testFilenameTrunc.c_str();
		try {
			PlinkInput plinkInput( tft );
			assert_fail( "Exception should be raised for unknown individual" );
		} catch ( Exception e ) {
		}

		tearDown( testDirname, true, true );
	}

	void PlinkInputTest::testDuplicateIndividuals () {
		const string testDirname( setUp( false, true, true ) );
		const string testFilenameTrunc( testDirname + "/" + filenameTrunc );
		const string covFilename( testDirname + "/" + filenameTrunc + covariateMatrixExtension );
		ofstream covStream( covFilename.c_str(), ios::app );
		covStream << "Skala\tEsra\t1.5e-1\t-1.3" << endl;
		covStream.close();

		const char * const tft = testFilenameTrunc.c_str();
		try {
			PlinkInput plinkInput( tft );
			assert_fail( "Exception should be raised for duplicate individual" );
		} catch ( Exception e ) {
		}

		tearDown( testDirname, true, true );
	}

	void PlinkInputTest::testIncompleteCovariates () {
		const string testDirname( setUp( false, true, false ) );
		const string testFilenameTrunc( testDirname + "/" + filenameTrunc );
		const string covFilename( testDirname + "/" + filenameTrunc + covariateMatrixExtension );
		ofstream covStream( covFilename.c_str(), ios::app );
		covStream << "Fohliks Flo\t-10" << endl;
		covStream.close();

		const char * const tft = testFilenameTrunc.c_str();
		try {
			PlinkInput plinkInput( tft );
			assert_fail( "Exception should be raised for individual with too few covariate values" );
		} catch ( Exception e ) {
		}

		tearDown( testDirname, true, false );
	}

}
