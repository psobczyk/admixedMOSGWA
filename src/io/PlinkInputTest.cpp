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

#include "PlinkInput.hpp"
#include "../TestSuite.hpp"
#include "../Exception.hpp"
#include <unistd.h>	// for rmdir(char[])
#include <cstdlib>	// for mkdtemp(char[])
#include <cstring>
#include <cmath>	// for isnan(...)
#include <string>
#include <vector>
#include <bitset>
#include <fstream>

using namespace std;
using namespace io;
using namespace linalg;
using namespace unitpp;

namespace test {

	/** Tests the class {@link io::PlinkInput}. */
	class PlinkInputTest : public TestSuite {

		/** Name template for temporary test data directory and files. */
		static const char
			* const tmpDirnameTemplate,
			* const filenameTrunc,
			* const snpListFileExtension,
			* const individualListFileExtension,
			* const genotypeMatrixFileExtension,
			* const covariateMatrixFileExtension;

		/** Prepare test setup.
		* @param snpMajour determines whether SNP majour mode should be used. Else Individual majour mode.
		* @returns name of the test directory.
		*/
		string setUp ( const bool snpMajour, const bool withCovariates );

		/** Remove test setup. */
		void tearDown ( const string& testDirname, const bool withCovariates );

		/** Test {@link io::Input} interface methods. */
		void testRead ();

		/** Test {@link io::InputCo} interface extension methods. */
		void testCovariates ();

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
		* const PlinkInputTest::filenameTrunc = "data",
		* const PlinkInputTest::snpListFileExtension = "bim",
		* const PlinkInputTest::individualListFileExtension = "fam",
		* const PlinkInputTest::genotypeMatrixFileExtension = "bed",
		* const PlinkInputTest::covariateMatrixFileExtension = "cov";

	PlinkInputTest::PlinkInputTest () : TestSuite( "PlinkInputTest" ) {
		addTestMethod( "PlinkInputTest::testRead", this, &PlinkInputTest::testRead );
		addTestMethod( "PlinkInputTest::testCovariates", this, &PlinkInputTest::testCovariates );
		addTestMethod( "PlinkInputTest::testUnknownIndividuals", this, &PlinkInputTest::testUnknownIndividuals );
		addTestMethod( "PlinkInputTest::testDuplicateIndividuals", this, &PlinkInputTest::testDuplicateIndividuals );
		addTestMethod( "PlinkInputTest::testIncompleteCovariates", this, &PlinkInputTest::testIncompleteCovariates );
	}

	string PlinkInputTest::setUp ( const bool snpMajour, const bool withCovariates ) {
		vector<char> tmpDirname( strlen( tmpDirnameTemplate ) + 1 );	// +1 for trailing 0
		strcpy( tmpDirname.data(), tmpDirnameTemplate );
		assert_true( "Create test directory failed", NULL != mkdtemp( tmpDirname.data() ) );
		const string testDirname( tmpDirname.data() );

		{
			const string snpFilename( testDirname + "/" + filenameTrunc + "." + snpListFileExtension );
			ofstream bim( snpFilename.c_str() );
			bim << "0\tAdenin\t0\t0\tA\tT" << endl;
			bim << "4\tCytosin\t1.7e1\t27\tC\tG" << endl;
			bim << "4\tGuanin\t-7e-1 271717\tG C" << endl;
			bim << "7 Thymin 7.0 1 T A" << endl;	// space instead tab
			bim.close();
		}

		{
			const string idvFilename( testDirname + "/" + filenameTrunc + "." + individualListFileExtension );
			ofstream fam( idvFilename.c_str() );
			fam << "Fohliks\tFlo\tStefan\tLudmilla\t1\t-1.0e-1" << endl;
			fam << "Grün\tGeorg\tKarl\tWaltraud\t1\t2e2" << endl;
			fam << "Skala Esra Reginald Rabia 2 7" << endl;		// space instead tab
			fam.close();
		}

		{
			const string genFilename( testDirname + "/" + filenameTrunc + "." + genotypeMatrixFileExtension );
			ofstream bed( genFilename.c_str(), ofstream::binary );
			bed << static_cast<char>( 0x6c );
			bed << static_cast<char>( 0x1b );
			if ( snpMajour ) {
				bed << static_cast<char>( 1 );
				{
					const bitset<8> adeninGenome( string( "001111" ) );
					bed << static_cast<char>( adeninGenome.to_ulong() );
				}
				{
					const bitset<8> cytosinGenome( string( "001101" ) );
					bed << static_cast<char>( cytosinGenome.to_ulong() );
				}
				{
					const bitset<8> guaninGenome( string( "001110" ) );
					bed << static_cast<char>( guaninGenome.to_ulong() );
				}
				{
					const bitset<8> thyminGenome( string( "001100" ) );
					bed << static_cast<char>( thyminGenome.to_ulong() );
				}
			} else {
				bed << static_cast<char>( 0 );
				{
					const bitset<8> floGenome( string( "00100111" ) );
					bed << static_cast<char>( floGenome.to_ulong() );
				}
				{
					const bitset<8> georgGenome( string( "11111111" ) );
					bed << static_cast<char>( georgGenome.to_ulong() );
				}
				{
					const bitset<8> esraGenome( string( "00000000" ) );
					bed << static_cast<char>( esraGenome.to_ulong() );
				}
			}
			bed.close();
		}

		if ( withCovariates ) {
			const string covFilename( testDirname + "/" + filenameTrunc + "." + covariateMatrixFileExtension );
			ofstream covStream( covFilename.c_str() );
			covStream << "FID IID\tAge Salary" << endl;		// space+tab mix
			covStream << "Skala Esra 1.5e-1 -1.3" << endl;	// space instead tab
			covStream << "Grün\tGeorg\t10\t2e7" << endl;
			// Flo missing
			covStream.close();
		}

		return testDirname;
	}

	void PlinkInputTest::tearDown ( const string& testDirname, const bool withCovariates ) {
		assert_eq( "Remove test SNP file failed", 0, unlink( ( testDirname + "/" + filenameTrunc + "." + snpListFileExtension ).c_str() ) );
		assert_eq( "Remove test Individual file failed", 0, unlink( ( testDirname + "/" + filenameTrunc + "." + individualListFileExtension ).c_str() ) );
		assert_eq( "Remove test genome file failed", 0, unlink( ( testDirname + "/" + filenameTrunc + "." + genotypeMatrixFileExtension ).c_str() ) );
		if ( withCovariates ) {
			assert_eq( "Remove test covariate file failed", 0, unlink( ( testDirname + "/" + filenameTrunc + "." + covariateMatrixFileExtension ).c_str() ) );
		}
		assert_eq( "Remove test directory failed", 0, rmdir( testDirname.c_str() ) );
	}

	void PlinkInputTest::testRead () {
		// Test both arrangements for PLink genotype data.
		for ( int snpMajour = 0; snpMajour < 2; ++snpMajour ) {
			const string testDirname( setUp( 0 < snpMajour, false ) );
			const string testFilenameTrunc( testDirname + "/" + filenameTrunc );
			const char * const tft = testFilenameTrunc.c_str();
			PlinkInput plinkInput( tft );
			AutoVector vector( plinkInput.countIndividuals() );

			{
				const SNP * snp = plinkInput.getSnps();
				assert_eq( "countSnps", 4, plinkInput.countSnps() );

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
				assert_eq( "countIndividuals", 3, plinkInput.countIndividuals() );

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

			{
				plinkInput.retrieveGenotypeVector( 0, vector );
				assert_eq( "genotypeMatrixNontransposed[0,0]", 1.0, vector.get( 0 ) );
				assert_eq( "genotypeMatrixNontransposed[1,0]", 1.0, vector.get( 1 ) );
				assert_eq( "genotypeMatrixNontransposed[2,0]", -1.0, vector.get( 2 ) );
			}

			{
				plinkInput.retrieveGenotypeVector( 1, vector );
				assert_true( "genotypeMatrixNontransposed[0,1]", ::isnan( vector.get( 0 ) ) );
				assert_eq( "genotypeMatrixNontransposed[1,1]", 1.0, vector.get( 1 ) );
				assert_eq( "genotypeMatrixNontransposed[2,1]", -1.0, vector.get( 2 ) );
			}

			{
				plinkInput.retrieveGenotypeVector( 2, vector );
				assert_eq( "genotypeMatrixNontransposed[0,2]", 0.0, vector.get( 0 ) );
				assert_eq( "genotypeMatrixNontransposed[1,2]", 1.0, vector.get( 1 ) );
				assert_eq( "genotypeMatrixNontransposed[2,2]", -1.0, vector.get( 2 ) );
			}

			{
				plinkInput.retrieveGenotypeVector( 3, vector );
				assert_eq( "genotypeMatrixNontransposed[0,3]", -1.0, vector.get( 0 ) );
				assert_eq( "genotypeMatrixNontransposed[1,3]", 1.0, vector.get( 1 ) );
				assert_eq( "genotypeMatrixNontransposed[2,3]", -1.0, vector.get( 2 ) );
			}

			{
				plinkInput.retrievePhenotypeVector( vector );
				assert_eq( "phenotypeVector[0]", -0.1, vector.get( 0 ) );
				assert_eq( "phenotypeVector[1]", 200.0, vector.get( 1 ) );
				assert_eq( "phenotypeVector[2]", 7.0, vector.get( 2 ) );
			}

			assert_eq( "No covariates", 0, plinkInput.countCovariates() );

			tearDown( testDirname, false );
		}
	}

	void PlinkInputTest::testCovariates () {
		const string testDirname( setUp( false, true ) );
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

		tearDown( testDirname, true );
	}

	void PlinkInputTest::testUnknownIndividuals () {
		const string testDirname( setUp( false, true ) );
		const string testFilenameTrunc( testDirname + "/" + filenameTrunc );
		const string covFilename( testDirname + "/" + filenameTrunc + "." + covariateMatrixFileExtension );
		ofstream covStream( covFilename.c_str(), ios::app );
		covStream << "Grau\tGustav\t1\t0" << endl;
		covStream.close();

		const char * const tft = testFilenameTrunc.c_str();
		try {
			PlinkInput plinkInput( tft );
			assert_fail( "Exception should be raised for unknown individual" );
		} catch ( Exception e ) {
		}

		tearDown( testDirname, true );
	}

	void PlinkInputTest::testDuplicateIndividuals () {
		const string testDirname( setUp( false, true ) );
		const string testFilenameTrunc( testDirname + "/" + filenameTrunc );
		const string covFilename( testDirname + "/" + filenameTrunc + "." + covariateMatrixFileExtension );
		ofstream covStream( covFilename.c_str(), ios::app );
		covStream << "Skala\tEsra\t1.5e-1\t-1.3" << endl;
		covStream.close();

		const char * const tft = testFilenameTrunc.c_str();
		try {
			PlinkInput plinkInput( tft );
			assert_fail( "Exception should be raised for duplicate individual" );
		} catch ( Exception e ) {
		}

		tearDown( testDirname, true );
	}

	void PlinkInputTest::testIncompleteCovariates () {
		const string testDirname( setUp( false, true ) );
		const string testFilenameTrunc( testDirname + "/" + filenameTrunc );
		const string covFilename( testDirname + "/" + filenameTrunc + "." + covariateMatrixFileExtension );
		ofstream covStream( covFilename.c_str(), ios::app );
		covStream << "Fohliks Flo\t-10" << endl;
		covStream.close();

		const char * const tft = testFilenameTrunc.c_str();
		try {
			PlinkInput plinkInput( tft );
			assert_fail( "Exception should be raised for individual with too few covariate values" );
		} catch ( Exception e ) {
		}

		tearDown( testDirname, true );
	}

}
