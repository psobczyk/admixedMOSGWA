#include "PlinkInput.hpp"
#include "../TestSuite.hpp"
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

	/** Tests the class {@link PlinkInput}. */
	class PlinkInputTest : public TestSuite {

		/** Name template for temporary test data directpry and files. */
		static const char
			* const tmpDirnameTemplate,
			* const filenameTrunc,
			* const snpListFileExtension,
			* const individualListFileExtension,
			* const genotypeMatrixFileExtension;

		/** Prepare test setup.
		* @param snpMajour determines whether SNP majour mode should be used. Else Individual majour mode.
		* @returns name of the test directory.
		*/
		string setUp ( const bool snpMajour = true );

		/** Remove test setup. */
		void tearDown ( const string& testDirname );

		/** Test interface methods. */
		void testRead ();

		public:

		/** Construct the test object. */
		PlinkInputTest ();

	} * plinkInputTest = new PlinkInputTest();	// automatically freed by unit++

	const char
		* const PlinkInputTest::tmpDirnameTemplate = "plinkinput-test.XXXXXX",
		* const PlinkInputTest::filenameTrunc = "data",
		* const PlinkInputTest::snpListFileExtension = "bim",
		* const PlinkInputTest::individualListFileExtension = "fam",
		* const PlinkInputTest::genotypeMatrixFileExtension = "bed";

	PlinkInputTest::PlinkInputTest () : TestSuite( "PlinkInputTest" ) {
		addTestMethod( "PlinkInputTest::testRead", this, &PlinkInputTest::testRead );
	}

	string PlinkInputTest::setUp ( const bool snpMajour ) {
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
			fam << "Skala Esra Reginald Rabia 2 7" << endl; // space instead tab
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

		return testDirname;
	}

	void PlinkInputTest::tearDown ( const string& testDirname ) {
		assert_eq( "Remove test SNP file failed", 0, unlink( ( testDirname + "/" + filenameTrunc + "." + snpListFileExtension ).c_str() ) );
		assert_eq( "Remove test Individual file failed", 0, unlink( ( testDirname + "/" + filenameTrunc + "." + individualListFileExtension ).c_str() ) );
		assert_eq( "Remove test genome file failed", 0, unlink( ( testDirname + "/" + filenameTrunc + "." + genotypeMatrixFileExtension ).c_str() ) );
		assert_eq( "Remove test directory failed", 0, rmdir( testDirname.c_str() ) );
	}

	void PlinkInputTest::testRead () {
		// Test both arrangements for PLink genotype data.
		for ( int snpMajour = 0; snpMajour < 2; ++snpMajour ) {
			const string testDirname( setUp( 0 < snpMajour ) );
			const string testFilenameTrunc( testDirname + "/" + filenameTrunc );
			PlinkInput plinkInput( testFilenameTrunc.c_str() );

			{
				const SNP snp0 = plinkInput.getSnp( 0 );
				assert_eq( "SNP[0].chromosome", 0, snp0.getChromosome() );
				assert_eq( "SNP[0].id", string( "Adenin" ), snp0.getSnpId() );
				assert_eq( "SNP[0].distance", 0.0, snp0.getGeneticDistance() );
				assert_eq( "SNP[0].position", 0, snp0.getBasePairPosition() );
				assert_eq( "SNP[0].allele1", 'A', snp0.getAllele1() );
				assert_eq( "SNP[0].allele2", 'T', snp0.getAllele2() );
			}

			{
				const SNP snp1 = plinkInput.getSnp( 1 );
				assert_eq( "SNP[1].chromosome", 4, snp1.getChromosome() );
				assert_eq( "SNP[1].id", string( "Cytosin" ), snp1.getSnpId() );
				assert_eq( "SNP[1].distance", 17.0, snp1.getGeneticDistance() );
				assert_eq( "SNP[1].position", 27, snp1.getBasePairPosition() );
				assert_eq( "SNP[1].allele1", 'C', snp1.getAllele1() );
				assert_eq( "SNP[1].allele2", 'G', snp1.getAllele2() );
			}

			{
				const SNP snp2 = plinkInput.getSnp( 2 );
				assert_eq( "SNP[2].chromosome", 4, snp2.getChromosome() );
				assert_eq( "SNP[2].id", string( "Guanin" ), snp2.getSnpId() );
				assert_eq( "SNP[2].distance", -0.7, snp2.getGeneticDistance() );
				assert_eq( "SNP[2].position", 271717, snp2.getBasePairPosition() );
				assert_eq( "SNP[2].allele1", 'G', snp2.getAllele1() );
				assert_eq( "SNP[2].allele2", 'C', snp2.getAllele2() );
			}

			{
				const SNP snp3 = plinkInput.getSnp( 3 );
				assert_eq( "SNP[3].chromosome", 7, snp3.getChromosome() );
				assert_eq( "SNP[3].id", string( "Thymin" ), snp3.getSnpId() );
				assert_eq( "SNP[3].distance", 7.0, snp3.getGeneticDistance() );
				assert_eq( "SNP[3].position", 1, snp3.getBasePairPosition() );
				assert_eq( "SNP[3].allele1", 'T', snp3.getAllele1() );
				assert_eq( "SNP[3].allele2", 'A', snp3.getAllele2() );
			}

			{
				const Individual idv0 = plinkInput.getIndividual( 0 );
				assert_eq( "Individual[0].familyId", string( "Fohliks" ), idv0.getFamilyID() );
				assert_eq( "Individual[0].id", string( "Flo" ), idv0.getIndividualID() );
				assert_eq( "Individual[0].paternalId", string( "Stefan" ), idv0.getPaternalID() );
				assert_eq( "Individual[0].maternalId", string( "Ludmilla" ), idv0.getMaternalID() );
				assert_eq( "Individual[0].sex", Individual::MALE, idv0.getSexCode() );
				assert_eq( "Individual[0].phenotype", -0.1, idv0.getPhenotype() );
			}

			{
				const Individual idv1 = plinkInput.getIndividual( 1 );
				assert_eq( "Individual[1].familyId", string( "Grün" ), idv1.getFamilyID() );
				assert_eq( "Individual[1].id", string( "Georg" ), idv1.getIndividualID() );
				assert_eq( "Individual[1].paternalId", string( "Karl" ), idv1.getPaternalID() );
				assert_eq( "Individual[1].maternalId", string( "Waltraud" ), idv1.getMaternalID() );
				assert_eq( "Individual[1].sex", Individual::MALE, idv1.getSexCode() );
				assert_eq( "Individual[1].phenotype", 200.0, idv1.getPhenotype() );
			}

			{
				const Individual idv2 = plinkInput.getIndividual( 2 );
				assert_eq( "Individual[2].familyId", string( "Skala" ), idv2.getFamilyID() );
				assert_eq( "Individual[2].id", string( "Esra" ), idv2.getIndividualID() );
				assert_eq( "Individual[2].paternalId", string( "Reginald" ), idv2.getPaternalID() );
				assert_eq( "Individual[2].maternalId", string( "Rabia" ), idv2.getMaternalID() );
				assert_eq( "Individual[2].sex", Individual::FEMALE, idv2.getSexCode() );
				assert_eq( "Individual[2].phenotype", 7.0, idv2.getPhenotype() );
			}

			{
				const Vector gv0 = plinkInput.getGenotypeVector( 0 );
				assert_eq( "genotypeMatrixNontransposed[0,0]", 1.0, gv0.get( 0 ) );
				assert_eq( "genotypeMatrixNontransposed[1,0]", 1.0, gv0.get( 1 ) );
				assert_eq( "genotypeMatrixNontransposed[2,0]", -1.0, gv0.get( 2 ) );
			}

			{
				const Vector gv1 = plinkInput.getGenotypeVector( 1 );
				assert_true( "genotypeMatrixNontransposed[0,1]", isnan( gv1.get( 0 ) ) );
				assert_eq( "genotypeMatrixNontransposed[1,1]", 1.0, gv1.get( 1 ) );
				assert_eq( "genotypeMatrixNontransposed[2,1]", -1.0, gv1.get( 2 ) );
			}

			{
				const Vector gv2 = plinkInput.getGenotypeVector( 2 );
				assert_eq( "genotypeMatrixNontransposed[0,2]", 0.0, gv2.get( 0 ) );
				assert_eq( "genotypeMatrixNontransposed[1,2]", 1.0, gv2.get( 1 ) );
				assert_eq( "genotypeMatrixNontransposed[2,2]", -1.0, gv2.get( 2 ) );
			}

			{
				const Vector gv3 = plinkInput.getGenotypeVector( 3 );
				assert_eq( "genotypeMatrixNontransposed[0,3]", -1.0, gv3.get( 0 ) );
				assert_eq( "genotypeMatrixNontransposed[1,3]", 1.0, gv3.get( 1 ) );
				assert_eq( "genotypeMatrixNontransposed[2,3]", -1.0, gv3.get( 2 ) );
			}

			{
				const Vector pv = plinkInput.getPhenotypeVector();
				// Mind that in the test data, phenotype has been stored from a float array.
				assert_eq( "phenotypeVector[0]", -0.1, pv.get( 0 ) );
				assert_eq( "phenotypeVector[1]", 200.0, pv.get( 1 ) );
				assert_eq( "phenotypeVector[2]", 7.0, pv.get( 2 ) );
			}

			tearDown( testDirname );
		}
	}

}
