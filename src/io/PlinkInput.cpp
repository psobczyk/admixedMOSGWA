#include "PlinkInput.hpp"
#include <string>
#include <fstream>
#include <iostream>
#include <cstdlib>
#include <cassert>
#include <cmath>	// for nan(...)

using namespace std;
using namespace linalg;

namespace io {

	const char
		* const PlinkInput::snpListExtension = ".bim",
		* const PlinkInput::individualListExtension = ".fam",
		* const PlinkInput::genotypeMatrixExtension = ".bed";

	/** Mind http://pngu.mgh.harvard.edu/~purcell/plink/binary.shtml specifies in the long translation example that the bits are swapped.
	* E.g. 10 is stored with LSB 1! Therefore <code>genotypeTranslation[1]</code> yields <code>NaN</code>, not <code>genotypeTranslation[2]</code>!
	*/
	const double PlinkInput::genotypeTranslation[] = {
		-1.0,	// 00 homocygote
		nan( "missing" ),	// 10 missing
		0.0,	// 01 heterocygote
		+1.0	// 11 homocygote
	};


	PlinkInput::PlinkInput ( const char* const filenameTrunc )
		:
		genotypeMatrixTransposed( 0, 0 ),
		phenotypeVector( 0 )
	{
		string line;
		{
			// Read SNP-list file
			string snpFilename( filenameTrunc );
			snpFilename += snpListExtension;
			ifstream bim( snpFilename.c_str() );
			while ( !bim.eof() )  {
				getline( bim, line );
				if ( bim.eof() ) break;
				const char
					* text = line.c_str(),	// refers to start of text portion
					* cursor = text;	// refers to current position in text portion

				// Parse chromosome number
				unsigned long chromosomeId = 0;
				switch ( *cursor ) {
					case 'X':
					case 'x':
						++cursor;
						if ( 'Y' == *cursor || 'y' == *cursor ) {
							chromosomeId = 25;
							++cursor;
						} else {
							chromosomeId = 23;
						}
						break;
					case 'Y':
					case 'y':
						++cursor;
						chromosomeId = 24;
						break;
					case 'M':
					case 'm':
						++cursor;
						if ( 'T' == *cursor || 't' == *cursor ) {
							++cursor;
						}
						chromosomeId = 26;
						break;
					default:
						assert( '0' <= *cursor && *cursor <= '9' );
						chromosomeId = strtoul( text, const_cast<char**>( &cursor ), 10 );
				}
				assert( ' ' == *cursor || '\t' == *cursor );
				++cursor;

				// Parse SNP id
				text = cursor;
				while ( *cursor && ' ' != *cursor && '\t' != *cursor ) ++cursor;
				const string snpId( text, cursor - text );
				assert( ' ' == *cursor || '\t' == *cursor );
				++cursor;

				// Parse genetic distance
				text = cursor;
				assert( '-' == *cursor || '.' == *cursor || '0' <= *cursor && *cursor <= '9' );
				const double geneticDistance = strtod( text, const_cast<char**>( &cursor ) );
				assert( ' ' == *cursor || '\t' == *cursor );
				++cursor;

				// Parse gene position
				text = cursor;
				assert( '0' <= *cursor && *cursor <= '9' );
				const unsigned long basePairPosition = strtoul( text, const_cast<char**>( &cursor ), 10 );
				assert( ' ' == *cursor || '\t' == *cursor );
				++cursor;

				// Parse allele 1
				text = cursor;
				assert( *cursor );
				const char allele1 = *cursor++;
				assert( ' ' == *cursor || '\t' == *cursor );
				++cursor;

				// Parse allele 2
				text = cursor;
				assert( *cursor );
				const char allele2 = *cursor++;
				assert( 0 == *cursor );

				const SNP snp(
					chromosomeId,
					snpId,
					geneticDistance,
					basePairPosition,
					allele1,
					allele2
				);
				snpList.push_back( snp );
			}
			bim.close();
		}

		{
			// Read Individual-list file
			string idvFilename( filenameTrunc );
			idvFilename += individualListExtension;
			ifstream fam( idvFilename.c_str() );
			while ( !fam.eof() )  {
				getline( fam, line );
				if ( fam.eof() ) break;
				const char
					* text = line.c_str(),	// refers to start of text portion
					* cursor = text;	// refers to current position in text portion

				// Parse family id
				while ( *cursor && ' ' != *cursor && '\t' != *cursor ) ++cursor;
				const string familyId( text, cursor - text );
				assert( ' ' == *cursor || '\t' == *cursor );
				++cursor;

				// Parse individual id
				text = cursor;
				while ( *cursor && ' ' != *cursor && '\t' != *cursor ) ++cursor;
				const string individualId( text, cursor - text );
				assert( ' ' == *cursor || '\t' == *cursor );
				++cursor;

				// Parse paternal id
				text = cursor;
				while ( *cursor && ' ' != *cursor && '\t' != *cursor ) ++cursor;
				const string paternalId( text, cursor - text );
				assert( ' ' == *cursor || '\t' == *cursor );
				++cursor;

				// Parse maternal id
				text = cursor;
				while ( *cursor && ' ' != *cursor && '\t' != *cursor ) ++cursor;
				const string maternalId( text, cursor - text );
				assert( ' ' == *cursor || '\t' == *cursor );
				++cursor;

				// Parse sex code
				text = cursor;
				Individual::Sex sex = Individual::MISSING;
				switch ( *cursor++ ) {
					case '1': sex = Individual::MALE; break;
					case '2': sex = Individual::FEMALE; break;
				}
				assert( ' ' == *cursor || '\t' == *cursor );
				++cursor;

				// Parse phenotype
				text = cursor;
				assert( '-' == *cursor || '.' == *cursor || '0' <= *cursor && *cursor <= '9' );
				const double phenotype = strtod( text, const_cast<char**>( &cursor ) );
				assert( 0 == *cursor );

				const Individual individual(
					familyId,
					individualId,
					paternalId,
					maternalId,
					sex,
					phenotype
				);
				individualList.push_back( individual );
			}
			fam.close();
		}

		{
			// Read genotype matrix file
			string genFilename( filenameTrunc );
			genFilename += genotypeMatrixExtension;
			ifstream bed( genFilename.c_str(), ifstream::binary );
			char byte;
			// check PLink magic number
			bed.read( &byte, 1 );
			assert( 0x6c == byte );		// prevent non-PLink or out-of-date PLink
			bed.read( &byte, 1 );
			assert( 0x1b == byte );		// prevent non-PLink or out-of-date PLink
			// read SNP<->Individual transpose flag
			bed.read( &byte, 1 );
			assert( 0 == byte || 1 == byte );
			// Loop, if flag==1 by SNP, Individual; else by Individual, SNP.
			const size_t
				snps = countSnps(),
				idvs = countIndividuals();
			size_t snp, idv;
			const size_t
				&outerLoopLimit = byte ? snps : idvs,
				&innerLoopLimit = byte ? idvs : snps;
			size_t
				&outerLoopVariable = byte ? snp : idv,
				&innerLoopVariable = byte ? idv : snp;
			genotypeMatrixTransposed.exactSize( snps, idvs );
			for ( outerLoopVariable = 0; outerLoopVariable < outerLoopLimit; ++outerLoopVariable ) {
				// Plink wastes bits a end of "line"
				// see http://pngu.mgh.harvard.edu/~purcell/plink/binary.shtml.html
				int nextBit = 8;	// position of the next bit to be processed; read next byte if >= 8
				for ( innerLoopVariable = 0; innerLoopVariable < innerLoopLimit; ++innerLoopVariable ) {
					if ( 8 <= nextBit ) {
						bed.read( &byte, 1 );
						if ( bed.eof() ) {
							assert( snps == snp + 1 && idvs == idv + 1 );
						}
						nextBit = 0;
					}
					// retrieve next 2 bits; algorithm relies on a genome of chromosome PAIRS.
					const unsigned int geneticPattern = ( byte >> nextBit ) & 0x3;
					nextBit += 2;
					const double value = genotypeTranslation[geneticPattern];
					genotypeMatrixTransposed.set( snp, idv, value );
				}
				assert( 0 == ( byte & 0xff ) >> nextBit );		// Expect null bits trailer, if any
			}
			bed.read( &byte, 1 );
			assert( bed.eof() );
			bed.close();
		}

		{
			const size_t idvs = countIndividuals();

			// Cache phenotype vector
			phenotypeVector.exactSize( countIndividuals() );
			for ( size_t idv = 0; idv < idvs; ++idv ) {
				const Individual individual( individualList.at( idv ) );
				const double phenotype = individual.getPhenotype();
				phenotypeVector.set( idv, phenotype );
			}
		}
	}

	size_t PlinkInput::countSnps () const {
		return snpList.size();
	}

	size_t PlinkInput::countIndividuals () const {
		return individualList.size();
	}

	SNP PlinkInput::getSnp ( const size_t snpIndex ) {
		return snpList.at( snpIndex );
	}

	Individual PlinkInput::getIndividual ( const size_t individualIndex ) {
		return individualList.at( individualIndex );
	}

	Vector PlinkInput::getGenotypeVector ( const size_t snpIndex ) {
		return genotypeMatrixTransposed.rowVector( snpIndex );
	}

	Vector PlinkInput::getPhenotypeVector () {
		return phenotypeVector;
	}

	PlinkInput::~PlinkInput () {
	}

}
