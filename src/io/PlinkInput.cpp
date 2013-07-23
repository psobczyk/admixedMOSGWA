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
#include "../Exception.hpp"
#include <string>
#include <fstream>
#include <iostream>
#include <cstdlib>
#include <cassert>
#include <cmath>	// for nan(...)
#include <map>

using namespace std;
using namespace linalg;

namespace io {

	const char
		* const PlinkInput::snpListExtension = ".bim",
		* const PlinkInput::individualListExtension = ".fam",
		* const PlinkInput::genotypeMatrixExtension = ".bed",
		* const PlinkInput::covariateMatrixExtension = ".cov";

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
		covariateMatrixTransposed( 0, 0 ),
		phenotypeVector( 0 )
	{
		map<const string,size_t> idvIndex;
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

		// Read Individual-list file
		string idvFilename( filenameTrunc );
		idvFilename += individualListExtension;
		{
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
				idvIndex[ familyId + '\t' + individualId ] = individualList.size();
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
			// Cache phenotype vector

			const size_t idvs = countIndividuals();
			phenotypeVector.exactSize( countIndividuals() );
			for ( size_t idv = 0; idv < idvs; ++idv ) {
				const Individual individual( individualList.at( idv ) );
				const double phenotype = individual.getPhenotype();
				phenotypeVector.set( idv, phenotype );
			}
		}

		{
			// Read covariate file
			string covFilename( filenameTrunc );
			covFilename += covariateMatrixExtension;
			// Note: C++ has no standard way to distinguish non-existence from other failures
			// http://bytes.com/topic/c/answers/129448-file-stream-open-failure-reason
			ifstream covStream( covFilename.c_str() );
			// First line: FID IID cov1 cov2 …
			getline( covStream, line );
			if ( covStream ) {
				const char
					* text = line.c_str(),	// refers to start of text portion
					* cursor = text;	// refers to current position in text portion
				// Parse "FID"
				bool firstLineOk
				=
					'F' == *cursor++
					&&
					'I' == *cursor++
					&&
					'D' == *cursor++
					&& (
						' ' == *cursor
						||
						'\t' == *cursor
					);
				while ( ' ' == *cursor || '\t' == *cursor ) ++cursor;
				firstLineOk &=
					'I' == *cursor++
					&&
					'I' == *cursor++
					&&
					'D' == *cursor++
					&& (
						0 == *cursor	// possibly 0 covariates
						||
						' ' == *cursor
						||
						'\t' == *cursor
					);
				if ( !firstLineOk ) {
					throw Exception(
						"First line of covariate file %s"
						" should start with \"FID\tIID\","
						" but starts with \"%s\".\n",
						covFilename.c_str(),
						text
					);
				}
				do {
					while ( ' ' == *cursor || '\t' == *cursor ) ++cursor;	// slurp whitespace
					if ( *( text = cursor ) ) {
						while ( *cursor && ' ' != *cursor && '\t' != *cursor ) ++cursor;
						const string covId( text, cursor - text );
						covariateList.push_back( covId );
					}
				} while ( *( text = cursor ) );

				// Lines after the first
				const size_t
					idvs = countIndividuals(),
					covs = countCovariateVectors();
				covariateMatrixTransposed.exactSize( covs, idvs );
				covariateMatrixTransposed.fill( ::nan( "missing" ) );
				vector<bool> idvEncountered( idvs );
				for (
					vector<bool>::iterator i = idvEncountered.begin();
					idvEncountered.end() != i;
					++i
				) {
					*i = false;	// initialise all flags to false
				}
				while ( ! covStream.eof() ) {
					getline( covStream, line );
					if ( covStream.eof() ) {
						break;
					}
					const char
						* text = line.c_str(),	// refers to start of text portion
						* cursor = text;	// refers to current position in text portion

					// Parse family id
					while ( *cursor && ' ' != *cursor && '\t' != *cursor ) ++cursor;
					const string familyId( text, cursor - text );
					assert( ' ' == *cursor || '\t' == *cursor );
					// TODO: Use Exception instead of assert and permit more than one whitespace. Like in the first line: FID IID
					++cursor;

					// Parse individual id
					text = cursor;
					while ( *cursor && ' ' != *cursor && '\t' != *cursor ) ++cursor;
					const string individualId( text, cursor - text );
					assert( ' ' == *cursor || '\t' == *cursor );
					++cursor;

					const map<string,size_t>::const_iterator idvIndexIterator = idvIndex.find( familyId + '\t' + individualId );
					if ( idvIndex.end() == idvIndexIterator ) {
						throw Exception(
							"Covariate file \"%s\" contains row for \"%s %s\","
							" who is not a known individual from \"%s\".",
							covFilename.c_str(),
							familyId.c_str(),
							individualId.c_str(),
							idvFilename.c_str()
						);
					}
					const size_t idv = idvIndexIterator->second;
					if ( idvEncountered.at( idv ) ) {
						throw Exception(
							"Covariate file \"%s\" contains more than one row for \"%s %s\".",
							covFilename.c_str(),
							familyId.c_str(),
							individualId.c_str()
						);
					} else {
						idvEncountered.at( idv ) = true;
					}

					// Parse covariate values
					Vector covariatesForIdv = covariateMatrixTransposed.columnVector( idv );
					for ( size_t cov = 0; cov < covs; ++cov ) {
						while ( ' ' == *cursor || '\t' == *cursor ) ++cursor;
						text = cursor;
						if ( ! (
							'-' == *cursor
							||
							'.' == *cursor
							||
							'0' <= *cursor && *cursor <= '9'
						) ) {
							throw Exception(
								"Expecting number"
								" in covariate file \"%s\""
								" for individual \"%s %s\""
								" and covariate \"%s\","
								" but encountered \"%s\".",
								covFilename.c_str(),
								familyId.c_str(),
								individualId.c_str(),
								covariateList.at( cov ).c_str(),
								cursor
							);
						}
						const double covValue = strtod( text, const_cast<char**>( &cursor ) );
						covariatesForIdv.set( cov, covValue );
					}
					while ( ' ' == *cursor || '\t' == *cursor ) ++cursor;
					if ( *( text = cursor ) ) {
						throw Exception(
							"Covariate file \"%s\" contains more than %u"
							" covariate values for \"%s %s\".",
							covFilename.c_str(),
							covs,
							familyId.c_str(),
							individualId.c_str()
						);
					}
				}

				// Remaining individuals have covariate values all missing
				for ( size_t idv = 0; idv < idvs; ++idv ) {
					if ( ! idvEncountered.at( idv ) ) {
						Vector covariatesForIdv = covariateMatrixTransposed.columnVector( idv );
						covariatesForIdv.fill( ::nan( "missing" ) );
					}
				}
			}
			covStream.close();
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

	void PlinkInput::retrieveGenotypesIntoVector ( const size_t snpIndex, Vector& vector ) {
		vector.copy( genotypeMatrixTransposed.rowVector( snpIndex ) );
	}

	void PlinkInput::retrievePhenotypesIntoVector ( Vector& vector ) {
		vector.copy( phenotypeVector );
	}

	size_t PlinkInput::countCovariateVectors () const {
		return covariateList.size();
	}

	string PlinkInput::getCovariateName ( const size_t covIndex ) const {
		return covariateList.at( covIndex );
	}

	void PlinkInput::retrieveCovariatesIntoVector ( const size_t covIndex, linalg::Vector& vector ) {
		vector.copy( covariateMatrixTransposed.rowVector( covIndex ) );
	}

	PlinkInput::~PlinkInput () {
	}

}
