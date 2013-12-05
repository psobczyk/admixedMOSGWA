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

	/** Parse a number from text with the understanding of dash as missing value. */
	double parseNumber( const char *text, const char **cursor ) {
		if ( '-' == *text ) {
			*cursor = text+1;
			if (
				'.' != text[1]
				&&
				! ( '0' <= **cursor && **cursor <= '9' )
			) {
				return ::nan( "missing" );
			}
		}
		const double value = strtod( text, const_cast<char**>( cursor ) );
		return value;
	}

	const char
		* const PlinkInput::snpListExtension = ".bim",
		* const PlinkInput::individualListExtension = ".fam",
		* const PlinkInput::genotypeMatrixExtension = ".bed",
		* const PlinkInput::covariateMatrixExtension = ".cov",
		* const PlinkInput::phenotypeMatrixExtension = ".yvm";

	/** Mind http://pngu.mgh.harvard.edu/~purcell/plink/binary.shtml specifies in the long translation example that the bits are swapped.
	* E.g. 10 is stored with LSB 1! Therefore <code>genotypeTranslation[1]</code> yields <code>NaN</code>, not <code>genotypeTranslation[2]</code>!
	*/
	const double PlinkInput::genotypeTranslation[] = {
		-1.0,	// 00 homocygote
		::nan( "missing" ),	// 10 missing
		0.0,	// 01 heterocygote
		+1.0	// 11 homocygote
	};


	PlinkInput::PlinkInput ( const char* const filenameTrunc )
		:
		genotypeMatrixTransposed( 0, 0 ),
		covariateMatrixTransposed( 0, 0 ),
		phenotypeMatrixTransposed( 0, 0 )
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
				const double geneticDistance = parseNumber( text, &cursor );
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
				snps.push_back( snp );
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

				const Individual individual(
					familyId,
					individualId,
					paternalId,
					maternalId,
					sex
				);

				const size_t idv = countIndividuals();
				idvIndex[ familyId + '\t' + individualId ] = idv;
				individuals.push_back( individual );

				// Parse phenotypes
				phenotypeMatrixTransposed.upSize( phenotypeMatrixTransposed.countRows(), idv+1 );
				phenotypeMatrixTransposed.columnVector( idv ).fill( ::nan( "missing" ) );
				for ( size_t traitIndex = 0; 0 != *cursor; ++traitIndex ) {
					const size_t traitCount = phenotypeMatrixTransposed.countRows();
					// if necessary, add auto-generated trait name(s)
					for ( size_t i = traits.size(); i <= traitIndex; ++i ) {
						traits.push_back( "" );
					}
					// if necessary, extend phenotype matrix
					if ( traitCount <= traitIndex ) {
						phenotypeMatrixTransposed.upSize(
							traitIndex+1,
							idv+1
						);
						phenotypeMatrixTransposed.subMatrix( traitCount, 0, traitIndex-traitCount+1, idv+1 ).fill( ::nan( "missing" ) );
					}
					text = cursor;
					if (
						'-' == *cursor
						||
						'.' == *cursor
						||
						'0' <= *cursor && *cursor <= '9'
					)  {
						const double phenotype = parseNumber( text, &cursor );
						phenotypeMatrixTransposed.set( traitIndex, idv, phenotype );
					} else {
						throw Exception(
							"Individuals file \"%s\""
							" contains unexpected character '%c'"
							" in a phenotype column.",
							idvFilename.c_str(),
							*cursor
						);
					}
				}
				assert( 0 == *cursor );
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
			// Read covariate file
			string covFilename( filenameTrunc );
			covFilename += covariateMatrixExtension;
			readExtraFile(
				"covariate",
				covFilename,
				idvIndex,
				covariates,
				covariateMatrixTransposed
			);
		}

		{
			// Read extra phenotypes file
			string yvmFilename( filenameTrunc );
			yvmFilename += phenotypeMatrixExtension;
			readExtraFile(
				"additional phenotype",
				yvmFilename,
				idvIndex,
				traits,
				phenotypeMatrixTransposed
			);
		}
	}

	void PlinkInput::readExtraFile (
		const char * const topic,
		const string& filename,
		const map<const string,size_t>& idvIndex,
		vector<string>& names,
		AutoMatrix& matrixTransposed
	) {
		const size_t
			idvs = countIndividuals(),
			offset = matrixTransposed.countRows();
		assert( 0 == offset || idvs == matrixTransposed.countColumns() );
		assert( offset == names.size() );

		// Note: C++ has no standard way to distinguish non-existence from other failures
		// http://bytes.com/topic/c/answers/129448-file-stream-open-failure-reason
		ifstream stream( filename.c_str() );
		// First line: FID IID cov1 cov2 …
		string line;
		getline( stream, line );
		if ( stream ) {
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
					"First line of %s file %s"
					" should start with \"FID\tIID\","
					" but starts with \"%s\".\n",
					topic,
					filename.c_str(),
					text
				);
			}
			size_t extras = 0;	// count the numeric data columns
			do {
				while ( ' ' == *cursor || '\t' == *cursor ) ++cursor;	// slurp whitespace
				if ( *( text = cursor ) ) {
					while ( *cursor && ' ' != *cursor && '\t' != *cursor ) ++cursor;
					const string name( text, cursor - text );
						names.push_back( name );
					++extras;
				}
			} while ( *( text = cursor ) );

			// Lines after the first
			matrixTransposed.exactSize( offset + extras, idvs );
			matrixTransposed.subMatrix( offset, 0, extras, idvs ).fill( ::nan( "missing" ) );
			vector<bool> idvEncountered( idvs );
			for (
				vector<bool>::iterator i = idvEncountered.begin();
				idvEncountered.end() != i;
				++i
			) {
				*i = false;	// initialise all flags to false
			}
			while ( ! stream.eof() ) {
				getline( stream, line );
				if ( stream.eof() ) {
					break;
				}
				const char
					* cursor = line.c_str(),	// current position in the line
					* text = cursor;	// start of an interesting text portion

				// Parse family id
				while ( *cursor && ' ' != *cursor && '\t' != *cursor ) ++cursor;
				const string familyId( text, cursor - text );

				while ( ' ' == *cursor || '\t' == *cursor ) ++cursor;	// slurp whitespace

				// Parse individual id
				text = cursor;
				while ( *cursor && ' ' != *cursor && '\t' != *cursor ) ++cursor;
				const string individualId( text, cursor - text );

				const map<string,size_t>::const_iterator idvIndexIterator = idvIndex.find( familyId + '\t' + individualId );
				if ( idvIndex.end() == idvIndexIterator ) {
					throw Exception(
						"The %s file \"%s\" contains row for \"%s %s\","
						" who is not a known individual.",
						topic,
						filename.c_str(),
						familyId.c_str(),
						individualId.c_str()
					);
				}
				const size_t idv = idvIndexIterator->second;
				if ( idvEncountered.at( idv ) ) {
					throw Exception(
						"The %s file \"%s\" contains more than one row for \"%s %s\".",
						topic,
						filename.c_str(),
						familyId.c_str(),
						individualId.c_str()
					);
				} else {
					idvEncountered.at( idv ) = true;
				}

				// Parse numeric values
				Vector valuesForIdv = matrixTransposed.columnVector( idv ).subVector( offset, extras );
				for ( size_t index = 0; index < extras; ++index ) {
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
							" in %s file \"%s\""
							" for individual \"%s %s\""
							" and %s \"%s\","
							" but encountered \"%s\".",
							topic,
							filename.c_str(),
							familyId.c_str(),
							individualId.c_str(),
							topic,
							names.at( index ).c_str(),
							cursor
						);
					}
					const double value = parseNumber( text, &cursor );
					valuesForIdv.set( index, value );
				}
				while ( ' ' == *cursor || '\t' == *cursor ) ++cursor;
				if ( *( text = cursor ) ) {
					throw Exception(
						"The %s file \"%s\" contains more than %u"
						" values for \"%s %s\": \"%s\".",
						topic,
						filename.c_str(),
						extras,
						familyId.c_str(),
						individualId.c_str(),
						text
					);
				}
			}
		}
		stream.close();
	}

	void PlinkInput::retrieveGenotypeVector ( const size_t snpIndex, Vector& vector ) {
		const Vector genotypeVector = genotypeMatrixTransposed.rowVector( snpIndex );
		vector.copy( genotypeVector );
	}

	void PlinkInput::retrievePhenotypeVector ( const size_t traitIndex, Vector& vector ) {
		const Vector phenotypeVector = phenotypeMatrixTransposed.rowVector( traitIndex );
		vector.copy( phenotypeVector );
	}

	void PlinkInput::retrieveCovariateVector ( const size_t covIndex, linalg::Vector& vector ) {
		const Vector covariateVector = covariateMatrixTransposed.rowVector( covIndex );
		vector.copy( covariateVector );
	}

	PlinkInput::~PlinkInput () {
	}

}
