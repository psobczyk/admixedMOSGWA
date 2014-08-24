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
#include "../Exception.hpp"
#include <iostream>
#include <cstdlib>
#include <cassert>
#include <cmath>	// for nan(...)
#include <cstring>
#include <cerrno>

using namespace std;
using namespace linalg;
using namespace io::PlinkConstants;

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

	PlinkInput::PlinkInput ( const char* const filenameTrunc )
		:
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
			if ( ! bim.good() ) {
				throw Exception(
					"SNP file \"%s\""
					" open failed: \"%s\".",
					snpFilename.c_str(),
					strerror( errno )	// not perfectly threadsafe
				);
			}
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
			if ( ! fam.good() ) {
				throw Exception(
					"Individuals file \"%s\""
					" open failed: \"%s\".",
					idvFilename.c_str(),
					strerror( errno )	// not perfectly threadsafe
				);
			}
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
			// Prepare reading genotype matrix file
			genotypeFilename = filenameTrunc;
			genotypeFilename += genotypeMatrixExtension;
			genotypeFile.open( genotypeFilename.c_str(), ifstream::binary );
			if ( ! genotypeFile.good() ) {
				throw Exception(
					"Genotype file \"%s\""
					" open failed: \"%s\".",
					genotypeFilename.c_str(),
					strerror( errno )	// not perfectly threadsafe
				);
			}
			// check PLink magic number
			{
				assert( 0 < sizeof( bedFileMagic ) );
				char header[sizeof( bedFileMagic )];
				genotypeFile.read( header, sizeof( header ) );
				if ( genotypeFile.fail() ) {
					throw Exception(
						"Genotype file \"%s\""
						" read first %u characters failed: \"%s\".",
						genotypeFilename.c_str(),
						sizeof( header ),
						strerror( errno )	// not perfectly threadsafe
					);
				}
				const size_t orderType = sizeof( header ) - 1;
				if ( memcmp( bedFileMagic, header, orderType ) ) {
					throw Exception(
						"Genotype file \"%s\""
						" is not recognised PLINK format.",
						genotypeFilename.c_str()
					);
				}
				switch ( header[orderType] ) {
					case 0: genotypeArrayTransposition = IDV_MAJOUR; break;
					case 1: genotypeArrayTransposition = SNP_MAJOUR; break;
					default:
						throw Exception(
							"Genotype file \"%s\""
							" has unrecognised data transposition flag %u.",
							genotypeFilename.c_str(),
							header[orderType]
						);
				}
			}

			genotypeArrayStart = genotypeFile.tellg();
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

	/** Shared file state breaks thread-safety. */
	void PlinkInput::retrieveGenotypeVector ( const size_t snpIndex, Vector& vector ) {
		const size_t
			snps = countSnps(),
			idvs = countIndividuals();
		streampos position( genotypeArrayStart );
		char byte;
		switch ( genotypeArrayTransposition ) {
			case IDV_MAJOUR: {
				const size_t snpsRoundUp = snps + ( 0x3 & -snps );
				assert( 0 == ( 0x3 & snpsRoundUp ) );
				const size_t dimStep = snpsRoundUp >> 2;
				position += snpIndex >> 2;
				const unsigned int patternLowBit = ( 0x3 & snpIndex ) << 1;
				for (
					size_t idv = 0;
					idv < idvs;
					++idv,
					position += dimStep
				) {
					genotypeFile.seekg( position );
					genotypeFile.read( &byte, sizeof( byte ) );
					if ( genotypeFile.fail() ) {
						throw Exception(
							"Genotype file \"%s\""
							" read failed: %s.",
							genotypeFilename.c_str(),
							strerror( errno )	// not perfectly threadsafe
						);
					}
					const unsigned int geneticPattern = ( byte >> patternLowBit ) & 0x3;
					const double value = genotypeTranslation[geneticPattern];
					vector.set( idv, value );
				}
			} break;
			case SNP_MAJOUR: {
				const size_t idvsRoundUp = idvs + ( 0x3 & -idvs );
				assert( 0 == ( 0x3 & idvsRoundUp ) );
				position += snpIndex * ( idvsRoundUp >> 2 );
				genotypeFile.seekg( position );
				int nextBit = 8;	// position of the next bit to be processed; read next byte if >= 8
				for (
					size_t idv = 0;
					idv < idvs;
					++idv,
					nextBit += 2
				) {
					if ( 8 <= nextBit ) {
						genotypeFile.read( &byte, sizeof( byte ) );
						if ( genotypeFile.fail() ) {
							throw Exception(
								"Genotype file \"%s\""
								" read failed: %s.",
								genotypeFilename.c_str(),
								strerror( errno )	// not perfectly threadsafe
							);
						}
						nextBit = 0;
					}
					// retrieve next 2 bits; algorithm relies on a genome of chromosome PAIRS.
					const unsigned int geneticPattern = ( byte >> nextBit ) & 0x3;
					const double value = genotypeTranslation[geneticPattern];
					vector.set( idv, value );
				}
				if ( ( byte & 0xff ) >> nextBit ) {	// Expect null bits trailer, if any
					throw Exception(
						"Genotype file \"%s\""
						" violates PLINK format for SNP index %u.",
						genotypeFilename.c_str(),
						snpIndex
					);
				}
			} break;
			default:
				assert( 0 == "Invalid genotypeArrayTransposition" );
		}
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
		genotypeFile.close();
	}

}
