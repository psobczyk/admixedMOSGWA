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

#include "RandomInput.hpp"
#include <cstdlib>
#include <sstream>
#include <cassert>

using namespace std;
using namespace linalg;

namespace io {

	string RandomInput::prefixedString ( const char * prefix, const size_t i ) {
		stringstream s;
		s << prefix << i;
		return s.str();
	}

	RandomInput::RandomInput (
		const size_t individualCount,
		const size_t snpCount,
		const size_t covariateCount,
		const size_t traitCount
	)
	{
		const char
			*idvPrefix = "I",
			*famPrefix = "F",
			*snpPrefix = "S",
			*covPrefix = "C",
			*traitPrefix = "T";

		// Generate individuals
		individuals.reserve( individualCount );
		for ( size_t idvIndex = 0; idvIndex < individualCount; ++idvIndex ) {
			const Individual individual(
				prefixedString( famPrefix, idvIndex ),
				prefixedString( idvPrefix, idvIndex ),
				"Dad",
				"Mom",
				Individual::MISSING
			);
			individuals.push_back( individual );
		}

		// Generate SNPs
		snps.reserve( snpCount );
		for ( size_t snpIndex = 0; snpIndex < snpCount; ++snpIndex ) {
			const SNP snp(
				"1",
				prefixedString( snpPrefix, snpIndex ),
				snpIndex,
				snpIndex,
				'A',
				'T'
			);
			snps.push_back( snp );
		}

		// Generate covariates
		covariates.reserve( covariateCount );
		for ( size_t covIndex = 0; covIndex < covariateCount; ++covIndex ) {
			const string covariate(
				prefixedString( covPrefix, covIndex )
			);
			covariates.push_back( covariate );
		}

		// Generate traits
		traits.reserve( traitCount );
		for ( size_t traitIndex = 0; traitIndex < traitCount; ++traitIndex ) {
			const string trait(
				prefixedString( traitPrefix, traitIndex )
			);
			traits.push_back( trait );
		}
	}

	void RandomInput::retrieveVector ( const size_t indexSeed, Vector& v, const int min, const int cases ) {
		const size_t individualCount = countIndividuals();
		assert( v.countDimensions() == individualCount );
		srand( indexSeed );
		for ( int idv = 0; idv < individualCount; ++idv ) {
			v.set( idv, min + rand() % cases );
		}
	}

	void RandomInput::retrieveGenotypeVector ( const size_t snpIndex, Vector& v ) {
		const size_t offset = 0;
		retrieveVector( offset + snpIndex, v, -1, 3 );
	}

	void RandomInput::retrieveCovariateVector ( const size_t covIndex, Vector& v ) {
		const size_t offset = countSnps();
		retrieveVector( offset + covIndex, v, 0, 2 );
	}

	void RandomInput::retrievePhenotypeVector ( const size_t traitIndex, Vector& v ) {
		const size_t offset = countSnps() + countCovariates();
		retrieveVector( offset + traitIndex, v, 0, 2 );
	}

}
