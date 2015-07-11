/********************************************************************************
 *	This file is part of the MOSGWA program code.				*
 *	Copyright ©2015, Bernhard Bodenstorfer.				*
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

#include "CachedInput.hpp"

using namespace std;
using namespace linalg;

namespace io {

	CachedInput::CachedInput ( Input& input, const size_t limit ) : input( input ), genotypeCache( input.countIndividuals(), limit ) {}

	size_t CachedInput::countIndividuals () const {
		return input.countIndividuals();
	}

	const Individual* CachedInput::getIndividuals () const {
		return input.getIndividuals();
	}

	size_t CachedInput::countSnps () const {
		return input.countSnps();
	}

	const SNP* CachedInput::getSnps () const {
		return input.getSnps();
	}

	size_t CachedInput::countCovariates () const {
		return input.countCovariates();
	}

	const string* CachedInput::getCovariates () const {
		return input.getCovariates();
	}

	size_t CachedInput::countTraits () const {
		return input.countTraits();
	}

	const string* CachedInput::getTraits () const {
		return input.getTraits();
	}

	void CachedInput::retrieveGenotypeVector ( const size_t snpIndex, Vector& vector ) {
		if ( !genotypeCache.retrieve( snpIndex, vector ) ) {
			input.retrieveGenotypeVector( snpIndex, vector );
			genotypeCache.store( snpIndex, vector );
		}
	}

	void CachedInput::retrieveCovariateVector ( const size_t covIndex, Vector& vector ) {
		input.retrieveCovariateVector( covIndex, vector );
	}

	void CachedInput::retrievePhenotypeVector ( const size_t traitIndex, Vector& vector ) {
		input.retrievePhenotypeVector( traitIndex, vector );
	}

}
