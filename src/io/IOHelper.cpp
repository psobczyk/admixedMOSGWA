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

#include "IOHelper.hpp"
#include <cmath>

using namespace std;

namespace io {

	IOHelper::IOHelper () {}

	IOHelper::IOHelper (
		const size_t individualCount,
		const size_t snpCount,
		const size_t covariateCount,
		const size_t traitCount
	)
	:
		individuals(
			individualCount,
			Individual( "", "", "", "", Individual::MISSING )
		),
		snps(
			snpCount,
			SNP( 0, "", ::nan( "missing" ), 0, 0, 0 )
		),
		covariates( covariateCount ),
		traits( traitCount )
	{}

	IOHelper::IOHelper (
		const vector<Individual>& individuals,
		const vector<SNP>& snps,
		const vector<string>& covariates,
		const vector<string>& traits
	)
	:
		individuals( individuals ),
		snps( snps ),
		covariates( covariates ),
		traits( traits )
	{}

	template<class T> void IOHelper::set ( const T * sources, vector<T>& targets ) {
		const size_t size = targets.size();
		for ( int index = 0; index < size; ++index ) {
			targets[index] = sources[index];
		}
	}

	size_t IOHelper::countIndividuals () const {
		return individuals.size();
	}

	const Individual * IOHelper::getIndividuals () const {
		return individuals.data();
	}

	void IOHelper::setIndividuals ( const Individual * individuals ) {
		set(
			individuals,
			this->individuals
		);
	}

	size_t IOHelper::countSnps () const {
		return snps.size();
	}

	const SNP * IOHelper::getSnps () const {
		return snps.data();
	}

	void IOHelper::setSnps ( const SNP * snps ) {
		set(
			snps,
			this->snps
		);
	}

	size_t IOHelper::countCovariates () const {
		return covariates.size();
	}

	const string * IOHelper::getCovariates () const {
		return covariates.data();
	}

	void IOHelper::setCovariates ( const std::string * covariates ) {
		set(
			covariates,
			this->covariates
		);
	}

	size_t IOHelper::countTraits () const {
		return traits.size();
	}

	const string * IOHelper::getTraits () const {
		return traits.data();
	}

	void IOHelper::setTraits ( const std::string * traits ) {
		set(
			traits,
			this->traits
		);
	}

	IOHelper::~IOHelper () {}

}
