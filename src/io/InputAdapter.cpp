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

#include "InputAdapter.hpp"

using namespace std;

namespace io {

	InputAdapter::InputAdapter () {}

	InputAdapter::InputAdapter (
		const vector<Individual>& individuals,
		const vector<SNP>& snps,
		const vector<string>& covariates
	)
	:
	individuals( individuals ),
	snps( snps ),
	covariates( covariates )
	{}

	size_t InputAdapter::countIndividuals () const {
		return individuals.size();
	}

	const Individual * InputAdapter::getIndividuals () const {
		return individuals.data();
	}

	size_t InputAdapter::countSnps () const {
		return snps.size();
	}

	const SNP * InputAdapter::getSnps () const {
		return snps.data();
	}

	size_t InputAdapter::countCovariates () const {
		return covariates.size();
	}

	const string * InputAdapter::getCovariates () const {
		return covariates.data();
	}

	InputAdapter::~InputAdapter () {}

}
