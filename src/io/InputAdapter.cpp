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
		const vector<string>& covariates,
		const vector<string>& traits
	)
	:
	IOHelper( individuals, snps, covariates, traits )
	{}

	InputAdapter::~InputAdapter () {}

}
