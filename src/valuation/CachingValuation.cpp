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

#include "CachingValuation.hpp"

using namespace std;
using namespace lookup;
using namespace util;

namespace valuation {

	CachingValuation::Retriever::Retriever ( Valuation& valuation )
		: valuation( valuation ) {
	}

	double CachingValuation::Retriever::retrieve ( const ModelIndex& modelIndex ) {
		const double value = valuation.valuate( modelIndex );
		return value;
	}

	CachingValuation::CachingValuation ( Valuation& decorated )
		: retriever( decorated ) {
	}

	double CachingValuation::valuate ( const lookup::ModelIndex& modelIndex ) {
		const double value = get( modelIndex, retriever );
		return value;
	}

}
