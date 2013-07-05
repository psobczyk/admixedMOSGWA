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

#include "TestValuation.hpp"
#include <cmath>
#include <algorithm>

using namespace std;
using namespace lookup;

namespace test {

	TestValuation::TestValuation () {}

	TestValuation::TestValuation ( const std::map<ModelIndex,double>& mappings )
		: mappings( mappings ) {
	}

	void TestValuation::put ( const lookup::ModelIndex& modelIndex, const double value ) {
		mappings[ modelIndex ] = value;
	}

	double TestValuation::valuate ( const lookup::ModelIndex& modelIndex ) {
		map<ModelIndex,double>::iterator iterator( mappings.lower_bound( modelIndex ) );
		if ( mappings.end() == iterator || modelIndex != iterator->first ) {
			return nan( "unset" );
		}
		return iterator->second;
	}

}
