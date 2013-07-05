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

#include "TestMinimizable.hpp"

//#define MINIMIZATION_TESTMINIMIZABLE_DEBUG
#ifdef MINIMIZATION_TESTMINIMIZABLE_DEBUG
#include <iostream>
using namespace std;
#endif

using namespace linalg;
using namespace minimization;

namespace test {

	TestMinimizable::TestMinimizable ( const Vector& minVec )
	: Minimizable(), minVec( minVec.countDimensions() ) {
		this->minVec.copy( minVec );
	}

	size_t TestMinimizable::countDimensions () const {
		return minVec.countDimensions();
	}

	void TestMinimizable::calculateFunction (
		const Vector& x,
		double& functionResult
	) {
		AutoVector v( minVec );
		v.axpy( -1.0, x );
		functionResult = v.innerProduct( v );
		#ifdef MINIMIZATION_TESTMINIMIZABLE_DEBUG
		cout << "f(" << x << ") = " << functionResult << endl;
		#endif
	}

	void TestMinimizable::calculateDerivative (
		const Vector& x,
		Vector& derivativeResult
	) {
		derivativeResult.fill( 0.0 );
		derivativeResult.axpy( 2.0, x );
		derivativeResult.axpy( -2.0, minVec );
		#ifdef MINIMIZATION_TESTMINIMIZABLE_DEBUG
		cout << "df/dx(" << x << ") = (" << derivativeResult << ")" << endl;
		#endif
	}

}
