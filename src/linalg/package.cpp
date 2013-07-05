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

#include "package.hpp"
#include <iostream>

using namespace std;

namespace linalg {

	/** If <code>required</code>==0, 0 is returned. */
	size_t upperPowerOf2 ( const size_t required ) {
		if ( 0 >= required ) return 0;
		size_t s = required - 1, p = 1;
		while ( s ) {
			s >>= 1;
			p <<= 1;
		}
		return p;
	}


	void printVector ( const Vector& v ) {
		cout << v;
	}

	void printMatrix ( const Matrix& m ) {
		cout << m;
	}

	void printPermutation ( const Permutation& p ) {
		cout << p;
	}

}
