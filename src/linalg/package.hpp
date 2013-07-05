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

#ifndef PACKAGE_LINALG_HPP
#define PACKAGE_LINALG_HPP

#include "Vector.hpp"
#include "Matrix.hpp"
#include "Permutation.hpp"

/** A simple C++ wrapper for GSL linear algebra.
* If you ask: "Why GSL and not BLAS oder LAPACK?"
* The decision was not strategic, but pragmatic.
* In the medium term, it is recommended
* to switch to a more efficient library for linear algebra.
* @author Bernhard Bodenstorfer
*/
namespace linalg {

	/** Calculate a power greater or equal <code>required</code>. */
	size_t upperPowerOf2 ( const size_t required );

	/** Output to <code>cout</code> for testing and debugging. */
	void printVector ( const Vector& v );

	/** Output to <code>cout</code> for testing and debugging. */
	void printMatrix ( const Matrix& m );

	/** Output to <code>cout</code> for testing and debugging. */
	void printPermutation ( const Permutation& p );

}

#endif	/* PACKAGE_LINALG_HPP */
