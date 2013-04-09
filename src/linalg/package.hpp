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
