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

#ifndef LINALG_VECTOR_HPP
#define LINALG_VECTOR_HPP

#include <string>
#include <ostream>
#include <gsl/gsl_vector.h>

#include "Matrix.hpp"
#include "Permutation.hpp"

// for friend declaration
namespace minimization {
	class Minimizer;
}

namespace linalg {

	class Matrix;

	/** A cheap C++ wrapper for a fixed-size GSL vector view.
	* @see http://www.gnu.org/software/gsl/manual/html_node/Vector-views.html
	*/
	class Vector : protected gsl_vector_view {

		/** Used to obtain 0-dimensional vectors. */
		friend class Matrix;

		/** Used to allow GSL minimizers to access the protected base gsl_vector_view. */
		friend class minimization::Minimizer;

		// Note: copy constructor and assignment operator are allowed
		// and equal default ones for a plain Vector, which actually is a view.

		protected:

		/** Set vector from a given array.
		* The current pointer will be overwritten, so take care to avoid memory leaks!
		* If the Vector is an {@link AutoVector},
		* this should be called only with base pointing to the current <code>vector.base</code>
		* or when proper freeing of the current is guaranteed.
		* Similarly for the destructor to work,
		* an {@link AutoVector} should only get a <code>base</code> which has been malloc'ed.
		* It will be free'd by {@link AutoVector::~AutoVector()},
		* but not by {@link Vector::~Vector()}!
		* @returns whether any existing views might be invalidated by the resize
		*/
		bool gslInit ( const size_t dims, double *base, const size_t stride = 1 );

		/** Construct a vector of given dimensions from an array. */
		Vector ( const size_t dims, double *base, const size_t stride );

		public:

		/** Construct from a GSL vector view, which is the logical C equivalent.
		* Data will be shared.
		*/
		Vector ( const gsl_vector_view& view );

		/** Construct from a GSL vector, but do not take ownership.
		* Data will be shared.
		* Hence, the gsl_vector must continue to exist and eventually be properly freed.
		*/
		Vector ( const gsl_vector& vector );

		/** Compare whether this equals that. */
		bool operator== ( const Vector& that ) const;

		/** Fill all entries of the vector with a value. */
		void fill ( const double value );

		/** Fill the vector with values from an array.
		* The fill affects the logical portion of the vector,
		* i.e. the array must be at least of size {@link countDimensions()}.
		*/
		void fill ( const double *array );

		/** Copy the values from another vector.
		* Both vectors must have the same number of dimensions.
		*/
		void copy ( const Vector& that );

		/** Get the number of logical dimensions. */
		size_t countDimensions () const;

		/** Determine whether <code>this</code> is the null vector. */
		bool isNull () const;

		/** Get the element for a given dimension. */
		double get ( const size_t dim ) const;

		/** Set the element for a given dimension. */
		void set ( const size_t dim, const double value );

		/** View part of a Vector as Vector.
		* @param dim specifies the start of the sub-vector in the vector
		* @param dims specifies the number of dimensions of the sub-vector
		*/
		Vector subVector ( const size_t dim, const size_t dims );

		/** Calculate the sum of squares. */
		double sumSquares () const;

		/** Calculate the scalar product of this vector and another. */
		double innerProduct ( const Vector& that ) const;

		/** Multiply all entries with a given factor. */
		double scale ( const double scalar );

		/** Add alpha times <code>that</code> other vector to <code>this</code>. */
		void axpy ( const double alpha, const Vector& that );

		/** Convert vector into its corresponding householder vector.
		* @returns \f[\tau=2/||v||^2\f]
		* @see http://www.gnu.org/software/gsl/manual/html_node/Householder-Transformations.html
		*/
		double householderize ();

		/** Apply a householder transform from the left side to <code>this</code> vector. */
		void householderTransform ( const double tau, const Vector& householder );

		/** Performs a gemv operation:
		* add alpha A * v to beta times this;
		* or with A transposed,
		* depending on transposeA being false or true.
		*/
		void gemv (
			const double alpha,
			const Matrix& a, const bool transposeA,
			const Vector& v,
			const double beta
		);

		/** Solves the system \f[ Rx = b \f],
		* where R is the right upper triangular part of the Matrix r.
		* The result goes to <code>this</code> Vector.
		*/
		void solveR ( const Matrix& r, const Vector& b );

		/** Permute elements according to given permutation or its inverse.
		* The result goes to <code>this</code> Vector.
		* @param inverse specifies whether the inverse of <code>p</code> should be applied,
		* otherwise <code>p</code>.
		*/
		void permute ( const Permutation& p, const bool inverse );

		/** For a QR-decomposed matrix, applies Q or Q^T from the left to <code>this</code> vector.
		* Be aware that the full square matrix Q is used even when QR is not square.
		* @param tau stores the top elements of the householder vectors, as calculated by {@link Matrix::factorizeQR}.
		* @param qr stores the upper right matrix and most of top elements of the householder vectors, as calculated by {@link Matrix::factorizeQR}.
		* @param transposeQ indicates whether Q^T should be used, otherwise Q.
		*/
		void multQ ( const Vector& tau, const Matrix& qr, const bool transposeQ );

		/** Create a string representation. Mainly for debugging. */
		std::string toString () const;

		/** Placeholder destructor for {@link AutoVector::~AutoVector}. */
		virtual ~Vector ();
	};

	/** Ouput a {@link Vector}. */
	std::ostream& operator<< ( std::ostream& s, const Vector& v );

}

#endif	/* LINALG_VECTOR_HPP */
