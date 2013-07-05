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

#ifndef LINALG_MATRIX_HPP
#define LINALG_MATRIX_HPP

#include <string>
#include <ostream>
#include <gsl/gsl_matrix.h>

#include "Vector.hpp"
#include "Permutation.hpp"

namespace linalg {

	class Vector;

	/** A cheap C++ wrapper for a fixed-size GSL matrix view.
	* @see http://www.gnu.org/software/gsl/manual/html_node/Matrix-views.html
	* @see http://stackoverflow.com/questions/5924140/how-to-multiply-to-matrix-subview-in-gsl-gnu-scientific-library
	*/
	class Matrix : protected gsl_matrix_view {

		/** Necessary to apply GSL RU solver function. */
		friend class Vector;

		// Note: copy constructor and assignment operator are allowed
		// and equal default ones for a plain Matrix, which actually is a view.

		protected:

		/** Set matrix from a given array.
		* The current pointer will be overwritten, so take care to avoid memory leaks!
		* If the Matrix is an {@link AutoMatrix},
		* this should be called only with base pointing to the current <code>matrix.base</code>
		* or when proper freeing of the current is guaranteed.
		* Similarly for the destructor to work,
		* an {@link AutoMatrix} should only get a <code>base</code> which has been malloc'ed.
		* It will be free'd by {@link AutoMatrix::~AutoMatrix()},
		* but not by {@link Matrix::~Matrix()}!
		* @returns whether any existing views might be invalidated by the resize
		*/
		bool gslInit ( const size_t rows, const size_t cols, double *base, const size_t tda );

		/** Construct a matrix of given dimensions from an array.
		* The array must be pre-allocated and freed after the Matrix ceases to exist.
		*/
		Matrix ( const size_t rows, const size_t cols, double *base, const size_t tda );

		public:

		/** Construct from a GSL matrix view, which is the logical C equivalent.
		* Data will be shared.
		*/
		Matrix ( const gsl_matrix_view& view );

		/** Construct from a GSL matrix, but do not take ownership.
		* Data will be shared.
		* Hence, the gsl_matrix must continue to exist and eventually be properly freed.
		*/
		Matrix ( const gsl_matrix& matrix );

		/** Fill all entries of the matrix with a value. */
		void fill ( const double value );

		/** Fill the matrix with values from an array.
		* The fill affects the logical portion of the matrix,
		* i.e. the array must be at least of size {@link countRows()} times {@link countColumns()}.
		* @param array contains data in C convention format, "row-major",
		* i.e. first row complete; opposite of FORTRAN "column-major" convention.
		* This is done, because GSL does it like that, see
		* <a href="http://www.gnu.org/software/gsl/manual/html_node/Matrices.html">GSL Matrix</a>.
		* Once GSL will be dropped, the class Matrix might store in "column-major"
		* for convenient communication with BLAS.
		*/
		void fill ( const double *array );

		/** Get the number of logical rows. */
		size_t countRows () const;

		/** Get the number of logical columns. */
		size_t countColumns () const;

		/** Access the element at a given position. */
		double get ( const size_t row, const size_t column ) const;

		/** Set the element for a given dimension. */
		void set ( const size_t row, const size_t col, const double value );

		/** View part of a Matrix as Matrix.
		* @param row specifies the start row of the sub-matrix in the vector
		* @param col specifies the start column of the sub-matrix in the vector
		* @param rows specifies the number of rows of the sub-matrix
		* @param cols specifies the number of columns of the sub-matrix
		*/
		Matrix subMatrix ( const size_t row, const size_t col, const size_t rows, const size_t cols );

		/** Const version of {@link subMatrix( size_t, size_t, size_t, size_t )}. */
		const Matrix subMatrix ( const size_t row, const size_t col, const size_t rows, const size_t cols ) const;

		/** Obtain access to a row as a vector view.
		* Warning: since the view directly accesses the matrix data,
		* it must not be used if the matrix has significantly changed
		* in its dimensions or data array.
		* When in doubt,
		* mind the return values of size-changing operations of {@link AutoMatrix}!
		* @see AutoMatrix::exactSize( size_t, size_t )
		* @see AutoMatrix::upSize( size_t, size_t )
		*/
		Vector rowVector ( const size_t row );

		/** Const version of {@link rowVector( size_t )}. */
		const Vector rowVector ( const size_t row ) const;

		/** Obtain access to a column as a vector view.
		* Warning: since the view directly accesses the matrix data,
		* it must not be used if the matrix has significantly changed
		* in its dimensions or data array.
		* When in doubt,
		* mind the return values of size-changing operations of {@link AutoMatrix}!
		* @see AutoMatrix::exactSize( size_t, size_t )
		* @see AutoMatrix::upSize( size_t, size_t )
		*/
		Vector columnVector ( const size_t col );

		/** Const version of {@link columnVector( size_t )}. */
		const Vector columnVector ( const size_t col ) const;

		/** Obtain access to subdiagonal as a vector view.
		* Warning: since the view directly accesses the matrix data,
		* it must not be used if the matrix has significantly changed
		* in its dimensions or data array.
		* When in doubt,
		* mind the return values of size-changing operations of {@link AutoMatrix}!
		* @param offset of subdiagonal from diagonal (0 meaning main diagonal)
		* @see AutoMatrix::exactSize( size_t, size_t )
		* @see AutoMatrix::upSize( size_t, size_t )
		*/
		Vector subdiagonalVector ( const size_t offset );

		/** Const version of {@link subdiagonalVector( size_t )}. */
		const Vector subdiagonalVector ( const size_t offset ) const;

		/** Obtain access to diagonal as a vector view.
		* Warning: since the view directly accesses the matrix data,
		* it must not be used if the matrix has significantly changed
		* in its dimensions or data array.
		* When in doubt,
		* mind the return values of size-changing operations of {@link AutoMatrix}!
		* @see AutoMatrix::exactSize( size_t, size_t )
		* @see AutoMatrix::upSize( size_t, size_t )
		*/
		Vector diagonalVector ();

		/** Const version of {@link diagonalVector()}. */
		const Vector diagonalVector () const;

		/** Obtain access to superdiagonal as a vector view.
		* Warning: since the view directly accesses the matrix data,
		* it must not be used if the matrix has significantly changed
		* in its dimensions or data array.
		* When in doubt,
		* mind the return values of size-changing operations of {@link AutoMatrix}!
		* @param offset of superdiagonal from diagonal (0 meaning main diagonal)
		* @see AutoMatrix::exactSize( size_t, size_t )
		* @see AutoMatrix::upSize( size_t, size_t )
		*/
		Vector superdiagonalVector ( const size_t offset );

		/** Const version of {@link superdiagonalVector( size_t )}. */
		const Vector superdiagonalVector ( const size_t offset ) const;

		/** Apply a householder transform from the left side to <code>this</code> matrix. */
		void householderTransform ( const double tau, const Vector& householder );

		/** Performs a gemm operation:
		* add alpha A * B to beta times this;
		* or with A and/or B transposed,
		* depending on whether transposeA and transposeB are false or true.
		*/
		void gemm (
			const double alpha,
			const Matrix& a, const bool transposeA,
			const Matrix& b, const bool transposeB,
			const double beta
		);

		/** QR-decomposition.
		* The matrix is decomposed into $$ Q\cdot R $$ with orthogonal Q and right upper triangular R.
		* The factorisation result is stored in <code>this</code> matrix and <code>tau</code>:
		* R is the right upper part of <code>this</code>; and <code>tau</code> and the left lower
		* remainder of <code>this</code> encodes the Householder vectors representing Q.
		*/
		void factorizeQR ( Vector& tau );

		/** QR-decomposition with column pivoting.
		* The matrix is decomposed into $$ Q\cdot R\cdot P^T $$
		* with orthogonal Q, right upper triangular R and permutation P.
		* The factorisation result is stored in <code>this</code> matrix,
		* <code>tau</code> and <code>permutation</code>:
		* R is the right upper part of <code>this</code>; and <code>tau</code> and the left lower
		* remainder of <code>this</code> encodes the Householder vectors representing Q.
		*/
		void factorizeQRP ( Vector& tau, Permutation& permutation );

		/** Extract Q from a given QR-decomposed matrix.
		* The orthogonal matrix Q is extracted from a given QR-decomposition,
		* in which Q is stored as Householder vectors.
		* The columns of Q will all have norm one, i.e. Q represents an isometry.
		* The number of columns of <code>this</code> matrix must equal
		* either the number of columns of qr or the number of rows of qr,
		* depending on whether only Q spanning the range of R should be generated
		* or an extension of such Q to an orthogonal square matrix.
		* @see #factorizeQR
		*/
		void extractQ ( const Vector& tau, const Matrix& qr );

		/** Create a string representation. Mainly for debugging. */
		std::string toString () const;

		/** Placeholder destructor for {@link AutoMatrix::~AutoMatrix}. */
		virtual ~Matrix ();
	};

	/** Ouput an {@link AutoMatrix}. */
	std::ostream& operator<< ( std::ostream& s, const Matrix& m );

}

#endif	/* LINALG_MATRIX_HPP */
