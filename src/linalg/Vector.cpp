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

#include "Vector.hpp"
#include <iostream>
#include <sstream>
#include <string.h>
#include <math.h>
#include <assert.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>

using namespace std;

namespace linalg {

	bool Vector::gslInit ( const size_t dims, double *base, const size_t stride ) {
		bool invalidates
			= base != vector.data
			|| dims < countDimensions()
			|| ( stride != vector.stride && 1 < countDimensions() );
		// Circumvent GSL limitation to allow 0-dimensional vectors
		if ( 0 == dims ) {
			assert( 0 <= stride );
			vector.data = base;
			vector.size = dims;
			vector.stride = stride;
		} else {
			*static_cast<gsl_vector_view*>( this ) = gsl_vector_view_array_with_stride( base, stride, dims );
		}
		return invalidates;
	}

	Vector::Vector ( const size_t dims, double *base, const size_t stride ) {
		gslInit( dims, base, stride );
	}

	Vector::Vector ( const gsl_vector_view& view ) : gsl_vector_view( view ) {}

	Vector::Vector ( const gsl_vector& vector ) {
		this->vector = vector;
	}

	bool Vector::operator== ( const Vector& that ) const {
		// later GSL:
		// return 0 != gsl_vector_equal( &vector, &that.vector );
		// earlier GSL: (not performance optimised, but anyway used only for testing)
		// TODO: make it more efficient
		const size_t dims = countDimensions();
		if ( dims != that.countDimensions() ) return false;
		for ( size_t dim = 0; dim < dims; ++dim ) {
			if ( get( dim ) !=  that.get( dim ) ) return false;
		}
		return true;
	}

	void Vector::fill ( const double value ) {
		gsl_vector_set_all( &vector, value );
	}

	void Vector::fill ( const double *array ) {
		// A harmless snappy const_cast circumvents a complete class "ConstVector"
		Vector that( countDimensions(), const_cast<double*>( array ), 1u );
		copy( that );
	}

	void Vector::copy ( const Vector& that ) {
		gsl_vector_memcpy( &vector, &that.vector );
	}

	size_t Vector::countDimensions () const {
		return vector.size;
	}

	bool Vector::isNull () const {
		return gsl_vector_isnull( &vector );
	}

	double Vector::get ( const size_t dim ) const {
		return gsl_vector_get( &vector, dim );
	}

	void Vector::set ( const size_t dim, const double value ) {
		gsl_vector_set( &vector, dim, value );
	}

	Vector Vector::subVector ( const size_t dim, const size_t dims ) {
		// Circumvent GSL limitation to allow 0-dimensional vectors
		if ( 0 == dims ) {
			return Vector( dims, vector.data + dim, 1u );
		} else {
			return Vector( gsl_vector_subvector( &vector, dim, dims ) );
		}
	}

	const Vector Vector::subVector ( const size_t dim, const size_t dims ) const {
		return const_cast<Vector*>( this )->subVector( dim, dims );
	}

	double Vector::sumSquares () const {
		return innerProduct( *this );
	}

	double Vector::innerProduct ( const Vector& that ) const {
		double sum;
		gsl_blas_ddot( &vector, &that.vector, &sum );
		return sum;
	}

	double Vector::scale ( const double scalar ) {
		gsl_blas_dscal( scalar, &vector );
	}

	void Vector::axpy ( const double alpha, const Vector& that ) {
		gsl_blas_daxpy( alpha, &that.vector, &vector );
	}

	double Vector::householderize () {
		const size_t dims = countDimensions();
		return 0 < dims ? gsl_linalg_householder_transform( &vector ) : nan( "unspecified" );
	}

	void Vector::householderTransform ( const double tau, const Vector& householder ) {
		gsl_linalg_householder_hv( tau, &householder.vector, &vector );
	}

	void Vector::gemv (
		const double alpha,
		const Matrix& a, const bool transposeA,
		const Vector& v,
		const double beta
	) {
		const size_t
			effRows = transposeA ? a.countColumns() : a.countRows(),
			effCols = transposeA ? a.countRows() : a.countColumns();
		if ( 0 == effRows || 0 == effCols ) {
			assert( effRows == countDimensions() );
			assert( effCols == v.countDimensions() );
			if ( 0 == effCols ) {
				scale( beta );	// adding no vectors yields beta times this vector (for case of 0-dim v)
			}
		} else {
			gsl_blas_dgemv(
				transposeA ? CblasTrans : CblasNoTrans,
				alpha, &a.matrix, &v.vector,
				beta, &vector
			);
		}
	}

	void Vector::solveR ( const Matrix& r, const Vector& b ) {
		const size_t cols = r.countColumns();
		assert( cols <= r.countRows() );
		if ( 0 == cols ) {
			assert( 0 == countDimensions() );
			// no operation in 0-dimensional space
		} else {
			Matrix ru = r.subMatrix( 0, 0, cols, cols );
			gsl_linalg_QR_Rsolve( &ru.matrix, &b.vector, &vector );
			// TODO<BB>: what is the difference to gsl_linalg_R_solve? Which should be used?
		}
	}

	void Vector::multQ ( const Vector& tau, const Matrix& qr, const bool transposeQ ) {
		const size_t
			rows = qr.countRows(),
			cols = qr.countColumns();
		assert( cols <= rows );
		if ( 0 == cols ) {
			assert( rows == countDimensions() );
			// no operation for 0-dimensional space for matrix R means identity
		} else if ( transposeQ ) {
			gsl_linalg_QR_QTvec ( &qr.matrix, &tau.vector, &vector );
		} else {
			gsl_linalg_QR_Qvec ( &qr.matrix, &tau.vector, &vector );
		}
	}

	string Vector::toString () const {
		ostringstream s;
		s << *this;
		return s.str();
	}

	Vector::~Vector () {}

	std::ostream& operator<< ( std::ostream& s, const Vector& v ) {
		for ( size_t dim = 0; dim < v.countDimensions(); ++dim ) {
			if ( 0 < dim ) s << '\t';
			s << v.get( dim );
		}
		s << endl;
		return s;
	}

}
