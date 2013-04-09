#include "Matrix.hpp"
#include "AutoVector.hpp"
#include <iostream>
#include <sstream>
#include <string.h>
#include <assert.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>

using namespace std;

namespace linalg {

	/** Note that GSL has rows and columns in C convention ordering. */
	bool Matrix::gslInit ( const size_t rows, const size_t cols, double *base, const size_t tda ) {
		bool invalidates
			= base != matrix.data
			|| rows < countRows()
			|| cols < countColumns()
			|| ( tda != matrix.tda && 1 < countRows() );
		// Circumvent GSL limitation to allow 0-row and/or 0-column matrices
		if ( 0 == rows || 0 == cols ) {
			assert( 0 <= rows );
			assert( 0 <= cols );
			assert( cols <= tda );
			matrix.data = base;
			matrix.size1 = rows;
			matrix.size2 = cols;
			matrix.tda = tda;
		} else {
			*static_cast<gsl_matrix_view*>( this ) = gsl_matrix_view_array_with_tda( base, rows, cols, tda );
		}
		return invalidates;
	}

	Matrix::Matrix ( const size_t rows, const size_t cols, double *base, const size_t tda ) {
		gslInit( rows, cols, base, tda );
	}

	Matrix::Matrix ( const gsl_matrix_view& view ) : gsl_matrix_view( view ) {}

	Matrix::Matrix ( const gsl_matrix& matrix ) {
		this->matrix = matrix;
	}

	void Matrix::fill ( const double value ) {
		gsl_matrix_set_all( &matrix, value );
	}

	void Matrix::fill ( const double *array ) {
		const size_t cols = countColumns();
		const size_t rows = countRows();
		if ( matrix.tda == cols || 1 >= rows ) {	// packed is copied as one lump
			const size_t blockSize = rows * cols * sizeof( double );
			memcpy( matrix.data, array, blockSize );
		} else {
			const size_t blockSize = cols * sizeof( double );
			double *target = matrix.data;
			for ( size_t row = 0; row < rows; ++row ) {
				memcpy( target, array, blockSize );
				array += cols;
				target += matrix.tda;
			}
		}
	}

	size_t Matrix::countRows () const {
		return matrix.size1;
	}

	size_t Matrix::countColumns () const {
		return matrix.size2;
	}

	double Matrix::get ( const size_t row, const size_t column ) const {
		return gsl_matrix_get( &matrix, row, column );
	}

	void Matrix::set ( const size_t row, const size_t col, const double value ) {
		gsl_matrix_set( &matrix, row, col, value );
	}

	Matrix Matrix::subMatrix ( const size_t row, const size_t col, const size_t rows, const size_t cols ) {
		// Circumvent GSL limitation to allow 0-row and/or 0-column matrices
		if ( 0 == rows || 0 == cols ) {
			assert( row + rows <= countRows() );
			assert( col + cols <= countColumns() );
			return Matrix( rows, cols, matrix.data + row * matrix.tda + col, matrix.tda );
		} else {
			return Matrix( gsl_matrix_submatrix( &matrix, row, col, rows, cols ) );
		}
	}

	const Matrix Matrix::subMatrix ( const size_t row, const size_t col, const size_t rows, const size_t cols ) const {
		return const_cast<Matrix*>( this )->subMatrix( row, col, rows, cols );
	}

	Vector Matrix::rowVector ( const size_t row ) {
		return Vector( gsl_matrix_row( &matrix, row ) );
	}

	const Vector Matrix::rowVector ( const size_t row ) const {
		return const_cast<Matrix*>( this )->rowVector( row );
	}

	Vector Matrix::columnVector ( const size_t column ) {
		return Vector( gsl_matrix_column( &matrix, column ) );
	}

	const Vector Matrix::columnVector ( const size_t column ) const {
		return const_cast<Matrix*>( this )->columnVector( column );
	}

	Vector Matrix::subdiagonalVector ( const size_t offset ) {
		// To circumvent GSL limitation to allow 0-dimensional vectors
		// This use of static data relies on the immutability of 0-dim Vectors.
		static Vector nullDimVector( 0, NULL, 0 );
		if ( 0 == countColumns() ) return nullDimVector;
		return Vector( gsl_matrix_subdiagonal( &matrix, offset ) );
	}

	const Vector Matrix::subdiagonalVector ( const size_t offset ) const {
		return const_cast<Matrix*>( this )->subdiagonalVector( offset );
	}

	Vector Matrix::diagonalVector () {
		return Vector( gsl_matrix_diagonal( &matrix ) );
	}

	const Vector Matrix::diagonalVector () const {
		return const_cast<Matrix*>( this )->diagonalVector();
	}

	Vector Matrix::superdiagonalVector ( const size_t offset ) {
		// To circumvent GSL limitation to allow 0-dimensional vectors
		// This use of static data relies on the immutability of 0-dim Vectors.
		static Vector nullDimVector( 0, NULL, 0 );
		if ( 0 == countRows() ) return nullDimVector;
		return Vector( gsl_matrix_superdiagonal( &matrix, offset ) );
	}

	const Vector Matrix::superdiagonalVector ( const size_t offset ) const {
		return const_cast<Matrix*>( this )->superdiagonalVector( offset );
	}

	void Matrix::householderTransform ( const double tau, const Vector& householder ) {
		gsl_linalg_householder_hm( tau, &householder.vector, &matrix );
	}

	void Matrix::gemm (
		const double alpha,
		const Matrix& a, const bool transposeA,
		const Matrix& b, const bool transposeB,
		const double beta
	) {
		gsl_blas_dgemm(
			transposeA ? CblasTrans : CblasNoTrans,
			transposeB ? CblasTrans : CblasNoTrans,
			alpha, &a.matrix, &b.matrix,
			beta, &matrix
		);
	}

	void Matrix::factorizeQR ( Vector& tau ) {
		gsl_linalg_QR_decomp( &matrix, &tau.vector );
	}

	void Matrix::factorizeQRP ( Vector& tau, Permutation& permutation ) {
		int sign;
		AutoVector workspace( countColumns() );
		gsl_linalg_QRPT_decomp( &matrix, &tau.vector, &permutation, &sign, &workspace.vector );
	}

	void Matrix::extractQ ( const Vector& tau, const Matrix& qr ) {
		const size_t
			rows = countRows(),
			cols = countColumns(),
			qrRows = qr.countRows(),
			qrCols = qr.countColumns(),
			tauDims = tau.countDimensions();
		assert( qrRows >= qrCols );
		assert( qrCols == tauDims );
		assert( rows == qrRows );
		assert( cols == qrCols || cols == qrRows );
		fill( 0.0 );
		Vector diagonal = diagonalVector();
		diagonal.fill( 1.0 );
		for ( size_t col = qrCols; 0 < col--; ) {
			Matrix affected = subMatrix( col, col, rows - col, cols - col );
			affected.householderTransform(
				tau.get( col ),
				qr.columnVector( col ).subVector( col, rows - col )
			);
		}
	}

	string Matrix::toString () const {
		ostringstream s;
		s << *this;
		return s.str();
	}

	Matrix::~Matrix () {}

	std::ostream& operator<< ( std::ostream& s, const Matrix& m ) {
		const size_t
			rows = m.countRows(),
			cols = m.countColumns();
		for ( size_t row = 0; row < rows; ++row ) {
			for ( size_t col = 0; col < cols; ++col ) {
				if ( 0 < col ) s << '\t';
				s << m.get( row, col );
			}
			s << endl;
		}
		return s;
	}

}
