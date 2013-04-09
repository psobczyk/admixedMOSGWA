#include "QRuncher.hpp"
#include <assert.h>

using namespace linalg;

Vector QRuncher::getHouseholderVector ( const size_t col ) {
	return hrMat.columnVector( col ).subVector( col, hrMat.countRows() - col );
}

QRuncher::QRuncher ( const Vector& yVec ) : yMat( yVec.countDimensions(), 1 ), hrMat( yVec.countDimensions(), 0 ), tauVec( 0 ) {
	yMat.columnVector( 0 ).copy( yVec );
}

double QRuncher::pushColumn ( const Vector& xVec ) {
	// dimensions sanity checks
	const size_t rows = hrMat.countRows();
	assert( rows == xVec.countDimensions() );
	const size_t cols = hrMat.countColumns();
	assert( cols == tauVec.countDimensions() );
	// TODO<BB>: What should happen in unintended case cols > rows?
	assert( cols < rows );

	// append the new column
	hrMat.upSize( rows, cols + 1 );
	Vector v = hrMat.columnVector( cols );
	v.copy( xVec );

	// Apply Householder transformations for the older columns
	const Matrix oldHR = hrMat.subMatrix( 0, 0, rows, cols );
	v.multQ( tauVec, oldHR, true );

	// Prepare new column as HH transformation
	Vector w = v.subVector( cols, rows - cols );
	tauVec.upSize( cols + 1 );
	tauVec.set( cols, w.householderize() );

	// Apply the new transformation to extend yMat by one column
	// note that yMat was already 1 wider than hrMat was
	assert( cols + 1 == yMat.countColumns() );
	assert( rows == yMat.countRows() );
	yMat.upSize( yMat.countRows(), cols + 2 );
	Vector yNew = yMat.columnVector( cols + 1 );
	yNew.copy( yMat.columnVector( cols ) );
	Vector subYnew = yNew.subVector( cols, rows - cols );
	subYnew.householderTransform( tauVec.get( cols ), getHouseholderVector( cols ) );

	// the diagonal element of R-matrix measures linear independence
	return w.get( 0 );
}

bool QRuncher::popColumn () {
	const size_t cols = hrMat.countColumns();
	assert( cols == tauVec.countDimensions() );
	assert( cols + 1 == yMat.countColumns() );
	if ( 0 < cols ) {
		hrMat.upSize( hrMat.countRows(), cols - 1 );
		tauVec.upSize( cols - 1 );
		yMat.upSize( yMat.countRows(), cols );
		return true;
	} else {
		return false;
	}
}

/** The algorithm assumes that all pushed columns have been linearly independent.
* Use the return value of {@link QRuncher::pushColumn} to ensure this.
*/
double QRuncher::calculateRSS () {
	// yMat has always at least 1 column
	const size_t col = yMat.countColumns() - 1;
	// full rank of R is assumed so that the residue is easy to calculate
	const Vector residue = yMat.columnVector( col ).subVector( col, yMat.countRows() - col );
	return residue.sumSquares();
}

/** The algorithm assumes that all pushed columns have been linearly independent.
* Use the return value of {@link QRuncher::pushColumn} to ensure this.
* Author's note:
* Once the C++ standard will firmly include rvalue references,
* the cumbersome copying of the coefficient vector should be elegantly avoided
* using a move constructor for {@link linalg::AutoMatrix}.
*/
AutoVector QRuncher::calculateCoefficients () {
	const size_t cols = hrMat.countColumns();
	const Matrix r = hrMat.subMatrix( 0, 0, cols, cols );
	const Vector qTy = yMat.columnVector( cols ).subVector( 0, cols );
	AutoVector x( cols );
	x.solveR( r, qTy );
	return x;
}

/** The algorithm exploits the fact that the R-Matrix with a column deleted is almost triangular.
* Only the columns right of the deletion zone each have one potential nonzero element below the diagonal.
* These elements can be efficiently removed by a sequence of only 2-dimensional Householder transformations.
* The corresponding transformations applied to y then yield the RSS.
* To avoid accumulation of error, the modified matrix and vector are not used for further calculations, but discarded.
* The calculated RSS is meant to guide the selection of SNPs for backward steps based on
* {@link popColumn()} and {@link pushColumn( Vector& )} operations.
* Like {@link calculateRSS}, the precondition is that the pushed vectors have been linearly independent.
*/
double QRuncher::calculateSkipColumnRSS ( const size_t col ) {
	const size_t
		rows = hrMat.countRows(),
		cols = hrMat.countColumns();
	assert( col < cols );

	if ( col + 1 == cols ) {
		/** If only the last column is affected, shortcut is possible: look into yMat. */
		Vector v = yMat.columnVector( col );
		Vector w = v.subVector( col, rows - col );
		return w.sumSquares();
	} else {
		// Copy relevant parts of y
		AutoVector yTemp( rows );
		yTemp.subVector( col, rows - col ).copy( yMat.columnVector( cols ).subVector( col, rows - col ) );

		// The 2-dim Householder transform needs only a little memory
		AutoMatrix hrLocal( 2, 2 );
		Vector householderColumn = hrLocal.columnVector( 0 );
		Vector nextColumn = hrLocal.columnVector( 1 );

		// In below loop: column index k in the virtual after-column-col-delete-matrix
		// corresponds to column k+1 in hrMatrix

		// initialise "next diagonal element"
		householderColumn.set( 0, hrMat.get( col, col+1 ) );
		// loop virtual columns col .. cols - 2
		for ( size_t k = col; k + 1 < cols && k + 1 < rows; ++k ) {
			// set "below diagonal element" from hrMatrix
			householderColumn.set( 1, hrMat.get( k+1, k+1 ) );
			const double tau = householderColumn.householderize();
			yTemp.subVector( k, 2 ).householderTransform( tau, householderColumn );
			// If next column exists, prepare
			if ( k + 2 < cols ) {
				nextColumn.set( 0, hrMat.get( k, k+2 ) );
				nextColumn.set( 1, hrMat.get( k+1, k+2 ) );
				nextColumn.householderTransform( tau, householderColumn );
				// prepare diagonal element for next iteration
				householderColumn.set( 0, nextColumn.get( 1 ) );
			}
		}

		return yTemp.subVector( cols - 1, rows - cols + 1 ).sumSquares();
	}
}
