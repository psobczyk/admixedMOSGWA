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

#include "QRuncher.hpp"
#include "linalg/AutoPermutation.hpp"
#include <assert.h>

using namespace std;
using namespace linalg;

bool QRuncher::isLinearlyIndependent ( const size_t col ) const {
	const size_t rank = ranks.at( col );
	return 0 < col
		? ranks.at( col - 1 ) < rank
		: 0 < rank;
}

size_t QRuncher::getRank () const {
	return 0 < ranks.size() ? ranks.back() : 0;
}

Vector QRuncher::getHouseholderVector ( const size_t col ) {
	const size_t
		rows = hrMat.countRows(),
		rank = ranks.at( col );
	assert( 0 < rank );	// otherwise this is a linearly dependent column and does not have a HH vector
	assert( rank <= rows );
	Vector v = hrMat.columnVector( col );
	Vector w = v.subVector( rank - 1, rows - ( rank - 1 ) );
	return w;
}

QRuncher::QRuncher ( const Vector& yVec ) : yMat( yVec.countDimensions(), 1 ), hrMat( yVec.countDimensions(), 0 ), tauVec( 0 ) {
	yMat.columnVector( 0 ).copy( yVec );
}

bool QRuncher::pushColumn ( const Vector& xVec ) {
	// dimensions sanity checks
	const size_t
		rows = hrMat.countRows(),
		cols = hrMat.countColumns(),
		rank = getRank();
	assert( rows == xVec.countDimensions() );
	assert( cols == tauVec.countDimensions() );
	assert( cols == ranks.size() );
	assert( rank <= cols );

	// append the new column
	hrMat.upSize( rows, cols + 1 );
	Vector v = hrMat.columnVector( cols );
	v.copy( xVec );

	// Apply Householder transforms for the older columns
	for ( size_t col = 0; col < cols; ++col ) {
		if ( isLinearlyIndependent( col ) ) {
			// apply its householder transform
			// Note: tauVec holds unused entries for the other,
			// i.e. linearly dependent columns.
			const size_t r = ranks.at( col );
			assert( 0 < r );	// otherwise not linearly independent and does not have HH vector
			const double tau = tauVec.get( col );
			const Vector hh = getHouseholderVector( col );
			Vector w = v.subVector( r - 1, rows - ( r - 1 ) );
			w.householderTransform( tau, hh );
		}
	}

	// Prepare to apply the new column's Householder transform.
	Vector w = v.subVector( rank, rows - rank );

	// Extend tauVec by one dimension.
	tauVec.upSize( cols + 1 );

	// Extend yMat by one column.
	// Note that yMat was already 1 wider than hrMat was
	assert( cols + 1 == yMat.countColumns() );
	assert( rows == yMat.countRows() );
	yMat.upSize( yMat.countRows(), cols + 2 );
	Vector yNew = yMat.columnVector( cols + 1 );
	yNew.copy( yMat.columnVector( cols ) );
	// Slightly suboptimal is that tauVec and yMat are extended in the linearly dependent case, too.

	// The below comparison looks numerically unstable, and yes, it is.
	// This reflects the instability of matrix rank.
	if ( w.isNull() ) {
		// If the new column was linearly dependent,
		// rank has not changed by adding the new column.
		ranks.push_back( rank );
		return false;
	} else {
		ranks.push_back( rank + 1 );
		// If the new column was linearly independent,
		// prepare it as Householder transform
		const double tau = w.householderize();
		tauVec.set( cols, tau );
		// and apply apply that.
		Vector subYnew = yNew.subVector( rank, rows - rank );
		subYnew.householderTransform( tauVec.get( cols ), w );
		return true;
	}
}

bool QRuncher::popColumn () {
	const size_t cols = hrMat.countColumns();
	assert( cols == ranks.size() );
	assert( cols == tauVec.countDimensions() );
	assert( cols + 1 == yMat.countColumns() );
	if ( 0 < cols ) {
		hrMat.upSize( hrMat.countRows(), cols - 1 );
		tauVec.upSize( cols - 1 );
		ranks.pop_back();
		yMat.upSize( yMat.countRows(), cols );
		return true;
	} else {
		return false;
	}
}

double QRuncher::calculateRSS () const {
	// yMat has always at least 1 column
	const size_t
		rows = yMat.countRows(),
		cols = yMat.countColumns(),
		rank = getRank();
	assert( rank <= rows );
	Vector v = const_cast<AutoMatrix&>( yMat ).columnVector( cols - 1 );
	const Vector w = v.subVector( rank, rows - rank );
	const double rss = w.sumSquares();
	return rss;
}

/** The algorithm assumes that all pushed columns have been linearly independent.
* Use the return value of {@link QRuncher::pushColumn} to ensure this.
* Author's note:
* Once the C++ standard will firmly include rvalue references,
* the cumbersome copying of the coefficient vector should be elegantly avoided
* using a move constructor for {@link linalg::AutoMatrix}.
*/
AutoVector QRuncher::calculateCoefficients () const {
	const size_t
		cols = hrMat.countColumns(),
		rank = getRank();
	Matrix
		r0 = const_cast<AutoMatrix&>( hrMat ).subMatrix( 0, 0, rank, cols ),
		*rp = &r0;	// to efficiently deal with linear in/dependent cases
	const Vector qTy = const_cast<AutoMatrix&>( yMat ).columnVector( cols ).subVector( 0, rank );
	bool permuteMatrix = false;
	Permutation *pp = NULL;

	if ( rank < cols ) {
		pp = new AutoPermutation( cols );
		size_t
			i = 0,
			j = rank;
		for ( size_t col = 0; col < cols; ++col ) {
			if ( isLinearlyIndependent( col ) ) {
				permuteMatrix |= col != i;
				pp->set( col, i++ );
			} else {
				permuteMatrix |= col != j;
				pp->set( col, j++ );
			}
		}
		if ( permuteMatrix ) {
			rp = new AutoMatrix( rank, rank );
			for ( size_t col = 0, r = 1; col < cols; ++col ) {
				const size_t targetCol = pp->get( col );
				if ( targetCol < rank ) {
					const Vector v = r0.columnVector( col ).subVector( 0, r );
					Vector w = rp->columnVector( targetCol ).subVector( 0, r );
					w.copy( v );
					++r;
				}
			}
		} else {
			// cutting trailing linearly dependent columns suffices
			r0 = r0.subMatrix( 0, 0, rank, rank );
		}
	}

	AutoVector x( cols );
	Vector
		px = x.subVector( 0, rank ),
		nx = x.subVector( rank, cols - rank );
	px.solveR( *rp, qTy );
	nx.fill( 0.0 );

	// in case of linearly dependent columns
	// re-permute x
	// and clean up extra resources
	if ( pp != NULL ) {
		if ( permuteMatrix ) {
			x.permute( *pp, false );		// backpermute x
			delete rp;	// using the virtual destructor -> AutoMatrix
		}
		delete pp;	// using the virtual destructor -> AutoPermutation
	}

	return x;
}

/** The algorithm exploits the fact that the R-Matrix with a column deleted still is almost triangular.
* Only the columns right of the deletion zone each have one potential nonzero element below the diagonal.
* These elements can be efficiently removed by a sequence of only 2-dimensional Householder transformations.
* The corresponding transformations applied to y then yield the RSS.
* To avoid accumulation of error, the modified matrix and vector are not used for further calculations, but discarded.
* The calculated RSS is meant to guide the selection of SNPs for backward steps based on
* {@link popColumn()} and {@link pushColumn( Vector& )} operations.
*/
double QRuncher::calculateSkipColumnRSS ( const size_t col ) const {
	const size_t
		rows = hrMat.countRows(),
		cols = hrMat.countColumns();
	assert( yMat.countRows() == rows );
	assert( yMat.countColumns() == 1 + cols );
	assert( col < cols );

	if ( isLinearlyIndependent( col ) ) {
		const size_t rank = getRank();
		assert( rank <= rows );
		assert( rank <= cols );
		// the "so far" unannihilated part of y (i.e. up to before col) starts at rank:
		// ranks[col] - 1 is the same as 0 < col ? ranks[ col - 1 ] : 0
		// due to the linear independence assured above.
		const size_t waterlineRank = ranks.at( col ) - 1;
		assert( waterlineRank <= rank );
		assert( waterlineRank <= col );

		// Make copies of the relevant remainders of hrMat (as of col+1) and yMat's last column
		// into a common memory-reduced space.
		const size_t
			rowsR = rank - waterlineRank,
			colsR = cols - col - 1;		// nonnegative due to col < cols
		AutoMatrix ryR( rowsR, colsR + 1 );
		Matrix rR = ryR.subMatrix( 0, 0, rowsR, colsR );
		for ( size_t c = 0; c < colsR; ++c ) {
			const size_t r = ranks.at( col + 1 + c ) - waterlineRank;
			Vector u = const_cast<AutoMatrix&>( hrMat ).columnVector( col + 1 + c );
			const Vector v = u.subVector( waterlineRank, r );
			Vector w = rR.columnVector( c ).subVector( 0, r );
			w.copy( v );
		}
		Vector yR = ryR.columnVector( colsR );
		{
			Vector u = const_cast<AutoMatrix&>( yMat ).columnVector( cols );
			const Vector v = u.subVector( waterlineRank, rowsR );
			yR.copy( v );
		}

		size_t waterlineR = 0;	// waterline in reduced matrix
		for ( size_t c = 0; c < colsR; ++c ) {
			assert( waterlineRank + waterlineR <= ranks.at( col + 1 + c ) );
			const size_t r = ranks.at( col + 1 + c ) - waterlineRank - waterlineR;
			if ( 0 < r ) {	// only columns with sufficiently ranked coefficients can annihilate parts of y
				Vector hh = rR.columnVector( c ).subVector( waterlineR, r );
				if ( ! hh.isNull() ) {
					const double tau = hh.householderize();
					Matrix transformable = ryR.subMatrix( waterlineR, c + 1, r, colsR - c );
					transformable.householderTransform( tau, hh );
					++waterlineR;	// dealt with one more row.
				}
			}
		}
		const Vector yResidue = yR.subVector( waterlineR, rowsR - waterlineR );
		const double
			rssAdditional = yResidue.sumSquares(),
			rssOld = calculateRSS();
		return rssAdditional + rssOld;
	} else {
		// if only a linearly dependent column is removed
		// then the residue is unchanged
		return calculateRSS();
	}
}
