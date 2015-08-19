/********************************************************************************
 *	This file is part of the MOSGWA program code.				*
 *	Copyright ©2011–2015, Erich Dolejsi, Bernhard Bodenstorfer.		*
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

#include "Helpfull.hpp"
#include "linalg/AutoVector.hpp"
#include "linalg/AutoMatrix.hpp"
#include "linalg/AutoPermutation.hpp"
#include <iostream>
#include <cassert>
#include <limits>

using namespace std;
using namespace linalg;

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// helpful functions
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

/** computes 3 Vectors for the Regression accounting for Population 1,2,3
* @param vecPos: which of the 3 Vectors
* @param individum: the Indivium < Sample number
*/
int populationVector ( const int vecPos, const int individum ) {
	switch ( vecPos ) {
		case 1:
			if ( individum < 45 ) return 0; //CHB
			if ( individum < 90 ) return 0; //JPT
			if ( individum < 150 ) return -1;	//CEU
			if ( individum < 210 ) return 1;	//YRI
			break;
		case 2:
			if ( individum < 45 ) return 1; //CHB
			if ( individum < 90 ) return 0; //JPT
			if ( individum < 150 ) return -1;	//CEU
			if ( individum < 210 ) return 0;	//YRI
			break;
		case 3:
			if ( individum < 45 ) return 0; //CHB
			if ( individum < 90 ) return 1; //JPT
			if ( individum < 150 ) return -1;	//CEU
			if ( individum < 210 ) return 0;	//YRI
			break;
		default: return 0;	// error
	}
	return 0;	// error
}

/** convert string to interger */
int str2int ( const string &str ) {
	stringstream ss(str);
	int n;
	ss >> n;	// TODO<BB>: here and in similar spots: this is not fault-tolerant.
	return n;
}


/** convert interger to string */
string int2str ( const int n ) {
	stringstream ss;
	ss << n;
	return ss.str();
}


/** convert integer to string */
string double2str ( const double x ) {
	stringstream ss;
	ss << x;
	return ss.str();
}


/** convert interger to a formated string (of a given width with filling in the character fillwith) */
string int2strPadWith ( const int n, const int width, const char fillwith = ' ' ) {
	stringstream ss;
	ss << setfill(fillwith) << setw (width) << n;
	return ss.str();
}

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// mathematical functions
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



int factorial ( int i ) {	// TODO<BB> is integer calculation enough i.e. i less than about 15?
	if ( 0 > i ) throw;	// TODO<BB> handle exceptions
	int prod;
	for ( prod = 1; 1 < i; --i ) {
		prod *= i;	// TODO<BB> make faster using precalculated table
	}
	return prod;
}

double log_factorial ( int i ) {
	if ( 0 > i ) throw;	// TODO<BB> handle exceptions
	double sum;
	for ( sum = 0.0; 1 < i; --i ) {
		sum += log( i );	// TODO<BB> make faster using precalculated table
	}
	return sum;
}

void calculate_p ( const Matrix& X, const Vector& beta, Vector& p ) {
	const size_t rows = p.countDimensions();
	p.gemv( -1.0, X, false, beta, 0.0 );
	for ( size_t i = 0; i < rows; ++i ) {
		p.set( i, 1.0 / ( 1.0  + exp( p.get( i ) ) ) );
	}
}

void calculate_loglik ( const Vector& p, const Vector& y, double* loglik ) {
	const size_t rows = p.countDimensions();
	assert( y.countDimensions() == rows );
	*loglik = 0.0;
	for ( size_t i = 0; i < rows; ++i ) {
		const double
			yi = y.get( i ),
			pi = p.get( i );
		*loglik += yi * log( pi ) + ( 1.0 - yi ) * log( 1.0 - pi );
	}
}

void calculate_XW2 ( const Matrix& X, const Vector& p, Matrix& XW2 ) {
	const size_t rows = X.countRows();
	assert( p.countDimensions() == rows );
	assert( XW2.countColumns() == rows );
	// XW2 = X^T W^(1/2) =	crossprod(x, diag(pi * (1 - pi))^0.5)
	for ( size_t i = 0; i < rows; ++i ) {
		Vector xw2col = XW2.columnVector( i );
		const double pi = p.get( i );
		xw2col.copy( const_cast<Matrix&>( X ).rowVector( i ) );
		xw2col.scale( sqrt( pi * ( 1.0 - pi ) ) );
	}
}

void calculate_Fisher_LU ( const Matrix& XW2, Matrix& Fisher, Permutation& perm_k ) {
	// Fisher matrix: X^T * W * X
	Fisher.gemm( 1.0, XW2, false, XW2, true, 0.0 );

	// to compute determinant or inverse, we need an LU-decompositon
	Fisher.factorizeLUP( perm_k );
}

/** computes logistic firth-regression for design matrix x and target 
 plain C code: logistf package for R from Georg Heinze
 All Parameters are needed to be set!!! */
bool logistffit (
	// Input
	const Matrix&		X,
	const Vector&		y,
	// In+Output
	Vector&			beta,		// coefficients (initial and then result)
	// Output
	double*			out_loglik,	// log-likelihood of the model
	// Control Parameters
	const int		maxit,
	const int		maxhs,
	const int		maxstep,
	const double	lconv,
	const double	gconv,
	const double	xconv
) {
	const size_t
		rows = X.countRows(),	// n … number of individuals
		cols = X.countColumns();	// k … 1 + number of variables
	assert( y.countDimensions() == rows );
	assert( beta.countDimensions() == cols );

	AutoMatrix
		XW2( cols, rows ),
		Fisher( cols, cols ),
		covs( cols, cols );
	AutoVector
		H( rows ),	// diag of H matrix vector
		p( rows ),	// probabilities for target
		helpVn( rows ),
		helpVk( cols ),
		UStar( cols ),
		delta( cols );
	AutoPermutation perm_k( cols );	// needed for LU_decomp
	double loglik;

	calculate_p( X, beta, p );
	calculate_loglik( p, y, &loglik );
	calculate_XW2( X, p, XW2 );
	calculate_Fisher_LU( XW2, Fisher, perm_k );

	loglik += 0.5 * Fisher.lnAbsDetLU();

	//Loop;
	double
		lchange = 5.0,
		loglik_old = numeric_limits<double>::quiet_NaN(),
		mx;
	int iter = 0;
	bool stop = false;

	while ( !stop ) {
		loglik_old = loglik;
		calculate_XW2( X, p, XW2 );
		calculate_Fisher_LU( XW2, Fisher, perm_k );

		// check if matrix is invertible, that is, all diagonal elements are non-zero
		for ( size_t i = 0; i < cols; ++i ) {
			if ( fabs( Fisher.get( i, i ) ) < 1e-15 ) {
				cerr << "Fisher Matrix is not invertible" << endl;
				return false;
			}
		}
		// covs = (X^T W X)^-1
		covs.invertLUP( Fisher, perm_k );

		// Compute the diagonal-Elements of (XW2)^T * covs * XW2
		for ( size_t i = 0; i < rows; ++i ) {
			const Vector xw2i = XW2.columnVector( i );
			helpVk.gemv( 1.0, covs, false, xw2i, 0.0 );
			const double hi = helpVk.innerProduct( xw2i );
			H.set( i, hi );
		}

		for ( size_t i = 0; i < rows; ++i ) {
			const double
				yi = y.get( i ),
				pi = p.get( i ),
				hi = H.get( i );
			helpVn.set( i, yi - pi + hi * ( 0.5 - pi ) );
		}
		UStar.gemv( 1.0, X, true, helpVn, 0.0 );
		delta.gemv( 1.0, covs, false, UStar, 0.0 );

		// set mx = max(abs(delta))/maxstep
		mx = 0.0;
		for ( size_t i = 0; i < cols; ++i ) {
			const double mmx = fabs( delta.get( i ) );
			if ( mx < mmx ) {
				mx = mmx;
			}
		}
		mx /= maxstep;

		if ( 1.0 < mx ) {
			delta.scale( 1.0/mx );
		}

		if ( 0 < maxit ) {
			++iter;

			// beta= beta + delta
			beta.axpy( 1.0, delta );

			// Half-Steps
			for ( int halfs = 1; halfs <= maxhs; ++halfs ) {
				calculate_p( X, beta, p );
				calculate_loglik( p, y, &loglik );
				calculate_XW2( X, p, XW2 );
				calculate_Fisher_LU( XW2, Fisher, perm_k );

				loglik += 0.5 * Fisher.lnAbsDetLU();

				lchange = loglik - loglik_old;

				if ( loglik > loglik_old ) {
					break; //	end Half-Step-loop
				}

				beta.axpy( -pow( 2.0, -halfs ), delta );
				// TODO<BB>: replace pow( 2, . ) by quick table lookup
				// and better re-start from original betas than back-correcting beta+t*delta
			}
		}

		//////////////////////////////////////////////////////////////////////////////
		// Stop-Condition:
		if ( iter >= maxit ) {
			stop = true;
		} else if ( lchange < lconv ) {
			// or convergence criteria are met
			stop = true;
			for ( size_t i = 0; stop && i < cols; ++i ) {
				if (
					fabs( delta.get( i ) ) >= xconv
					||
					fabs( UStar.get( i ) ) >= gconv
				) {
					stop = false;
				}
			}
		}
	}

	*out_loglik = loglik;
	return true;
}
