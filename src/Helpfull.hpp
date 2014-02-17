/********************************************************************************
 *	This file is part of the MOSGWA program code.				*
 *	Copyright ©2011–2013, Erich Dolejsi, Bernhard Bodenstorfer.		*
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

#ifndef HELPFULL_HPP
#define HELPFULL_HPP

#include <iomanip>
#include <sstream>
#include <math.h>

#include <string.h>
#include <math.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>

#include <bitset> 

////++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// helpful functions
////++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

const int bitsSize = 2000000;      // AG: maximal number of SNPs
typedef std::bitset<bitsSize> TBitset;  // AG: for an initial population of GA


/** computes */
int populationVector ( const int VecPos, const int Individum );

/** convert string to integer */
int str2int ( const std::string &str );

/** convert integer to string */
std::string int2str ( const int n );

/** convert double to string */
std::string double2str ( const double x );

/** convert integer to a formated string (of a given width with filling in the character fillwith) */
std::string int2strPadWith ( const int n, const int width, const char fillwith );

////++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// mathematical functions
////++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

/** computes i! */
int factorial( int i );

/** computes \f$\ln(i!)\f$ */
double log_factorial( int i );

/** computes logistic firth-regression for design matrix x and target y.
 plain C code: logistf package for R from Georg Heinze (with slight modifications: error control in case of singular matrices)
 All Parameters are needed to be set!!! */
bool logistffit(
	// Input
	const int		*in_k,		// # of rows of x (# of varibles) BB: I'd rather call that "columns"
	const int		*in_n,		// # of columns of x (# of samples) BB: I'd rather call that "rows"
	const gsl_matrix *	x,		// design matrix ( n * k )
	const int*		y,		// target vector ( n )
	// Output
	double*			beta,		// regression coefficients (vector k)
	double*			var,		// var matrix ( k * k )
	double*			pi,		// probabilities for target ( n )
	double*			H,		// diag of H matrix vector (n)
	double*			out_loglik,	// log-likelihood of the model
	int*			iter,		// # of main iterations;
	int*			evals,		// # of evalations of evaluatios
	double*			lchange,	// the change in the log-likelihood between steps
	double*			ret_max_U_star,
	double*			ret_max_delta,
	// Optional Input
	const double*		weight,		// the weights
	const double*		offset,		// offset
	const int*		firth,		// use firth-regression
	const int*		col_fit,	// a "boolean" vector indication which colums to use
	const double*		init,		// initial values for beta
	const int*		eval_llh,	// only evaluate likelihood
	// Control Parameters
	const int*		maxit,
	const int*		maxhs,
	const int*		maxstep,
	const double*	lconv,
	const double*	gconv,
	const double*	xconv
);

#endif
