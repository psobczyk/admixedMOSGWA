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

#include "linalg/Vector.hpp"
#include "linalg/Matrix.hpp"

////++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// helpful functions
////++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

const int bitsSize = 2000000;      // AG: maximal number of SNPs
typedef std::bitset<bitsSize> TBitset;  // AG: for an initial population of GA


/** computes */
int populationVector ( const int VecPos, const int Individum );

/** convert integer to string */
std::string int2str ( const int n );

/** convert integer to a formated string (of a given width with filling in the character fillwith) */
std::string int2strPadWith ( const int n, const int width, const char fillwith );

////++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// mathematical functions
////++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

/** computes \f$\ln(i!)\f$ */
double log_factorial( int i );

/** computes logistic firth-regression for design matrix x and target y.
 plain C code: logistf package for R from Georg Heinze (with slight modifications: error control in case of singular matrices)
 All Parameters are needed to be set!!! */
bool logistffit(
	// Input
	const linalg::Matrix&	X,		// design matrix ( n * k )
	const linalg::Vector&	y,		// target vector ( n )
	// In+Output
	linalg::Vector&		beta,		// regression coefficients (vector k)
	// Output
	double*			out_loglik,	// log-likelihood of the model
	// Control Parameters
	const int		maxit,
	const int 		maxhs,
	const int 		maxstep,
	const double		lconv,
	const double		gconv,
	const double		xconv
);

#endif
