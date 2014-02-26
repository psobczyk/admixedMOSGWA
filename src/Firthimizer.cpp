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

#include "Firthimizer.hpp"
#include <assert.h>
#include <math.h>
#include <iostream>

using namespace std;
using namespace linalg;
using namespace minimization;

Firthimizer::Firthimizer ( const Vector& yVec ) : status( INITIAL ), yVec( yVec.countDimensions() ), xMat( yVec.countDimensions(), 0 ), coefficientVec( 0 ), priorVec( yVec.countDimensions() ), swDiag( yVec.countDimensions() ), swxMat( yVec.countDimensions(), 0 ), tauVec( 0 ), logLikelihood( 0.0 ) {
	this->yVec.copy( yVec );
}

size_t Firthimizer::countDimensions () const {
	size_t cols = xMat.countColumns();
	return cols;
}

void Firthimizer::insertColumn ( const size_t col, const Vector& xVec, const double guess ) {
	// dimensions sanity checks
	const size_t rows = xMat.countRows();
	assert( rows == xVec.countDimensions() );
	size_t cols = xMat.countColumns();
	// TODO<BB>: What should happen in unintended case cols > rows?
	assert( cols < rows );
	assert( cols == coefficientVec.countDimensions() );

	// after insertion, internal data will need recalculation
	status = INITIAL;

	// Get space for insertion
	xMat.upSize( rows, ++cols );
	coefficientVec.upSize( cols );

	// shift forward all columns and regression coefficients after new column position
	for ( size_t i = cols - 1; i > col; --i ) {
		const Vector v = xMat.columnVector( i - 1 );
		Vector w = xMat.columnVector( i );
		w.copy( v );
		const double b = coefficientVec.get( i - 1 );
		coefficientVec.set( i, b );
	}

	// insert the new column.
	Vector v = xMat.columnVector( col );
	v.copy( xVec );
	coefficientVec.set( col, guess );

	// prepare new space in calculation matrix
	swxMat.upSize( rows, cols );

	// prepare space in QR-factorisation diagonal
	tauVec.upSize( cols );
}

void Firthimizer::replaceColumn ( const size_t col, const Vector& xVec, const double guess ) {
	// dimensions sanity checks
	const size_t rows = xMat.countRows();
	assert( rows == xVec.countDimensions() );
	size_t cols = xMat.countColumns();
	// TODO<BB>: What should happen in unintended case cols > rows?
	assert( cols <= rows );

	// after insertion, internal data will need recalculation
	status = INITIAL;

	// replace the column.
	Vector v = xMat.columnVector( col );
	v.copy( xVec );
	coefficientVec.set( col, guess );
}

void Firthimizer::removeColumn ( const size_t col ) {
	const size_t rows = xMat.countRows();
	size_t cols = xMat.countColumns();
	assert( cols == coefficientVec.countDimensions() );

	if ( col < cols ) {

		// after removal, internal data will need recalculation
		status = INITIAL;

		// shift backward all columns and regression coefficients after new column position
		for ( size_t i = col + 1; i < cols; ++i ) {
			const Vector v = xMat.columnVector( i );
			Vector w = xMat.columnVector( i - 1 );
			w.copy( v );
			const double b = coefficientVec.get( i );
			coefficientVec.set( i - 1, b );
		}

		// downsize all depending internal variables
		xMat.upSize( xMat.countRows(), --cols );
		coefficientVec.upSize( cols );
		swxMat.upSize( rows, cols );
		tauVec.upSize( cols );
	} else {
		// TODO<BB> implement exception
	}
}

void Firthimizer::setCoefficients ( const Vector& coefficients ) {
	const size_t
		rows = xMat.countRows(),
		cols = xMat.countColumns();

	// dimensions sanity check
	assert( cols == coefficientVec.countDimensions() );
	assert( cols == coefficients.countDimensions() );

	// after setting regression coefficients, all is calculated up to and including log likelihood
	status = FUNCTION;

	// copy correlation coefficients
	coefficientVec.copy( coefficients );

	// prepare to calculate priorVec distribution p, sqrt(W)X and log-likelihood
	priorVec.gemv( -1.0, xMat, false, coefficientVec, 0.0 );	// -1 due to exp(-Xt) in logistic regression in loop below

	// prepare to calculate swxMat = sqrt(W)X
	swxMat.fill( 0.0 );

	// prepare to calculate log-likelihood
	logLikelihood = 0.0;

	// calculation loop
	for ( size_t i = 0; i < rows; ++i ) {
		const double p = 1.0 / ( 1.0 + exp( priorVec.get( i ) ) );
		priorVec.set( i, p );	// only now priorVec[i] is the approximate priorVec value

		const double swii = sqrt( p * ( 1.0 - p ) );
		swDiag.set( i, swii );

		const Vector xi = xMat.rowVector( i );
		Vector swxi = swxMat.rowVector( i );
		swxi.axpy( swii, xi );

		const double yi = yVec.get( i );
		// possible optimisation would be: calculate log of products of e.g. groups of 10 terms, but danger of float under/overflow
		if ( 0.0 == yi ) {
			logLikelihood += log( 1.0 - p );
		} else if ( 1.0 == yi ) {
			logLikelihood += log( p );
		} else {
			// cout << "error: invalid y vector entry " << yi << endl;
			logLikelihood += yi * log( p ) + ( 1.0 - yi ) * log( 1.0 - p );
		}
	}

	// Write D=swxMat and note that the matrix P=D(D^TD)^{-1}D^T is an orthogonal projection.
	// Using QR decomposition D=QR, where Q=(Q_1,Q_0) meaningfully,
	// it is an easy exercise to see P=Q_1 Q_1^T, i.e. projects onto range(Q_1)=range(D).
	// Hence, perform QR-factorisation of swxMat into (tauVec,swxMat).
	swxMat.factorizeQR( tauVec );

	// multiply log L with log sqrt det I
	// = 1/2 log( det swx^T * swx ) = 1/2 log( det R^TQ^TQR ) = log( det R * signum Q )
	// = log( prod( diag( R ) ) * signum Q ) = log abs prod diag R
	// = log prod abs diag R = sum_i log abs R_ii
	const Vector diag = swxMat.diagonalVector();	// since swxMat is already QR-factorized, its diagonal is diag R
	for ( size_t i = 0; i < cols; ++i ) {
		logLikelihood += log( fabs( diag.get( i ) ) );
	}
}

void Firthimizer::calculateFunction ( const Vector& coefficients, double& negativeLogLikelihood ) {
	setCoefficients( coefficients );
	negativeLogLikelihood = -logLikelihood;
}

void Firthimizer::calculateDerivative ( const Vector& coefficients, Vector& negativeLogLikelihoodDerivative ) {
	const size_t
		rows = xMat.countRows(),
		cols = xMat.countColumns();

	setCoefficients( coefficients );
	AutoMatrix q( rows, cols );
	AutoVector ypf( rows );
	q.extractQ( tauVec, swxMat );
	// calculate derivative as u* = -f + y - priors
	for ( size_t i = 0; i < rows; ++i ) {
		const Vector qRow = q.rowVector( i );
		const double hi = qRow.sumSquares();
		const double firthCorrection = hi * ( 0.5 - priorVec.get( i ) );
		ypf.set( i, firthCorrection );
	}
	ypf.axpy( 1.0, yVec );
	ypf.axpy( -1.0, priorVec );
	negativeLogLikelihoodDerivative.gemv( -1.0, xMat, true, ypf, 0.0 );
}

void Firthimizer::calculateFunctionAndDerivative ( const Vector& coefficients, double& negativeLogLikelihood, Vector& negativeLogLikelihoodDerivative ) {
	Firthimizer::calculateDerivative( coefficients, negativeLogLikelihoodDerivative );
	negativeLogLikelihood = -logLikelihood; // this has been calculated on the way to here in setCoefficients
}

double Firthimizer::calculateLogLikelihood ( Minimizer& minimizer ) {
	if ( status < FUNCTION ) {
		// The assignment is moderately redundant, except if iteration fails and nan is returned.
		// The minus sign comes from the formulation of the regression as minimization.
		logLikelihood = -minimizer.minimize( *this, coefficientVec );
	}
	return logLikelihood;
}

const Vector& Firthimizer::calculateCoefficients ( Minimizer& minimizer ) {
	calculateLogLikelihood( minimizer );
	return coefficientVec;
}

double Firthimizer::scoreTest ( const linalg::Vector& zVec ) const {
	const size_t
		rows = xMat.countRows(),
		cols = xMat.countColumns();

	// dimensions sanity checks
	assert( rows == zVec.countDimensions() );

	AutoVector uVec( rows );
	uVec.copy( yVec );
	uVec.axpy( -1, priorVec );
	const double u = uVec.innerProduct( zVec );

	// sqrt(W) z
	AutoVector swzVec( rows );
	for ( size_t i = 0; i < rows; ++i ) {
		swzVec.set( i, swDiag.get( i ) * zVec.get( i ) );
	}

	// Now mixing notations from Firth and He+Lin
	// I_gamma,gamma = sum z*W_ii*z = norm( sqrt(W)z )^2.
	// I_gamma,eta * I_eta,eta^{-1} * I_gamma,eta^T = [Bernhard's Denkbuch anno 2010, p. 19.Jan. f] = norm(P sqrt(W)z)^2,
	// where P = QQ^T is the Euclidean-orthogonal projection onto range(sqrt(W)X).
	// Hence, V = norm(sqrt(W)z)^2 - norm(P sqrt(W)z)^2 = norm( (1-P) sqrt(W)z )^2.
	AutoVector pswzVec( swzVec );
	pswzVec.multQ( tauVec, swxMat, true );
	{
		Vector pswzTail = pswzVec.subVector( cols, rows - cols );
		pswzTail.fill( 0.0 );	// project onto range(R), assuming invertibility of R
	}
	pswzVec.multQ( tauVec, swxMat, false ); // now this is QQ^T sqrt(W)z
	pswzVec.axpy( -1.0, swzVec );
	const double v = pswzVec.innerProduct( pswzVec );
	if ( 0.0 == v ) {
		return nan( "linearly dependent" );
	} else {
		const double scoreTest = u / sqrt( v );
		return fabs( scoreTest );
	}
}
