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

#ifndef FIRTHIMIZER_HPP
#define FIRTHIMIZER_HPP

#include "linalg/AutoMatrix.hpp"
#include "linalg/AutoVector.hpp"
#include "minimization/Minimizer.hpp"

/** Incremental-interfaced regression solver for logistic Firth regression.
* The incremental interface is purposefully similar to that of {@link QRuncher}.
* The idea is to be able to re-use information calculated earlier
* when considering a subset of the regression variables.
* @author Bernhard Bodenstorfer
* @see QRuncher
*/
class Firthimizer : private minimization::Minimizable {

	protected:

	/** Keeps track of the calculation status of internal state. */
	enum {
		/** Initial state, not even {@link coefficientVec} is up-to-date. */
		INITIAL = 0,
		/** Intermediate results up to and including {@link logLikelihood} are up-to-date. */
		FUNCTION = 1
	} status;

	/** The sample result vector. */
	linalg::AutoVector yVec;

	/** The logistic regression matrix for the variables added so far. */
	linalg::AutoMatrix xMat;

	/** The (approximated) regression coefficients. */
	linalg::AutoVector coefficientVec;

	/** The probabilities corresponding to the regression coefficients. */
	linalg::AutoVector priorVec;

	/** The diagonal of matrix sqrt(W) with W = diag( prior_i * ( 1 - prior_i ) ) as described in Firth's algorithm. */
	linalg::AutoVector swDiag;

	/** The matrix sqrt(W)*X as described in Firth's algorithm.
	* Mostly this is stored QR-factorised.
	* @see #tauVec
	*/
	linalg::AutoMatrix swxMat;

	/** The diagonal vector from the QR-factorisation of
	* the matrix sqrt(W)*X as described in Firth's algorithm.
	* @see #swxMat
	*/
	linalg::AutoVector tauVec;

	/** The logarithm of the likelihood density of the model corresponding to the regression coefficients. */
	double logLikelihood;

	/** Provide regression coefficients
	* and calculate those internal variables which depend on them
	* up to and including {@link logLikelihood}.
	* No minimization loop is run to improve coefficients.
	*/
	void setCoefficients ( const linalg::Vector& coefficients );

	/** Retrieve the current number of columns in the regression matrix.
	* @overrides Minimizable::countDimensions
	*/
	virtual size_t countDimensions () const;

	/** Calculate the negative log likelihood of the Firth regression, given the vector of regression parameters.
	* Negative value is needed to formulate the regression as minimization rather than maximization.
	* @overrides Minimizable::calculateFunction
	*/
	virtual void calculateFunction ( const linalg::Vector& coefficients, double& negativeLogLikelihood );

	/** Calculate the gradient of the negative log likelihood of the Firth regression, given the vector of regression parameters.
	* Negative value is needed to formulate the regression as minimization rather than maximization.
	* @overrides Minimizable::calculateDerivative
	*/
	virtual void calculateDerivative ( const linalg::Vector& coefficients, linalg::Vector& negativeLogLikelihoodDerivative );

	/** Calculate the negative log likelihood of the Firth regression, and its derivative,
	* given the vector of regression parameters.
	* Negative value is needed to formulate the regression as minimization rather than maximization.
	* @overrides Minimizable::calculateFunctionAndDerivative
	*/
	virtual void calculateFunctionAndDerivative (
		const linalg::Vector& coefficients,
		double& negativeLogLikelihood,
		linalg::Vector& negativeLogLikelihoodDerivative
	);

	public:

	/** Set the result vector.
	* This also determines the number of rows.
	*/
	Firthimizer ( const linalg::Vector& yVec );

	/** Extend the regression matrix by inserting one column at the specified position.
	* Not more columns may be added than what is the dimension of the space as given by
	* the number of dimensions of the vector <code>yVec</code>
	* with which the <code>Firthimizer</code> has been constructed.
	* @param col where the new column should be inserted. Insertion after the currently last column is fastest.
	* @param xVec must have the same number of dimensions as the <code>yVec</code>
	* with which the <code>Firthimizer</code> has been constructed.
	* @param guess specifies a starting point for the regression coefficient in the minimization algorithm.
	* @see #Firthimizer( const Vector& )
	* @see #countDimensions()
	*/
	void insertColumn ( const size_t col, const linalg::Vector& xVec, const double guess = 0.0 );

	/** Exchange the specified column in the regression matrix.
	* @param col which column should be replaced.
	* @param xVec must have the same number of dimensions as the <code>yVec</code>
	* with which the <code>Firthimizer</code> has been constructed.
	* @param guess specifies a starting point for the regression coefficient in the minimization algorithm.
	* @see #Firthimizer( const Vector& )
	* @see #countDimensions()
	*/
	void replaceColumn ( const size_t col, const linalg::Vector& xVec, const double guess = 0.0 );

	/** Diminish the regression matrix by removing one column at the specified position.
	* @param col where the new column should be removed. Removal of the currently last column is fastest.
	* @see #countDimensions()
	*/
	void removeColumn ( const size_t col );

	/** Calculate the log(likelihood).
	* @param minimizer provides minimization for the regression.
	*/
	double calculateLogLikelihood ( minimization::Minimizer& minimizer );

	/** Calculate the logistic regression coefficients.
	* @param minimizer provides minimization for the regression.
	* @param coefficients vector (output)
	* to store the optimal coefficients
	* for the currently pushed variables.
	* @returns a reference to the internally calculated and stored coefficient vector.
	* The referenced object is not guaranteed to persist after calling any further methods on <code>this</code> Firthimizer.
	*/
	const linalg::Vector& calculateCoefficients ( minimization::Minimizer& minimizer );

	/** Calculate the score test value to quickly assess eligibility of a SNP.
	* @returns <code>NaN</code> if <code>zVec</code>
	* is in the range of the internal regression matrix sqrt(W(coefficientVec))X
	* for the currently set regression coefficients {@link coefficientVec}.
	*/
	double scoreTest ( const linalg::Vector& zVec ) const;

};

#endif	/* FIRTHIMIZER_HPP */
