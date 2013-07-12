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

#ifndef QRUNCHER_HPP
#define QRUNCHER_HPP

#include "linalg/AutoMatrix.hpp"
#include "linalg/AutoVector.hpp"
#include <vector>

/** Incremental QR solver for least squares regression.
* @author Bernhard Bodenstorfer
*/
class QRuncher {

	protected:

	/** The sample result vector and its householder transform images defined inductively.
	* The leftmost column vector (column 0) is the result vector determined by the sample,
	* customarily called y.
	* If 0 < col, then the vector in column col
	* is the image of the vector in column col-1 transformed by the householder transformation
	* in hrMat and tauVec in column resp. dimension col-1.
	* Hence, yMat always has one column more than hrMat.
	*/
	linalg::AutoMatrix yMat;

	/** Holds the QR decomposition of a number of columns in the order they were added.
	* It is actually a triangular upper right matrix (the "R" from QR)
	* mingled with the Householder vectors below the diagonal.
	* @see <a href="http://www.gnu.org/software/gsl/manual/html_node/QR-Decomposition.html">GSL QR decomposition</a>
	*/
	linalg::AutoMatrix hrMat;

	/** Holds the Householder vector tau values for the respective columns.
	* @see {@link QRuncher::factorizeQR}
	*/
	linalg::AutoVector tauVec;

	/** Tracks the rank for the columns as they have been added to {@link QRuncher::hrMat}.
	* Element <code>m</code> specifies the dimension of the linear hull of
	* the <code>m+1</code> pushed vectors from <code>0</code> up to <code>m</code>.
	*/
	std::vector<size_t> ranks;

	/** Return whether the given column was linearly independent from its predecessors.
	*/
	bool isLinearlyIndependent ( const size_t col ) const;

	/** Get the rank of the current regression matrix.
	* That is the number of linearly independent vectors pushed.
	*/
	size_t getRank () const;

	/** Retrieve the householder vector for the given column.
	* It is assumed that the pushed vector for the column was not linearly dependent
	* from its predecessors.
	* Then the sought householder vector is the column in {@link QRuncher::hrMat},
	* all elements downward from the rank-element to the number of rows of the matrix.
	*/
	linalg::Vector getHouseholderVector ( const size_t col );

	public:

	/** Set the results vector.
	* This also determines the number of rows.
	*/
	QRuncher ( const linalg::Vector& yVec );

	/** Extend the regression matrix at its right end by one column.
	* Since the residual square sums algorithm relies on all columns to be linearly independent,
	* not more columns may be added than what is the dimension of the space as given by
	* the number of dimensions of the vector <code>yVec</code>
	* with which the <code>QRuncher</code> has been constructed.
	* @param xVec must have the same number of dimensions as the <code>yVec</code>
	* with which the <code>QRuncher</code> has been constructed.
	* @returns whether the added column was linearly independent from the previous columns.
	* @see {@link QRuncher::QRuncher( const Vector& )}
	*/
	bool pushColumn ( const linalg::Vector& xVec );

	/** Diminish the regression matrix at its right end by one column.
	* @returns whether there was still a column left to chop away.
	*/
	bool popColumn ();

	/** Calculate the residual sum of squares. */
	double calculateRSS () const;

	/** Calculate the regression coefficients.
	* In case of linearly dependent column vectors,
	* one particular solution is found,
	* where the coefficients for linearly dependent columns are set 0.
	* @returns a coefficient vector for the currently pushed column vectors
	*/
	linalg::AutoVector calculateCoefficients () const;

	/** Quickly calculate the residual sum of squares (RSS)
	* if the given column in the push hierarchy (0 being the oldest push) had been skipped.
	* @param col must be in the range
	* from 0 to the number of calls to {@link pushColumn} minus calls to {@link popColumn} minus 1.
	* @returns the RSS of the current model minus the specified pushed column,
	* enumerated by age of the push.
	*/
	double calculateSkipColumnRSS ( const size_t col ) const;
};

#endif	/* QRUNCHER_HPP */
