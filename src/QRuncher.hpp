#ifndef QRUNCHER_HPP
#define QRUNCHER_HPP

#include "linalg/AutoMatrix.hpp"
#include "linalg/AutoVector.hpp"

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

	/** Retrieve a householder vector from below the diagonal of {@link QRuncher::hrMat}. */
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
	* @returns a measure of linear independence,
	* (close to) zero if the added column was (almost) a linear combination of the previous columns.
	* @see {@link QRuncher::QRuncher( const Vector& )}
	*/
	double pushColumn ( const linalg::Vector& xVec );

	/** Diminish the regression matrix at its right end by one column.
	* @returns whether there was still a column left to chop away.
	*/
	bool popColumn ();

	/** Calculate the residual sum of squares. */
	double calculateRSS ();

	/** Calculate the regression coefficients.
	* Aborts execution in case of error,
	* i.e. linear dependence of regression matrix columns,
	* in particular, if there are more columns than rows,
	* i.e. more regression-variables than equations (individuals).
	* @returns a coefficient vector for the currently pushed variables
	*/
	linalg::AutoVector calculateCoefficients ();

	/** Quickly calculate the residual sum of squares (RSS)
	* if the given column in the push hierarchy (0 being the oldest push) had been skipped.
	* @param col must be in the range
	* from 0 to the number of calls to {@link pushColumn} minus calls to {@link popColumn} minus 1.
	* @returns the RSS of the current model minus the specified pushed column,
	* enumerated by age of the push.
	*/
	double calculateSkipColumnRSS ( const size_t col );
};

#endif	/* QRUNCHER_HPP */
