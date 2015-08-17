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
#include "TestSuite.hpp"
#include "linalg/AutoPermutation.hpp"
#include <cmath>
#include <cassert>

using namespace std;
using namespace linalg;
using namespace unitpp;

namespace test {

	/** Reference implementation using Gram-Schmidt algorithm.
	* This looks somewhat like a {@link QRuncher},
	* but uses a completely different * (numerically less accurate) algorithm,
	* which is, however, much simpler.
	* Also it does not work for linear dependence.
	* This is used to get reasonable comparison values for testing {@link QRuncher}.
	*/
	class TestQRuncher {
		AutoMatrix
			/** Stores Gram-Schmidt orthogonalised pushed columns. */
			xMat,
			/** Stores y-Vector remainders corresponding to the evolution of xMat. */
			yMat;

		/** Perform one Gram-Schmidt step on <code>w</code>
		* That is, change <code>w</code>
		* to its remainder orthogonal to <code>v</code>
		* and return the square sum of that.
		* This is used to test the {@link QRuncher::pushColumn} return values.
		*/
		static void stepGramSchmidt ( const Vector& v, Vector& w );

		public:
		TestQRuncher ( const Vector& yVec );
		double pushColumn ( const Vector& xVec );
		void popColumn ();
		double calculateRSS () const;
	};

	TestQRuncher::TestQRuncher ( const Vector& yVec )
	:
		xMat( yVec.countDimensions(), 0 ),
		yMat( yVec.countDimensions(), 1 )
	{
		Vector yCol = yMat.columnVector( 0 );
		yCol.copy( yVec );
	}

	void TestQRuncher::stepGramSchmidt ( const Vector& v, Vector& w ) {
		w.axpy( - w.innerProduct( v ) / v.sumSquares(), v );
	}

	double TestQRuncher::pushColumn ( const Vector& xVec ) {
		const size_t
			rows = xMat.countRows(),
			cols = xMat.countColumns();
		xMat.upSize( rows, cols + 1 );
		yMat.upSize( rows, cols + 2 );
		Vector
			xCol = xMat.columnVector( cols ),
			yCol = yMat.columnVector( cols + 1 );
		xCol.copy( xVec );
		yCol.copy( yMat.columnVector( cols ) );
		for ( size_t col = 0; col < cols; ++col ) {
			const Vector oldCol = xMat.columnVector( col );
			stepGramSchmidt( oldCol, xCol );
		}
		stepGramSchmidt( xCol, yCol );
		return sqrt( xCol.sumSquares() );
	}

	void TestQRuncher::popColumn () {
		const size_t
			rows = xMat.countRows(),
			cols = xMat.countColumns();
		assert( 0 < cols );
		xMat.upSize( rows, cols - 1 );
		yMat.upSize( rows, cols );
	}

	double TestQRuncher::calculateRSS () const {
		const size_t cols = xMat.countColumns();
		const Vector yCol = const_cast<AutoMatrix&>( yMat ).columnVector( cols );
		return yCol.sumSquares();
	}

	/** Tests the class {@link QRuncher}. */
	class QRuncherTest : public TestSuite {

		/** Test with a sequence of actions. */
		void testSequence ();

		/** Test calculation of regression coefficients. */
		void testCoefficients ();

		/** Test the (default) copy constructor. */
		void testCopyConstructor ();

		/** Test the (default) assignment operator. */
		void testAssignment ();

		/** Tests {@link QRuncher::calculateSkipColumnRSS}. */
		void testSkipColumn ();

		/** Test with linearly dependent regression matrix columns. */
		void testLinearlyDependent ();

		/** Tests {@link QRuncher::calculateSkipColumnRSS}
		* in the case of skipping the last added linearly independent column.
		*/
		void testSkipLastLinearlyIndependentColumn ();

		/** Tests with linearly dependent regression matrix columns,
		* but in a way which reduces rounding errors.
		*/
		void testLimitedRounding ();

		/** Tests the mathematical equality of (X^T X)^-1 right lower corner
		* equal square of inverse of last pushColumn result,
		* which is the absolute value of the R matrix right lower corner.
		*/
		void testRightLowerCorner ();

		public:

		QRuncherTest () : TestSuite( "QRuncherTest" ) {
			addTestMethod( "QRuncherTest::testSequence", this, &QRuncherTest::testSequence );
			addTestMethod( "QRuncherTest::testCoefficients", this, &QRuncherTest::testCoefficients );
			addTestMethod( "QRuncherTest::testCopyConstructor", this, &QRuncherTest::testCopyConstructor );
			addTestMethod( "QRuncherTest::testAssignment", this, &QRuncherTest::testAssignment );
			addTestMethod( "QRuncherTest::testSkipColumn", this, &QRuncherTest::testSkipColumn );
			addTestMethod( "QRuncherTest::testSkipLastLinearlyIndependentColumn", this, &QRuncherTest::testSkipLastLinearlyIndependentColumn );
			addTestMethod( "QRuncherTest::testLinearlyDependent", this, &QRuncherTest::testLinearlyDependent );
			addTestMethod( "QRuncherTest::testLimitedRounding", this, &QRuncherTest::testLimitedRounding );
			addTestMethod( "QRuncherTest::testRightLowerCorner", this, &QRuncherTest::testRightLowerCorner );
		}
	} * qrUncherTest = new QRuncherTest();	// automatically freed by unit++

	void QRuncherTest::testSequence () {
		// col 1 norm 9, col 2 norm 7, col 3 = col 2 - col 1
		const double data[] = {
			1., 2.,		1.,
			4., 3.,		-1.,
			8., 6.,		-2.
		};

		AutoMatrix dataMat( 3, 3 );
		dataMat.fill( data );
		const Vector
			x0 = dataMat.columnVector( 0 ),
			x1 = dataMat.columnVector( 1 ),
			x2 = dataMat.columnVector( 2 ),
			y = dataMat.columnVector( 0 );

		TestQRuncher tqr( y );
		QRuncher qr( y );

		// 0 columns
		assert_eq( "Zero variables tqr", 81.0, tqr.calculateRSS() );
		assert_eq( "Zero variables qr", 81.0, qr.calculateRSS() );
		AutoVector coefficients0a( 0 );
		qr.calculateCoefficients( coefficients0a );

		// 1 column
		assert_eq( "No linear dep 1 tqr", 9.0, tqr.pushColumn( x0 ) );
		assert_eq( "No linear dep 1 qr", 9.0, qr.pushColumn( x0 ) );
		assert_close( "Exact match tqr", 0.0, qr.calculateRSS() );
		assert_close( "Exact match qr", 0.0, qr.calculateRSS() );
		AutoVector coefficients1a( 1 );
		qr.calculateCoefficients( coefficients1a );
		assert_close( "Exact match coefficients", 1.0, coefficients1a.get( 0 ) );

		// 2 columns
		assert_close( "No linear dep 2", tqr.pushColumn( x1 ), qr.pushColumn( x1 ) );
		assert_close( "Still exact match", 0.0, qr.calculateRSS() );
		AutoVector coefficients2a( 2 );
		qr.calculateCoefficients( coefficients2a );
		assert_close( "Exact match coefficient 0", 1.0, coefficients2a.get( 0 ) );
		assert_close( "Exact match coefficient 1", 0.0, coefficients2a.get( 1 ) );

		// 3 columns (linear dependence, though probably obscured by rounding errors)
		assert_close( "Full linear dep", 0.0, qr.pushColumn( x2 ) );
		assert_close( "And still exact match", 0.0, qr.calculateRSS() );
		AutoVector coefficients3( 3 );
		qr.calculateCoefficients( coefficients3 );
		assert_close( "Exact match coefficient 0", 1.0, coefficients3.get( 0 ) );
		assert_close( "Exact match coefficient 1", 0.0, coefficients3.get( 1 ) );
		assert_close( "Exact match coefficient 2", 0.0, coefficients3.get( 2 ) );

		// back to 2 columns
		qr.popColumn();
		assert_close( "Back to still exact match", 0.0, qr.calculateRSS() );
		AutoVector coefficients2b( 2 );
		qr.calculateCoefficients( coefficients2b );
		assert_eq( "Back to two-dim regression coefficients", coefficients2a, coefficients2b );

		// back to 1 column
		qr.popColumn();
		assert_close( "Back to exact match", 0.0, qr.calculateRSS() );
		AutoVector coefficients1b( 1 );
		qr.calculateCoefficients( coefficients1b );
		assert_eq( "Back to one-dim regression coefficients", coefficients1a, coefficients1b );

		// back to 0 columns
		qr.popColumn();
		assert_eq( "Back to zero variables", 81.0, qr.calculateRSS() );
		AutoVector coefficients0b( 0 );
		qr.calculateCoefficients( coefficients0b );
		assert_eq( "Back to zero-dim regression coefficients", coefficients0a, coefficients0b );
	}

	void QRuncherTest::testCoefficients () {
		const double data[] = {
			3., 4.,		8.,
			4., 3.,		8.
		};

		AutoMatrix dataMat( 2, 3 );
		dataMat.fill( data );
		const Vector
			x0 = dataMat.columnVector( 0 ),
			x1 = dataMat.columnVector( 1 ),
			y = dataMat.columnVector( 2 );

		TestQRuncher tqr( y );
		QRuncher qr( y );
		AutoVector coefficients( 0 );

		// no variables yet; must run for dim 0, too
		qr.calculateCoefficients( coefficients );

		assert_eq( "No linear depe 0 qr", 5.0, qr.pushColumn( x0 ) );
		assert_eq( "No linear depe 0 tqr", 5.0, tqr.pushColumn( x0 ) );
		// Minimise || (8,8) - t * (3,4) ||^2 = (8-3t)^2+(8-4t)^2.
		// Derivative proportional 3(8-3t) + 4(8-4t) = 24-9t + 32-16t = 56-25t
		coefficients.upSize( 1 );
		qr.calculateCoefficients( coefficients );
		assert_close( "Approximate match coefficient 0", 56./25., coefficients.get( 0 ) );

		assert_close( "No linear dep 1", tqr.pushColumn( x1 ), qr.pushColumn( x1 ) );
		assert_close( "Exact match", 0, qr.calculateRSS() );
		coefficients.upSize( 2 );
		qr.calculateCoefficients( coefficients );
		assert_close( "Exact match coefficients 0", 8./7., coefficients.get( 0 ) );
		assert_close( "Exact match coefficients 1", 8./7., coefficients.get( 1 ) );

		qr.popColumn();
		coefficients.upSize( 1 );
		qr.calculateCoefficients( coefficients );
		assert_close( "Back to approximate match coefficient 0", 56./25., coefficients.get( 0 ) );
	}

	void QRuncherTest::testCopyConstructor () {
		// test vectors with norms: 15, 13, 25
		const double data[] = {
			3., 12.,	2.,
			4., 15.,	5.,
			12., 16.,	14.
		};

		AutoMatrix dataMat( 3, 3 );
		dataMat.fill( data );
		const Vector
			x0 = dataMat.columnVector( 0 ),
			x1 = dataMat.columnVector( 1 ),
			y = dataMat.columnVector( 2 );

		QRuncher qr( y );
		QRuncher qrCopy0( qr );

		assert_eq(
			"Zero variables same RSS",
			qr.calculateRSS(),
			qrCopy0.calculateRSS()
		);
		AutoVector coefficients( 0 );
		qr.calculateCoefficients( coefficients );
		AutoVector coefficientsCopy0( 0 );
		qrCopy0.calculateCoefficients( coefficientsCopy0 );
		assert_eq(
			"Zero variables same coefficients",
			coefficients,
			coefficientsCopy0
		);

		assert_close(
			"Original push linearly independent",
			13.0,
			qr.pushColumn( x0 ),
			1e-14
		);
		assert_eq(
			"Copy still zero variables RSS",
			225.0,
			qrCopy0.calculateRSS()
		);
		// copy has still zero coefficients
		qrCopy0.calculateCoefficients( coefficientsCopy0 );
		assert_close(
			"Copy push same linearly independent",
			13.0,
			qrCopy0.pushColumn( x0 ),
			1e-14
		);
		assert_eq(
			"Zero plus 1 variables same RSS 0",
			qr.calculateRSS(),
			qrCopy0.calculateRSS()
		);
		coefficients.upSize( 1 );
		qr.calculateCoefficients( coefficients );
		coefficientsCopy0.upSize( 1 );
		qrCopy0.calculateCoefficients( coefficientsCopy0 );
		assert_eq(
			"Zero plus 1 variables same coefficients 1",
			coefficients,
			coefficientsCopy0
		);

		QRuncher qrCopy1( qr );
		assert_eq(
			"One variable same RSS 1",
			qr.calculateRSS(),
			qrCopy1.calculateRSS()
		);
		AutoVector coefficientsCopy1( 1 );
		qrCopy1.calculateCoefficients( coefficientsCopy1 );
		assert_eq(
			"One variable same coefficients 1",
			coefficients,
			coefficientsCopy1
		);

		TestQRuncher tqr( y );
		tqr.pushColumn( x0 );
		const double qs1 = tqr.pushColumn( x1 );

		assert_close(
			"Original push second linearly independent",
			qs1,
			qr.pushColumn( x1 ),
			1e-14
		);
		assert_eq(
			"Copies still same one variable RSS",
			qrCopy0.calculateRSS(),
			qrCopy1.calculateRSS()
		);
		coefficients.upSize( 2 );
		qr.calculateCoefficients( coefficients );
		qrCopy0.calculateCoefficients( coefficientsCopy0 );
		assert_eq(
			"Copies still same one variable coefficients",
			coefficientsCopy0,
			coefficientsCopy1
		);
		qrCopy1.calculateCoefficients( coefficientsCopy1 );
		assert_eq(
			"Copies again still same one variable coefficients",
			coefficientsCopy0,
			coefficientsCopy1
		);
		assert_close(
			"Copy push second linearly independent",
			qs1,
			qrCopy1.pushColumn( x1 ),
			1e-14
		);
		assert_eq(
			"Zero plus 1 variables same RSS 2",
			qr.calculateRSS(),
			qrCopy1.calculateRSS()
		);
		coefficientsCopy1.upSize( 2 );
		qrCopy1.calculateCoefficients( coefficientsCopy1 );
		assert_eq(
			"Zero plus 1 variables same coefficients 2",
			coefficients,
			coefficientsCopy1
		);

		QRuncher qrCopy2( qr );
		assert_eq(
			"Two variables same RSS",
			qr.calculateRSS(),
			qrCopy2.calculateRSS()
		);
		AutoVector coefficientsCopy2( 2 );
		qrCopy2.calculateCoefficients( coefficientsCopy2 );
		assert_eq(
			"Two variables same coefficients",
			coefficients,
			coefficientsCopy2
		);

		qrCopy2.popColumn();
		assert_eq(
			"Same one var RSS after pop",
			qrCopy0.calculateRSS(),
			qrCopy2.calculateRSS()
		);
		coefficientsCopy2.upSize( 1 );
		qrCopy2.calculateCoefficients( coefficientsCopy2 );
		assert_eq(
			"Same one var coefficients after pop",
			coefficientsCopy0,
			coefficientsCopy2
		);
		qrCopy2.popColumn();
		assert_eq(
			"Initial RSS after second pop",
			225.0,
			qrCopy2.calculateRSS()
		);
	}

	void QRuncherTest::testAssignment () {
		// test vectors with norms: 17, 25, 29
		const double data[] = {
			7., 20.,	8.,
			24., 21.,	15.
		};

		AutoMatrix dataMat( 2, 3 );
		dataMat.fill( data );
		const Vector
			x0 = dataMat.columnVector( 0 ),
			x1 = dataMat.columnVector( 1 ),
			y = dataMat.columnVector( 2 );

		TestQRuncher tqr( y );
		QRuncher qr( y );
		QRuncher qrCopy0( AutoVector( 0 ) );
		assert_eq( "Trivial RSS", 0, qrCopy0.calculateRSS() );
		AutoVector
			coefficients( 0 ),
			coefficientsCopy0( 0 );

		qrCopy0 = qr;
		assert_eq(
			"Zero variables same RSS",
			qr.calculateRSS(),
			qrCopy0.calculateRSS()
		);
		qr.calculateCoefficients( coefficients );
		qrCopy0.calculateCoefficients( coefficientsCopy0 );
		assert_eq(
			"Zero variables same coefficients",
			coefficients,
			coefficientsCopy0
		);

		const double qs1 = tqr.pushColumn( x0 );
		assert_close(
			"Original push linearly independent 1",
			qs1,
			qr.pushColumn( x0 )
		);
		assert_eq(
			"Copy still zero variables RSS",
			289.0,
			qrCopy0.calculateRSS()
		);
		coefficients.upSize( 1 );
		qr.calculateCoefficients( coefficients );
		qrCopy0.calculateCoefficients( coefficientsCopy0 );
		assert_close(
			"Copy push linearly independent 1",
			qs1,
			qrCopy0.pushColumn( x0 )
		);
		assert_eq(
			"Zero plus 1 variables same RSS 1",
			qr.calculateRSS(),
			qrCopy0.calculateRSS()
		);
		coefficientsCopy0.upSize( 1 );
		qrCopy0.calculateCoefficients( coefficientsCopy0 );
		assert_eq(
			"Zero plus 1 variables same coefficients 1",
			coefficients,
			coefficientsCopy0
		);

		QRuncher qrCopy1( x0 );
		qrCopy1 = qr;
		assert_eq(
			"One variable same RSS",
			qr.calculateRSS(),
			qrCopy1.calculateRSS()
		);
		AutoVector coefficientsCopy1( 1 );
		qrCopy1.calculateCoefficients( coefficientsCopy1 );
		assert_eq(
			"One variable same coefficients",
			coefficients,
			coefficientsCopy1
		);

		const double qs2 = tqr.pushColumn( x1 );
		assert_close(
			"Original push second linearly independent",
			qs2,
			qr.pushColumn( x1 )
		);
		assert_eq(
			"Copies still same one variable RSS",
			qrCopy0.calculateRSS(),
			qrCopy1.calculateRSS()
		);
		qrCopy0.calculateCoefficients( coefficientsCopy0 );
		qrCopy1.calculateCoefficients( coefficientsCopy1 );
		assert_eq(
			"Copies still same one variable coefficients",
			coefficientsCopy0,
			coefficientsCopy1
		);
		assert_close(
			"Copy push second linearly independent",
			qs2,
			qrCopy1.pushColumn( x1 )
		);
		assert_eq(
			"Zero plus 1 variables same RSS 2",
			qr.calculateRSS(),
			qrCopy1.calculateRSS()
		);
		coefficients.upSize( 2 );
		coefficientsCopy1.upSize( 2 );
		qr.calculateCoefficients( coefficients );
		qrCopy1.calculateCoefficients( coefficientsCopy1 );
		assert_eq(
			"Zero plus 1 variables same coefficients 2",
			coefficients,
			coefficientsCopy1
		);

		QRuncher qrCopy2( y );
		qrCopy2 = qr;
		assert_eq(
			"Two variables same RSS",
			qr.calculateRSS(),
			qrCopy2.calculateRSS()
		);
		qr.calculateCoefficients( coefficients );
		AutoVector coefficientsCopy2( 2 );
		qrCopy2.calculateCoefficients( coefficientsCopy2 );
		assert_eq(
			"Two variables same coefficients",
			coefficients,
			coefficientsCopy2
		);

		qrCopy2.popColumn();
		assert_eq(
			"Same one var RSS after pop",
			qrCopy0.calculateRSS(),
			qrCopy2.calculateRSS()
		);
		qrCopy0.calculateCoefficients( coefficientsCopy0 );
		coefficientsCopy2.upSize( 1 );
		qrCopy2.calculateCoefficients( coefficientsCopy2 );
		assert_eq(
			"Same one var coefficients after pop",
			coefficientsCopy0,
			coefficientsCopy2
		);
		qrCopy2.popColumn();
		assert_eq(
			"Initial RSS after second pop",
			289.0,
			qrCopy2.calculateRSS()
		);
	}

	void QRuncherTest::testSkipColumn () {
		// col norms: 11, 15, 18, 21, 18, letzte ist minus zweite
		const double data[] = {
			1.,	2.,	3.,	-4.,	-5.,	-2.,
			-2.,	6.,	5.,	-5.,	7.,	-6.,
			4.,	8.,	-11.,	12.,	-9.,	-8.,
			-10.,	11.,	-13.,	16.,	13.,	-11.
		};

		AutoMatrix dataMat( 4, 6 );
		dataMat.fill( data );
		const Vector
			x0 = dataMat.columnVector( 0 ),
			x1 = dataMat.columnVector( 1 ),
			x2 = dataMat.columnVector( 2 ),
			x3 = dataMat.columnVector( 3 ),
			x4 = dataMat.columnVector( 4 ),
			y = dataMat.columnVector( 5 );

		TestQRuncher tqr( y );
		QRuncher qr( y );

		// Prepare comparison model RSS values
		tqr.pushColumn( x0 );
		assert_close(
			"[0] push linearly independent",
			11.0,
			qr.pushColumn( x0 ),
			1e-14
		);
		double rss0 = qr.calculateRSS();
		assert_close(
			"[0] partial match, significant remainder",
			tqr.calculateRSS(),
			rss0,
			1e-13
		);
		assert_close(
			"[02] push linearly independent",
			tqr.pushColumn( x2 ),
			qr.pushColumn( x2 )
		);
		double rss02 = qr.calculateRSS();
		assert_close(
			"[02] partial match",
			tqr.calculateRSS(),
			rss02
		);
		assert_close(
			"[023] push linearly independent",
			tqr.pushColumn( x3 ),
			qr.pushColumn( x3 ),
			1e-14
		);
		const double rss023 = qr.calculateRSS();
		assert_close(
			"[023] partial match",
			tqr.calculateRSS(),
			rss023,
			1e-12
		);
		assert_true( "[023] pop", qr.popColumn() );
		assert_eq( "[02]b RSS", rss02, qr.calculateRSS() );
		assert_true( "[02]b pop", qr.popColumn() );
		assert_eq( "[0]b RSS", rss0, qr.calculateRSS() );
		assert_true( "[0]b pop", qr.popColumn() );

		tqr.popColumn();
		tqr.popColumn();
		tqr.popColumn();
		assert_eq(
			"[] zero variables, full RSS tqr check",
			225.0,
			tqr.calculateRSS()
		);
		assert_eq(
			"[] zero variables, full RSS",
			225.0,
			qr.calculateRSS()
		);

		assert_close(
			"[0]c push linearly independent",
			tqr.pushColumn( x0 ),
			qr.pushColumn( x0 ),
			1e-14
		);
		assert_eq( "[0]c RSS", rss0, qr.calculateRSS() );
		assert_close(
			"[0]c back to zero variables",
			225.0,
			qr.calculateSkipColumnRSS( 0 ),
			1e-13
		);
		assert_eq( "[0]c RSS unchanged", rss0, qr.calculateRSS() );

		assert_close(
			"[01] push linearly independent",
			tqr.pushColumn( x1 ),
			qr.pushColumn( x1 ),
			1e-14
		);
		const double rss01 = qr.calculateRSS();
		assert_close( "[01] exact match", 0.0, rss01 );
		assert_close( "[01] back to exact match", 0.0, qr.calculateSkipColumnRSS( 0 ) );
		assert_close(
			"[01] back to one variable, from y-stack",
			rss0,
			qr.calculateSkipColumnRSS( 1 ),
			1e-12
		);
		assert_eq( "[01] RSS unchanged", rss01, qr.calculateRSS() );

		assert_close(
			"[012] push linearly independent",
			tqr.pushColumn( x2 ),
			qr.pushColumn( x2 ),
			1e-14
		);
		const double rss012 = qr.calculateRSS();
		assert_close( "[012] exact match", 0.0, rss012 );
		assert_close(
			"[012] back to still exact match",
			0.0,
			qr.calculateSkipColumnRSS( 0 )
		);
		assert_close(
			"[012] back to inexact match",
			rss02,
			qr.calculateSkipColumnRSS( 1 ),
			1e-12
		);
		assert_close(
			"[012] back to one variable, from y-stack",
			rss01,
			qr.calculateSkipColumnRSS( 2 ),
			1e-12
		);
		assert_eq( "[012] RSS unchanged", rss012, qr.calculateRSS() );

		assert_close(
			"[0123] push linearly independent",
			tqr.pushColumn( x3 ),
			qr.pushColumn( x3 ),
			1e-14
		);
		const double rss0123 = qr.calculateRSS();
		assert_close( "[0123] exact match", 0.0, rss0123 );
		assert_close(
			"[0123] Back to still exact match",
			0.0,
			qr.calculateSkipColumnRSS( 0 )
		);
		assert_close(
			"[0123] Back to inexact match",
			rss023,
			qr.calculateSkipColumnRSS( 1 ),
			1e-12
		);
		assert_close(
			"[0123] other back to still exact match",
			0.0,
			qr.calculateSkipColumnRSS( 2 )
		);
		assert_close(
			"[0123] back to one variable, from y-stack",
			rss012,
			qr.calculateSkipColumnRSS( 3 ),
			1e-12
		);
		assert_eq( "[0123] RSS unchanged", rss0123, qr.calculateRSS() );

		assert_true( "Can pop", qr.popColumn() );
		tqr.popColumn();

		assert_close(
			"[0124] push linearly independent",
			tqr.pushColumn( x4 ),
			qr.pushColumn( x4 ),
			1e-14
		);
		const double rss0124 = qr.calculateRSS();
		assert_close( "[0124] exact match", 0, rss0124 );
		assert_close(
			"[0124] Back to still exact match",
			0.0,
			qr.calculateSkipColumnRSS( 0 )
		);
		assert_close(
			"[0124] again exact match",
			0.0,
			qr.calculateSkipColumnRSS( 2 )
		);
		assert_close(
			"[0124] back to one variable, from y-stack",
			rss012,
			qr.calculateSkipColumnRSS( 3 ),
			1e-12
		);
		assert_eq( "[0124] RSS unchanged", rss0124, qr.calculateRSS() );
	}

	void QRuncherTest::testSkipLastLinearlyIndependentColumn () {
		const double data[] = {
			1., 0.,		3.,
			0., 1.,		4.
		};

		AutoMatrix dataMat( 2, 3 );
		dataMat.fill( data );
		const Vector
			x0 = dataMat.columnVector( 0 ),
			x1 = dataMat.columnVector( 1 ),
			y = dataMat.columnVector( 2 );

		QRuncher qr( y );
		assert_eq( "[] RSS", 25.0, qr.calculateRSS() );
		assert_eq( "[0] linearly independent", 1.0, qr.pushColumn( x0 ) );
		assert_eq( "[0] RSS", 16.0, qr.calculateRSS() );
		assert_eq( "[0]-0 RSS", 25.0, qr.calculateSkipColumnRSS( 0 ) );
		assert_eq( "[01] push linearly independent", 1.0, qr.pushColumn( x1 ) );
		assert_eq( "[01] RSS", 0.0, qr.calculateRSS() );
		assert_eq( "[01]-0 RSS", 9.0, qr.calculateSkipColumnRSS( 0 ) );
		assert_eq( "[01]-1 RSS", 16.0, qr.calculateSkipColumnRSS( 1 ) );
		assert_true( "pop 1", qr.popColumn() );
		assert_eq( "[0]b RSS", 16.0, qr.calculateRSS() );
		assert_eq( "[0]b-0 RSS", 25.0, qr.calculateSkipColumnRSS( 0 ) );
		assert_true( "pop 0", qr.popColumn() );
		assert_eq( "[]b RSS", 25.0, qr.calculateRSS() );
		assert_false( "pop too many", qr.popColumn() );
		assert_eq( "postpoptoomany RSS", 25.0, qr.calculateRSS() );
	}

	void QRuncherTest::testLinearlyDependent () {
		// col 0 norm 0, col 1 norm 2, col 2 norm 5, col 3 = col 2 - col 1, col 4 = col 2 + col 1
		const double data[] = {
			0., 1., 1., 0., 2.,	2.,
			0., 1., 2., 1., 3.,	4.,
			0., 1., 2., 1., 3.,	5.,
			0., 1., 4., 3., 5.,	6.
		};
		// y norm 9

		AutoMatrix dataMat( 4, 6 );
		dataMat.fill( data );
		const Vector
			x0 = dataMat.columnVector( 0 ),
			x1 = dataMat.columnVector( 1 ),
			x2 = dataMat.columnVector( 2 ),
			x3 = dataMat.columnVector( 3 ),
			x4 = dataMat.columnVector( 4 ),
			y = dataMat.columnVector( 5 );

		// Calculate master comparison values by Gram-Schmidt algorithm
		const double
			squareSum1 = x1.sumSquares(),
			innerProductY1 = y.innerProduct( x1 ),
			mean = innerProductY1 / squareSum1;	// in case of column[1], that's the mean indeed

		TestQRuncher tqr( y );
		QRuncher qr( y );

		assert_eq( "[] RSS", 81.0, qr.calculateRSS() );
		// 0 coefficients are enough
		AutoVector coefficients( 0 );
		qr.calculateCoefficients( coefficients );

		// push column 0
		assert_eq( "[0] zero vector is linearly dependent", 0.0, qr.pushColumn( x0 ) );
		assert_close( "[0] full remainder", 81.0, qr.calculateRSS() );
		AutoVector coefficients0( 1 );
		qr.calculateCoefficients( coefficients0 );
		assert_eq( "[0] no match attempt coefficients", 0.0, coefficients0.get( 0 ) );

		// push column 1
		assert_close(
			"[01] push linearly independent",
			tqr.pushColumn( x1 ),
			qr.pushColumn( x1 )
		);
		const double rss1 = tqr.calculateRSS();
		assert_close( "[01] reduced remainder", rss1, qr.calculateRSS(), 1e-14 );
		AutoVector coefficients01( 2 );
		qr.calculateCoefficients( coefficients01 );
		assert_eq( "[01] no match attempt coefficient 0", 0.0, coefficients01.get( 0 ) );
		assert_close(
			"[01] match attempt coefficient 1",
			mean,
			coefficients01.get( 1 ),
			1e-14
		);

		// push column 1 again
		assert_eq( "[011] push linearly dependent", 0.0, qr.pushColumn( x1 ) );
		assert_close( "[011] reduced remainder remains", rss1, qr.calculateRSS() , 1e-14 );
		AutoVector coefficients011( 3 );
		qr.calculateCoefficients( coefficients011 );
		assert_eq( "[011] no match attempt coefficient 0", 0.0, coefficients011.get( 0 ) );
		assert_close(
			"[011] match attempt coefficient 1",
			mean,
			coefficients011.get( 1 ),
			1e-14
		);
		assert_eq( "[011] no match attempt coefficient 2", 0.0, coefficients011.get( 2 ) );

		// push column 2
		assert_close(
			"[0112] push linearly independent",
			tqr.pushColumn( x2 ),
			qr.pushColumn( x2 )
		);
		const double rss12 = tqr.calculateRSS();
		assert_close( "[0112] still no exact match", rss12, qr.calculateRSS() );
		AutoVector coefficients0112( 4 );
		qr.calculateCoefficients( coefficients0112 );
		assert_eq( "[0112] no match attempt coefficient 0", 0.0, coefficients0112.get( 0 ) );
		assert_eq( "[0112] no match attempt coefficient 2", 0.0, coefficients0112.get( 2 ) );
		AutoMatrix r( 4, 4 );	// the regression matrix reflecting current push sequence
		r.columnVector( 0 ).copy( x0 );
		r.columnVector( 1 ).copy( x1 );
		r.columnVector( 2 ).copy( x1 );
		r.columnVector( 3 ).copy( x2 );
		AutoVector c( 4 );
		c.copy( y );
		c.gemv( -1.0, r, false, coefficients0112, 1.0 );	// remaincer c = y - r*coefficients
		// remainder is minimal if and only if it is orthogonal the spanning vectors
		assert_close( "[0112] remainder orthogonal 1", 0.0, c.innerProduct( x1 ), 1e-12 );
		assert_close( "[0112] remainder orthogonal 2", 0.0, c.innerProduct( x2 ), 1e-12 );

		// push column 2 again
		// Due to rounding error, pushColumn does cannot detect linear dependence:
		// assert_false( "[01122] push linearly dependent", qr.pushColumn( x2 ) );

		// Also the following would fail,
		// because the push implies approximation by a 3-dimensional space
		// where one direction is virtually random
		// qr.pushColumn( x2 );
		// assert_close( "[01122] still no exact match", rss12, qr.calculateRSS() );

		// This is the reason why a separate test testLimitedRounding() with no rounding necessary has been added.
	}

	void QRuncherTest::testLimitedRounding () {
		const double data[] = {
			1., 0., 0., 0., 0., 0.,		84.,
			0., 1., 4., 0., 12., 0.,	12.,
			0., 0., 0., 4., 4., 0.,		4.,
			0., 0., 0., 3., 3., 2.,		3.
		};
		// y norm is 85

		AutoMatrix dataMat( 4, 7 );
		dataMat.fill( data );
		const Vector
			x0 = dataMat.columnVector( 0 ),
			x1 = dataMat.columnVector( 1 ),
			x2 = dataMat.columnVector( 2 ),
			x3 = dataMat.columnVector( 3 ),
			x4 = dataMat.columnVector( 4 ),
			x5 = dataMat.columnVector( 5 ),
			y = dataMat.columnVector( 6 );

		QRuncher qr( y );
		assert_eq( "[] rss", 85.0*85.0, qr.calculateRSS() );
		AutoVector coefficients( 0 );
		qr.calculateCoefficients( coefficients );	// asserting that 0 dim is enough

		assert_eq( "[0] push linearly independent", 1.0, qr.pushColumn( x0 ) );
		assert_eq( "[0] rss", 169.0, qr.calculateRSS() );
		assert_eq( "[0]-0 skip rss", 85.0*85.0, qr.calculateSkipColumnRSS( 0 ) );
		AutoVector coefficients0( 1 );
		qr.calculateCoefficients( coefficients0 );
		assert_eq( "[0] coefficient 0", 84.0, coefficients0.get( 0 ) );

		assert_eq( "[01] push linearly independent", 1.0, qr.pushColumn( x1 ) );
		assert_eq( "[01] rss", 25.0, qr.calculateRSS() );
		assert_eq( "[01]-0 skip rss", 84.0*84.0 + 25.0, qr.calculateSkipColumnRSS( 0 ) );
		assert_eq( "[01]-1 skip rss", 169.0, qr.calculateSkipColumnRSS( 1 ) );
		AutoVector coefficients01( 2 );
		qr.calculateCoefficients( coefficients01 );
		assert_eq( "[01] coefficient 0", 84.0, coefficients01.get( 0 ) );
		assert_eq( "[01] coefficient 1", 12.0, coefficients01.get( 1 ) );

		assert_eq( "[012] push linearly dependent", 0.0, qr.pushColumn( x2 ) );
		assert_eq( "[012] rss", 25.0, qr.calculateRSS() );
		assert_eq( "[012]-0 skip rss", 84.0*84.0 + 25.0, qr.calculateSkipColumnRSS( 0 ) );
		assert_eq( "[012]-1 skip rss, 2 replacing 1", 25.0, qr.calculateSkipColumnRSS( 1 ) );
		assert_eq( "[012]-2 skip rss, 2 superfluous", 25.0, qr.calculateSkipColumnRSS( 2 ) );
		AutoVector coefficients012( 3 );
		qr.calculateCoefficients( coefficients012 );
		assert_eq( "[012] coefficient 0", 84.0, coefficients012.get( 0 ) );
		assert_eq( "[012] coefficient 1", 12.0, coefficients012.get( 1 ) );
		assert_eq( "[012] coefficient 2", 0.0, coefficients012.get( 2 ) );

		assert_eq( "[0123] push linearly independent", 5.0, qr.pushColumn( x3 ) );
		assert_eq( "[0123] rss", 0.0, qr.calculateRSS() );
		assert_eq( "[0123]-0 skip rss", 84.0*84.0, qr.calculateSkipColumnRSS( 0 ) );
		assert_eq( "[0123]-1 skip rss, 2 replacing 1", 0.0, qr.calculateSkipColumnRSS( 1 ) );
		assert_eq( "[0123]-2 skip rss, 2 superfluous", 0.0, qr.calculateSkipColumnRSS( 2 ) );
		assert_eq( "[0123]-3 skip rss", 25.0, qr.calculateSkipColumnRSS( 3 ) );
		AutoVector coefficients0123( 4 );
		qr.calculateCoefficients( coefficients0123 );
		assert_eq( "[0123] coefficient 0", 84.0, coefficients0123.get( 0 ) );
		assert_eq( "[0123] coefficient 1", 12.0, coefficients0123.get( 1 ) );
		assert_eq( "[0123] coefficient 2", 0.0, coefficients0123.get( 2 ) );
		assert_eq( "[0123] coefficient 3", 1.0, coefficients0123.get( 3 ) );

		assert_eq( "[01234] push linearly dependent", 0.0, qr.pushColumn( x4 ) );
		assert_eq( "[01234] rss", 0.0, qr.calculateRSS() );
		assert_eq( "[01234]-0 skip rss", 84.0*84.0, qr.calculateSkipColumnRSS( 0 ) );
		assert_eq( "[01234]-1 skip rss, 2 replacing 1", 0.0, qr.calculateSkipColumnRSS( 1 ) );
		assert_eq( "[01234]-2 skip rss, 2 superfluous", 0.0, qr.calculateSkipColumnRSS( 2 ) );
		assert_eq( "[01234]-3 skip rss, 4+1 replacing 3", 0.0, qr.calculateSkipColumnRSS( 3 ) );
		assert_eq( "[01234]-4 skip rss, 4 superfluous", 0.0, qr.calculateSkipColumnRSS( 4 ) );
		AutoVector coefficients01234( 5 );
		qr.calculateCoefficients( coefficients01234 );
		assert_eq( "[01234] coefficient 0", 84.0, coefficients01234.get( 0 ) );
		assert_eq( "[01234] coefficient 1", 12.0, coefficients01234.get( 1 ) );
		assert_eq( "[01234] coefficient 2", 0.0, coefficients01234.get( 2 ) );
		assert_eq( "[01234] coefficient 3", 1.0, coefficients01234.get( 3 ) );
		assert_eq( "[01234] coefficient 4", 0.0, coefficients01234.get( 4 ) );

		assert_true( "[0123]b pop", qr.popColumn() );
		assert_eq( "[0123]b rss", 0.0, qr.calculateRSS() );
		AutoVector coefficients0123b( 4 );
		qr.calculateCoefficients( coefficients0123b );
		assert_eq( "[0123]b coefficients", coefficients0123, coefficients0123b );

		assert_true( "[012]b pop", qr.popColumn() );
		assert_eq( "[012]b rss", 25.0, qr.calculateRSS() );
		AutoVector coefficients012b( 3 );
		qr.calculateCoefficients( coefficients012b );
		assert_eq( "[012]b coefficients", coefficients012, coefficients012b );

		assert_true( "[01]b pop", qr.popColumn() );
		assert_eq( "[01]b rss", 25.0, qr.calculateRSS() );
		AutoVector coefficients01b( 2 );
		qr.calculateCoefficients( coefficients01b );
		assert_eq( "[01]b coefficients", coefficients01, coefficients01b );

		assert_true( "[0]b pop", qr.popColumn() );
		assert_eq( "[0]b rss", 169.0, qr.calculateRSS() );
		AutoVector coefficients0b( 1 );
		qr.calculateCoefficients( coefficients0b );
		assert_eq( "[0]b coefficients", coefficients0, coefficients0b );

		assert_true( "[04] push linearly independent", 1.0 < qr.pushColumn( x4 ) );
		assert_eq( "[04] rss", 0.0, qr.calculateRSS() );
		assert_eq( "[04]-0 skip rss", 84.0*84.0, qr.calculateSkipColumnRSS( 0 ) );
		assert_eq( "[04]-4 skip rss", 169.0, qr.calculateSkipColumnRSS( 1 ) );
		AutoVector coefficients04( 2 );
		qr.calculateCoefficients( coefficients04 );
		assert_eq( "[04] coefficient 0", 84.0, coefficients04.get( 0 ) );
		assert_eq( "[04] coefficient 4", 1.0, coefficients04.get( 1 ) );

		assert_true( "[043] push linearly independent", 1.0 < qr.pushColumn( x3 ) );
		assert_eq( "[043] rss", 0.0, qr.calculateRSS() );
		assert_eq( "[043]-0 skip rss", 84.0*84.0, qr.calculateSkipColumnRSS( 0 ) );
		assert_close( "[043]-4 skip rss", 144.0, qr.calculateSkipColumnRSS( 1 ), 1e-12 );
		assert_eq( "[043]-3 skip rss", 0.0, qr.calculateSkipColumnRSS( 2 ) );
		AutoVector coefficients043( 3 );
		qr.calculateCoefficients( coefficients043 );
		assert_eq( "[043] coefficient 0", 84.0, coefficients043.get( 0 ) );
		assert_eq( "[043] coefficient 4", 1.0, coefficients043.get( 1 ) );
		assert_eq( "[043] coefficient 3", 0.0, coefficients043.get( 2 ) );

		assert_eq( "[0432] push linearly dependent", 0.0, qr.pushColumn( x4 ) );
		assert_eq( "[0432] rss", 0.0, qr.calculateRSS() );
		assert_eq( "[0432]-0 skip rss", 84.0*84.0, qr.calculateSkipColumnRSS( 0 ) );
		assert_eq( "[0432]-4 skip rss, 3+2 replacing 4", 0.0, qr.calculateSkipColumnRSS( 1 ) );
		assert_eq( "[0432]-3 skip rss, 3 not needed for y match", 0.0, qr.calculateSkipColumnRSS( 2 ) );
		assert_eq( "[0432]-2 skip rss, 2 superfluous", 0.0, qr.calculateSkipColumnRSS( 3 ) );
		AutoVector coefficients0432( 4 );
		qr.calculateCoefficients( coefficients0432 );
		assert_eq( "[0432] coefficient 0", 84.0, coefficients0432.get( 0 ) );
		assert_eq( "[0432] coefficient 4", 1.0, coefficients0432.get( 1 ) );
		assert_eq( "[0432] coefficient 3", 0.0, coefficients0432.get( 2 ) );
		assert_eq( "[0432] coefficient 2", 0.0, coefficients0432.get( 3 ) );

		assert_true( "[04325] push linearly independent", 1.0 < qr.pushColumn( x5 ) );
		assert_eq( "[04325] rss", 0.0, qr.calculateRSS() );
		assert_eq( "[04325]-0 skip rss", 84.0*84.0, qr.calculateSkipColumnRSS( 0 ) );
		assert_eq( "[04325]-4 skip rss, 3+2 replacing 4", 0.0, qr.calculateSkipColumnRSS( 1 ) );
		assert_eq( "[04325]-3 skip rss, 3 not needed for y match", 0.0, qr.calculateSkipColumnRSS( 2 ) );
		assert_eq( "[04325]-2 skip rss, 2 superfluous", 0.0, qr.calculateSkipColumnRSS( 3 ) );
		assert_eq( "[04325]-5 skip rss, 5 not needed for y match", 0.0, qr.calculateSkipColumnRSS( 4 ) );
		AutoVector coefficients04325( 5 );
		qr.calculateCoefficients( coefficients04325 );
		assert_eq( "[04325] coefficient 0", 84.0, coefficients04325.get( 0 ) );
		assert_eq( "[04325] coefficient 4", 1.0, coefficients04325.get( 1 ) );
		assert_eq( "[04325] coefficient 3", 0.0, coefficients04325.get( 2 ) );
		assert_eq( "[04325] coefficient 2", 0.0, coefficients04325.get( 3 ) );
		assert_eq( "[04325] coefficient 5", 0.0, coefficients04325.get( 4 ) );
	}

	void QRuncherTest::testRightLowerCorner () {
		const double data[] = {
			9., 2., -1.,
			4., -8., 5.,
			3., -5., 7.
		};
		AutoMatrix dataMat( 3, 3 );
		dataMat.fill( data );

		// test even dimension 2x2
		{
			const Matrix x = dataMat.subMatrix( 0, 0, 2, 2 );

			const Vector
				x0 = const_cast<Matrix&>( x ).columnVector( 0 ),
				x1 = const_cast<Matrix&>( x ).columnVector( 1 );
			QRuncher qr( x1 );	// y is irrelevant for this test case
			qr.pushColumn( x0 );
			const double rii = qr.pushColumn( x1 );

			AutoMatrix
				xTx( 2, 2 ),
				inverse( 2, 2 );
			AutoPermutation p( 2 );
			xTx.gemm( 1., x, true, x, false, 0. );
			xTx.factorizeLUP( p );
			inverse.invertLUP( xTx, p );

			assert_close(
				"2x2 corner close",
				inverse.get( 1, 1 ),
				1. / ( rii * rii )
			);
		}

		// test mixed dimension 3x2
		{
			const Matrix x = dataMat.subMatrix( 0, 0, 3, 2 );

			const Vector
				x0 = const_cast<Matrix&>( x ).columnVector( 0 ),
				x1 = const_cast<Matrix&>( x ).columnVector( 1 );
			QRuncher qr( x0 );	// y is irrelevant for this test case
			qr.pushColumn( x0 );
			const double rii = qr.pushColumn( x1 );

			AutoMatrix
				xTx( 2, 2 ),
				inverse( 2, 2 );
			AutoPermutation p( 2 );
			xTx.gemm( 1., x, true, x, false, 0. );
			xTx.factorizeLUP( p );
			inverse.invertLUP( xTx, p );

			assert_close(
				"3x2 corner close",
				inverse.get( 1, 1 ),
				1. / ( rii * rii )
			);
		}

		// test odd dimension 3x3
		{
			const Matrix x = dataMat;

			const Vector
				x0 = const_cast<Matrix&>( x ).columnVector( 0 ),
				x1 = const_cast<Matrix&>( x ).columnVector( 1 ),
				x2 = const_cast<Matrix&>( x ).columnVector( 2 );
			QRuncher qr( x2 );	// y is irrelevant for this test case
			qr.pushColumn( x0 );
			qr.pushColumn( x1 );
			const double rii = qr.pushColumn( x2 );

			AutoMatrix
				xTx( 3, 3 ),
				inverse( 3, 3 );
			AutoPermutation p( 3 );
			xTx.gemm( 1., x, true, x, false, 0. );
			xTx.factorizeLUP( p );
			inverse.invertLUP( xTx, p );

			assert_close(
				"3x3 corner close",
				inverse.get( 2, 2 ),
				1. / ( rii * rii )
			);
		}
	}

}
