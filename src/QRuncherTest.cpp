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
#include <math.h>

using namespace std;
using namespace linalg;
using namespace unitpp;

namespace test {

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
		}
	} * qrUncherTest = new QRuncherTest();	// automatically freed by unit++

	void QRuncherTest::testSequence () {
		const double
			// col 1 norm 9, col 2 norm 7, col 3 = col 2 - col 1
			xData[] = {
				1, 2, 1,
				4, 3, -1,
				8, 6, -2
			},
			yData[] = { 1, 4, 8 };

		AutoMatrix x( 3, 3 );
		x.fill( xData );
		AutoVector y( 3 );
		y.fill( yData );

		QRuncher qr( y );
		// 0 columns
		assert_eq( "Zero variables", 81.0, qr.calculateRSS() );
		const AutoVector coefficients0a = qr.calculateCoefficients();
		assert_eq( "Zero-dim regression coefficients", 0, coefficients0a.countDimensions() );

		// 1 column
		assert_true( "No linear dep", qr.pushColumn( x.columnVector( 0 ) ) );
		assert_close( "Exact match", 0.0, qr.calculateRSS() );
		const AutoVector coefficients1a = qr.calculateCoefficients();
		assert_eq( "One-dim regression coefficients", 1, coefficients1a.countDimensions() );
		assert_close( "Exact match coefficients", 1.0, coefficients1a.get( 0 ) );

		// 2 columns
		assert_true( "No linear dep", qr.pushColumn( x.columnVector( 1 ) ) );
		assert_close( "Still exact match", 0.0, qr.calculateRSS() );
		const AutoVector coefficients2a = qr.calculateCoefficients();
		assert_eq( "Two-dim regression coefficients", 2, coefficients2a.countDimensions() );
		assert_close( "Exact match coefficient 0", 1.0, coefficients2a.get( 0 ) );
		assert_close( "Exact match coefficient 1", 0.0, coefficients2a.get( 1 ) );

		// 3 columns (linear dependence)
		assert_false( "Full linear dep", qr.pushColumn( x.columnVector( 2 ) ) );
		assert_close( "And still exact match", 0.0, qr.calculateRSS() );
		const AutoVector coefficients3 = qr.calculateCoefficients();
		assert_eq( "Three-dim regression coefficients", 3, coefficients3.countDimensions() );
		assert_close( "Exact match coefficient 0", 1.0, coefficients3.get( 0 ) );
		assert_close( "Exact match coefficient 1", 0.0, coefficients3.get( 1 ) );
		assert_close( "Exact match coefficient 2", 0.0, coefficients3.get( 2 ) );

		// back to 2 columns
		qr.popColumn();
		assert_close( "Back to still exact match", 0.0, qr.calculateRSS() );
		const AutoVector coefficients2b = qr.calculateCoefficients();
		assert_eq( "Back to two-dim regression coefficients", coefficients2a, coefficients2b );

		// back to 1 column
		qr.popColumn();
		assert_close( "Back to exact match", 0.0, qr.calculateRSS() );
		const AutoVector coefficients1b = qr.calculateCoefficients();
		assert_eq( "Back to one-dim regression coefficients", coefficients1a, coefficients1b );

		// back to 0 columns
		qr.popColumn();
		assert_eq( "Back to zero variables", 81.0, qr.calculateRSS() );
		const AutoVector coefficients0b = qr.calculateCoefficients();
		assert_eq( "Back to zero-dim regression coefficients", coefficients0a, coefficients0b );
	}

	void QRuncherTest::testCoefficients () {
		const double
			xData[] = {
				3, 4,
				4, 3
			},
			yData[] = { 8, 8 };

		AutoMatrix x( 2, 2 );
		x.fill( xData );
		AutoVector y( 2 );
		y.fill( yData );

		QRuncher qr( y );
		assert_eq( "No variables yet", 0, qr.calculateCoefficients().countDimensions() );

		assert_true( "No linear dep 0", qr.pushColumn( x.columnVector( 0 ) ) );
		// Minimise || (8,8) - t * (3,4) ||^2 = (8-3t)^2+(8-4t)^2.
		// Derivative proportional 3(8-3t) + 4(8-4t) = 24-9t + 32-16t = 56-25t
		assert_close( "Approximate match coefficient 0", 56./25., qr.calculateCoefficients().get( 0 ) );

		assert_true( "No linear dep 1", qr.pushColumn( x.columnVector( 1 ) ) );
		assert_close( "Exact match", 0, qr.calculateRSS() );
		assert_close( "Exact match coefficients 0", 8./7., qr.calculateCoefficients().get( 0 ) );
		assert_close( "Exact match coefficients 1", 8./7., qr.calculateCoefficients().get( 1 ) );

		qr.popColumn();
		assert_close( "Back to approximate match coefficient 0", 56./25., qr.calculateCoefficients().get( 0 ) );
	}

	void QRuncherTest::testCopyConstructor () {

		// test vectors with norms: 15, 13, 25
		const double data[][3] = {
			{ 2, 5, 14 },
			{ 3, 4, 12 },
			{ 12, 15, 16 }
		};

		AutoVector u( 3 ), v( 3 ), w( 3 );
		u.fill( data[0] );
		v.fill( data[1] );
		w.fill( data[2] );

		QRuncher qr( u );
		QRuncher qrc0( qr );
		assert_eq( "Zero variables same RSS", qr.calculateRSS(), qrc0.calculateRSS() );
		assert_eq( "Zero variables same coefficients", qr.calculateCoefficients(), qrc0.calculateCoefficients() );

		assert_true( "Original push linearly independent", qr.pushColumn( v ) );
		assert_eq( "Copy still zero variables RSS", 225, qrc0.calculateRSS() );
		assert_eq( "Copy still zero variables coefficients", 0, qrc0.calculateCoefficients().countDimensions() );
		assert_true( "Copy push same linearly independent", qrc0.pushColumn( v ) );
		assert_eq( "Zero plus 1 variables same RSS", qr.calculateRSS(), qrc0.calculateRSS() );
		assert_eq( "Zero plus 1 variables same coefficients", qr.calculateCoefficients(), qrc0.calculateCoefficients() );

		QRuncher qrc1( qr );
		assert_eq( "One variable same RSS", qr.calculateRSS(), qrc1.calculateRSS() );
		assert_eq( "One variable same coefficients", qr.calculateCoefficients(), qrc1.calculateCoefficients() );

		assert_true( "Original push second linearly independent", qr.pushColumn( w ) );
		assert_eq( "Copies still same one variable RSS", qrc0.calculateRSS(), qrc1.calculateRSS() );
		assert_eq( "Copies still same one variable coefficients", qrc0.calculateCoefficients(), qrc1.calculateCoefficients() );
		assert_true( "Copy push second linearly independent", qrc1.pushColumn( w ) );
		assert_eq( "Zero plus 1 variables same RSS", qr.calculateRSS(), qrc1.calculateRSS() );
		assert_eq( "Zero plus 1 variables same coefficients", qr.calculateCoefficients(), qrc1.calculateCoefficients() );

		QRuncher qrc2( qr );
		assert_eq( "Two variables same RSS", qr.calculateRSS(), qrc2.calculateRSS() );
		assert_eq( "Two variables same coefficients", qr.calculateCoefficients(), qrc2.calculateCoefficients() );

		qrc2.popColumn();
		assert_eq( "Same one var RSS after pop", qrc0.calculateRSS(), qrc2.calculateRSS() );
		assert_eq( "Same one var coefficients after pop", qrc0.calculateCoefficients(), qrc2.calculateCoefficients() );
		qrc2.popColumn();
		assert_eq( "Initial RSS after second pop", 225, qrc2.calculateRSS() );
	}

	void QRuncherTest::testAssignment () {

		// test vectors with norms: 17, 25, 29
		const double data[][2] = {
			{ 8, 15 },
			{ 7, 24 },
			{ 20, 21 }
		};

		AutoVector u( 2 ), v( 2 ), w( 2 );
		u.fill( data[0] );
		v.fill( data[1] );
		w.fill( data[2] );

		QRuncher qr( u );
		QRuncher qrc0( AutoVector( 0 ) );
		assert_eq( "Trivial RSS", 0, qrc0.calculateRSS() );
		assert_eq( "Trivial coefficients", 0, qrc0.calculateCoefficients().countDimensions() );

		qrc0 = qr;
		assert_eq( "Zero variables same RSS", qr.calculateRSS(), qrc0.calculateRSS() );
		assert_eq( "Zero variables same coefficients", qr.calculateCoefficients(), qrc0.calculateCoefficients() );

		assert_true( "Original push linearly independent", qr.pushColumn( v ) );
		assert_eq( "Copy still zero variables RSS", 289, qrc0.calculateRSS() );
		assert_eq( "Copy still zero variables coefficients", 0, qrc0.calculateCoefficients().countDimensions() );
		assert_true( "Copy push linearly independent", qrc0.pushColumn( v ) );
		assert_eq( "Zero plus 1 variables same RSS", qr.calculateRSS(), qrc0.calculateRSS() );
		assert_eq( "Zero plus 1 variables same coefficients", qr.calculateCoefficients(), qrc0.calculateCoefficients() );

		QRuncher qrc1( v );
		qrc1 = qr;
		assert_eq( "One variable same RSS", qr.calculateRSS(), qrc1.calculateRSS() );
		assert_eq( "One variable same coefficients", qr.calculateCoefficients(), qrc1.calculateCoefficients() );

		assert_true( "Original push second linearly independent", qr.pushColumn( w ) );
		assert_eq( "Copies still same one variable RSS", qrc0.calculateRSS(), qrc1.calculateRSS() );
		assert_eq( "Copies still same one variable coefficients", qrc0.calculateCoefficients(), qrc1.calculateCoefficients() );
		assert_true( "Copy push second linearly independent", qrc1.pushColumn( w ) );
		assert_eq( "Zero plus 1 variables same RSS", qr.calculateRSS(), qrc1.calculateRSS() );
		assert_eq( "Zero plus 1 variables same coefficients", qr.calculateCoefficients(), qrc1.calculateCoefficients() );

		QRuncher qrc2( u );
		qrc2 = qr;
		assert_eq( "Two variables same RSS", qr.calculateRSS(), qrc2.calculateRSS() );
		assert_eq( "Two variables same coefficients", qr.calculateCoefficients(), qrc2.calculateCoefficients() );

		qrc2.popColumn();
		assert_eq( "Same one var RSS after pop", qrc0.calculateRSS(), qrc2.calculateRSS() );
		assert_eq( "Same one var coefficients after pop", qrc0.calculateCoefficients(), qrc2.calculateCoefficients() );
		qrc2.popColumn();
		assert_eq( "Initial RSS after second pop", 289, qrc2.calculateRSS() );
	}

	void QRuncherTest::testSkipColumn () {
		const double
			// col norms: 11, 15, 18, 21, 18
			xData[] = {
				1,	2,	3,	-4,	-5,
				-2,	6,	5,	-5,	7,
				4,	8,	-11,	12,	-9,
				-10,	11,	-13,	16,	13
			},
			yData[] = { -2, -6, -8, -11 };

		AutoMatrix x( 4, 5 );
		x.fill( xData );
		AutoVector y( 4 );
		y.fill( yData );

		QRuncher qr( y );
		// Prepare comparison model RSS values
		assert_true( "[0] push linearly independent", qr.pushColumn( x.columnVector( 0 ) ) );
		double rss0 = qr.calculateRSS();
		assert_true( "[0] partial match, significant remainder", 0.5 < rss0 );
		assert_true( "[02] push linearly independent", qr.pushColumn( x.columnVector( 2 ) ) );
		double rss02 = qr.calculateRSS();
		assert_true( "[023] push linearly independent", qr.pushColumn( x.columnVector( 3 ) ) );
		double rss023 = qr.calculateRSS();
		assert_true( "[023] pop", qr.popColumn() );
		assert_eq( "[02]b RSS", rss02, qr.calculateRSS() );
		assert_true( "[02]b pop", qr.popColumn() );
		assert_eq( "[0]b RSS", rss0, qr.calculateRSS() );
		assert_true( "[0]b pop", qr.popColumn() );

		assert_eq( "[] zero variables, full RSS", 225, qr.calculateRSS() );

		assert_true( "[0]c push linearly independent", qr.pushColumn( x.columnVector( 0 ) ) );
		assert_eq( "[0]c RSS", rss0, qr.calculateRSS() );
		assert_eq( "[0]c back to zero variables", 225, qr.calculateSkipColumnRSS( 0 ) );
		assert_eq( "[0]c RSS unchanged", rss0, qr.calculateRSS() );

		assert_true( "[01] push linearly independent", qr.pushColumn( x.columnVector( 1 ) ) );
		const double rss01 = qr.calculateRSS();
		assert_close( "[01] exact match", 0, rss01 );
		assert_close( "[01] back to exact match", 0, qr.calculateSkipColumnRSS( 0 ) );
		assert_close( "[01] back to one variable, from y-stack", rss0, qr.calculateSkipColumnRSS( 1 ), 1e-12 );
		assert_eq( "[01] RSS unchanged", rss01, qr.calculateRSS() );

		assert_true( "[012] push linearly independent", qr.pushColumn( x.columnVector( 2 ) ) );
		const double rss012 = qr.calculateRSS();
		assert_close( "[012] exact match", 0, rss012 );
		assert_close( "[012] back to still exact match", 0, qr.calculateSkipColumnRSS( 0 ) );
		assert_close( "[012] back to inexact match", rss02, qr.calculateSkipColumnRSS( 1 ), 1e-12 );
		assert_close( "[012] back to one variable, from y-stack", rss01, qr.calculateSkipColumnRSS( 2 ), 1e-12 );
		assert_eq( "[012] RSS unchanged", rss012, qr.calculateRSS() );

		assert_true( "[0123] push linearly independent", qr.pushColumn( x.columnVector( 3 ) ) );
		const double rss0123 = qr.calculateRSS();
		assert_close( "[0123] exact match", 0, rss0123 );
		assert_close( "[0123] Back to still exact match", 0, qr.calculateSkipColumnRSS( 0 ) );
		assert_close( "[0123] Back to inexact match", rss023, qr.calculateSkipColumnRSS( 1 ), 1e-12 );
		assert_close( "[0123] other back to still exact match", 0, qr.calculateSkipColumnRSS( 2 ) );
		assert_close( "[0123] back to one variable, from y-stack", rss012, qr.calculateSkipColumnRSS( 3 ), 1e-12 );
		assert_eq( "[0123] RSS unchanged", rss0123, qr.calculateRSS() );

		assert_true( "Can pop", qr.popColumn() );

		assert_true( "[0124] push linearly independent", qr.pushColumn( x.columnVector( 4 ) ) );
		const double rss0124 = qr.calculateRSS();
		assert_close( "[0124] exact match", 0, rss0124 );
		assert_close( "[0124] Back to still exact match", 0, qr.calculateSkipColumnRSS( 0 ) );
		assert_close( "[0124] again exact match", 0, qr.calculateSkipColumnRSS( 2 ) );
		assert_close( "[0124] back to one variable, from y-stack", rss012, qr.calculateSkipColumnRSS( 3 ), 1e-12 );
		assert_eq( "[0124] RSS unchanged", rss0124, qr.calculateRSS() );
	}

	void QRuncherTest::testSkipLastLinearlyIndependentColumn () {
		const double
			xData[] = {
				1,	0,
				0,	1
			},
			yData[] = { 3, 4 };

		AutoMatrix x( 2, 2 );
		x.fill( xData );
		AutoVector y( 2 );
		y.fill( yData );

		QRuncher qr( y );
		assert_eq( "[] RSS", 25.0, qr.calculateRSS() );
		assert_true( "[0] linearly independent", qr.pushColumn( x.columnVector( 0 ) ) );
		assert_eq( "[0] RSS", 16.0, qr.calculateRSS() );
		assert_eq( "[0]-0 RSS", 25.0, qr.calculateSkipColumnRSS( 0 ) );
		assert_true( "[01] push linearly independent", qr.pushColumn( x.columnVector( 1 ) ) );
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
		const double
			// col 0  norm 0, col 1 norm 2, col 2 norm 5, col 3 = col 2 - col 1, col 4 = col 2 + col 1
			xData[] = {
				0.0, 1.0, 1.0, 0.0, 2.0,
				0.0, 1.0, 2.0, 1.0, 3.0,
				0.0, 1.0, 2.0, 1.0, 3.0,
				0.0, 1.0, 4.0, 3.0, 5.0
			},
			// y norm 9
			yData[] = { 2.0, 4.0, 5.0, 6.0 };

		AutoMatrix x( 4, 5 );
		x.fill( xData );
		AutoVector y( 4 );
		y.fill( yData );

		// Calculate master comparison values by Gram-Schmidt algorithm
		const double
			squareSum1 = x.columnVector( 1 ).sumSquares(),
			squareSum2 = x.columnVector( 2 ).sumSquares(),
			innerProduct12 = x.columnVector( 1 ).innerProduct( x.columnVector( 2 ) );
		AutoVector a( 4 ), b( 4 );
		a.copy( x.columnVector( 1 ) ),	// for the orthognal part of x[*,1]
		b.copy( x.columnVector( 2 ) );	// for the orthognal part of x[*,2]
		a.axpy( - innerProduct12 / squareSum2, x.columnVector( 2 ) );
		b.axpy( - innerProduct12 / squareSum1, x.columnVector( 1 ) );
		const double
			squareSumA = a.sumSquares(),
			squareSumB = b.sumSquares(),
			innerProductY1 = y.innerProduct( x.columnVector( 1 ) ),
			mean = innerProductY1 / squareSum1;	// in case of column[1], that's the mean indeed
		AutoVector remainder( y );
		remainder.axpy( -mean, x.columnVector( 1 ) );
		const double
			rss1 = remainder.sumSquares(),
			innerProductY1B = remainder.innerProduct( b ),
			moment = innerProductY1B / squareSumB;
		remainder.axpy( -moment, b );
		const double
			rss12 = remainder.sumSquares(),
			innerProductY2 = y.innerProduct( x.columnVector( 2 ) ),
			restant12 = innerProductY2 / squareSum2;
		remainder.copy( y );
		remainder.axpy( -restant12, x.columnVector( 2 ) );
		const double
			rss2 = remainder.sumSquares(),
			innerProductY2A = remainder.innerProduct( a ),
			restant21 = innerProductY2A / squareSumA;
		remainder.axpy( -restant21, a );
		const double rss21 = remainder.sumSquares();
		assert_close( "Gram-Schmidt implementation test", rss12, rss21 );

		QRuncher qr( y );
		assert_eq( "[] RSS", 81.0, qr.calculateRSS() );
		assert_eq( "[] coefficients", 0, qr.calculateCoefficients().countDimensions() );

		// push column 0
		assert_false( "[0] zero vector is linearly dependent", qr.pushColumn( x.columnVector( 0 ) ) );
		assert_close( "[0] full remainder", 81.0, qr.calculateRSS() );
		const AutoVector coefficients0 = qr.calculateCoefficients();
		assert_eq( "[0] coefficients", 1, coefficients0.countDimensions() );
		assert_eq( "[0] no match attempt coefficients", 0.0, coefficients0.get( 0 ) );

		// push column 1
		assert_true( "[01] push linearly independent", qr.pushColumn( x.columnVector( 1 ) ) );
		assert_close( "[01] reduced remainder", rss1, qr.calculateRSS() );
		const AutoVector coefficients01 = qr.calculateCoefficients();
		assert_eq( "[01] coefficients", 2, coefficients01.countDimensions() );
		assert_eq( "[01] no match attempt coefficient 0", 0.0, coefficients01.get( 0 ) );
		assert_close( "[01] match attempt coefficient 1", mean, coefficients01.get( 1 ) );

		// push column 1 again
		assert_false( "[011] push linearly dependent", qr.pushColumn( x.columnVector( 1 ) ) );
		assert_close( "[011] reduced remainder remains", rss1, qr.calculateRSS() );
		const AutoVector coefficients011 = qr.calculateCoefficients();
		assert_eq( "[011] coefficients", 3, coefficients011.countDimensions() );
		assert_eq( "[011] no match attempt coefficient 0", 0.0, coefficients011.get( 0 ) );
		assert_close( "[011] match attempt coefficient 1", mean, coefficients011.get( 1 ) );
		assert_eq( "[011] no match attempt coefficient 2", 0.0, coefficients011.get( 2 ) );

		// push column 2
		assert_true( "[0112] push linearly independent", qr.pushColumn( x.columnVector( 2 ) ) );
		assert_close( "[0112] still no exact match", rss12, qr.calculateRSS() );
		const AutoVector coefficients0112 = qr.calculateCoefficients();
		assert_eq( "[0112] coefficients", 4, coefficients0112.countDimensions() );
		assert_eq( "[0112] no match attempt coefficient 0", 0.0, coefficients0112.get( 0 ) );
		assert_eq( "[0112] no match attempt coefficient 2", 0.0, coefficients0112.get( 2 ) );
		AutoMatrix r( 4, 4 );	// the regression matrix reflecting current push sequence
		r.columnVector( 0 ).copy( x.columnVector( 0 ) );
		r.columnVector( 1 ).copy( x.columnVector( 1 ) );
		r.columnVector( 2 ).copy( x.columnVector( 1 ) );
		r.columnVector( 3 ).copy( x.columnVector( 2 ) );
		AutoVector c( 4 );
		c.copy( y );
		c.gemv( -1.0, r, false, coefficients0112, 1.0 );	// remaincer c = y - r*coefficients
		// remainder is minimal if and only if it is orthogonal the spanning vectors
		assert_close( "[0112] remainder orthogonal 1", 0.0, c.innerProduct( x.columnVector( 1 ) ), 1e-12 );
		assert_close( "[0112] remainder orthogonal 2", 0.0, c.innerProduct( x.columnVector( 2 ) ), 1e-12 );

		// push column 2 again
		// Due to rounding error, pushColumn does cannot detect linear dependence:
		// assert_false( "[01122] push linearly dependent", qr.pushColumn( x.columnVector( 2 ) ) );

		// Also the following would fail,
		// because the push implies approximation by a 3-dimensional space
		// where one direction is virtually random
		// qr.pushColumn( x.columnVector( 2 ) );
		// assert_close( "[01122] still no exact match", rss12, qr.calculateRSS() );

		// This is the reason why a separate test testLimitedRounding() with no rounding necessary has been added.
	}

	void QRuncherTest::testLimitedRounding () {
		const double
			xData[] = {
				1.0, 0.0, 0.0, 0.0, 0.0, 0.0,
				0.0, 1.0, 4.0, 0.0, 12.0, 0.0,
				0.0, 0.0, 0.0, 4.0, 4.0, 0.0,
				0.0, 0.0, 0.0, 3.0, 3.0, 2.0
			},
			yData[] = { 84.0, 12.0, 4.0, 3.0 };	// Norm is 85

		AutoMatrix x( 4, 6 );
		x.fill( xData );
		AutoVector y( 4 );
		y.fill( yData );

		QRuncher qr( y );
		assert_eq( "[] rss", 85.0*85.0, qr.calculateRSS() );
		const AutoVector coefficients = qr.calculateCoefficients();
		assert_eq( "[] coefficients dim", 0, coefficients.countDimensions() );

		assert_true( "[0] push linearly independent", qr.pushColumn( x.columnVector( 0 ) ) );
		assert_eq( "[0] rss", 169.0, qr.calculateRSS() );
		assert_eq( "[0]-0 skip rss", 85.0*85.0, qr.calculateSkipColumnRSS( 0 ) );
		const AutoVector coefficients0 = qr.calculateCoefficients();
		assert_eq( "[0] coefficients dim", 1, coefficients0.countDimensions() );
		assert_eq( "[0] coefficient 0", 84.0, coefficients0.get( 0 ) );

		assert_true( "[01] push linearly independent", qr.pushColumn( x.columnVector( 1 ) ) );
		assert_eq( "[01] rss", 25.0, qr.calculateRSS() );
		assert_eq( "[01]-0 skip rss", 84.0*84.0 + 25.0, qr.calculateSkipColumnRSS( 0 ) );
		assert_eq( "[01]-1 skip rss", 169.0, qr.calculateSkipColumnRSS( 1 ) );
		const AutoVector coefficients01 = qr.calculateCoefficients();
		assert_eq( "[01] coefficients dim", 2, coefficients01.countDimensions() );
		assert_eq( "[01] coefficient 0", 84.0, coefficients01.get( 0 ) );
		assert_eq( "[01] coefficient 1", 12.0, coefficients01.get( 1 ) );

		assert_false( "[012] push linearly dependent", qr.pushColumn( x.columnVector( 2 ) ) );
		assert_eq( "[012] rss", 25.0, qr.calculateRSS() );
		assert_eq( "[012]-0 skip rss", 84.0*84.0 + 25.0, qr.calculateSkipColumnRSS( 0 ) );
		assert_eq( "[012]-1 skip rss, 2 replacing 1", 25.0, qr.calculateSkipColumnRSS( 1 ) );
		assert_eq( "[012]-2 skip rss, 2 superfluous", 25.0, qr.calculateSkipColumnRSS( 2 ) );
		const AutoVector coefficients012 = qr.calculateCoefficients();
		assert_eq( "[012] coefficients dim", 3, coefficients012.countDimensions() );
		assert_eq( "[012] coefficient 0", 84.0, coefficients012.get( 0 ) );
		assert_eq( "[012] coefficient 1", 12.0, coefficients012.get( 1 ) );
		assert_eq( "[012] coefficient 2", 0.0, coefficients012.get( 2 ) );

		assert_true( "[0123] push linearly independent", qr.pushColumn( x.columnVector( 3 ) ) );
		assert_eq( "[0123] rss", 0.0, qr.calculateRSS() );
		assert_eq( "[0123]-0 skip rss", 84.0*84.0, qr.calculateSkipColumnRSS( 0 ) );
		assert_eq( "[0123]-1 skip rss, 2 replacing 1", 0.0, qr.calculateSkipColumnRSS( 1 ) );
		assert_eq( "[0123]-2 skip rss, 2 superfluous", 0.0, qr.calculateSkipColumnRSS( 2 ) );
		assert_eq( "[0123]-3 skip rss", 25.0, qr.calculateSkipColumnRSS( 3 ) );
		const AutoVector coefficients0123 = qr.calculateCoefficients();
		assert_eq( "[0123] coefficients dim", 4, coefficients0123.countDimensions() );
		assert_eq( "[0123] coefficient 0", 84.0, coefficients0123.get( 0 ) );
		assert_eq( "[0123] coefficient 1", 12.0, coefficients0123.get( 1 ) );
		assert_eq( "[0123] coefficient 2", 0.0, coefficients0123.get( 2 ) );
		assert_eq( "[0123] coefficient 3", 1.0, coefficients0123.get( 3 ) );

		assert_false( "[01234] push linearly dependent", qr.pushColumn( x.columnVector( 4 ) ) );
		assert_eq( "[01234] rss", 0.0, qr.calculateRSS() );
		assert_eq( "[01234]-0 skip rss", 84.0*84.0, qr.calculateSkipColumnRSS( 0 ) );
		assert_eq( "[01234]-1 skip rss, 2 replacing 1", 0.0, qr.calculateSkipColumnRSS( 1 ) );
		assert_eq( "[01234]-2 skip rss, 2 superfluous", 0.0, qr.calculateSkipColumnRSS( 2 ) );
		assert_eq( "[01234]-3 skip rss, 4+1 replacing 3", 0.0, qr.calculateSkipColumnRSS( 3 ) );
		assert_eq( "[01234]-4 skip rss, 4 superfluous", 0.0, qr.calculateSkipColumnRSS( 4 ) );
		const AutoVector coefficients01234 = qr.calculateCoefficients();
		assert_eq( "[01234] coefficients dim", 5, coefficients01234.countDimensions() );
		assert_eq( "[01234] coefficient 0", 84.0, coefficients01234.get( 0 ) );
		assert_eq( "[01234] coefficient 1", 12.0, coefficients01234.get( 1 ) );
		assert_eq( "[01234] coefficient 2", 0.0, coefficients01234.get( 2 ) );
		assert_eq( "[01234] coefficient 3", 1.0, coefficients01234.get( 3 ) );
		assert_eq( "[01234] coefficient 4", 0.0, coefficients01234.get( 4 ) );

		assert_true( "[0123]b pop", qr.popColumn() );
		assert_eq( "[0123]b rss", 0.0, qr.calculateRSS() );
		assert_eq( "[0123]b coefficients", coefficients0123, qr.calculateCoefficients() );

		assert_true( "[012]b pop", qr.popColumn() );
		assert_eq( "[012]b rss", 25.0, qr.calculateRSS() );
		assert_eq( "[012]b coefficients", coefficients012, qr.calculateCoefficients() );

		assert_true( "[01]b pop", qr.popColumn() );
		assert_eq( "[01]b rss", 25.0, qr.calculateRSS() );
		assert_eq( "[01]b coefficients", coefficients01, qr.calculateCoefficients() );

		assert_true( "[0]b pop", qr.popColumn() );
		assert_eq( "[0]b rss", 169.0, qr.calculateRSS() );
		assert_eq( "[0]b coefficients", coefficients0, qr.calculateCoefficients() );

		assert_true( "[04] push linearly independent", qr.pushColumn( x.columnVector( 4 ) ) );
		assert_eq( "[04] rss", 0.0, qr.calculateRSS() );
		assert_eq( "[04]-0 skip rss", 84.0*84.0, qr.calculateSkipColumnRSS( 0 ) );
		assert_eq( "[04]-4 skip rss", 169.0, qr.calculateSkipColumnRSS( 1 ) );
		const AutoVector coefficients04 = qr.calculateCoefficients();
		assert_eq( "[04] coefficients dim", 2, coefficients04.countDimensions() );
		assert_eq( "[04] coefficient 0", 84.0, coefficients04.get( 0 ) );
		assert_eq( "[04] coefficient 4", 1.0, coefficients04.get( 1 ) );

		assert_true( "[043] push linearly independent", qr.pushColumn( x.columnVector( 3 ) ) );
		assert_eq( "[043] rss", 0.0, qr.calculateRSS() );
		assert_eq( "[043]-0 skip rss", 84.0*84.0, qr.calculateSkipColumnRSS( 0 ) );
		assert_close( "[043]-4 skip rss", 144.0, qr.calculateSkipColumnRSS( 1 ), 1e-12 );
		assert_eq( "[043]-3 skip rss", 0.0, qr.calculateSkipColumnRSS( 2 ) );
		const AutoVector coefficients043 = qr.calculateCoefficients();
		assert_eq( "[043] coefficients dim", 3, coefficients043.countDimensions() );
		assert_eq( "[043] coefficient 0", 84.0, coefficients043.get( 0 ) );
		assert_eq( "[043] coefficient 4", 1.0, coefficients043.get( 1 ) );
		assert_eq( "[043] coefficient 3", 0.0, coefficients043.get( 2 ) );

		assert_false( "[0432] push linearly dependent", qr.pushColumn( x.columnVector( 4 ) ) );
		assert_eq( "[0432] rss", 0.0, qr.calculateRSS() );
		assert_eq( "[0432]-0 skip rss", 84.0*84.0, qr.calculateSkipColumnRSS( 0 ) );
		assert_eq( "[0432]-4 skip rss, 3+2 replacing 4", 0.0, qr.calculateSkipColumnRSS( 1 ) );
		assert_eq( "[0432]-3 skip rss, 3 not needed for y match", 0.0, qr.calculateSkipColumnRSS( 2 ) );
		assert_eq( "[0432]-2 skip rss, 2 superfluous", 0.0, qr.calculateSkipColumnRSS( 3 ) );
		const AutoVector coefficients0432 = qr.calculateCoefficients();
		assert_eq( "[0432] coefficients dim", 4, coefficients0432.countDimensions() );
		assert_eq( "[0432] coefficient 0", 84.0, coefficients0432.get( 0 ) );
		assert_eq( "[0432] coefficient 4", 1.0, coefficients0432.get( 1 ) );
		assert_eq( "[0432] coefficient 3", 0.0, coefficients0432.get( 2 ) );
		assert_eq( "[0432] coefficient 2", 0.0, coefficients0432.get( 3 ) );

		assert_true( "[04325] push linearly independent", qr.pushColumn( x.columnVector( 5 ) ) );
		assert_eq( "[04325] rss", 0.0, qr.calculateRSS() );
		assert_eq( "[04325]-0 skip rss", 84.0*84.0, qr.calculateSkipColumnRSS( 0 ) );
		assert_eq( "[04325]-4 skip rss, 3+2 replacing 4", 0.0, qr.calculateSkipColumnRSS( 1 ) );
		assert_eq( "[04325]-3 skip rss, 3 not needed for y match", 0.0, qr.calculateSkipColumnRSS( 2 ) );
		assert_eq( "[04325]-2 skip rss, 2 superfluous", 0.0, qr.calculateSkipColumnRSS( 3 ) );
		assert_eq( "[04325]-5 skip rss, 5 not needed for y match", 0.0, qr.calculateSkipColumnRSS( 4 ) );
		const AutoVector coefficients04325 = qr.calculateCoefficients();
		assert_eq( "[04325] coefficients dim", 5, coefficients04325.countDimensions() );
		assert_eq( "[04325] coefficient 0", 84.0, coefficients04325.get( 0 ) );
		assert_eq( "[04325] coefficient 4", 1.0, coefficients04325.get( 1 ) );
		assert_eq( "[04325] coefficient 3", 0.0, coefficients04325.get( 2 ) );
		assert_eq( "[04325] coefficient 2", 0.0, coefficients04325.get( 3 ) );
		assert_eq( "[04325] coefficient 5", 0.0, coefficients04325.get( 4 ) );
	}

}
