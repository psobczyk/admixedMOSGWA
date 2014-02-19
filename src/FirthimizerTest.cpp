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
#include "linalg/package.hpp"
#include "TestSuite.hpp"
#include <math.h>
#include <stdio.h>

using namespace std;
using namespace linalg;
using namespace minimization;
using namespace unitpp;

namespace test {

	/** Trivial subclass to make protected methods of {@link Firthimizer} public for testing. */
	struct TestFirthimizer : public Firthimizer {

		/** Construct an instance. */
		TestFirthimizer ( const Vector& yVec );

		/** Publicly visible {@link Firthimizer::scoreTest} for testing. */
		using Firthimizer::setCoefficients;

		/** Publicly visible {@link Firthimizer::countDimensions} for testing. */
		using Firthimizer::countDimensions;

		/** Publicly visible {@link Firthimizer::calculateFunction} for testing. */
		using Firthimizer::calculateFunction;

		/** Publicly visible {@link Firthimizer::calculateDerivative} for testing. */
		using Firthimizer::calculateDerivative;

		/** Publicly visible {@link Firthimizer::calculateFunctionAndDerivative} for testing. */
		using Firthimizer::calculateFunctionAndDerivative;
	};

	TestFirthimizer::TestFirthimizer ( const linalg::Vector& yVec ) : Firthimizer( yVec ) {}

	/** Tests the class {@link Firthimizer}. */
	class FirthimizerTest : public TestSuite {

		/** Test {@link minimization::Minimizable} interface methods. */
		void testMinimizable ();

		/** Test {@link Firthimizer::calculateDerivative}.
		* This is also done in {@link FirthimizerTest::testMinimizable},
		* but here a large regression matrix with random numbers ist tested.
		*/
		void testDerivative ();

		/** Test {@link Firthimizer::calculateLogLikelihood} and {@link Firthimizer::calculateCoefficients}. */
		void testRegression ();

		/** Test {@link Firthimizer::scoreTest}. */
		void testScoreTest ();

		public:

		FirthimizerTest () : TestSuite( "FirthimizerTest" ) {
			addTestMethod( "FirthimizerTest::testMinimizable", this, &FirthimizerTest::testMinimizable );
			addTestMethod( "FirthimizerTest::testDerivative", this, &FirthimizerTest::testDerivative );
			addTestMethod( "FirthimizerTest::testRegression", this, &FirthimizerTest::testRegression );
			addTestMethod( "FirthimizerTest::testScoreTest", this, &FirthimizerTest::testScoreTest );
		}
	} * firthimizerTest = new FirthimizerTest();	// automatically freed by unit++

	void FirthimizerTest::testMinimizable () {
		// Test zero-dim Firthimizer
		{
			const AutoVector y( 0 );
			TestFirthimizer firthimizer( y );
			assert_eq( "0-dim count", 0, firthimizer.countDimensions() );

			const AutoVector coefficients( 0 );
			double f = nan( "undefined" );
			AutoVector df( 0 );
			firthimizer.calculateFunction( coefficients, f );
			assert_eq( "0-dim total log probability", 0.0, -f );
			firthimizer.calculateDerivative( coefficients, df );

			f = nan( "undefined" );
			firthimizer.calculateFunctionAndDerivative( coefficients, f, df );
			assert_eq( "0-dim total log probability combined", 0.0, -f );
		}

		// Test one-dim Firthimizer
		{
			AutoVector y( 1 );
			y.set( 0, 1.0 );

			TestFirthimizer firthimizer( y );
			assert_eq( "1-dim-0 count", 0, firthimizer.countDimensions() );

			AutoVector coefficients( 0 );
			double f = nan( "undefined" );
			AutoVector df( 0 );

			firthimizer.calculateFunction( coefficients, f );
			assert_eq( "1-dim-0 p=1/2 log probability", log( 0.5 ), -f );
			firthimizer.calculateDerivative( coefficients, df );

			f = nan( "undefined" );
			firthimizer.calculateFunctionAndDerivative( coefficients, f, df );
			assert_eq( "1-dim-0 p=1/2 log probability combined", log( 0.5 ), -f );

			// Test with one-dim coefficient space at beta=0
			AutoVector x( 1 );
			x.set( 0, 1.0 );
			firthimizer.insertColumn( 0, x, 0.0 );
			assert_eq( "1-dim-1 count", 1, firthimizer.countDimensions() );

			coefficients.upSize( 1 );
			coefficients.set( 0, 0.0 );

			f = nan( "undefined" );
			df.upSize( 1 );
			df.set( 0, nan( "undefined" ) );
			firthimizer.calculateFunction( coefficients, f );
			assert_eq( "1-dim-1 p=1/2 log probability", log( 0.25 ), -f );
			firthimizer.calculateDerivative( coefficients, df );
			assert_eq( "1-dim-1 p=1/2 log probability derivative", 0.5, -df.get( 0 ) );

			f = nan( "undefined" );
			df.set( 0, nan( "undefined" ) );
			firthimizer.calculateFunctionAndDerivative( coefficients, f, df );
			assert_eq( "1-dim-1 p=1/2 log probability combined", log( 0.25 ), -f );
			assert_eq( "1-dim-1 p=1/2 log probability derivative", 0.5, -df.get( 0 ) );

			const double
				delta = 0.0001,
				epsilon = 0.000001;
			double fMinusDelta, fPlusDelta;
			coefficients.set( 0, -delta );
			firthimizer.calculateFunction( coefficients, fMinusDelta );
			coefficients.set( 0, +delta );
			firthimizer.calculateFunction( coefficients, fPlusDelta );
			assert_true( "1-dim-1 numerical differentiation", fabs( df.get( 0 ) - ( fPlusDelta - fMinusDelta ) / ( 2.0 * delta ) ) < epsilon );

			coefficients.set( 0, 10.0 );
			firthimizer.calculateDerivative( coefficients, df );
			coefficients.set( 0, 10.0-delta );
			firthimizer.calculateFunction( coefficients, fMinusDelta );
			coefficients.set( 0, 10.0+delta );
			firthimizer.calculateFunction( coefficients, fPlusDelta );
			assert_true( "1-dim-1 numerical differentiation at 10", fabs( df.get( 0 ) - ( fPlusDelta - fMinusDelta ) / ( 2.0 * delta ) ) < epsilon );

			// Test with one-dim coefficient space at beta close to infinity
			coefficients.set( 0, 1000.0 );

			f = nan( "undefined" );
			df.set( 0, nan( "undefined" ) );

			firthimizer.calculateFunctionAndDerivative( coefficients, f, df );
			assert_true( "1-dim-1 p=1 log probability combined infinite", isinf( -f ) );
			assert_true( "1-dim-1 p=1 log probability combined negative", -f < 0 );
			assert_eq( "1-dim-1 p=1 log probability derivative", -0.5, -df.get( 0 ) );

			// Test removeColumn
			firthimizer.removeColumn( 0 );
			coefficients.upSize( 0 );
			df.upSize( 0 );
			firthimizer.calculateFunction( coefficients, f );
			assert_eq( "1-dim-0 p=1/2 log probability after remove", log( 0.5 ), -f );
			firthimizer.calculateDerivative( coefficients, df );
			f = nan( "undefined" );
			firthimizer.calculateFunctionAndDerivative( coefficients, f, df );
			assert_eq( "1-dim-0 p=1/2 log probability combined after remove", log( 0.5 ), -f );
		}

		// Numerically test 3-dim-2 Firthimizer.
		{
			AutoVector y( 3 );
			y.set( 0, 0.0 );
			y.set( 1, 1.0 );
			y.set( 2, 1.0 );

			AutoMatrix x( 3, 2 );
			x.set( 0, 0, 1.0 );
			x.set( 1, 0, 1.0 );
			x.set( 2, 0, 0.0 );
			x.set( 0, 1, -1.0 );
			x.set( 1, 1, 0.0 );
			x.set( 2, 1, 1.0 );

			TestFirthimizer firthimizer( y );
			assert_eq( "3-dim-0 count", 0, firthimizer.countDimensions() );

			firthimizer.insertColumn( 0, x.columnVector( 0 ) );
			assert_eq( "3-dim-1 count", 1, firthimizer.countDimensions() );

			firthimizer.insertColumn( 1, x.columnVector( 1 ) );
			assert_eq( "3-dim-2 count", 2, firthimizer.countDimensions() );

			double f = nan( "undefined" );
			AutoVector df( 2 );
			df.set( 0, f );
			df.set( 1, f );

			const double
				delta = 0.0001,
				epsilon = 0.000001,
				log2 = log( 2.0 );
			double
				fMinusDelta0, fPlusDelta0,
				fMinusDelta1, fPlusDelta1;

			AutoVector coefficients( 2 );

			coefficients.set( 1, -log2 );
			coefficients.set( 0, log2-delta );
			firthimizer.calculateFunction( coefficients, fMinusDelta0 );
			coefficients.set( 0, log2+delta );
			firthimizer.calculateFunction( coefficients, fPlusDelta0 );
			coefficients.set( 0, log2 );

			coefficients.set( 1, -log2-delta );
			firthimizer.calculateFunction( coefficients, fMinusDelta1 );
			coefficients.set( 1, -log2+delta );
			firthimizer.calculateFunction( coefficients, fPlusDelta1 );
			coefficients.set( 1, -log2 );

			firthimizer.calculateFunctionAndDerivative( coefficients, f, df );

			assert_true( "3-dim-2 numerical differentiation 0", fabs( df.get( 0 ) - ( fPlusDelta0 - fMinusDelta0 ) / ( 2.0 * delta ) ) < epsilon );
			assert_true( "3-dim-2 numerical differentiation 1", fabs( df.get( 1 ) - ( fPlusDelta1 - fMinusDelta1 ) / ( 2.0 * delta ) ) < epsilon );

			// Test replaceColumn by exchanging columns.
			const AutoVector coefficients0( coefficients );
			const double f0 = f;
			const AutoVector df0( df );
			firthimizer.replaceColumn( 0, x.columnVector( 1 ) );
			firthimizer.replaceColumn( 1, x.columnVector( 0 ) );
			coefficients.set( 0, -log2 );
			coefficients.set( 1, log2 );

			firthimizer.calculateFunctionAndDerivative( coefficients, f, df );

			assert_true( "Numerical after exchange same min", fabs( f0 - f ) < epsilon );
			assert_true( "Numerical after exchange df[0]", fabs( df0.get( 0 ) - df.get( 1 ) ) < epsilon );
			assert_true( "Numerical after exchange df[1]", fabs( df0.get( 1 ) - df.get( 0 ) ) < epsilon );
		}
	}

	void FirthimizerTest::testDerivative () {
		const int
			columns = 100,
			rows = 1000;
		const double
			delta = 0.0001,
			epsilon = 0.0001;
		srand48( -2083648846L );		// some seed value to make it reproducible
		AutoMatrix x( rows, columns );
		AutoVector
			y( rows ),
			beta( columns ),
			betaDelta( columns ),
			gradient( columns ),
			gradientApprox( columns );
		TestFirthimizer firthimizer( y );

		for ( int row = 0; row < rows; ++row ) {
			for ( int column = 0; column < columns; ++column ) {
				x.set( row, column, 0.2 * drand48() );
			}
			y.set( row, 0.6 + 0.0 * drand48() );	// eigenartig: kaum Einfluss
		}
		for ( int column = 0; column < columns; ++column ) {
			beta.set( column, 1.0 - 0.01 * column + 0.0 * drand48() );	// eigenartig: kaum Einfluss
			firthimizer.insertColumn( column, x.columnVector( column ) );
		}

		firthimizer.calculateDerivative( beta, gradient );

		for ( int column = 0; column < columns; ++column ) {
			double fMinus, fPlus;
			betaDelta.copy( beta );
			betaDelta.set( column, beta.get( column ) - delta );
			firthimizer.calculateFunction( betaDelta, fMinus );
			betaDelta.set( column, beta.get( column ) + delta );
			firthimizer.calculateFunction( betaDelta, fPlus );
			gradientApprox.set( column, ( fPlus - fMinus ) / ( 2.0 * delta ) );
		}

		printf( "Gradient:\n" );
		printVector( gradient );
		printf( "GradientApprox:\n" );
		printVector( gradientApprox );

		for ( int column = 0; column < columns; ++column ) {
			char buffer[255];
			sprintf( buffer, "Gradient[%u]", column );
			assert_close(
				buffer,
				gradientApprox.get( column ),
				gradient.get( column ),
				1e-6 * fabs( gradientApprox.get( column ) )
			);
		}
	}

	void FirthimizerTest::testRegression () {
		Minimizer minimizer( 0.01, 0.0001, 100 );

		// Test zero-dim regression
		{
			const AutoVector y( 0 );
			Firthimizer firthimizer( y );
			assert_eq( "Zero remainder in 0-dim space", 0.0, firthimizer.calculateLogLikelihood( minimizer ) );
			assert_eq( "Zero-dim coefficients in 0-dim space", 0, firthimizer.calculateCoefficients( minimizer ).countDimensions() );
		}

		// Test one-dim regression
		{
			AutoVector y( 1 );
			y.set( 0, 1.0 );

			Firthimizer firthimizer( y );
			assert_eq( "Full remainder in 1-dim space", log( 0.5 ), firthimizer.calculateLogLikelihood( minimizer ) );
			assert_eq( "Zero-dim coefficients in 1-dim space", 0, firthimizer.calculateCoefficients( minimizer ).countDimensions() );

			// Test with one-dim coefficient space at beta=0
			AutoVector x( 1 );
			x.set( 0, 1.0 );
			firthimizer.insertColumn( 0, x, 0.0 );
			// See calculation of exact solution in Bernhard's Denkbuch 2010, S. 29-30 April: p=3/4.
			assert_true(
				"Remainder in 1-dim space",
				fabs (
					0.5 * log( 0.75*0.75*0.75*(1.0-0.75) )
					-
					firthimizer.calculateLogLikelihood( minimizer )
				) < 0.0001
			);
			const Vector& coefficients( firthimizer.calculateCoefficients( minimizer ) );
			assert_eq( "1-dim coefficients in 1-dim space", 1, coefficients.countDimensions() );
			assert_true(
				"Coefficient in 1-dim space",
				fabs (
					log( 3.0 )
					-
					coefficients.get( 0 )
				) < 0.0001
			);
		}
	}

	void FirthimizerTest::testScoreTest () {
		// Test zero-dim Firthimizer
		{
			const AutoVector y( 0 );
			const AutoVector beta( 0 );
			const AutoVector z( 0 );

			TestFirthimizer firthimizer( y );
			firthimizer.setCoefficients( beta );
			assert_true( "Total regression in 0 dim space", isnan( firthimizer.scoreTest( z ) ) );
		}

		// Test empty Firthimizer
		{
			AutoVector y( 1 );
			y.set( 0, 1.0 );
			const AutoVector beta( 0 );
			AutoVector z( 1 );
			z.set( 0, 1.0 );

			TestFirthimizer firthimizer( y );
			firthimizer.setCoefficients( beta );
			assert_eq( "Initially full gain in 1 dim space", 1.0, firthimizer.scoreTest( z ) );
		}

		// Test one-column Firthimizer in 2-dim space
		{
			AutoMatrix x( 2, 2 );
			x.set( 0, 0, 1.0 );
			x.set( 1, 0, 0.0 );
			x.set( 0, 1, 0.0 );
			x.set( 1, 1, 1.0 );
			const Vector
				y = x.columnVector( 0 ),
				z = x.columnVector( 1 );
			const AutoVector betaEmpty( 0 );
			AutoVector
				betaZero( 1 ),
				betaBig( 1 );
			betaZero.set( 0, 0.0 );
			betaBig.set( 0, 1000.0 );	// exp( 1000 ) is roughly infinite.

			TestFirthimizer firthimizer( y );

			firthimizer.setCoefficients( betaEmpty );
			assert_eq( "Initially full gain in 2 dim space", 1.0, firthimizer.scoreTest( z ) );

			firthimizer.insertColumn( 0, y );
			firthimizer.setCoefficients( betaBig ); // Fit first component
			// Calculated the following by hand [Bernhard's Denkbuch anno 2010, p. 23.Jan.]
			assert_eq(
				"After fitting 1st dim, remaining gain from 2nd dim in 2 dim space",
				1.0,
				firthimizer.scoreTest( z )
			);
			assert_true(
				"After fitting push linearly dependent in 2 dim space",
				isnan( firthimizer.scoreTest( y ) )
			);

			firthimizer.removeColumn( 0 );

			firthimizer.insertColumn( 0, z );
			firthimizer.setCoefficients( betaZero );
			assert_eq( "After non-fitting full gain in 2 dim space", 1.0, firthimizer.scoreTest( y ) );
			assert_true(
				"After non-fitting push linearly dependent in 2 dim space",
				isnan( firthimizer.scoreTest( z ) )
			);

			firthimizer.removeColumn( 0 );

			firthimizer.setCoefficients( betaEmpty );
			assert_eq( "After pop full gain in 2 dim space", 1.0, firthimizer.scoreTest( z ) );
		}

		// Test one-column Firthimizer in 3-dim space
		{
			const double data[] = {
				1, 0, 0,
				0, 1, 1,
				0, 0, 1
			};
			AutoMatrix x( 3, 3 );
			x.fill( data );
			const Vector
				y = x.columnVector( 1 ),
				z = x.columnVector( 2 );
			const AutoVector betaEmpty( 0 );
			AutoVector
				betaZero( 1 ),
				betaSmall( 1 );
			betaZero.set( 0, 0.0 );
			betaSmall.set( 0, -1000.0 );	// exp( -1000 ) is roughly zero.

			TestFirthimizer firthimizer( y );

			firthimizer.setCoefficients( betaEmpty );
			assert_eq( "Initially full gain in 3 dim space", 1.0, firthimizer.scoreTest( y ) );
			assert_eq( "Initially useless gain in 3 dim space", 0.0, firthimizer.scoreTest( z ) );

			firthimizer.insertColumn( 0, x.columnVector( 0 ) );
			firthimizer.setCoefficients( betaSmall );
			assert_eq(
				"After fitting 1st dim, remaining gain from 2nd dim in 3 dim space",
				1.0,
				firthimizer.scoreTest( y )
			);
			assert_eq(
				"After fitting 1st dim, no gain from useless 2nd dim in 3 dim space",
				0.0,
				firthimizer.scoreTest( z )
			);
		}
	}
}
