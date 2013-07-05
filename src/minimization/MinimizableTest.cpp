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

#include "TestMinimizable.hpp"
#include "../linalg/AutoVector.hpp"
#include "../TestSuite.hpp"

using namespace unitpp;
using namespace linalg;
using namespace minimization;

namespace test {

	/** Tests the class {@link minimization::Minimizer} and thereby {@link minimization::Minimizable}. */
	class MinimizableTest : public TestSuite {

		void test0dim ();
		void test1dim ();
		void test2dim ();

		public:

		MinimizableTest () : TestSuite( "minimization::MinimizableTest" ) {
			addTestMethod( "MinimizableTest::test0dim", this, &MinimizableTest::test0dim );
			addTestMethod( "MinimizableTest::test1dim", this, &MinimizableTest::test1dim );
			addTestMethod( "MinimizableTest::test2dim", this, &MinimizableTest::test2dim );
		}
	} * minimizableTest = new MinimizableTest();	// automatically freed by unit++

	/** Test {@link Minimizable::calculateFunctionAndDerivative} in 0-dimensional case. */
	void MinimizableTest::test0dim () {
		const AutoVector minVec( 0 );
		TestMinimizable minimizable( minVec );
		const AutoVector x( 0 );
		double f;
		AutoVector df( 0 );
		minimizable.calculateFunction( x, f );
		assert_eq( "function value", 0.0, f );

		minimizable.calculateDerivative( x, df );
		assert_eq( "derivative dimension", 0, df.countDimensions() );

		minimizable.calculateFunctionAndDerivative( x, f, df );
		assert_eq( "again function value", 0.0, f );
		assert_eq( "again derivative dimension", 0, df.countDimensions() );
	}

	/** Test {@link Minimizable::calculateFunctionAndDerivative} in 1-dimensional case. */
	void MinimizableTest::test1dim () {
		AutoVector minVec( 1 );
		minVec.set( 0, 0.0 );
		TestMinimizable minimizable( minVec );
		AutoVector x( 1 );
		double f;
		AutoVector df( 1 );

		x.set( 0, 1.0 );
		minimizable.calculateFunction( x, f );
		assert_eq( "function value", 1.0, f );

		minimizable.calculateDerivative( x, df );
		assert_eq( "derivative", 2.0, df.get( 0 ) );

		minimizable.calculateFunctionAndDerivative( x, f, df );
		assert_eq( "again function value", 1.0, f );
		assert_eq( "again derivative", 2.0, df.get( 0 ) );
	}

	/** Test {@link Minimizable::calculateFunctionAndDerivative} in 2-dimensional case. */
	void MinimizableTest::test2dim () {
		AutoVector minVec( 2 );
		minVec.set( 0, 0.0 );
		minVec.set( 1, -1.0 );

		TestMinimizable minimizable( minVec );
		AutoVector x( 2 );
		double f;
		AutoVector df( 2 );

		x.set( 0, 1.0 );
		x.set( 1, 0.0 );
		minimizable.calculateFunction( x, f );
		assert_eq( "function value", 2.0, f );

		minimizable.calculateDerivative( x, df );
		assert_eq( "derivative[0]", 2.0, df.get( 0 ) );
		assert_eq( "derivative[1]", 2.0, df.get( 1 ) );

		minimizable.calculateFunctionAndDerivative( x, f, df );
		assert_eq( "again function value", 2.0, f );
		assert_eq( "again derivative[0]", 2.0, df.get( 0 ) );
		assert_eq( "again derivative[1]", 2.0, df.get( 1 ) );
	}

}
