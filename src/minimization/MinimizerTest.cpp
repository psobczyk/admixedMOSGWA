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

#include "Minimizer.hpp"
#include "TestMinimizable.hpp"
#include "../linalg/AutoVector.hpp"
#include "../TestSuite.hpp"
#include <math.h>

using namespace unitpp;
using namespace linalg;
using namespace minimization;

namespace test {

	/** Tests the class {@link minimization::Minimizer} and thereby {@link minimization::Minimizable}. */
	class MinimizerTest : public TestSuite {

		void test0dim ();
		void test1dim ();
		void test2dim ();
		void test0iterations ();
		void test1iteration ();
		void test2iterations ();

		public:

		MinimizerTest () : TestSuite( "minimization::MinimizerTest" ) {
			addTestMethod( "MinimizerTest::test0dim", this, &MinimizerTest::test0dim );
			addTestMethod( "MinimizerTest::test1dim", this, &MinimizerTest::test1dim );
			addTestMethod( "MinimizerTest::test2dim", this, &MinimizerTest::test2dim );
			addTestMethod( "MinimizerTest::test0iterations", this, &MinimizerTest::test0iterations );
			addTestMethod( "MinimizerTest::test1iteration", this, &MinimizerTest::test1iteration );
			addTestMethod( "MinimizerTest::test2iterations", this, &MinimizerTest::test2iterations );
		}
	} * minimizerTest = new MinimizerTest();	// automatically freed by unit++

	/** Test {@link Minimizer::minimize} in 0-dimensional space. */
	void MinimizerTest::test0dim () {
		const AutoVector minVec( 0 );
		TestMinimizable minimizable( minVec );
		AutoVector x( 0 );

		Minimizer minimizer( 0.01, 0.0001, 100 );
		const double minimum = minimizer.minimize( minimizable, x );
		assert_eq( "minimum value", 0.0, minimum );
	}

	/** Test {@link Minimizer::minimize} in 1-dimensional space. */
	void MinimizerTest::test1dim () {
		AutoVector minVec( 1 );
		minVec.set( 0, 2.0 );
		TestMinimizable minimizable( minVec );

		AutoVector x( 1 );
		x.set( 0, 1.0 );

		Minimizer minimizer( 0.01, 0.0001, 100 );
		const double minimum = minimizer.minimize( minimizable, x );
		assert_close( "minimum value", 0.0, minimum, 1e-10 );
		assert_close( "minimum x[0]", 2.0, x.get( 0 ), 1e-10 );
	}

	/** Test {@link Minimizer::minimize} in 2-dimensional space. */
	void MinimizerTest::test2dim () {
		AutoVector minVec( 2 );
		minVec.set( 0, -1.0 );
		minVec.set( 1, -5.0 );
		TestMinimizable minimizable( minVec );

		AutoVector x( 2 );
		x.set( 0, 200.0 );
		x.set( 1, -160.0 );

		Minimizer minimizer( 0.01, 0.0001, 100 );
		const double minimum = minimizer.minimize( minimizable, x );
		assert_close( "minimum value", 0.0, minimum, 1e-10 );
		assert_close( "minimum x[0]", -1.0, x.get( 0 ), 1e-10 );
		assert_close( "minimum x[1]", -5.0, x.get( 1 ), 1e-10 );
	}

	/** Test {@link Minimizer::minimize} with zero iterations. */
	void MinimizerTest::test0iterations () {
		Minimizer minimizer( 0.01, 0.0001, 0 );
		AutoVector minVec( 0 );
		AutoVector x( 0 );

		{
			TestMinimizable minimizable( minVec );
			const double minimum = minimizer.minimize( minimizable, x );
			assert_eq( "no iterations necessary in 0-dim", 0.0, minimum );
		}

		minVec.upSize( 1 );
		x.upSize( 1 );

		{
			TestMinimizable minimizable( minVec );
			x.set( 0, 0.0 );
			const double minFail = minimizer.minimize( minimizable, x );
			assert_true( "zero iterations insufficient in 1-dim", isnan( minFail ) );
			assert_eq( "no iterations endpoint in 1-dim is startpoint", 0.0, x.get( 0 ) );
		}
	}

	/** Test {@link Minimizer::minimize} with one iteration. */
	void MinimizerTest::test1iteration () {
		Minimizer minimizer( 0.01, 0.0001, 1 );
		AutoVector minVec( 0 );
		AutoVector x( 0 );

		{
			TestMinimizable minimizable( minVec );
			const double minimum = minimizer.minimize( minimizable, x );
			assert_eq( "one iteration not bad in 0-dim", 0.0, minimum );
		}

		minVec.upSize( 1 );
		minVec.set( 0, 1.0 );
		x.upSize( 1 );

		{
			TestMinimizable minimizable( minVec );
			x.set( 0, 0.0 );
			const double minimum = minimizer.minimize( minimizable, x );
			assert_eq( "one iteration with start=end in 1-dim", 0.0, minimum );
			assert_eq( "one iteration with start=end in 1-dim finds point", 1.0, x.get( 0 ) );
		}
	}

	/** Test {@link Minimizer::minimize} with two iterations. */
	void MinimizerTest::test2iterations () {
		Minimizer minimizer( 0.01, 0.0001, 2 );
		AutoVector minVec( 1 );
		AutoVector x( 1 );
		x.set( 0, 1.0 );

		{
			TestMinimizable minimizable( minVec );
			x.set( 0, 0.01 );
			const double minimum = minimizer.minimize( minimizable, x );
			assert_close( "two iterations with start!=end in 1-dim", 0.0, minimum, 1e-10 );
			assert_close( "two iterations with start!=end in 1-dim finds point", 0.0, x.get( 0 ), 1e-10 );
		}
	}

}
