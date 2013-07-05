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

#include "AutoPermutation.hpp"
#include "../TestSuite.hpp"
#include <math.h>

using namespace std;
using namespace unitpp;
using namespace linalg;

namespace test {

	/** Tests the class {@link linalg::AutoPermutation} and thereby {@link linalg::Permutation}. */
	class AutoPermutationTest : public TestSuite {

		void testCopyConstructor ();
		void testAssignment ();
		void testOperators ();
		void testFillGetSwap ();
		void testResize ();
		void testUpSize ();

		public:

		AutoPermutationTest () : TestSuite( "linalg::AutoPermutationTest" ) {
			addTestMethod( "AutoPermutationTest::testCopyConstructor", this, &AutoPermutationTest::testCopyConstructor );
			addTestMethod( "AutoPermutationTest::testAssignment", this, &AutoPermutationTest::testAssignment );
			addTestMethod( "AutoPermutationTest::testOperators", this, &AutoPermutationTest::testOperators );
			addTestMethod( "AutoPermutationTest::testFillGetSwap", this, &AutoPermutationTest::testFillGetSwap );
			addTestMethod( "AutoPermutationTest::testResize", this, &AutoPermutationTest::testResize );
			addTestMethod( "AutoPermutationTest::testUpSize", this, &AutoPermutationTest::testUpSize );
		}
	} * permutationTest = new AutoPermutationTest();	// automatically freed by unit++

	/** Test {@link AutoPermutation::AutoPermutation( AutoPermutation& )}. */
	void AutoPermutationTest::testCopyConstructor () {
		const size_t
			data1[] = { 1, 2, 3, 4, 5, 6, 7, 0 },
			data2[] = { 7, 0, 1, 2, 3, 4, 5, 6 };
		AutoPermutation p( 8, 10 );
		p.fill( data1 );

		AutoPermutation q( p );
		assert_eq( "dimension", 8, q.countDimensions() );
		assert_eq( "q[0]", 1, q.get( 0 ) );
		assert_eq( "q[3]", 4, q.get( 3 ) );

		p.fill( data2 );
		assert_eq( "detached q[1]", 2, q.get( 1 ) );
		assert_eq( "detached q[4]", 5, q.get( 4 ) );

		q.fill( data1 );
		assert_eq( "detached p[2]", 1, p.get( 2 ) );
		assert_eq( "detached p[5]", 4, p.get( 5 ) );
	}

	/** Test {@link AutoPermutation::operator=( AutoPermutation& )}. */
	void AutoPermutationTest::testAssignment () {
		const size_t
			data1[] = { 1, 2, 3, 4, 5, 6, 7, 0 },
			data2[] = { 7, 0, 1, 2, 3, 4, 5, 6 };
		AutoPermutation p( 8, 10 );
		p.fill( data1 );

		AutoPermutation q( 0 );
		q = p;
		assert_eq( "dimension", 8, q.countDimensions() );
		assert_eq( "w[0]", 1, q.get( 0 ) );
		assert_eq( "w[3]", 4, q.get( 3 ) );

		p.fill( data2 );
		assert_eq( "detached q[1]", 2, q.get( 1 ) );
		assert_eq( "detached q[4]", 5, q.get( 4 ) );

		q.fill( data1 );
		assert_eq( "detached p[2]", 1, p.get( 2 ) );
		assert_eq( "detached p[5]", 4, p.get( 5 ) );
	}

	/** Test {@link Permutation::operator==( Permutation& )}. */
	void AutoPermutationTest::testOperators () {
		const size_t
			data1[] = { 1, 2, 3, 4, 5, 6, 7, 0 },
			data2[] = { 7, 0, 1, 2, 3, 4, 5, 6 };
		AutoPermutation o( 8 ), p( 8, 10 ), q( 8 );
		o.fill( data1 );
		p.fill( data1 );
		q.fill( data2 );
		assert_true( "o=o", o == o );
		assert_true( "o=p", o == p );
		assert_true( "p=p", p == p );
		assert_false( "! o=q", o == q );
		assert_false( "! q=o", q == o );
		assert_false( "! p=q", p == q );
		assert_false( "! q=p", q == p );
		assert_true( "q=q", q == q );
	}

	/** Test {@link Permutation::get( size_t )}. */
	void AutoPermutationTest::testFillGetSwap () {
		const size_t data[] = { 4, 2, 3, 1, 0 };

		AutoPermutation p( 5 );
		p.fill( data );
		assert_eq( "p[0]", 4, p.get( 0 ) );
		assert_eq( "p[1]", 2, p.get( 1 ) );
		assert_eq( "p[2]", 3, p.get( 2 ) );
		assert_eq( "p[3]", 1, p.get( 3 ) );
		assert_eq( "p[4]", 0, p.get( 4 ) );

		AutoPermutation q( 5, 10 );
		q.fill( data );
		assert_true( "p == q", p == q );
		assert_true( "q == p", q == p );
		assert_eq( "q[0]", 4, q.get( 0 ) );
		assert_eq( "q[1]", 2, q.get( 1 ) );
		assert_eq( "q[2]", 3, q.get( 2 ) );
		assert_eq( "q[3]", 1, q.get( 3 ) );
		assert_eq( "q[4]", 0, q.get( 4 ) );

		p.swap( 0, 1 );
		p.swap( 2, 3 );
		assert_eq( "swap p[0]", 2, p.get( 0 ) );
		assert_eq( "swap p[1]", 4, p.get( 1 ) );
		assert_eq( "swap p[2]", 1, p.get( 2 ) );
		assert_eq( "swap p[3]", 3, p.get( 3 ) );
		assert_eq( "swap p[4]", 0, p.get( 4 ) );
	}

	/** Test {@link AutoPermutation::exactSize( size_t )}. */
	void AutoPermutationTest::testResize () {
		const size_t data[] = { 0, 1, 2, 3, 4, 5, 6, 7 };

		AutoPermutation p( 0 );
		assert_eq( "0 dims", 0, p.countDimensions() );
		p.fill( data );

		p.exactSize( 1 );
		assert_eq( "1 dim", 1, p.countDimensions() );
		p.fill( data );
		assert_eq( "1 0", 0, p.get( 0 ) );

		p.exactSize( 2 );
		assert_eq( "2 dims", 2, p.countDimensions() );
		p.fill( data );
		assert_eq( "2 0", 0, p.get( 0 ) );
		assert_eq( "2 1", 1, p.get( 1 ) );

		p.exactSize( 3 );
		assert_eq( "3 dims", 3, p.countDimensions() );
		p.fill( data );
		assert_eq( "3 0", 0, p.get( 0 ) );
		assert_eq( "3 1", 1, p.get( 1 ) );
		assert_eq( "3 2", 2, p.get( 2 ) );

		p.exactSize( 8 );
		assert_eq( "8 dims", 8, p.countDimensions() );
		p.fill( data );
		assert_eq( "8 0", 0, p.get( 0 ) );
		assert_eq( "8 6", 6, p.get( 6 ) );
		assert_eq( "8 7", 7, p.get( 7 ) );

		p.exactSize( 4 );
		assert_eq( "4 dims", 4, p.countDimensions() );
		p.fill( data );
		assert_eq( "4 0", 0, p.get( 0 ) );
		assert_eq( "4 3", 3, p.get( 3 ) );

		p.exactSize( 0 );
		assert_eq( "0 dims again", 0, p.countDimensions() );

		p.exactSize( 1 );
		assert_eq( "1 dim again", 1, p.countDimensions() );
		p.fill( data );
		assert_eq( "1 again", 0, p.get( 0 ) );
	}

	/** Test {@link AutoPermutation::upSize( size_t, size_t )}. */
	void AutoPermutationTest::testUpSize () {
		const size_t data[] = { 1, 2, 3, 4, 5, 6, 7, 8 };

		AutoPermutation p( 0 );
		assert_eq( "0 dims", 0, p.countDimensions() );

		p.upSize( 1 );
		assert_eq( "1 dim", 1, p.countDimensions() );
		p.fill( data );
		assert_eq( "1 -> 1", 1, p.get( 0 ) );

		p.upSize( 2 );
		assert_eq( "2 dims", 2, p.countDimensions() );
		p.fill( data );
		assert_eq( "2 -> 1", 1, p.get( 0 ) );
		assert_eq( "2 -> 2", 2, p.get( 1 ) );

		p.upSize( 3 );
		assert_eq( "3 dims", 3, p.countDimensions() );
		p.fill( data );
		assert_eq( "3 -> 1", 1, p.get( 0 ) );
		assert_eq( "3 -> 2", 2, p.get( 1 ) );
		assert_eq( "3 -> 3", 3, p.get( 2 ) );

		p.upSize( 8 );
		assert_eq( "8 dims", 8, p.countDimensions() );
		p.fill( data );
		assert_eq( "8 -> 1", 1, p.get( 0 ) );
		assert_eq( "8 -> 8", 8, p.get( 7 ) );

		p.upSize( 6 );
		assert_eq( "6 dims", 6, p.countDimensions() );
		p.fill( data );
		assert_eq( "6 -> 1", 1, p.get( 0 ) );
		assert_eq( "6 -> 6", 6, p.get( 5 ) );
	}

}
