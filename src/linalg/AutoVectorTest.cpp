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

#include "AutoVector.hpp"
#include "AutoMatrix.hpp"
#include "../TestSuite.hpp"
#include <math.h>

using namespace std;
using namespace unitpp;
using namespace linalg;

namespace test {

	/** Tests the class {@link linalg::AutoVector} and thereby {@link linalg::Vector}. */
	class AutoVectorTest : public TestSuite {

		void testCopyConstructor ();
		void testAssignment ();
		void testOperators ();
		void testGetSet ();
		void testFill ();
		void testIsNull ();
		void testResizePacked ();
		void testResizeNonPacked ();
		void testUpSize ();
		void testPackedSubVector ();
		void testStriddenSubVector ();
		void testSumSquares ();
		void testInnerProduct ();
		void testScale ();
		void testAxpy ();
		void testHouseholder ();
		void testGemvStraight ();
		void testGemvTransposed ();
		void testSolveRU ();
		void testMultQ ();

		public:

		AutoVectorTest () : TestSuite( "linalg::AutoVectorTest" ) {
			addTestMethod( "AutoVectorTest::testCopyConstructor", this, &AutoVectorTest::testCopyConstructor );
			addTestMethod( "AutoVectorTest::testAssignment", this, &AutoVectorTest::testAssignment );
			addTestMethod( "AutoVectorTest::testOperators", this, &AutoVectorTest::testOperators );
			addTestMethod( "AutoVectorTest::testGetSet", this, &AutoVectorTest::testGetSet );
			addTestMethod( "AutoVectorTest::testFill", this, &AutoVectorTest::testFill );
			addTestMethod( "AutoVectorTest::testIsNull", this, &AutoVectorTest::testIsNull );
			addTestMethod( "AutoVectorTest::testResizePacked", this, &AutoVectorTest::testResizePacked );
			addTestMethod( "AutoVectorTest::testResizeNonPacked", this, &AutoVectorTest::testResizeNonPacked );
			addTestMethod( "AutoVectorTest::testUpSize", this, &AutoVectorTest::testUpSize );
			addTestMethod( "AutoVectorTest::testPackedSubVector", this, &AutoVectorTest::testPackedSubVector );
			addTestMethod( "AutoVectorTest::testStriddenSubVector", this, &AutoVectorTest::testStriddenSubVector );
			addTestMethod( "AutoVectorTest::testSumSquares", this, &AutoVectorTest::testSumSquares );
			addTestMethod( "AutoVectorTest::testInnerProduct", this, &AutoVectorTest::testInnerProduct );
			addTestMethod( "AutoVectorTest::testScale", this, &AutoVectorTest::testScale );
			addTestMethod( "AutoVectorTest::testAxpy", this, &AutoVectorTest::testAxpy );
			addTestMethod( "AutoVectorTest::testHouseholder", this, &AutoVectorTest::testHouseholder );
			addTestMethod( "AutoVectorTest::testGemvStraight", this, &AutoVectorTest::testGemvStraight );
			addTestMethod( "AutoVectorTest::testGemvTransposed", this, &AutoVectorTest::testGemvTransposed );
			addTestMethod( "AutoVectorTest::testSolveRU", this, &AutoVectorTest::testSolveRU );
			addTestMethod( "AutoVectorTest::testMultQ", this, &AutoVectorTest::testMultQ );
		}
	} * vectorTest = new AutoVectorTest();	// automatically freed by unit++

	/** Test {@link AutoVector::AutoVector( AutoVector& )}. */
	void AutoVectorTest::testCopyConstructor () {
		const double
			data1[] = { 1, 2, 3, 4, 5, 6, 7, 8 },
			data2[] = { -1, -2, -3, -4, -5, -6, -7, -8 };
		AutoVector v( 8, 10 );
		v.fill( data1 );

		AutoVector w( v );
		assert_eq( "dimension", 8, w.countDimensions() );
		assert_eq( "w[0]", 1, w.get( 0 ) );
		assert_eq( "w[3]", 4, w.get( 3 ) );

		v.fill( data2 );
		assert_eq( "detached w[1]", 2, w.get( 1 ) );
		assert_eq( "detached w[4]", 5, w.get( 4 ) );

		w.fill( data1 );
		assert_eq( "detached v[2]", -3, v.get( 2 ) );
		assert_eq( "detached v[5]", -6, v.get( 5 ) );
	}

	/** Test {@link AutoVector::operator=( AutoVector& )}. */
	void AutoVectorTest::testAssignment () {
		const double
			data1[] = { 1, 2, 3, 4, 5, 6, 7, 8 },
			data2[] = { -1, -2, -3, -4, -5, -6, -7, -8 };
		AutoVector v( 8, 10 );
		v.fill( data1 );

		AutoVector w( 0 );
		w = v;
		assert_eq( "dimension", 8, w.countDimensions() );
		assert_eq( "w[0]", 1, w.get( 0 ) );
		assert_eq( "w[3]", 4, w.get( 3 ) );

		v.fill( data2 );
		assert_eq( "detached w[1]", 2, w.get( 1 ) );
		assert_eq( "detached w[4]", 5, w.get( 4 ) );

		w.fill( data1 );
		assert_eq( "detached v[2]", -3, v.get( 2 ) );
		assert_eq( "detached v[5]", -6, v.get( 5 ) );
	}

	/** Test {@link Vector::operator==( Vector& )}. */
	void AutoVectorTest::testOperators () {
		const double
			data1[] = { 1, 2, 3, 4, 5, 6, 7, 8 },
			data2[] = { -1, -2, -3, -4, -5, -6, -7, -8 };
		AutoVector u( 8 ), v( 8, 10 ), w( 8 );
		u.fill( data1 );
		v.fill( data1 );
		w.fill( data2 );
		assert_true( "u=u", u == u );
		assert_true( "u=v", u == v );
		assert_true( "v=v", v == v );
		assert_false( "! u=w", u == w );
		assert_false( "! w=u", w == u );
		assert_false( "! v=w", v == w );
		assert_false( "! w=v", w == v );
		assert_true( "w=w", w == w );
	}

	/** Test {@link Vector::get( size_t )} and {@link Vector::set( size_t, double )}. */
	void AutoVectorTest::testGetSet () {
		AutoVector v( 8 );
		v.set( 0, 1 );
		assert_eq( "v[0]", 1, v.get( 0 ) );
		v.set( 3, 2 );
		assert_eq( "v[3]", 2, v.get( 3 ) );
		v.set( 7, -10.5 );
		assert_eq( "v[7]", -10.5, v.get( 7 ) );
	}

	/** Test {@link Vector::fill( double )} and {@link Vector::fill( double* )}. */
	void AutoVectorTest::testFill () {
		const double data[] = { -1, -2, -3, -4 };
		AutoVector v( 4 );

		v.fill( 4 );
		assert_eq( "v[0]", 4, v.get( 0 ) );
		assert_eq( "v[1]", 4, v.get( 1 ) );
		assert_eq( "v[2]", 4, v.get( 2 ) );
		assert_eq( "v[3]", 4, v.get( 3 ) );

		v.fill( data );
		assert_eq( "v'[0]", -1, v.get( 0 ) );
		assert_eq( "v'[1]", -2, v.get( 1 ) );
		assert_eq( "v'[2]", -3, v.get( 2 ) );
		assert_eq( "v'[3]", -4, v.get( 3 ) );
	}

	/** Test {@link Vector::isNull()}. */
	void AutoVectorTest::testIsNull () {
		const AutoVector o( 0 );
		assert_true( "0-dim vector is always null", o.isNull() );

		const double data[] = { 1, 2, 3 };
		AutoVector v( 3 );
		v.fill( data );
		assert_false( "123 not null", v.isNull() );
		v.fill( 0.0 );
		assert_true( "000 is null", v.isNull() );
		v.set( 0, 1.0 );
		assert_false( "100 not null", v.isNull() );
		v.set( 0, 0.0 );
		v.set( 1, 1.0 );
		assert_false( "010 not null", v.isNull() );
		v.set( 1, 0.0 );
		v.set( 2, 1.0 );
		assert_false( "001 not null", v.isNull() );
	}

	/** Test {@link AutoVector::exactSize( size_t )}. */
	void AutoVectorTest::testResizePacked () {
		const double data[] = { 1, 2, 3, 4, 5, 6, 7, 8 };

		AutoVector v( 0 );
		assert_eq( "0 dims", 0, v.countDimensions() );
		v.fill( data );

		v.exactSize( 1 );
		assert_eq( "1 dim", 1, v.countDimensions() );
		v.fill( data );
		assert_eq( "1", 1, v.get( 0 ) );

		v.exactSize( 2 );
		assert_eq( "2 dims", 2, v.countDimensions() );
		v.fill( data );
		assert_eq( "2", 1, v.get( 0 ) );
		assert_eq( "2", 2, v.get( 1 ) );

		v.exactSize( 3 );
		assert_eq( "3 dims", 3, v.countDimensions() );
		v.fill( data );
		assert_eq( "3", 1, v.get( 0 ) );
		assert_eq( "3", 2, v.get( 1 ) );
		assert_eq( "3", 3, v.get( 2 ) );

		v.exactSize( 8 );
		assert_eq( "8 dims", 8, v.countDimensions() );
		v.fill( data );
		assert_eq( "8", 1, v.get( 0 ) );
		assert_eq( "8", 7, v.get( 6 ) );
		assert_eq( "8", 8, v.get( 7 ) );

		v.exactSize( 4 );
		assert_eq( "4 dims", 4, v.countDimensions() );
		v.fill( data );
		assert_eq( "4", 1, v.get( 0 ) );
		assert_eq( "4", 4, v.get( 3 ) );

		v.exactSize( 0 );
		assert_eq( "0 dims again", 0, v.countDimensions() );

		v.exactSize( 1 );
		assert_eq( "1 dim again", 1, v.countDimensions() );
		v.fill( data );
		assert_eq( "1 again", 1, v.get( 0 ) );
	}

	/** Test resizing a vector with extra space. */
	void AutoVectorTest::testResizeNonPacked () {
		const double data[] = { 1, 2, 3, 4, 5, 6, 7, 8 };
		const size_t dims = 4;
		for ( size_t tda = dims; tda < 6; ++tda ) {
			char buf[256];
			sprintf( buf, "with %u dims %u tda\t", dims, tda );
			string msg( buf );
			AutoVector v( 4, tda );
			assert_eq( msg + "dims", dims, v.countDimensions() );
			v.fill( data );
			assert_eq( msg + "(0)", 1, v.get( 0 ) );
			assert_eq( msg + "(1)", 2, v.get( 1 ) );
			assert_eq( msg + "(2)", 3, v.get( 2 ) );
			assert_eq( msg + "(3)", 4, v.get( 3 ) );
			v.exactSize( 2 );
			assert_eq( "2 dims", 2, v.countDimensions() );
			assert_eq( msg + "resized (0)", 1, v.get( 0 ) );
			assert_eq( msg + "resized (1)", 2, v.get( 1 ) );
			v.exactSize( 3 );
			assert_eq( "3 dims", 3, v.countDimensions() );
			assert_eq( msg + "resized3 keeps (0)", 1, v.get( 0 ) );
			assert_eq( msg + "resized3 keeps (1)", 2, v.get( 1 ) );
		}
	}

	/** Test {@link AutoVector::upSize( size_t, size_t )}. */
	void AutoVectorTest::testUpSize () {
		const double data[] = { 1, 2, 3, 4, 5, 6, 7, 8 };

		AutoVector v( 0 );
		assert_eq( "0 dims", 0, v.countDimensions() );

		v.upSize( 1 );
		assert_eq( "1 dim", 1, v.countDimensions() );
		v.fill( data );
		assert_eq( "1 -> 1", 1, v.get( 0 ) );

		v.upSize( 2 );
		assert_eq( "2 dims", 2, v.countDimensions() );
		v.fill( data );
		assert_eq( "2 -> 1", 1, v.get( 0 ) );
		assert_eq( "2 -> 2", 2, v.get( 1 ) );

		v.upSize( 3 );
		assert_eq( "3 dims", 3, v.countDimensions() );
		v.fill( data );
		assert_eq( "3 -> 1", 1, v.get( 0 ) );
		assert_eq( "3 -> 2", 2, v.get( 1 ) );
		assert_eq( "3 -> 3", 3, v.get( 2 ) );

		v.upSize( 8 );
		assert_eq( "8 dims", 8, v.countDimensions() );
		v.fill( data );
		assert_eq( "8 -> 1", 1, v.get( 0 ) );
		assert_eq( "8 -> 8", 8, v.get( 7 ) );

		v.upSize( 6 );
		assert_eq( "6 dims", 6, v.countDimensions() );
		v.fill( data );
		assert_eq( "6 -> 1", 1, v.get( 0 ) );
		assert_eq( "6 -> 6", 6, v.get( 5 ) );

		v.upSize( 0 );
		assert_eq( "0 dims again", 0, v.countDimensions() );
	}

	/** Test {@link Vector::subVector( size_t, size_t )}. */
	void AutoVectorTest::testPackedSubVector () {
		const double data[] = { 1, 2, 3, 4, 5, 6, 7, 8 };
		AutoVector v( 8 );
		v.fill( data );

		Vector s = v.subVector( 3, 3 );
		assert_eq( "sub dim", 3, s.countDimensions() );
		assert_eq( "s[0]", 4, s.get( 0 ) );
		assert_eq( "s[1]", 5, s.get( 1 ) );
		assert_eq( "s[2]", 6, s.get( 2 ) );

		s = s.subVector( 1, 2 );
		assert_eq( "ss dim", 2, s.countDimensions() );
		assert_eq( "ss[0]", 5, s.get( 0 ) );
		assert_eq( "ss[1]", 6, s.get( 1 ) );

		s = s.subVector( 0, 1 );
		assert_eq( "s1 dim", 1, s.countDimensions() );
		assert_eq( "s1[0]", 5, s.get( 0 ) );

		s = s.subVector( 1, 0 );
		assert_eq( "s0 dim", 0, s.countDimensions() );
	}

	/** Test {@link Vector::subVector( size_t, size_t )} with nontrivial stride. */
	void AutoVectorTest::testStriddenSubVector () {
		const double data[] = { 1, 2, 3, 4, 5, 6, 7, 8 };
		AutoMatrix m( 3, 2 );
		m.fill( data );
		Vector v = m.columnVector( 1 );
		assert_eq( "v dim", 3, v.countDimensions() );

		Vector s( v.subVector( 1, 2 ) );
		assert_eq( "s dim", 2, s.countDimensions() );
		assert_eq( "s[0]", 4, s.get( 0 ) );
		assert_eq( "s[1]", 6, s.get( 1 ) );
	}

	/** Test {@link Vector::sumSquares()}. */
	void AutoVectorTest::testSumSquares () {
		const double data[] = { 7, 9, 15, 16, 17 };
		AutoVector v( 5 );
		v.fill( data );
		assert_eq( "Norm square", 900, v.sumSquares() );
	}

	/** Test {@link Vector::innerProduct()}. */
	void AutoVectorTest::testInnerProduct () {
		const double
			dataV[] = { 1, 3, -2 },
			dataW[] = { 1, 2, 7 };
		AutoVector v( 3 ), w( 3 );
		v.fill( dataV );
		w.fill( dataW );
		assert_eq( "v*w", -7, v.innerProduct( w ) );
		assert_eq( "w*v", -7, w.innerProduct( v ) );
	}

	/** Test {@link Vector::scale()}. */
	void AutoVectorTest::testScale () {
		const double data[] = { 1, 0, -1 };
		AutoVector v( 3 );
		v.fill( data );
		v.scale( 2.0 );
		assert_eq( "v[0]", 2, v.get( 0 ) );
		assert_eq( "v[1]", 0, v.get( 1 ) );
		assert_eq( "v[2]", -2, v.get( 2 ) );
	}

	/** Test {@link Vector::axpy()}. */
	void AutoVectorTest::testAxpy () {
		const double
			dataV[] = { 1, 2, 3 },
			dataW[] = { 1, 5, 11 };
		AutoVector v( 3 ), w( 3 );
		v.fill( dataV );
		w.fill( dataW );
		v.axpy( 2, w );
		assert_eq( "v[0]", 3, v.get( 0 ) );
		assert_eq( "v[1]", 12, v.get( 1 ) );
		assert_eq( "v[2]", 25, v.get( 2 ) );
	}

	/** Test {@link Vector::householderize()}
	* and {@link Vector::householderTransform( double, Vector )}.
	*/
	void AutoVectorTest::testHouseholder () {
		AutoVector v( 0 );
		// the value of tau is undefined this just asserts that the method does not fail
		double tau = v.householderize();

		const double data1[] = { 3, 4 };
		v.upSize( 2 );
		v.fill( data1 );
		tau = v.householderize();

		AutoVector w( 2 );
		w.fill( data1 );
		w.householderTransform( tau, v );

		assert_eq( "w1[0]", 5, fabs( w.get( 0 ) ) );
		assert_eq( "w1[1]", 0, w.get( 1 ) );

		const double data2[] = { -12, 15, 16 };		// Vector of length 25
		v.upSize( 3 );
		v.fill( data2 );
		tau = v.householderize();

		w.upSize( 3 );
		w.fill( data2 );
		w.householderTransform( tau, v );

		assert_eq( "w2[0]", 25, fabs( w.get( 0 ) ) );
		assert_close( "w2[1]", 0, w.get( 1 ) );
		assert_close( "w2[2]", 0, w.get( 2 ) );

		const double data3[] = { 9, -13, 16, -17, 19 };		// Vector of length 34
		v.upSize( 5 );
		v.fill( data3 );
		tau = v.householderize();

		w.upSize( 5 );
		w.fill( data3 );
		w.householderTransform( tau, v );

		assert_eq( "w3[0]", 34, fabs( w.get( 0 ) ) );
		assert_close( "w3[1]", 0, w.get( 1 ) );
		assert_close( "w3[2]", 0, w.get( 2 ) );
		assert_close( "w3[3]", 0, w.get( 3 ) );
		assert_close( "w3[4]", 0, w.get( 4 ) );
	}

	/** Test {@link Vector::gemv} for non-transposed matrices. */
	void AutoVectorTest::testGemvStraight () {
		// Test 0-dim vector space automorphism
		{
			const AutoMatrix a( 0, 0 );
			const AutoVector v( 0 );
			AutoVector w( 0 );
			w.gemv( 1.0, a, false, v, 1.0 );
		}

		// Test 0-dim to 1-dim vector space homomorphism
		{
			const AutoMatrix a( 1, 0 );
			const AutoVector v( 0 );
			AutoVector w( 1 );
			w.set( 0, 13.6 );	// set to non-zero
			w.gemv( 1.0, a, false, v, -1.0 );
			assert_eq( "0-1 w[0]", -13.6, w.get( 0 ) );
		}

		// Test 1-dim to 0-dim vector space homomorphism
		{
			const AutoMatrix a( 0, 1 );
			AutoVector v( 1 );
			AutoVector w( 0 );
			v.set( 0, -13.6 );	// set to non-zero
			w.gemv( 1.0, a, false, v, 1.0 );
			assert_eq( "1-0 v[0]", -13.6, v.get( 0 ) );
		}

		// Test some less trivial setting
		AutoMatrix a( 3, 2 );
		AutoVector v( 2 ), w( 3 );
		const double dataA[] = {
			1, 8,
			2, 16,
			4, 32
		};
		a.fill( dataA );

		v.fill( 0.0 );
		w.gemv( 1.0, (const Matrix) a, false, (const Vector) v, 0.0 );
		assert_eq( "alpha=1, beta=0, w[0]", 0.0, w.get( 0 ) );
		assert_eq( "alpha=1, beta=0, w[1]", 0.0, w.get( 1 ) );
		assert_eq( "alpha=1, beta=0, w[2]", 0.0, w.get( 2 ) );

		v.set( 0, 1.0 );
		v.set( 1, -1.0 );
		w.gemv( 0.0, (const Matrix) a, false, (const Vector) v, 1.0 );
		assert_eq( "alpha=0, beta=1, w[0]", 0.0, w.get( 0 ) );
		assert_eq( "alpha=0, beta=1, w[1]", 0.0, w.get( 1 ) );
		assert_eq( "alpha=0, beta=1, w[2]", 0.0, w.get( 2 ) );

		w.gemv( 1.0, (const Matrix) a, false, (const Vector) v, 1.0 );
		assert_eq( "alpha=1, beta=1, w[0]", -7, w.get( 0 ) );
		assert_eq( "alpha=1, beta=1, w[1]", -14, w.get( 1 ) );
		assert_eq( "alpha=1, beta=1, w[2]", -28, w.get( 2 ) );
	}

	/** Test {@link Vector::gemv} for transposed matrices. */
	void AutoVectorTest::testGemvTransposed () {
		// Test 0-dim vector space automorphism
		{
			const AutoMatrix a( 0, 0 );
			const AutoVector v( 0 );
			AutoVector w( 0 );
			w.gemv( 1.0, a, true, v, 1.0 );
		}

		// Test 0-dim to 1-dim vector space homomorphism
		{
			const AutoMatrix a( 1, 0 );
			AutoVector v( 1 );
			AutoVector w( 0 );
			v.set( 0, -13.6 );	// set to non-zero
			w.gemv( 1.0, a, true, v, 1.0 );
			assert_eq( "1-0 v[0]", -13.6, v.get( 0 ) );
		}

		// Test 1-dim to 0-dim vector space homomorphism
		{
			const AutoMatrix a( 0, 1 );
			const AutoVector v( 0 );
			AutoVector w( 1 );
			w.set( 0, 13.6 );	// set to non-zero
			w.gemv( 1.0, a, true, v, -1.0 );
			assert_eq( "0-1 w[0]", -13.6, w.get( 0 ) );
		}

		// Test some less trivial setting
		AutoMatrix a( 3, 2 );
		AutoVector v( 3 ), w( 2 );
		const double dataA[] = {
			1, 8,
			2, 16,
			4, 32
		};
		a.fill( dataA );

		v.fill( 0.0 );
		w.gemv( 1.0, (const Matrix) a, true, (const Vector) v, 0.0 );
		assert_eq( "alpha=1, beta=0, w[0]", 0.0, w.get( 0 ) );
		assert_eq( "alpha=1, beta=0, w[1]", 0.0, w.get( 1 ) );

		v.set( 0, 1.0 );
		v.set( 1, -1.0 );
		v.set( 2, 100.0 );
		w.gemv( 0.0, (const Matrix) a, true, (const Vector) v, 1.0 );
		assert_eq( "alpha=0, beta=1, w[0]", 0.0, w.get( 0 ) );
		assert_eq( "alpha=0, beta=1, w[1]", 0.0, w.get( 1 ) );

		w.gemv( 1.0, (const Matrix) a, true, (const Vector) v, 1.0 );
		assert_eq( "alpha=1, beta=1, w[0]", 399, w.get( 0 ) );
		assert_eq( "alpha=1, beta=1, w[1]", 3192, w.get( 1 ) );
	}

	/** Test {@link Vector::solveR}. */
	void AutoVectorTest::testSolveRU () {
		const AutoMatrix om( 0, 0 );
		const AutoVector ov( 0 );
		AutoVector zv( 0 );
		zv.solveR( om, ov );
		assert_eq( "Zero dim system", 0, zv.countDimensions() );

		const double data[] = {
			1, 0,
			5.8, 1,
			7.3, 6.6	// ignored
		};

		AutoMatrix m( 3, 2 );
		m.fill( data );

		AutoVector b( 2 );
		b.set( 0, 5 );
		b.set( 1, -4);

		AutoVector x( 2 );

		{
			x.solveR( m, b );
			assert_eq( "x[0]", 5, x.get( 0 ) );
			assert_eq( "x[1]", -4, x.get( 1 ) );
		}

		m.set( 0, 1, -1 );
		{
			x.solveR( m, b );
			assert_eq( "x[0]", 1, x.get( 0 ) );
			assert_eq( "x[1]", -4, x.get( 1 ) );
		}

		m.set( 1, 1, 2 );
		{
			x.solveR( m, b );
			assert_eq( "x[0]", 3, x.get( 0 ) );
			assert_eq( "x[1]", -2, x.get( 1 ) );
		}
	}

	/** Test {@link Vector::multQ}. */
	void AutoVectorTest::testMultQ () {
		// Test 0-dim vector space automorphism
		{
			const AutoMatrix qr( 0, 0 );
			const AutoVector tau( 0 );
			AutoVector v( 0 );
			v.multQ( tau, qr, false );
			v.multQ( tau, qr, true );
		}

		// Test 0-dim to 1-dim vector space homomorphism (autoextended to 1-dim space identity)
		{
			const AutoMatrix qr( 1, 0 );
			const AutoVector tau( 0 );
			AutoVector v( 1 );

			v.set( 0, 13.6 );
			v.multQ( tau, qr, false );
			assert_eq( "0-1 Qv[0]", 13.6, v.get( 0 ) );

			v.set( 0, -13.6 );
			v.multQ( tau, qr, true );
			assert_eq( "0-1 QTv[0]", -13.6, v.get( 0 ) );
		}

		// Test 1-dim to 0-dim vector space homomorphism (autoextended to 1-dim space identity)
		{
			const AutoMatrix qr( 1, 0 );
			const AutoVector tau( 0 );
			AutoVector v( 1 );

			v.set( 0, 17.9 );
			v.multQ( tau, qr, false );
			assert_eq( "0-1 Qv[0]", 17.9, v.get( 0 ) );

			v.set( 0, -17.9 );
			v.multQ( tau, qr, true );
			assert_eq( "0-1 QTv[0]", -17.9, v.get( 0 ) );
		}

		// Test some less trivial setting
		const double data[] = {
			0, 3,
			5, 4,
			0, 0
		};
		AutoMatrix qr( 3, 2 );
		qr.fill( data );
		AutoVector tau( 2 );
		qr.factorizeQR( tau );

		AutoVector v( 3 );
		v.set( 0, 7 );
		v.set( 1, 13 );
		v.set( 2, 29 );
		v.multQ( tau, qr, false );
		assert_eq( "Qv[0]", 13.0, fabs( v.get( 0 ) ) );
		assert_eq( "Qv[1]", 7.0, fabs( v.get( 1 ) ) );
		assert_eq( "Qv[2]", 29.0, v.get( 2 ) );		// unchanged due to 3rd column of Q = (0,0,1)^T

		AutoVector w( 3 );
		w.set( 0, 17 );
		w.set( 1, 19 );
		w.set( 2, 23 );
		w.multQ( tau, qr, true );
		assert_eq( "Qw[0]", 19.0, fabs( w.get( 0 ) ) );
		assert_eq( "Qw[1]", 17.0, fabs( w.get( 1 ) ) );
		assert_eq( "Qw[2]", 23.0, w.get( 2 ) );		// unchanged due to 3rd column of Q = (0,0,1)^T
	}

}
