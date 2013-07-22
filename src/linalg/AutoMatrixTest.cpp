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

#include "AutoMatrix.hpp"
#include "AutoVector.hpp"
#include "AutoPermutation.hpp"
#include "../TestSuite.hpp"
#include <math.h>

using namespace std;
using namespace unitpp;
using namespace linalg;

namespace test {

	/** Tests the class {@link linalg::AutoMatrix} and thereby {@link linalg::Matrix}. */
	class AutoMatrixTest : public TestSuite {

		void testCopyConstructorEqualTda ();
		void testCopyConstructorExtraTda ();
		void testAssignmentEqualData ();
		void testAssignmentExtraData ();
		void testGetSet ();
		void testFill ();
		void testResizePacked ();
		void testResizeNonPacked ();
		void testUpSize ();
		void testRowColumnVector ();
		void testDiagonalVector ();
		void testNewRowColumn ();
		void testRemoveRowColumn ();
		void testHouseholder ();
		void testGemmStraight ();
		void testGemmTransposed ();
		void testFactorizeQR ();
		void testFactorizeQRP ();

		public:

		AutoMatrixTest () : TestSuite( "linalg::AutoMatrixTest" ) {
			addTestMethod( "AutoMatrixTest::testCopyConstructorEqualTda", this, &AutoMatrixTest::testCopyConstructorEqualTda );
			addTestMethod( "AutoMatrixTest::testCopyConstructorExtraTda", this, &AutoMatrixTest::testCopyConstructorExtraTda );
			addTestMethod( "AutoMatrixTest::testAssignmentEqualData", this, &AutoMatrixTest::testAssignmentEqualData );
			addTestMethod( "AutoMatrixTest::testAssignmentExtraData", this, &AutoMatrixTest::testAssignmentExtraData );
			addTestMethod( "AutoMatrixTest::testGetSet", this, &AutoMatrixTest::testGetSet );
			addTestMethod( "AutoMatrixTest::testFill", this, &AutoMatrixTest::testFill );
			addTestMethod( "AutoMatrixTest::testResizePacked", this, &AutoMatrixTest::testResizePacked );
			addTestMethod( "AutoMatrixTest::testResizeNonPacked", this, &AutoMatrixTest::testResizeNonPacked );
			addTestMethod( "AutoMatrixTest::testUpSize", this, &AutoMatrixTest::testUpSize );
			addTestMethod( "AutoMatrixTest::testRowColumnVector", this, &AutoMatrixTest::testRowColumnVector );
			addTestMethod( "AutoMatrixTest::testDiagonalVector", this, &AutoMatrixTest::testDiagonalVector );
			addTestMethod( "AutoMatrixTest::testNewRowColumn", this, &AutoMatrixTest::testNewRowColumn );
			addTestMethod( "AutoMatrixTest::testRemoveRowColumn", this, &AutoMatrixTest::testRemoveRowColumn );
			addTestMethod( "AutoMatrixTest::testHouseholder", this, &AutoMatrixTest::testHouseholder );
			addTestMethod( "AutoMatrixTest::testGemmStraight", this, &AutoMatrixTest::testGemmStraight );
			addTestMethod( "AutoMatrixTest::testGemmTransposed", this, &AutoMatrixTest::testGemmTransposed );
			addTestMethod( "AutoMatrixTest::testFactorizeQR", this, &AutoMatrixTest::testFactorizeQR );
			addTestMethod( "AutoMatrixTest::testFactorizeQRP", this, &AutoMatrixTest::testFactorizeQRP );
		}
	} * matrixTest = new AutoMatrixTest();	// automatically freed by unit++

	/** Test {@link AutoMatrix::AutoMatrix( AutoMatrix& )}. */
	void AutoMatrixTest::testCopyConstructorEqualTda () {
		const double
			data1[] = { 1, 2, 3, 4, 5, 6 },
			data2[] = { -1, -2, -3, -4, -5, -6 };
		AutoMatrix m( 3, 2 );
		m.fill( data1 );

		AutoMatrix n( m );
		assert_eq( "rows", 3, n.countRows() );
		assert_eq( "cols", 2, n.countColumns() );
		assert_eq( "copied 2", 2, n.get( 0, 1 ) );
		assert_eq( "copied 6", 6, n.get( 2, 1 ) );

		m.fill( data2 );
		assert_eq( "detached copy 1", 1, n.get( 0, 0 ) );
		assert_eq( "detached copy 5", 5, n.get( 2, 0 ) );

		n.fill( data1 );
		assert_eq( "detached original 3", -3, m.get( 1, 0 ) );
		assert_eq( "detached original 4", -4, m.get( 1, 1 ) );
	}

	/** Test {@link AutoMatrix::AutoMatrix( AutoMatrix& )} with supersize original. */
	void AutoMatrixTest::testCopyConstructorExtraTda () {
		const double
			data1[] = { 1, 2, 3, 4, 5, 6 },
			data2[] = { -1, -2, -3, -4, -5, -6 };
		AutoMatrix m( 3, 2, 4, 5 );
		m.fill( data1 );

		AutoMatrix n( m );
		assert_eq( "rows", 3, n.countRows() );
		assert_eq( "cols", 2, n.countColumns() );
		assert_eq( "copied 2", 2, n.get( 0, 1 ) );
		assert_eq( "copied 6", 6, n.get( 2, 1 ) );

		m.fill( data2 );
		assert_eq( "detached copy 1", 1, n.get( 0, 0 ) );
		assert_eq( "detached copy 5", 5, n.get( 2, 0 ) );

		n.fill( data1 );
		assert_eq( "detached original 3", -3, m.get( 1, 0 ) );
		assert_eq( "detached original 4", -4, m.get( 1, 1 ) );
	}

	/** Test {@link AutoMatrix::operator=( AutoMatrix& )}. */
	void AutoMatrixTest::testAssignmentEqualData () {
		const double
			data1[] = { 1, 2, 3, 4, 5, 6 },
			data2[] = { -1, -2, -3, -4, -5, -6 };
		AutoMatrix m( 3, 2 );
		m.fill( data1 );

		AutoMatrix n( 10, 10 );
		n = m;
		assert_eq( "rows", 3, n.countRows() );
		assert_eq( "cols", 2, n.countColumns() );
		assert_eq( "copied 2", 2, n.get( 0, 1 ) );
		assert_eq( "copied 6", 6, n.get( 2, 1 ) );

		m.fill( data2 );
		assert_eq( "detached copy 1", 1, n.get( 0, 0 ) );
		assert_eq( "detached copy 5", 5, n.get( 2, 0 ) );

		n.fill( data1 );
		assert_eq( "detached original 3", -3, m.get( 1, 0 ) );
		assert_eq( "detached original 4", -4, m.get( 1, 1 ) );
	}

	/** Test {@link AutoMatrix::operator=( AutoMatrix& )} with supersize original. */
	void AutoMatrixTest::testAssignmentExtraData () {
		const double
			data1[] = { 1, 2, 3, 4, 5, 6 },
			data2[] = { -1, -2, -3, -4, -5, -6 };
		AutoMatrix m( 3, 2, 4, 5 );
		m.fill( data1 );

		AutoMatrix n( 0, 0 );
		n = m;
		assert_eq( "rows", 3, n.countRows() );
		assert_eq( "cols", 2, n.countColumns() );
		assert_eq( "copied 2", 2, n.get( 0, 1 ) );
		assert_eq( "copied 6", 6, n.get( 2, 1 ) );

		m.fill( data2 );
		assert_eq( "detached copy 1", 1, n.get( 0, 0 ) );
		assert_eq( "detached copy 5", 5, n.get( 2, 0 ) );

		n.fill( data1 );
		assert_eq( "detached original 3", -3, m.get( 1, 0 ) );
		assert_eq( "detached original 4", -4, m.get( 1, 1 ) );
	}

	/** Test {@link AutoMatrix::get(size_t,size_t)} and {@link AutoMatrix::set(size_t,size_t,double)}. */
	void AutoMatrixTest::testGetSet () {
		const double data[] = {
			1, 2,
			3, 4,
			5, 6
		};

		AutoMatrix m( 3, 2 );
		m.fill( data );
		assert_eq( "get", 1, m.get( 0, 0 ) );
		assert_eq( "get", 2, m.get( 0, 1 ) );
		assert_eq( "get", 3, m.get( 1, 0 ) );
		assert_eq( "get", 4, m.get( 1, 1 ) );
		assert_eq( "get", 5, m.get( 2, 0 ) );
		assert_eq( "get", 6, m.get( 2, 1 ) );

		m.set( 0, 0, -1 );
		m.set( 1, 0, -3 );
		m.set( 2, 1, -6 );
		assert_eq( "get", -1, m.get( 0, 0 ) );
		assert_eq( "get", -3, m.get( 1, 0 ) );
		assert_eq( "get", -6, m.get( 2, 1 ) );
	}

	/** Test {@link Matrix::fill( double )} and {@link Matrix::fill( double* )}. */
	void AutoMatrixTest::testFill () {
		const double data[] = { -1.0, -2.0, -3.0, -4.0, -5.0, -6.0, -7.0, -8.0, -9.0 };
		AutoMatrix m( 2, 2 );

		m.fill( 4 );
		assert_eq( "m[0,0]", 4.0, m.get( 0, 0 ) );
		assert_eq( "m[0,1]", 4.0, m.get( 0, 1 ) );
		assert_eq( "m[1,0]", 4.0, m.get( 1, 0 ) );
		assert_eq( "m[1,1]", 4.0, m.get( 1, 1 ) );

		// below also documents: C-style: row major
		m.fill( data );
		assert_eq( "m22[0,0]", -1.0, m.get( 0, 0 ) );
		assert_eq( "m22[0,1]", -2.0, m.get( 0, 1 ) );
		assert_eq( "m22[1,0]", -3.0, m.get( 1, 0 ) );
		assert_eq( "m22[1,1]", -4.0, m.get( 1, 1 ) );

		AutoMatrix m23( 2, 3 );
		m23.fill( data );
		assert_eq( "m23[0,0]", -1.0, m23.get( 0, 0 ) );
		assert_eq( "m23[0,1]", -2.0, m23.get( 0, 1 ) );
		assert_eq( "m23[0,2]", -3.0, m23.get( 0, 2 ) );
		assert_eq( "m23[1,0]", -4.0, m23.get( 1, 0 ) );
		assert_eq( "m23[1,1]", -5.0, m23.get( 1, 1 ) );
		assert_eq( "m23[1,2]", -6.0, m23.get( 1, 2 ) );

		AutoMatrix m32( 3, 2 );
		m32.fill( data );
		assert_eq( "m32[0,0]", -1.0, m32.get( 0, 0 ) );
		assert_eq( "m32[0,1]", -2.0, m32.get( 0, 1 ) );
		assert_eq( "m32[1,0]", -3.0, m32.get( 1, 0 ) );
		assert_eq( "m32[1,1]", -4.0, m32.get( 1, 1 ) );
		assert_eq( "m32[2,0]", -5.0, m32.get( 2, 0 ) );
		assert_eq( "m32[2,1]", -6.0, m32.get( 2, 1 ) );
	}

	/** Test resizing a matrix with rows == tda. */
	void AutoMatrixTest::testResizePacked () {
		const double data[] = {
			1, 2,
			3, 4,
			5, 6,
			7, 8
		};

		AutoMatrix m( 0, 0 );
		assert_eq( "0x0 rows", 0, m.countRows() );
		assert_eq( "0x0 columns", 0, m.countColumns() );
		m.fill( data );

		m.exactSize( 0, 1 );
		assert_eq( "0x1 rows", 0, m.countRows() );
		assert_eq( "0x1 columns", 1, m.countColumns() );
		m.fill( data );

		m.exactSize( 1, 0 );
		assert_eq( "1x0 rows", 1, m.countRows() );
		assert_eq( "1x0 columns", 0, m.countColumns() );
		m.fill( data );

		m.exactSize( 4, 2 );
		assert_eq( "4x2 rows", 4, m.countRows() );
		assert_eq( "4x2 columns", 2, m.countColumns() );
		m.fill( data );
		assert_eq( "[1,0]", 3, m.get( 1, 0 ) );
		assert_eq( "[3,1]", 8, m.get( 3, 1 ) );

		m.exactSize( 2, 2 );
		assert_eq( "2x2 rows", 2, m.countRows() );
		assert_eq( "2x2 columns", 2, m.countColumns() );
		assert_eq( "2x2 sized [1,0]", 3, m.get( 1, 0 ) );
		assert_eq( "2x2 sized [1,1]", 4, m.get( 1, 1 ) );

		m.exactSize( 3, 3 );
		assert_eq( "3x3 rows", 3, m.countRows() );
		assert_eq( "3x3 columns", 3, m.countColumns() );
		assert_eq( "3x3 sized [0,0]", 1, m.get( 0, 0 ) );
		assert_eq( "3x3 sized [0,1]", 2, m.get( 0, 1 ) );
		assert_eq( "3x3 sized [1,0]", 3, m.get( 1, 0 ) );
		assert_eq( "3x3 sized [1,1]", 4, m.get( 1, 1 ) );

		m.exactSize( 2, 1 );
		assert_eq( "2x1 rows", 2, m.countRows() );
		assert_eq( "2x1 columns", 1, m.countColumns() );
		assert_eq( "2x1 sized [0,0]", 1, m.get( 0, 0 ) );
		assert_eq( "2x1 sized [1,0]", 3, m.get( 1, 0 ) );

		m.exactSize( 4, 7 );
		assert_eq( "4x7 rows", 4, m.countRows() );
		assert_eq( "4x7 columns", 7, m.countColumns() );

		m.exactSize( 0, 0 );
		assert_eq( "0x0 rows again", 0, m.countRows() );
		assert_eq( "0x0 columns again", 0, m.countColumns() );
		m.fill( data );

		m.exactSize( 1, 1 );
		assert_eq( "1x1 rows finally", 1, m.countRows() );
		assert_eq( "1x1 columns finally", 1, m.countColumns() );
		m.fill( data );
		assert_eq( "1x1 sized [0,0] finally", 1, m.get( 0, 0 ) );
	}

	/** Test resizing a matrix with the other constructor with rows <= tda. */
	void AutoMatrixTest::testResizeNonPacked () {
		const double data[] = {
			1, 2,
			3, 4,
			5, 6,
			7, 8
		};
		const size_t rows = 4;
		const size_t cols = 2;

		for ( size_t tda = rows; tda < 6; ++tda ) {
			char buf[256];
			sprintf( buf, "with %u rows %u cols %u tda\t", rows, cols, tda );
			string msg( buf );

			AutoMatrix m( 4, 2, tda, cols );
			assert_eq( msg + "rows", rows, m.countRows() );
			assert_eq( msg + "columns", cols, m.countColumns() );
			m.fill( data );
			assert_eq( msg + "[1,0]", 3, m.get( 1, 0 ) );
			assert_eq( msg + "[3,1]", 8, m.get( 3, 1 ) );

			m.exactSize( 2, 2 );
			assert_eq( "2x2 rows", 2, m.countRows() );
			assert_eq( "2x2 columns", 2, m.countColumns() );
			assert_eq( msg + "2x2 sized [1,0]", 3, m.get( 1, 0 ) );
			assert_eq( msg + "2x2 sized [1,1]", 4, m.get( 1, 1 ) );

			m.exactSize( 3, 3 );
			assert_eq( "3x3 rows", 3, m.countRows() );
			assert_eq( "3x3 columns", 3, m.countColumns() );
			assert_eq( msg + "3x3 sized [0,0]", 1, m.get( 0, 0 ) );
			assert_eq( msg + "3x3 sized [0,1]", 2, m.get( 0, 1 ) );
			assert_eq( msg + "3x3 sized [1,0]", 3, m.get( 1, 0 ) );
			assert_eq( msg + "3x3 sized [1,1]", 4, m.get( 1, 1 ) );
		}
	}

	/** Test {@link AutoMatrix::upSize( size_t, size_t )}. */
	void AutoMatrixTest::testUpSize () {
		const double data[] = { 1, 2, 3, 4, 5, 6, 7, 8 };

		// GSL does not like 0 as dimension extent
		AutoMatrix m( 1, 1 );
		m.fill( data );
		assert_eq( "1x1[0,0]", 1, m.get( 0, 0 ) );

		m.upSize( 2, 1 );
		assert_eq( "2x1 rows", 2, m.countRows() );
		assert_eq( "2x1 columns", 1, m.countColumns() );
		m.fill( data );
		assert_eq( "2x1[0,0]", 1, m.get( 0, 0 ) );
		assert_eq( "2x1[1,0]", 2, m.get( 1, 0 ) );

		m.upSize( 1, 2 );
		assert_eq( "1x2 rows", 1, m.countRows() );
		assert_eq( "1x2 columns", 2, m.countColumns() );
		m.fill( data );
		assert_eq( "1x2[0,0]", 1, m.get( 0, 0 ) );
		assert_eq( "1x2[0,1]", 2, m.get( 0, 1 ) );

		m.upSize( 2, 3 );
		assert_eq( "2x3 rows", 2, m.countRows() );
		assert_eq( "2x3 columns", 3, m.countColumns() );
		m.fill( data );
		assert_eq( "2x3[0,0]", 1, m.get( 0, 0 ) );
		assert_eq( "2x3[0,2]", 3, m.get( 0, 2 ) );
		assert_eq( "2x3[1,0]", 4, m.get( 1, 0 ) );
		assert_eq( "2x3[1,2]", 6, m.get( 1, 2 ) );

		m.upSize( 3, 2 );
		assert_eq( "3x2 rows", 3, m.countRows() );
		assert_eq( "3x2 columns", 2, m.countColumns() );
		m.fill( data );
		assert_eq( "3x2[0,0]", 1, m.get( 0, 0 ) );
		assert_eq( "3x2[0,1]", 2, m.get( 0, 1 ) );
		assert_eq( "3x2[2,0]", 5, m.get( 2, 0 ) );
		assert_eq( "3x2[2,1]", 6, m.get( 2, 1 ) );
	}

	/** Test {@link AutoMatrix::rowVector( size_t )} and {@link AutoMatrix::columnVector( size_t )}. */
	void AutoMatrixTest::testRowColumnVector () {
		const double data[] = {
			1, 2, 3, 4, 5,
			6, 7, 8, 9, 10,
			11, 12, 13, 14, 15,
			16, 17, 18, 19, 20
		};

		AutoMatrix m( 4, 5 );
		m.fill( data );

		const Vector
			vv0 = m.rowVector( 0 ),
			vv3 = m.rowVector( 3 ),
			v0 = m.columnVector( 0 ),
			v3 = m.columnVector( 3 ),
			v4 = m.columnVector( 4 );

		assert_eq( "vv0 dim", 5, vv0.countDimensions() );
		assert_eq( "vv0[0]", 1, vv0.get( 0 ) );
		assert_eq( "vv0[1]", 2, vv0.get( 1 ) );
		assert_eq( "vv0[2]", 3, vv0.get( 2 ) );
		assert_eq( "vv0[3]", 4, vv0.get( 3 ) );
		assert_eq( "vv0[4]", 5, vv0.get( 4 ) );
		assert_eq( "vv3 dim", 5, vv3.countDimensions() );
		assert_eq( "vv3[0]", 16, vv3.get( 0 ) );
		assert_eq( "vv3[1]", 17, vv3.get( 1 ) );
		assert_eq( "vv3[2]", 18, vv3.get( 2 ) );
		assert_eq( "vv3[3]", 19, vv3.get( 3 ) );
		assert_eq( "vv3[4]", 20, vv3.get( 4 ) );
		assert_eq( "v0 dim", 4, v0.countDimensions() );
		assert_eq( "v0[0]", 1, v0.get( 0 ) );
		assert_eq( "v0[1]", 6, v0.get( 1 ) );
		assert_eq( "v0[2]", 11, v0.get( 2 ) );
		assert_eq( "v0[3]", 16, v0.get( 3 ) );
		assert_eq( "v3 dim", 4, v3.countDimensions() );
		assert_eq( "v3[0]", 4, v3.get( 0 ) );
		assert_eq( "v3[1]", 9, v3.get( 1 ) );
		assert_eq( "v3[2]", 14, v3.get( 2 ) );
		assert_eq( "v3[3]", 19, v3.get( 3 ) );
		assert_eq( "v4 dim", 4, v4.countDimensions() );
		assert_eq( "v4[0]", 5, v4.get( 0 ) );
		assert_eq( "v4[1]", 10, v4.get( 1 ) );
		assert_eq( "v4[2]", 15, v4.get( 2 ) );
		assert_eq( "v4[3]", 20, v4.get( 3 ) );

		AutoMatrix n( 5, 4, 7, 9 );
		n.fill( data );

		const Vector
			ww0 = n.rowVector( 0 ),
			ww4 = n.rowVector( 4 ),
			w0 = n.columnVector( 0 ),
			w3 = n.columnVector( 3 );

		assert_eq( "ww0 dim", 4, ww0.countDimensions() );
		assert_eq( "ww0[0]", 1, ww0.get( 0 ) );
		assert_eq( "ww0[1]", 2, ww0.get( 1 ) );
		assert_eq( "ww0[2]", 3, ww0.get( 2 ) );
		assert_eq( "ww0[3]", 4, ww0.get( 3 ) );
		assert_eq( "ww4 dim", 4, ww4.countDimensions() );
		assert_eq( "ww4[0]", 17, ww4.get( 0 ) );
		assert_eq( "ww4[1]", 18, ww4.get( 1 ) );
		assert_eq( "ww4[2]", 19, ww4.get( 2 ) );
		assert_eq( "ww4[3]", 20, ww4.get( 3 ) );
		assert_eq( "w0 dim", 5, w0.countDimensions() );
		assert_eq( "w0[0]", 1, w0.get( 0 ) );
		assert_eq( "w0[1]", 5, w0.get( 1 ) );
		assert_eq( "w0[2]", 9, w0.get( 2 ) );
		assert_eq( "w0[3]", 13, w0.get( 3 ) );
		assert_eq( "w0[4]", 17, w0.get( 4 ) );
		assert_eq( "w3 dim", 5, w3.countDimensions() );
		assert_eq( "w3[0]", 4, w3.get( 0 ) );
		assert_eq( "w3[1]", 8, w3.get( 1 ) );
		assert_eq( "w3[2]", 12, w3.get( 2 ) );
		assert_eq( "w3[3]", 16, w3.get( 3 ) );
		assert_eq( "w3[4]", 20, w3.get( 4 ) );

		AutoMatrix o( 0, 1 );

		const Vector o0 = o.columnVector( 0 );
		assert_eq( "o0 dim", 0, o0.countDimensions() );
	}

	/** Test {@link AutoMatrix::diagonalVector( size_t )}. */
	void AutoMatrixTest::testDiagonalVector () {
		const double data[] = {
			1, 2, 3,
			4, 5, 6
		};

		AutoMatrix m( 2, 3, 4, 5 );
		m.fill( data );

		const Vector
			subdiag1 = m.subdiagonalVector( 1 ),
			subdiag0 = m.subdiagonalVector( 0 ),
			diag = m.diagonalVector(),
			superdiag0 = m.superdiagonalVector( 0 ),
			superdiag1 = m.superdiagonalVector( 1 ),
			superdiag2 = m.superdiagonalVector( 2 );

		assert_eq( "subdiag1 dim", 1, subdiag1.countDimensions() );
		assert_eq( "subdiag1[0]", 4, subdiag1.get( 0 ) );

		assert_eq( "diag dim", 2, diag.countDimensions() );
		assert_eq( "diag[0]", 1, diag.get( 0 ) );
		assert_eq( "diag[1]", 5, diag.get( 1 ) );

		assert_eq( "subdiag0", diag, subdiag0 );
		assert_eq( "superdiag0", diag, superdiag0 );

		assert_eq( "superdiag1 dim", 2, superdiag1.countDimensions() );
		assert_eq( "superdiag1[0]", 2, superdiag1.get( 0 ) );
		assert_eq( "superdiag1[1]", 6, superdiag1.get( 1 ) );

		assert_eq( "superdiag2 dim", 1, superdiag2.countDimensions() );
		assert_eq( "superdiag2[0]", 3, superdiag2.get( 0 ) );

		{
			// Test that no exception occurs for allowed values
			AutoMatrix oo( 1, 0 );
			const Vector db = oo.diagonalVector();
			assert_eq( "db dim", 0, db.countDimensions() );
			const Vector b0 = oo.subdiagonalVector( 1 );
			assert_eq( "b0 dim", 0, b0.countDimensions() );
			const Vector b1 = oo.subdiagonalVector( 1 );
			assert_eq( "b1 dim", 0, b1.countDimensions() );
		}

		{
			// Test that no exception occurs for allowed values
			AutoMatrix o( 0, 1 );
			const Vector dp = o.diagonalVector();
			assert_eq( "dp dim", 0, dp.countDimensions() );
			const Vector p0 = o.superdiagonalVector( 0 );
			assert_eq( "p0 dim", 0, p0.countDimensions() );
			const Vector p1 = o.superdiagonalVector( 1 );
			assert_eq( "p1 dim", 0, p1.countDimensions() );
		}
	}

	/** Test {@link AutoMatrix::newRow()} and {@link AutoMatrix::newColumn()}. */
	void AutoMatrixTest::testNewRowColumn () {
		const double data[] = {
			1, 2,
			3, 4,
			5, 6
		};

		AutoMatrix m( 2, 1 );
		m.fill( data );

		Vector v = m.newColumn();
		v.set( 0, -3 );
		v.set( 1, -4 );
		assert_eq( "m[0,0]", 1, m.get( 0, 0 ) );
		assert_eq( "m[1,0]", 2, m.get( 1, 0 ) );
		assert_eq( "m[0,1]", -3, m.get( 0, 1 ) );
		assert_eq( "m[1,1]", -4, m.get( 1, 1 ) );

		m.fill( data );
		assert_eq( "v[0]", 2, v.get( 0 ) );
		assert_eq( "v[1]", 4, v.get( 1 ) );

		Vector w = m.newRow();
		w.set( 0, -5 );
		w.set( 1, -6 );
		assert_eq( "m'[0,0]", 1, m.get( 0, 0 ) );
		assert_eq( "m'[0,1]", 2, m.get( 0, 1 ) );
		assert_eq( "m'[1,0]", 3, m.get( 1, 0 ) );
		assert_eq( "m'[1,1]", 4, m.get( 1, 1 ) );
		assert_eq( "m'[2,0]", -5, m.get( 2, 0 ) );
		assert_eq( "m'[2,1]", -6, m.get( 2, 1 ) );

		m.fill( data );
		assert_eq( "w[0]", 5, w.get( 0 ) );
		assert_eq( "w[1]", 6, w.get( 1 ) );
	}

	/** Test {@link AutoMatrix::removeRow()} and {@link AutoMatrix::removeColumn()}. */
	void AutoMatrixTest::testRemoveRowColumn () {
		const double data[] = {
			0.0, 0.1, 0.2, 0.3, 0.4,
			1.0, 1.1, 1.2, 1.3, 1.4,
			2.0, 2.1, 2.2, 2.3, 2.4,
			3.0, 3.1, 3.2, 3.3, 3.4,
			4.0, 4.1, 4.2, 4.3, 4.4
		};

		AutoMatrix m( 5, 5 );
		m.fill( data );

		m.removeRow( 4 );
		assert_eq( "First row-count", 4, m.countRows() );
		assert_eq( "First col-count", 5, m.countColumns() );
		assert_eq( "m1[0,0]", 0.0, m.get( 0, 0 ) );
		assert_eq( "m1[1,4]", 1.4, m.get( 1, 4 ) );
		assert_eq( "m1[2,0]", 2.0, m.get( 2, 0 ) );
		assert_eq( "m1[3,4]", 3.4, m.get( 3, 4 ) );

		m.removeColumn( 0 );
		assert_eq( "Second row-count", 4, m.countRows() );
		assert_eq( "Second col-count", 4, m.countColumns() );
		assert_eq( "m2[0,0]", 0.1, m.get( 0, 0 ) );
		assert_eq( "m2[1,3]", 1.4, m.get( 1, 3 ) );
		assert_eq( "m2[2,0]", 2.1, m.get( 2, 0 ) );
		assert_eq( "m2[3,3]", 3.4, m.get( 3, 3 ) );

		m.removeRow( 0 );
		assert_eq( "Third row-count", 3, m.countRows() );
		assert_eq( "Third col-count", 4, m.countColumns() );
		assert_eq( "m3[0,0]", 1.1, m.get( 0, 0 ) );
		assert_eq( "m3[1,3]", 2.4, m.get( 1, 3 ) );
		assert_eq( "m3[2,0]", 3.1, m.get( 2, 0 ) );

		m.removeColumn( 3 );
		assert_eq( "Fourth row-count", 3, m.countRows() );
		assert_eq( "Fourth col-count", 3, m.countColumns() );
		assert_eq( "m4[0,0]", 1.1, m.get( 0, 0 ) );
		assert_eq( "m4[1,2]", 2.3, m.get( 1, 2 ) );
		assert_eq( "m4[2,2]", 3.3, m.get( 2, 2 ) );

		m.removeRow( 1 );
		assert_eq( "Fifth row-count", 2, m.countRows() );
		assert_eq( "Fifth col-count", 3, m.countColumns() );
		assert_eq( "m5[0,0]", 1.1, m.get( 0, 0 ) );
		assert_eq( "m5[1,2]", 3.3, m.get( 1, 2 ) );

		m.removeColumn( 1 );
		assert_eq( "Sixth row-count", 2, m.countRows() );
		assert_eq( "Sixth col-count", 2, m.countColumns() );
		assert_eq( "m6[0,0]", 1.1, m.get( 0, 0 ) );
		assert_eq( "m6[0,1]", 1.3, m.get( 0, 1 ) );
		assert_eq( "m6[1,0]", 3.1, m.get( 1, 0 ) );
		assert_eq( "m6[1,1]", 3.3, m.get( 1, 1 ) );
	}

	/** Test {@link Matrix::householderTransform( double, Vector )}. */
	void AutoMatrixTest::testHouseholder () {
		const double
			dataV[] = { 3, 4 },
			dataA[] = {
				3, 5,
				4, 0
			};

		AutoVector v( 2 );
		v.fill( dataV );
		double tau = v.householderize();

		AutoMatrix a( 2, 2 );
		a.fill( dataA );
		a.householderTransform( tau, v );

		assert_true( "a1[0,0]", fabs( 5 - fabs( a.get( 0, 0 ) ) ) < 1e-8 );
		assert_true( "a1[1,0]", fabs( 0 - a.get( 1, 0 ) ) < 1e-8 );
		assert_true( "a1[0,1]", fabs( 3 - fabs( a.get( 0, 1 ) ) ) < 1e-8 );
		assert_true( "a1[1,1]", fabs( 4 - fabs( a.get( 1, 1 ) ) ) < 1e-8 );

		a.householderTransform( tau, v );

		assert_true( "a2[0,0]", fabs( 3 - a.get( 0, 0 ) ) < 1e-8 );
		assert_true( "a2[1,0]", fabs( 4 - a.get( 1, 0 ) ) < 1e-8 );
		assert_true( "a2[0,1]", fabs( 5 - a.get( 0, 1 ) ) < 1e-8 );
		assert_true( "a2[1,1]", fabs( 0 - a.get( 1, 1 ) ) < 1e-8 );
	}

	/** Test {@link Matrix::gemm} for non-transposed matrices. */
	void AutoMatrixTest::testGemmStraight () {
		AutoMatrix a( 3, 2 ), b( 2, 4 ), c( 3, 4 );
		const double dataA[] = {
			1, 8,
			2, 16,
			4, 32
		};
		a.fill( dataA );
		b.fill( 0.0 );
		c.gemm( 1.0, (const Matrix) a, false, (const Matrix) b, false, 0.0 );
		for ( int row = 0; row < b.countRows(); ++row ) {
			for ( int col = 0; col < b.countColumns(); ++col ) {
				char buf[256];
				sprintf( buf, "b[%d,%d]", row, col );
				assert_eq( buf, 0.0, b.get( row, col ) );
			}
		}
		for ( int row = 0; row < c.countRows(); ++row ) {
			for ( int col = 0; col < c.countColumns(); ++col ) {
				char buf[256];
				sprintf( buf, "alpha=1, beta=0, c[%d,%d]", row, col );
				assert_eq( buf, 0.0, c.get( row, col ) );
			}
		}
		b.set( 0, 0, 1.0 );
		b.set( 1, 0, -1.0 );
		c.gemm( 0.0, (const Matrix) a, false, (const Matrix) b, false, 1.0 );
		for ( int row = 0; row < c.countRows(); ++row ) {
			for ( int col = 0; col < c.countColumns(); ++col ) {
				char buf[256];
				sprintf( buf, "alpha=0, beta=1, c[%d,%d]", row, col );
				assert_eq( buf, 0.0, c.get( row, col ) );
			}
		}
		c.gemm( 1.0, (const Matrix) a, false, (const Matrix) b, false, 1.0 );
		assert_eq( "alpha=1, beta=1, c[0,0]", -7.0, c.get( 0, 0 ) );
		assert_eq( "alpha=1, beta=1, c[1,0]", -14.0, c.get( 1, 0 ) );
		assert_eq( "alpha=1, beta=1, c[2,0]", -28.0, c.get( 2, 0 ) );
		for ( int row = 0; row < c.countRows(); ++row ) {
			for ( int col = 1; col < c.countColumns(); ++col ) {
				char buf[256];
				sprintf( buf, "alpha=1, beta=1, c[%d,%d]", row, col );
				assert_eq( buf, 0.0, c.get( row, col ) );
			}
		}
		b.set( 1, 1, 2.0 );
		b.set( 0, 2, -1.0 );
		b.set( 1, 3, -1.0 );
		c.gemm( 101.0, (const Matrix) a, false, (const Matrix) b, false, -1.0 );
		assert_eq( "alpha=101, beta=-1, c[0,0]", -700.0, c.get( 0, 0 ) );
		assert_eq( "alpha=101, beta=-1, c[1,0]", -1400.0, c.get( 1, 0 ) );
		assert_eq( "alpha=101, beta=-1, c[2,0]", -2800.0, c.get( 2, 0 ) );
		assert_eq( "alpha=101, beta=-1, c[0,1]", 1616.0, c.get( 0, 1 ) );
		assert_eq( "alpha=101, beta=-1, c[1,1]", 3232.0, c.get( 1, 1 ) );
		assert_eq( "alpha=101, beta=-1, c[2,1]", 6464.0, c.get( 2, 1 ) );
		assert_eq( "alpha=101, beta=-1, c[0,2]", -101.0, c.get( 0, 2 ) );
		assert_eq( "alpha=101, beta=-1, c[1,2]", -202.0, c.get( 1, 2 ) );
		assert_eq( "alpha=101, beta=-1, c[2,2]", -404.0, c.get( 2, 2 ) );
		assert_eq( "alpha=101, beta=-1, c[0,3]", -808.0, c.get( 0, 3 ) );
		assert_eq( "alpha=101, beta=-1, c[1,3]", -1616.0, c.get( 1, 3 ) );
		assert_eq( "alpha=101, beta=-1, c[2,3]", -3232.0, c.get( 2, 3 ) );
	}

	/** Test {@link Matrix::gemm} for transposed matrices. */
	void AutoMatrixTest::testGemmTransposed () {
		AutoMatrix a( 3, 2 ), b( 3, 1 ), c( 1, 2 ), d( 2, 1 ), e( 3, 1 );
		const double dataA[] = {
			1, 8,
			2, 16,
			4, 32
		};
		const double dataB[] = {
			1,
			3,
			9
		};
		const double dataC[] = {
			2, -1
		};
		a.fill( dataA );
		b.fill( dataB );
		c.fill( dataC );
		d.gemm( 1.0, (const Matrix) a, true, (const Matrix) b, false, 0.0 );
		assert_eq( "alpha=1, beta=0, d[0,0]", 43.0, d.get( 0, 0 ) );
		assert_eq( "alpha=1, beta=0, d[1,0]", 344.0, d.get( 1, 0 ) );
		e.gemm( 1.0, (const Matrix) a, false, (const Matrix) c, true, 0.0 );
		assert_eq( "alpha=1, beta=0, e[0,0]", -6.0, e.get( 0, 0 ) );
		assert_eq( "alpha=1, beta=0, e[1,0]", -12.0, e.get( 1, 0 ) );
		assert_eq( "alpha=1, beta=0, e[2,0]", -24.0, e.get( 2, 0 ) );
	}

	/** Test {@link Matrix::factorizeQR} and {@link Matrix::extractQ}. */
	void AutoMatrixTest::testFactorizeQR () {
		const double dataA[] = {
			5, 3,
			0, 4,
			0, 0
		};
		AutoMatrix a( 3, 2 );
		AutoVector tau( 2 );
		AutoMatrix b( 3, 2 );

		a.fill( dataA );
		a.factorizeQR( tau );

		// right upper triangular part is R
		assert_eq( "a[0,0]", 5.0, a.get( 0, 0 ) );
		assert_eq( "a[0,1]", 3.0, a.get( 0, 1 ) );
		assert_eq( "a[1,1]", 4.0, a.get( 1, 1 ) );

		b.extractQ( tau, a );
		assert_eq( "b[0,0]", 1.0, b.get( 0, 0 ) );
		assert_eq( "b[1,0]", 0.0, b.get( 1, 0 ) );
		assert_eq( "b[2,0]", 0.0, b.get( 2, 0 ) );
		assert_eq( "b[0,1]", 0.0, b.get( 0, 1 ) );
		assert_eq( "b[1,1]", 1.0, b.get( 1, 1 ) );
		assert_eq( "b[2,1]", 0.0, b.get( 2, 1 ) );

		a.fill( dataA );
		a.set( 2, 1, 4.0 );
		a.factorizeQR( tau );

		// right upper triangular part is R
		assert_eq( "a[0,0]", 5.0, a.get( 0, 0 ) );
		assert_eq( "a[0,1]", 3.0, a.get( 0, 1 ) );
		assert_true( "a[1,1]", fabs( sqrt(32.0) - fabs( a.get( 1, 1 ) ) ) < 1e-8 );

		b.extractQ( tau, a );
		assert_eq( "b2[0,0]", 1.0, b.get( 0, 0 ) );
		assert_eq( "b2[1,0]", 0.0, b.get( 1, 0 ) );
		assert_eq( "b2[2,0]", 0.0, b.get( 2, 0 ) );
		assert_eq( "b2[0,1]", 0.0, b.get( 0, 1 ) );
		assert_true( "b2[1,1]", fabs( 1.0/sqrt(2.0) - fabs( b.get( 1, 1 ) ) ) < 1e-8 );
		assert_true( "b2[2,1]", fabs( 1.0/sqrt(2.0) - fabs( b.get( 2, 1 ) ) ) < 1e-8 );
		assert_true( "b2[1,1] - b[2,1]", fabs( b.get( 1, 1 ) - b.get( 2, 1 ) ) < 1e-8 );
	}

	/** Test {@link Matrix::factorizeQRP} and {@link Matrix::extractQ}. */
	void AutoMatrixTest::testFactorizeQRP () {
		const double dataA[] = {
			5, 3,
			0, 4,
			0, 0
		};
		AutoMatrix a( 3, 2 );
		AutoVector tau( 2 );
		AutoPermutation permutation( 2 );
		AutoMatrix b( 3, 2 );

		a.fill( dataA );
		a.factorizeQRP( tau, permutation );

		// permutation is identity
		assert_eq( "permutation[0]", 0, permutation.get( 0 ) );
		assert_eq( "permutation[1]", 1, permutation.get( 1 ) );

		// right upper triangular part is R
		assert_eq( "a[0,0]", 5.0, a.get( 0, 0 ) );
		assert_eq( "a[0,1]", 3.0, a.get( 0, 1 ) );
		assert_eq( "a[1,1]", 4.0, a.get( 1, 1 ) );

		b.extractQ( tau, a );
		assert_eq( "b[0,0]", 1.0, b.get( 0, 0 ) );
		assert_eq( "b[1,0]", 0.0, b.get( 1, 0 ) );
		assert_eq( "b[2,0]", 0.0, b.get( 2, 0 ) );
		assert_eq( "b[0,1]", 0.0, b.get( 0, 1 ) );
		assert_eq( "b[1,1]", 1.0, b.get( 1, 1 ) );
		assert_eq( "b[2,1]", 0.0, b.get( 2, 1 ) );

		a.fill( dataA );
		a.set( 2, 1, 12.0 );

		// original A columns needed for simulating Gram-Schmidt orthogonalisation
		AutoVector
			v( 3 ),
			w( 3 );
		v.copy( a.columnVector( 0 ) );
		w.copy( a.columnVector( 1 ) );

		a.factorizeQRP( tau, permutation );

		// right upper triangular part is R
		assert_eq( "a2[0,0]", 13.0, fabs( a.get( 0, 0 ) ) );
		// No easy expected values for a2[0,1] and a2[1,1]; simulating Gram-Schmidt orthogonalisation
		assert_true( "a2[0,1]", fabs( v.innerProduct( w ) / sqrt( w.sumSquares() ) - fabs( a.get( 0, 1 ) ) ) < 1e-8 );
		v.axpy( - v.innerProduct( w ) / w.sumSquares(), w );
		assert_true( "a2[1,1]", fabs( sqrt( v.sumSquares() ) - fabs( a.get( 1, 1 ) ) ) < 1e-8 );

		// permutation is exchange
		assert_eq( "permutation[0]", 1, permutation.get( 0 ) );
		assert_eq( "permutation[1]", 0, permutation.get( 1 ) );

		b.extractQ( tau, a );

		assert_true( "b2[0,0]", fabs( 3./13. - fabs( b.get( 0, 0 ) ) ) < 1e-8 );
		assert_true( "b2[1,0]", fabs( 4./13. - fabs( b.get( 1, 0 ) ) ) < 1e-8 );
		assert_true( "b2[2,0]", fabs( 12./13. - fabs( b.get( 2, 0 ) ) ) < 1e-8 );
		// No easy expected values for b2[*,1]; continue simulating Gram-Schmidt orthogonalisation
		v.scale( 1. / sqrt( v.sumSquares() ) );
		if ( v.get( 0 ) * b.get( 0, 1 ) < 0.0 ) {
			v.scale( -1. );
		}
		assert_true( "b2[0,1]", fabs( v.get( 0 ) - b.get( 0, 1 ) ) < 1e-8 );
		assert_true( "b2[1,1]", fabs( v.get( 1 ) - b.get( 1, 1 ) ) < 1e-8 );
		assert_true( "b2[2,1]", fabs( v.get( 2 ) - b.get( 2, 1 ) ) < 1e-8 );
	}

}
