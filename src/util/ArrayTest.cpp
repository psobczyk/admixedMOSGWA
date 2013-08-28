/********************************************************************************
 *	This file is part of the MOSGWA program code.				*
 *	Copyright ©2013, Bernhard Bodenstorfer.					*
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

#include "Array.hpp"
#include "../TestSuite.hpp"
#include <cstdarg>

using namespace unitpp;
using namespace util;

namespace test {

	/** Tests the template class {@link util::Array}. */
	struct ArrayTest : public TestSuite {

		ArrayTest ();

		void testElliptic ();
		void testVaList ();
		void testStruct ();

	} * arrayTest = new ArrayTest();	// automatically freed by unit++

	ArrayTest::ArrayTest () : TestSuite( "util::Array Test" ) {
		addTestMethod( "ArrayTest::testElliptic", this, &ArrayTest::testElliptic );
		addTestMethod( "ArrayTest::testVaList", this, &ArrayTest::testVaList );
		addTestMethod( "ArrayTest::testStruct", this, &ArrayTest::testStruct );
	}

	void ArrayTest::testElliptic () {
		{
			Array<int,0> array0;
			assert_eq( "zero size", 0, sizeof( array0 ) );
		}

		{
			Array<int,1> array1( 17 );
			assert_eq( "int size", sizeof( int ), sizeof( array1 ) );
			assert_eq( "int initialised", 17, array1[0] );
		}

		{
			Array<double,2> array2( 13.6, 11.2 );
			assert_eq( "double size", sizeof( double[2] ), sizeof( array2 ) );
			assert_eq( "double initialised [0]", 13.6, array2[0] );
			assert_eq( "double initialised [1]", 11.2, array2[1] );
		}
	}

	void ellipsor ( int value, ... ) {
		va_list arglist;
		va_start( arglist, value );
		Array<int,0> array0( arglist );
		Array<int,1> array1( arglist );
		va_end( arglist );
		assert_eq( "zero size", 0, sizeof( array0 ) );
		assert_eq( "int size", sizeof( int ), sizeof( array1 ) );
		assert_eq( "first from variable argument list", value, array1[0] );
	}

	void ArrayTest::testVaList () {
		ellipsor( 17, 17 );
	}

	struct DD { double a, b; };

	void ArrayTest::testStruct () {
		struct DD
			x = { 1.0, 2.0 },
			y = { 3.0, 4.0 },
			z = { 5.0, 6.0 };
		Array<DD,3> array3( x, y, z );
		assert_eq( "triple size", sizeof( DD[3] ), sizeof( array3 ) );
		assert_eq( "triple initialised [0].a", 1.0, array3[0].a );
		assert_eq( "triple initialised [0].b", 2.0, array3[0].b );
		assert_eq( "triple initialised [1].a", 3.0, array3[1].a );
		assert_eq( "triple initialised [1].b", 4.0, array3[1].b );
		assert_eq( "triple initialised [2].a", 5.0, array3[2].a );
		assert_eq( "triple initialised [2].b", 6.0, array3[2].b );
	}

}
