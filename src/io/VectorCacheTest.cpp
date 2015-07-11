/********************************************************************************
 *	This file is part of the MOSGWA program code.				*
 *	Copyright ©2015, Bernhard Bodenstorfer.				*
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

#include "VectorCache.hpp"
#include "../linalg/AutoVector.hpp"
#include "../TestSuite.hpp"

using namespace unitpp;
using namespace linalg;
using namespace io;

namespace test {

	/** Tests the class {@link linalg::VectorCache} and thereby {@link linalg::Matrix}. */
	class VectorCacheTest : public TestSuite {

		void testLimitZero ();
		void testLimitOne ();
		void testLimitTwo ();

		public:

		VectorCacheTest () : TestSuite( "linalg::VectorCacheTest" ) {
			addTestMethod( "VectorCacheTest::testLimitZero", this, &VectorCacheTest::testLimitZero );
			addTestMethod( "VectorCacheTest::testLimitOne", this, &VectorCacheTest::testLimitOne );
			addTestMethod( "VectorCacheTest::testLimitTwo", this, &VectorCacheTest::testLimitTwo );
		}
	} * vectorCacheTest = new VectorCacheTest();	// automatically freed by unit++

	/** Test {@link VectorCache} without storage, i.e. limit zero. */
	void VectorCacheTest::testLimitZero () {
		const double data[] = { 1.0, 2.0 };
		const size_t dim = sizeof( data ) / sizeof( data[0] );
		AutoVector v( dim ), w( dim );
		VectorCache cache( dim, 0 );

		v.fill( data );
		w.fill( -9.0 );

		assert_false( "Zero cache initially empty", cache.retrieve( 0, w ) );
		assert_eq( "Zero cache w[0]", -9.0, w.get( 0 ) );
		assert_eq( "Zero cache w[1]", -9.0, w.get( 1 ) );

		cache.store( 0, v );

		assert_false( "Zero cache full is empty", cache.retrieve( 0, w ) );
		assert_eq( "Zero cache still w[0]", -9.0, w.get( 0 ) );
		assert_eq( "Zero cache still w[1]", -9.0, w.get( 1 ) );
	}

	/** Test {@link VectorCache} with one slot storage, i.e. limit one. */
	void VectorCacheTest::testLimitOne () {
		const double data[] = { 1.0, 2.0 };
		const size_t dim = sizeof( data ) / sizeof( data[0] );
		AutoVector v( dim ), w( dim );
		VectorCache cache( dim, 1 );

		v.fill( data );
		w.fill( 0.0 );

		assert_false( "One cache initially empty", cache.retrieve( 0, w ) );
		assert_eq( "One cache w[0]", 0.0, w.get( 0 ) );
		assert_eq( "One cache w[1]", 0.0, w.get( 1 ) );

		cache.store( 2, v );

		assert_false( "One cache one miss 0", cache.retrieve( 0, w ) );
		assert_eq( "One cache miss 0 w[0]", 0.0, w.get( 0 ) );
		assert_eq( "One cache miss 0 w[1]", 0.0, w.get( 1 ) );

		assert_false( "One cache one miss 1", cache.retrieve( 1, w ) );
		assert_eq( "One cache miss 1 w[0]", 0.0, w.get( 0 ) );
		assert_eq( "One cache miss 1 w[1]", 0.0, w.get( 1 ) );

		assert_true( "One cache one hit 2", cache.retrieve( 2, w ) );
		assert_eq( "One cache hit 2 w[0]", 1.0, w.get( 0 ) );
		assert_eq( "One cache hit 2 w[1]", 2.0, w.get( 1 ) );

		v.set( 0, 5.0 );
		v.set( 1, 4.0 );
		cache.store( 4, v );

		assert_false( "One cache still miss 0", cache.retrieve( 0, w ) );
		assert_eq( "One cache still miss 0 w[0]", 1.0, w.get( 0 ) );
		assert_eq( "One cache still miss 0 w[1]", 2.0, w.get( 1 ) );

		assert_false( "One cache miss overwritten 2", cache.retrieve( 2, w ) );
		assert_eq( "One cache miss overwritten 2 w[0]", 1.0, w.get( 0 ) );
		assert_eq( "One cache miss overwritten 2 w[1]", 2.0, w.get( 1 ) );

		assert_true( "One cache now hit 4", cache.retrieve( 4, w ) );
		assert_eq( "One cache now hit 4 w[0]", 5.0, w.get( 0 ) );
		assert_eq( "One cache now hit 4 w[1]", 4.0, w.get( 1 ) );

		v.set( 0, 19.0 );
		v.set( 1, 13.0 );
		cache.store( 4, v );

		assert_true( "One cache hit changed 4", cache.retrieve( 4, w ) );
		assert_eq( "One cache now hit 4 w[0]", 19.0, w.get( 0 ) );
		assert_eq( "One cache now hit 4 w[1]", 13.0, w.get( 1 ) );
	}

	/** Test {@link VectorCache} with one slot storage, i.e. limit one. */
	void VectorCacheTest::testLimitTwo () {
		const size_t dim = 1;
		AutoVector v( dim ), w( dim );
		VectorCache cache( dim, 2 );

		v.fill( 17.0 );
		w.fill( -4.0 );

		assert_false( "Two cache initially empty", cache.retrieve( 0, w ) );
		assert_eq( "Two cache w[0]", -4.0, w.get( 0 ) );

		cache.store( 3, v );

		assert_false( "Two cache one miss 0", cache.retrieve( 0, w ) );
		assert_eq( "Two cache miss 0 w[0]", -4.0, w.get( 0 ) );

		assert_true( "Two cache one hit 3", cache.retrieve( 3, w ) );
		assert_eq( "Two cache hit 3 w[0]", 17.0, w.get( 0 ) );

		v.set( 0, 27.0 );
		cache.store( 1, v );
		w.fill( -4.0 );

		assert_false( "Two cache still miss 0", cache.retrieve( 0, w ) );
		assert_eq( "Two cache still miss 0 w[0]", -4.0, w.get( 0 ) );

		assert_true( "Two cache newly hit 1", cache.retrieve( 1, w ) );
		assert_eq( "Two cache newly hit 1 w[0]", 27.0, w.get( 0 ) );

		assert_true( "Two cache still hit 3", cache.retrieve( 3, w ) );
		assert_eq( "Two cache still hit 3 w[0]", 17.0, w.get( 0 ) );

		v.set( 0, 254.0 );
		cache.store( 1, v );
		w.fill( -22.0 );

		assert_false( "Two cache again still miss 0", cache.retrieve( 0, w ) );
		assert_eq( "Two cache again still miss 0 w[0]", -22.0, w.get( 0 ) );

		assert_true( "Two cache hit changed 1", cache.retrieve( 1, w ) );
		assert_eq( "Two cache hit changed 1 w[0]", 254.0, w.get( 0 ) );

		assert_true( "Two cache hit unchanged 3", cache.retrieve( 3, w ) );
		assert_eq( "Two cache hit unchanged 3 w[0]", 17.0, w.get( 0 ) );

		v.set( 0, 1593.0 );
		cache.store( 0, v );
		w.fill( -12.19 );

		assert_true( "Two cache now hit 0", cache.retrieve( 0, w ) );
		assert_eq( "Two cache now hit 0 w[0]", 1593.0, w.get( 0 ) );

		assert_false( "Two cache hit miss 1", cache.retrieve( 1, w ) );
		assert_eq( "Two cache hit miss 1 w[0]", 1593.0, w.get( 0 ) );

		assert_true( "Two cache still hit unchanged 3", cache.retrieve( 3, w ) );
		assert_eq( "Two cache still hit unchanged 3 w[0]", 17.0, w.get( 0 ) );

		v.set( 0, 192.0 );
		cache.store( 1, v );
		w.fill( -13.6 );

		assert_true( "Two cache now still hit 0", cache.retrieve( 0, w ) );
		assert_eq( "Two cache now still hit 0 w[0]", 1593.0, w.get( 0 ) );

		assert_true( "Two cache hit again 1", cache.retrieve( 1, w ) );
		assert_eq( "Two cache hit again 1 w[0]", 192.0, w.get( 0 ) );

		assert_false( "Two cache now miss 3", cache.retrieve( 3, w ) );
		assert_eq( "Two cache now miss 3 w[0]", 192.0, w.get( 0 ) );
	}

}
