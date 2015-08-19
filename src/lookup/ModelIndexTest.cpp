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

#include "ModelIndex.hpp"
#include "../TestSuite.hpp"
#include <stdio.h>
#include <set>
#include <iostream>

using namespace std;
using namespace unitpp;
using namespace lookup;

namespace test {

	/** Tests the class {@link lookup::ModelIndex}. */
	class ModelIndexTest : public TestSuite {

		private:

		/** Test constructor for empty model indices. */
		void testEmpty ();

		/** Test constructor. */
		void testConstruct ();

		/** Test convenience constructor. */
		void testConvenienceConstruct ();

		/** Test relational operators. */
		void testRelations ();

		/** Test containment-relations. */
		void testSetRelations ();

		/** Probe prototype construction in a somewhat more realistic setting. */
		void testComplicated ();

		public:

		ModelIndexTest () : TestSuite( "Test ModelIndex" ) {
			addTestMethod( "testEmpty", this, &ModelIndexTest::testEmpty );
			addTestMethod( "testConstruct", this, &ModelIndexTest::testConstruct );
			addTestMethod( "testConvenienceConstruct", this, &ModelIndexTest::testConvenienceConstruct );
			addTestMethod( "testRelations", this, &ModelIndexTest::testRelations );
			addTestMethod( "testSetRelations", this, &ModelIndexTest::testSetRelations );
			addTestMethod( "testComplicated", this, &ModelIndexTest::testComplicated );
		}
	} * modelIndexTest = new ModelIndexTest();	// automatically freed by unit++

	void ModelIndexTest::testEmpty () {
		set<size_t> s;
		const ModelIndex mi0, mi1( s );
		assert_eq( "empty equal each other", mi1, mi0 );
		assert_eq( "empty 0 size", 0, mi0.size() );
		assert_eq( "empty 1 size", 0, mi1.size() );
		for ( unsigned int i = 0; i < 65536; ++i ) {
			assert_true( "empty 0 does not contain anything", ! mi0.contains( i ) );
			assert_true( "empty 1 does not contain anything", ! mi1.contains( i ) );
		}
		assert_eq( "empty 0 iterators", mi0.end(), mi0.begin() );
		assert_eq( "empty 1 iterators", mi1.end(), mi1.begin() );
	}

	void ModelIndexTest::testConstruct () {
		set<size_t> s;
		size_t positions[] = { 2, 7, 9, 12, 3, 5, 9, 14 };
		for ( unsigned int i = 0; i < sizeof( positions ) / sizeof( size_t ); ++i ) {
			s.insert( positions[i] );
		}

		const ModelIndex mi1( s );
		assert_eq( "Size", 7, mi1.size() );

		ModelIndex::const_iterator iterator = mi1.begin();
		assert_eq( "Indices are sorted 0", 2, *iterator++ );
		assert_eq( "Indices are sorted 1", 3, *iterator++ );
		assert_eq( "Indices are sorted 2", 5, *iterator++ );
		assert_eq( "Indices are sorted 3", 7, *iterator++ );
		assert_eq( "Indices are sorted 4", 9, *iterator++ );
		assert_eq( "Indices are sorted 5", 12, *iterator++ );
		assert_eq( "Indices are sorted 6", 14, *iterator++ );
		assert_eq( "Indices are sorted end", mi1.end(), iterator );

		assert_false( "Does not contain 0", mi1.contains( 0 ) );
		assert_false( "Does not contain 1", mi1.contains( 1 ) );
		assert_true( "Contains 2", mi1.contains( 2 ) );
		assert_true( "Contains 3", mi1.contains( 3 ) );
		assert_false( "Does not contain 4", mi1.contains( 4 ) );
		assert_true( "Contains 5", mi1.contains( 5 ) );
		assert_false( "Does not contain 6", mi1.contains( 6 ) );
		assert_true( "Contains 7", mi1.contains( 7 ) );
		assert_false( "Does not contain 8", mi1.contains( 8 ) );
		assert_true( "Contains 9", mi1.contains( 9 ) );
		assert_false( "Does not contain 10", mi1.contains( 10 ) );
		assert_false( "Does not contain 11", mi1.contains( 11 ) );
		assert_true( "Contains 12", mi1.contains( 12 ) );
		assert_false( "Does not contain 13", mi1.contains( 13 ) );
		assert_true( "Contains 14", mi1.contains( 14 ) );
		assert_false( "Does not contain 15", mi1.contains( 15 ) );

		const ModelIndex mi1copy( mi1 );
		assert_eq( "Copy size", 7, mi1copy.size() );

		ModelIndex::const_iterator iterator1copy = mi1copy.begin();
		ModelIndex::const_iterator iterator1 = mi1.begin();

		for ( unsigned int i = 0; i < mi1.size(); ++i ) {
			char buf[256];
			sprintf( buf, "Element equality original and copy at %d", i );
			assert_eq( buf, *iterator1, *iterator1copy );
			++iterator1;
			++iterator1copy;
		}
		assert_eq( "Copy iterator at end", mi1copy.end(), iterator1copy );

		for ( size_t i = 0; i < 17; ++i ) {
			char buf[256];
			sprintf( buf, "Containment equality original and copy at %d", i );
			assert_eq( buf, mi1.contains( i ), mi1copy.contains( i ) );
		}

		assert_eq( "Copy equals original", (ModelIndex&)mi1, (ModelIndex&)mi1copy );
		// Since one of the two will have got the other's internal array
		assert_eq( "Copy still equals original", mi1, mi1copy );

		const ModelIndex mi1re( s );
		assert_eq( "Reconstructed equals original", (ModelIndex&)mi1, (ModelIndex&)mi1re );
		// Since one of the two will have got the other's internal array
		assert_eq( "Reconstructed still equals original", mi1, mi1re );
	}

	void ModelIndexTest::testConvenienceConstruct () {
		// Data with a few duplicates
		size_t positions[] = { 0, 18, 19, 24, 24, 23, 6, 23, 18 };
		vector<size_t> v;
		set<size_t> s;
		for ( unsigned int i = 0; i < sizeof( positions ) / sizeof( size_t ); ++i ) {
			s.insert( positions[i] );
			v.push_back( positions[i] );
		}
		const ModelIndex ms( s ), mv( v );
		assert_eq( "Equal from set as from vector", ms, mv );
		assert_eq( "Size s", 6, ms.size() );
		assert_eq( "Size v", 6, mv.size() );
	}

	void ModelIndexTest::testRelations () {
		set<size_t> s;

		const ModelIndex mi0( s );
		assert_true( "0 == 0", mi0 == mi0 );
		assert_false( "not 0 != 0", mi0 != mi0 );
		assert_false( "not 0 < 0", mi0 < mi0 );
		assert_false( "not 0 > 0", mi0 > mi0 );
		assert_true( "0 <= 0", mi0 <= mi0 );
		assert_true( "0 >= 0", mi0 >= mi0 );

		s.insert( 0 );

		const ModelIndex mi1( s );
		assert_false( "not 0 == 1", mi0 == mi1 );
		assert_false( "not 1 == 0", mi1 == mi0 );
		assert_true( "0 != 1", mi0 != mi1 );
		assert_true( "1 != 0", mi1 != mi0 );
		assert_true( "0 < 1", mi0 < mi1 );
		assert_false( "not 1 < 0", mi1 < mi0 );
		assert_false( "not 0 > 1", mi0 > mi1 );
		assert_true( "1 > 0", mi1 > mi0 );
		assert_true( "0 <= 1", mi0 <= mi1 );
		assert_false( "not 1 <= 0", mi1 <= mi0 );
		assert_false( "not 0 >= 1", mi0 >= mi1 );
		assert_true( "1 >= 0", mi1 >= mi0 );

		const ModelIndex mi0plus( mi0, 0 );
		assert_eq( "Size of 0+{0}", 1, mi0plus.size() );
		assert_eq( "0+{0} == {0}", mi1, mi0plus );

		const ModelIndex mi1minus( mi1, 0 );
		assert_eq( "Size of {0}-{0}", 0, mi1minus.size() );
		assert_eq( "{0}-{0} == 0", mi0, mi1minus );

		s.erase( 0 );
		s.insert( 1 );

		const ModelIndex mi2( s );;
		assert_false( "not 1 == 2", mi1 == mi2 );
		assert_false( "not 2 == 1", mi2 == mi1 );
		assert_true( "1 != 2", mi1 != mi2 );
		assert_true( "2 != 1", mi2 != mi1 );
		assert_true( "1 < 2", mi1 < mi2 );
		assert_false( "not 2 < 1", mi2 < mi1 );
		assert_false( "not 1 > 2", mi1 > mi2 );
		assert_true( "2 > 1", mi2 > mi1 );
		assert_true( "1 <= 2", mi1 <= mi2 );
		assert_false( "not 2 <= 1", mi2 <= mi1 );
		assert_false( "not 1 >= 2", mi1 >= mi2 );
		assert_true( "2 >= 1", mi2 >= mi1 );

		s.insert( 3 );
		const ModelIndex mi3( s );;
		assert_false( "not 2 == 3", mi2 == mi3 );
		assert_false( "not 3 == 2", mi3 == mi2 );
		assert_true( "2 != 3", mi2 != mi3 );
		assert_true( "3 != 2", mi3 != mi2 );
		assert_true( "2 < 3", mi2 < mi3 );
		assert_false( "not 3 < 2", mi3 < mi2 );
		assert_false( "not 2 > 3", mi2 > mi3 );
		assert_true( "3 > 2", mi3 > mi2 );
		assert_true( "2 <= 3", mi2 <= mi3 );
		assert_false( "not 3 <= 2", mi3 <= mi2 );
		assert_false( "not 2 >= 3", mi2 >= mi3 );
		assert_true( "3 >= 2", mi3 >= mi2 );

		const ModelIndex mi2plus( mi2, 3 );
		assert_eq( "2+{3}", 2, mi2plus.size() );
		assert_eq( "2+{3} == 3", mi3, mi2plus );

		const ModelIndex mi2minus( mi2, 1 );
		assert_eq( "2-{1}", 0, mi2minus.size() );
		assert_eq( "2-{1} == 0", mi0, mi2minus );
	}

	void ModelIndexTest::testSetRelations () {
		set<size_t> s;

		const ModelIndex mi0( s );
		assert_true( "0 subset 0", mi0.isSubsetOf( mi0 ) );
		assert_true( "0 superset 0", mi0.isSupersetOf( mi0 ) );

		s.insert( 0 );
		s.insert( 4 );
		const ModelIndex mi1( s );
		assert_true( "0 subset 1", mi0.isSubsetOf( mi1 ) );
		assert_false( "not 1 subset 0", mi1.isSubsetOf( mi0 ) );
		assert_false( "not 0 superset 1", mi0.isSupersetOf( mi1 ) );
		assert_true( "1 superset 0", mi1.isSupersetOf( mi0 ) );

		s.insert( 3 );
		const ModelIndex mi2( s );
		assert_true( "1 subset 2", mi1.isSubsetOf( mi2 ) );
		assert_false( "not 2 subset 1", mi2.isSubsetOf( mi1 ) );
		assert_false( "not 1 superset 2", mi1.isSupersetOf( mi2 ) );
		assert_true( "2 superset 1", mi2.isSupersetOf( mi1 ) );

		s.erase( 0 );
		const ModelIndex mi3( s );
		assert_true( "3 subset 2", mi3.isSubsetOf( mi2 ) );
		assert_false( "not 2 subset 3", mi2.isSubsetOf( mi3 ) );
		assert_false( "not 3 superset 2", mi3.isSupersetOf( mi2 ) );
		assert_true( "2 superset 3", mi2.isSupersetOf( mi3 ) );
		assert_false( "not 3 subset 1", mi3.isSubsetOf( mi1 ) );
		assert_false( "not 1 subset 3", mi1.isSubsetOf( mi3 ) );
		assert_false( "not 3 superset 1", mi3.isSupersetOf( mi1 ) );
		assert_false( "not 1 superset 3", mi1.isSupersetOf( mi3 ) );
	}

	void ModelIndexTest::testComplicated () {
		set<size_t> s;
		size_t positions[] = { 0, 4, 7, 8, 20, 21 };
		for ( unsigned int i = 0; i < sizeof( positions ) / sizeof( size_t ); ++i ) {
			s.insert( positions[i] );
		}

		const ModelIndex mi1( s );
		assert_eq( "Size", 6, mi1.size() );

		s.insert( 43 );
		s.insert( 44 );
		s.insert( 30 );
		s.insert( 19 );
		const ModelIndex mi2( s );
		const ModelIndex mi2prime(
			ModelIndex(
				ModelIndex(
					ModelIndex(
						mi1, 43
					), 44
				), 19
			), 30
		);
		for ( size_t i = 0; i < 50; ++i ) {
			char buf[256];
			sprintf( buf, "Containment equality 2 and 2' at %d", i );
			assert_eq( buf, mi2.contains( i ), mi2prime.contains( i ) );
		}

		s.erase( 7 );
		s.erase( 8 );
		const ModelIndex mi3( s );
		const ModelIndex mi3prime(
			ModelIndex(
				mi2, 8
			), 7
		);
		for ( size_t i = 0; i < 50; ++i ) {
			char buf[256];
			sprintf( buf, "Containment equality 3 and 3' at %d", i );
			assert_eq( buf, mi3.contains( i ), mi3prime.contains( i ) );
		}

		assert_eq( "2 == 2'", mi2, mi2prime );
		assert_eq( "3 == 3'", mi3, mi3prime );
	}

}
