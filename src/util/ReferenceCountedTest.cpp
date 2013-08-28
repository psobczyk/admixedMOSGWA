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

#include "ReferenceCounted.hpp"
#include "../TestSuite.hpp"

using namespace unitpp;
using namespace util;

namespace test {

	/** Tests the template class {@link util::ReferenceCounted}. */
	struct ReferenceCountedTest : public TestSuite {

		ReferenceCountedTest ();
		void testPrimitive ();
		void testDestructible ();
		void testCounting ();

	} * referenceCountedTest = new ReferenceCountedTest();	// automatically freed by unit++

	ReferenceCountedTest::ReferenceCountedTest () : TestSuite( "util::ReferenceCounted Test" ) {
		addTestMethod( "ReferenceCountedTest::testPrimitive", this, &ReferenceCountedTest::testPrimitive );
		addTestMethod( "ReferenceCountedTest::testDestructible", this, &ReferenceCountedTest::testDestructible );
		addTestMethod( "ReferenceCountedTest::testCounting", this, &ReferenceCountedTest::testCounting );
	}

	struct TestDestructible {

		/** Counts constructor calls. */
		static unsigned constructed;

		/** Counts destructor calls. */
		static unsigned destructed;

		/** Constructor indicates its calling via class counter. */
		TestDestructible ();

		/** Constructor indicates its calling via class counter. */
		TestDestructible ( const TestDestructible& original );

		/** Destructor indicates its calling via class counter. */
		~TestDestructible ();
	};

	/** C++ needs the declared static variable to be defined somewhere. */
	unsigned TestDestructible::constructed;
	unsigned TestDestructible::destructed;

	TestDestructible::TestDestructible () {
		++TestDestructible::constructed;
	}

	TestDestructible::TestDestructible ( const TestDestructible& original ) {
		++TestDestructible::constructed;
	}

	TestDestructible::~TestDestructible () {
		++TestDestructible::destructed;
	}

	/** Test {@link util::ReferenceCounted} for primitive types. */
	void ReferenceCountedTest::testPrimitive () {
		ReferenceCounted<int> reference( 100 );
		int& ref( reference );
		assert_eq( "initial value", 100, ref );
		ref = 1000;
		int& modifiedRef( reference );
		assert_eq( "modified value", 1000, modifiedRef );

		ReferenceCounted<const int> constReference( reference );
		const int& constRef( constReference );
		assert_eq( "copz-constructed value", 1000, constRef );
	}

	/** Test {@link util::ReferenceCounted} for types requiring destruction. */
	void ReferenceCountedTest::testDestructible () {
		TestDestructible::constructed = 0;
		TestDestructible::destructed = 0;

		TestDestructible
			testDestructible,
			otherDestructible;
		assert_eq( "Test objects have been constructed", 2, TestDestructible::constructed );

		{
			ReferenceCounted<TestDestructible> reference( testDestructible );
			assert_eq( "Internally referenced object has been constructed", 3, TestDestructible::constructed );

			TestDestructible& extractReference( reference );
			assert_eq( "Extract counted reference does not construct", 3, TestDestructible::constructed );

			{
				ReferenceCounted<TestDestructible> copiedReference( reference );
				assert_eq( "Copy counted reference does not construct", 3, TestDestructible::constructed );

				ReferenceCounted<TestDestructible> otherReference( otherDestructible );
				assert_eq( "One more internal constructed", 4, TestDestructible::constructed );

				ReferenceCounted<TestDestructible> otherCopy( otherReference );
				assert_eq( "Once more, copy does not construct", 4, TestDestructible::constructed );

				assert_eq( "None destructed", 0, TestDestructible::destructed );

				otherCopy = copiedReference;
				assert_eq( "Assignment did not construct", 4, TestDestructible::constructed );
				assert_eq( "Assignment did not destruct", 0, TestDestructible::destructed );

				otherReference = reference;
				assert_eq( "Second assignment did not construct", 4, TestDestructible::constructed );
				assert_eq( "Second assignment did destruct", 1, TestDestructible::destructed );
			}
			assert_eq( "Delete counted references does not construct", 4, TestDestructible::constructed );
			assert_eq( "Delete copied counted references did not destruct", 1, TestDestructible::destructed );
		}
		assert_eq( "No more object has been constructed", 4, TestDestructible::constructed );
		assert_eq( "Internally referenced one has been destructed", 2, TestDestructible::destructed );
	}

	/** Test {@link util::ReferenceCounted} reference counting. */
	void ReferenceCountedTest::testCounting () {
		TestDestructible::constructed = 0;
		TestDestructible::destructed = 0;

		TestDestructible testDestructible;
		assert_eq( "Test object has been constructed", 1, TestDestructible::constructed );

		ReferenceCounted<TestDestructible>
			reference0( testDestructible ),
			reference1( testDestructible ),
			reference2( testDestructible ),
			reference3( testDestructible ),
			reference4( testDestructible ),
			copy0( reference0 );
		assert_eq( "Internally referenced objects have been constructed", 6, TestDestructible::constructed );
		assert_eq( "Construction does not destruct", 0, TestDestructible::destructed );

		reference1 = reference0;
		assert_eq( "First assignment has destructed", 1, TestDestructible::destructed );

		reference0 = reference2;
		reference1 = reference2;
		assert_eq( "More assignments have not destructed due to copy", 1, TestDestructible::destructed );

		reference2 = reference3;
		reference1 = reference3;
		assert_eq( "Further assignments have not destructed due to reference0", 1, TestDestructible::destructed );
		reference0 = reference4;
		assert_eq( "Losing destructs original reference2", 2, TestDestructible::destructed );

		reference3 = reference3;
		reference2 = reference3;
		reference1 = reference3;
		assert_eq( "Assignments of same do not construct", 6, TestDestructible::constructed );
		assert_eq( "Assignments of same do not destructed", 2, TestDestructible::destructed );

		copy0 = reference0;
		assert_eq( "Losing destructs original reference0", 3, TestDestructible::destructed );

		TestDestructible
			&ref0( reference0 ),
			&ref1( reference1 ),
			&ref2( reference2 ),
			&ref3( reference3 ),
			&ref4( reference4 ),
			&cp0( copy0 );

		assert_eq( "Ref 1 eq 2", &ref1, &ref2 );
		assert_eq( "Ref 1 eq 3", &ref1, &ref3 );
		assert_eq( "Ref 0 eq 4", &ref0, &ref4 );
		assert_eq( "Ref 0 eq copy 0", &ref0, &cp0 );
		assert_true( "Ref 0 ne 1", &ref0 != &ref1 );

		reference1 = reference0;
		reference2 = reference0;
		assert_eq( "So far no more destructed", 3, TestDestructible::destructed );
		reference3 = reference0;
		assert_eq( "Now one but last internal destructed", 4, TestDestructible::destructed );

		TestDestructible
			&newRef1( reference1 ),
			&newRef2( reference2 ),
			&newRef3( reference3 );
		assert_eq( "Ref 0 eq 1", &ref0, &newRef1 );
		assert_eq( "Ref 0 eq 2", &ref0, &newRef2 );
		assert_eq( "Ref 0 eq 3", &ref0, &newRef3 );
	}

}
