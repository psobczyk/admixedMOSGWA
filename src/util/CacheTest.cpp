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

#include "Cache.hpp"
#include "../TestSuite.hpp"

using namespace unitpp;
using namespace util;

namespace test {

	/** Helps testing the template class {@link util::Cache}. */
	template <class Key, class Value> struct TestRetriever : public Cache<Key,Value>::Retriever {

		/** Store last queried key which was not in the cache.
		* Declaration does not use "Key restrieved;", because Key might be a const type.
		* Anyway, it's only needed for int values.
		*/
		int retrieved;

		/** Initialise the factory. */
		TestRetriever ();

		/** Create the asked-for object. */
		virtual Value retrieve ( const Key& key );
	};

	template <class Key, class Value> TestRetriever<Key,Value>::TestRetriever () : retrieved( 0 ) {}

	template <class Key, class Value> Value TestRetriever<Key,Value>::retrieve ( const Key& key ) {
		retrieved = key;
		return -key;
	}

	/** Tests the template class {@link util::Cache}. */
	struct CacheTest : public TestSuite {

		CacheTest ();
		void testConstCache ();
		void testVarCache ();

	} * cacheTest = new CacheTest();	// automatically freed by unit++

	CacheTest::CacheTest () : TestSuite( "util::Cache Test" ) {
		addTestMethod( "CacheTest::testConstCache", this, &CacheTest::testConstCache );
		addTestMethod( "CacheTest::testVarCache", this, &CacheTest::testVarCache );
	}

	/** Test {@link util::Cache} for const types. */
	void CacheTest::testConstCache () {
		TestRetriever<const int, const double> retriever;
		Cache<const int, const double> cache;

		// query from empty cache implies retrieve
		assert_eq( "initial test cache", 0, retriever.retrieved );
		const double * const me0 = &cache.get( 1, retriever );
		assert_eq( "retrieve[1]", -1, *me0 );
		assert_eq( "retrieved 1", 1, retriever.retrieved );

		// repeated query from cache implies no retrieve
		retriever.retrieved = 0;
		const double * const me1 = &cache.get( 1, retriever );
		assert_eq( "cached[1]", -1, *me1 );
		assert_eq( "same cached[1]", me0, me1 );
		assert_eq( "cached 1", 0, retriever.retrieved );

		// new query from cache implies retrieve
		const double * const mz0 = &cache.get( 20, retriever );
		assert_eq( "retrieve[20]", -20, *mz0 );
		assert_eq( "retrieved 20", 20, retriever.retrieved );

		// again repeated query from cache implies no retrieve
		retriever.retrieved = 0;
		const double * const me2 = &cache.get( 1, retriever );
		assert_eq( "again cached[1]", -1, *me2 );
		assert_eq( "again same cached[1]", me0, me2 );
		assert_eq( "again cached 1", 0, retriever.retrieved );

		// repeated new query from cache implies no retrieve
		retriever.retrieved = 0;
		const double * const mz1 = &cache.get( 20, retriever );
		assert_eq( "cached[20]", -20, *mz1 );
		assert_eq( "same cached[20]", mz0, mz1 );
		assert_eq( "cached 20", 0, retriever.retrieved );
	}

	/** Test {@link util::Cache} for variable types. */
	void CacheTest::testVarCache () {
		TestRetriever<int,double> retriever;
		Cache<int,double> cache;

		// query from empty cache implies retrieve
		assert_eq( "initial test cache", 0, retriever.retrieved );
		double * const me0 = &cache.get( 1, retriever );
		assert_eq( "retrieve[1]", -1, *me0 );
		assert_eq( "retrieved 1", 1, retriever.retrieved );

		// retrieved object is modifiable
		*me0 = 100;

		// repeated query from cache implies no retrieve
		retriever.retrieved = 0;
		const double * const me1 = &cache.get( 1, retriever );
		assert_eq( "modified cached[1]", 100, *me1 );
		assert_eq( "but same cached[1]", me0, me1 );
		assert_eq( "cached 1", 0, retriever.retrieved );

		// retrieved lvalue
		cache.get( 20, retriever ) = 200;
		assert_eq( "retrieved 20", 20, retriever.retrieved );

		// again repeated query from cache implies no retrieve
		retriever.retrieved = 0;
		const double * const me2 = &cache.get( 1, retriever );
		assert_eq( "again modified cached[1]", 100, *me2 );
		assert_eq( "again same cached[1]", me0, me2 );
		assert_eq( "again cached 1", 0, retriever.retrieved );

		// repeated new query from cache implies no retrieve
		retriever.retrieved = 0;
		assert_eq( "modified cached[20]", 200, cache.get( 20, retriever ) );
		assert_eq( "cached 20", 0, retriever.retrieved );
	}

}
