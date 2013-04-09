#include "CachingValuation.hpp"
#include "TestValuation.hpp"
#include "../TestSuite.hpp"
#include <cmath>

using namespace test;
using namespace lookup;
using namespace valuation;
using namespace unitpp;

namespace test {

	/** Tests the class {@link CachingValuation}. */
	class CachingValuationTest : public TestSuite {

		/** Test  methods. */
		void testCache ();

		public:

		CachingValuationTest () : TestSuite( "CachingValuationTest" ) {
			addTestMethod( "CachingValuationTest::testCache", this, &CachingValuationTest::testCache );
		}
	} * cachingValuationTest = new CachingValuationTest();	// automatically freed by unit++

	void CachingValuationTest::testCache () {
		const ModelIndex
			mi0,
			mi1( mi0, 1 ),
			mi2( mi0, 2 );

		TestValuation testValuation;
		CachingValuation cachingValuation( testValuation );

		{
			const double nan = cachingValuation.valuate( mi0 );
			assert_true( "Empty cache yields NaN.", isnan( nan ) );
			assert_true( "Empty cache does not yield infinity.", ! isinf( nan ) );
		}

		testValuation.put( mi0, 1.0 );
		{
			assert_eq( "One is on.", 1.0, testValuation.valuate( mi0 ) );
			const double nan0 = cachingValuation.valuate( mi0 );
			assert_true( "Second query uses cached NaN.", isnan( nan0 ) );
			assert_true( "Second query uses cached NaN, not infinity.", ! isinf( nan0 ) );
			const double nan1 = cachingValuation.valuate( mi1 );
			assert_true( "Not in cache yields NaN.", isnan( nan1 ) );
			assert_true( "Not in cache does not yield infinity.", ! isinf( nan1 ) );
		}

		testValuation.put( mi2, 2.0 );
		{
			assert_eq( "Two is on.", 2.0, testValuation.valuate( mi2 ) );
			assert_eq( "Query retrieves value.", 2.0, cachingValuation.valuate( mi2 ) );
		}

		testValuation.put( mi2, 3.0 );
		{
			assert_eq( "Three is on.", 3.0, testValuation.valuate( mi2 ) );
			assert_eq( "Query uses cached previous value.", 2.0, cachingValuation.valuate( mi2 ) );
		}
	}

}
