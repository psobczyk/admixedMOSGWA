#include "package.hpp"
#include "../TestSuite.hpp"

using namespace unitpp;
using namespace linalg;

namespace test {

	/** Tests the package commons. */
	class PackageTest : public TestSuite {

		private:

		/** Test function {@link upperPowerOf2}. */
		void testUpperPowerOf2 () {
			assert_eq( "0", 0, upperPowerOf2( 0 ) );
			assert_eq( "1", 1, upperPowerOf2( 1 ) );
			assert_eq( "2", 2, upperPowerOf2( 2 ) );
			assert_eq( "3", 4, upperPowerOf2( 3 ) );
			assert_eq( "4", 4, upperPowerOf2( 4 ) );
			assert_eq( "5", 8, upperPowerOf2( 5 ) );
			assert_eq( "6", 8, upperPowerOf2( 6 ) );
			assert_eq( "7", 8, upperPowerOf2( 7 ) );
			assert_eq( "8", 8, upperPowerOf2( 8 ) );
			assert_eq( "9", 16, upperPowerOf2( 9 ) );
			assert_eq( "15", 16, upperPowerOf2( 15 ) );
			assert_eq( "1023", 1024, upperPowerOf2( 1023 ) );
			assert_eq( "1024", 1024, upperPowerOf2( 1024 ) );
			assert_eq( "1025", 2048, upperPowerOf2( 1025 ) );
		}

		public:

		PackageTest () : TestSuite( "Test linalg package" ) {
			addTestMethod( "testUpperPowerOf2", this, &PackageTest::testUpperPowerOf2 );
		}
	} * packageTest = new PackageTest();	// automatically freed by unit++

}
