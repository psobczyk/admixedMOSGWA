#include "Exception.hpp"
#include "TestSuite.hpp"

using namespace std;
using namespace unitpp;

namespace test {

	/** Tests the class {@link Exception}. */
	class ExceptionTest : public TestSuite {

		/** Test {@link Exception::Exception} and {@link Exception::getMessage}. */
		void testMessaging ();

		public:

		ExceptionTest () : TestSuite( "ExceptionTest" ) {
			addTestMethod( "ExceptionTest::testMessaging", this, &ExceptionTest::testMessaging );
		}
	} * exceptionTest = new ExceptionTest();	// automatically freed by unit++

	void ExceptionTest::testMessaging () {
		Exception exception( "Want to report error %d in file \"%s\".", 42, "ExceptionTest" );
		string message = exception.what();
		assert_eq( "Exception message", "Want to report error 42 in file \"ExceptionTest\".", message );
	}
}
