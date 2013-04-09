#include "GslMinimizerHolder.hpp"
#include "../TestSuite.hpp"

using namespace unitpp;
using namespace minimization;

namespace test {

	/** Tests the class {@link minimization::GslMinimizerHolder}. */
	class GslMinimizerHolderTest : public TestSuite {

		void test1dim ();

		public:

		GslMinimizerHolderTest () : TestSuite( "minimization::GslMinimizerHolderTest" ) {
			addTestMethod( "GslMinimizerHolderTest::test1dim", this, &GslMinimizerHolderTest::test1dim );
		}
	} * gslMinimizerHolderTest = new GslMinimizerHolderTest();	// automatically freed by unit++

	/** Test {@link GslMinimizerHolder} in 1-dimensional space.
	* Allocation of a <code>gsl_multimin_fdfminimizer*</code> does not work for 0-dimensional space.
	*/
	void GslMinimizerHolderTest::test1dim () {
		GslMinimizerHolder( 1 );
	}

}
