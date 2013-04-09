#include "TestSuite.hpp"
#include <cmath>
#include <cstdio>

using namespace std;
using namespace unitpp;

namespace test {

	TestSuite::TestSuite ( const string &title ) : suite ( title ) {
		suite::main().add( title, this );
	}

	void TestSuite::assert_close(
		const string& msg,
		const double expected,
		const double actual,
		const double tolerance
	) {
		char buffer[1024];
		sprintf(
			buffer,
			" expected close to %f and got %f, difference was %f, tolerance %f",
			expected,
			actual,
			actual-expected,
			tolerance
		);
		string xmsg = msg + buffer;
		assert_true( xmsg, fabs( expected - actual ) < tolerance );
	}

}
