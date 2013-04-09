#ifndef _TEST_SUITE_TPL_
#define _TEST_SUITE_TPL_

#include "TestSuite.hpp"
#include <iostream>

namespace test {

	template<class C> void TestSuite::addTestMethod(
		const std::string& name,
		C* par,
		typename unitpp::test_mfun<C>::mfp fp,
		const char* file,
		unsigned int line
	) {
		add( name, unitpp::testcase( par, name, fp, file, line ) );
	}

	template <class T> void TestSuite::assertMatrixEqual_f (
		const char* file,
		const unsigned int line,
		const std::string &msg,
		const int x,
		const int y,
		const T * expected,
		const T * actual
	) {
		for ( int j = 0; j < y; ++j ) {
			for ( int i = 0; i < x; ++i ) {
				assert_eq(
					msg,
					expected[ x * j + i ],
					actual[ x * j + i ]
				);
			}
		}
	}

	template <class T> void TestSuite::printMatrix (
		const std::string &msg,
		const int x,
		const int y,
		const T * m
	) {
		if ( ! msg.empty() ) std::cout << msg << std::endl;
		for ( int j = 0; j < y; ++j ) {
			for ( int i = 0; i < x; ++i ) {
				std::cout << "\t" << m[ x * j + i ];
			}
			std::cout << std::endl;
		}
		std::cout << std::endl;
	}

}

#endif /* _TEST_SUITE_TPL_ */
