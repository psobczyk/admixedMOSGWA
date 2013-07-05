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
