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

#ifndef _TEST_SUITE_HPP_
#define _TEST_SUITE_HPP_

#include "TestSuite.hpp"
#include <string>
#include <unit++.h>

/** Contains unit tests.
* Install the package <a href="http://sourceforge.net/projects/unitpp/">unit++</a>
* to run unit tests.
*/
namespace test {

	/** Helpful extension for more convenient use of unit++. */
	class TestSuite : public unitpp::suite {

		public:

		/** Construct with a title for reports. */
		TestSuite ( const std::string &title );

		/** Add a method to the list of those to be run to perform the test. */
		template<class C> void addTestMethod(
			const std::string& name,
			C* par,
			typename unitpp::test_mfun<C>::mfp fp,
			const char* file = "",
			unsigned int line = 0
		);

		/** Assert proximity of double values. */
		void assert_close (
			const std::string& msg,
			const double expected,
			const double actual,
			const double tolerance = 1e-15
		);

		/** Assert equality of two rectangular arrays. */
		template <class T> void assertMatrixEqual_f (
			const char* file,
			const unsigned int line,
			const std::string &msg,
			const int x,
			const int y,
			const T * expected,
			const T * actual
		);

#define assertMatrixEqual(m, x, y, e, g) assertMatrixEqual_f(__FILE__, __LINE__, m, x, y, e, g)

		/** Print the elements of a rectangular array.
		* For debugging.
		*/
		template <class T> void printMatrix (
			const std::string &msg,
			const int x,
			const int y,
			const T * m
		);
	};

}

// C++ compilers want template functions source code available when the template functions are used
#include "TestSuite.tpl"

#endif /* _TEST_SUITE_HPP_ */
