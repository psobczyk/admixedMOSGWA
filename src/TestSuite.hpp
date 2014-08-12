/********************************************************************************
 *	This file is part of the MOSGWA program code.				*
 *	Copyright ©2012–2014, Bernhard Bodenstorfer.				*
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

#include <string>
#include <istream>
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

		/** Prepare temporary directory to help setting up a test scenario.
		* @param tmpDirNameTemplate indicates a wish for how the directory name should look.
		* It must contain the substring "XXXXXX" somewhere,
		* which will be replaced by some random characters.
		* @returns name of the newly created test directory.
		*/
		std::string createTmpDir ( const char* tmpDirNameTemplate );

		/** Remove temporary directory after testing is over.
		* @param tmpDirName names the directory to be deleted. It must be empty (again).
		*/
		void removeTmpDir ( const std::string &tmpDirName );

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

		/** Assert NaN, but not infinity. */
		void assert_nan (
			const std::string& msg,
			const double actual
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

		/** Assert equality of binary stream content with array of <code>char</code>. */
		void assertData (
			const std::string message,
			const size_t countBytes,
			const char * const data,
			std::istream& input
		);

		/** Assert equality of text stream content with array of lines. */
		void assertText (
			const std::string message,
			const size_t countLines,
			const char * const lines[],
			std::istream& input
		);
	};

}

// C++ compilers want template functions source code available when the template functions are used
#include "TestSuite.tpl"

#endif /* _TEST_SUITE_HPP_ */
