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
