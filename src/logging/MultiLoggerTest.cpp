/********************************************************************************
 *	This file is part of the MOSGWA program code.				*
 *	Copyright ©2015, Bernhard Bodenstorfer.					*
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

#include "MultiLogger.hpp"
#include "TestLogger.hpp"

using namespace std;
using namespace logging;
using namespace unitpp;

namespace test {

	/** Tests the class {@link logging::MultiLogger}. */
	class MultiLoggerTest : public TestSuite {

		/** Test logging. */
		void testLog ();

		public:

		/** Construct the test object. */
		MultiLoggerTest ();

	} * multiLoggerTest = new MultiLoggerTest();	// automatically freed by unit++

	MultiLoggerTest::MultiLoggerTest () : TestSuite( "MultiLoggerTest" ) {
		addTestMethod( "MultiLoggerTest::testLog", this, &MultiLoggerTest::testLog );
	}

	void MultiLoggerTest::testLog () {
		TestLogger la, lb;
		MultiLogger logger;

		assert_eq( "Initial logger threshold INFO", Logger::INFO, logger.getThreshold() );
		assert_eq( "Initial la threshold INFO", Logger::INFO, la.getThreshold() );
		assert_eq( "Initial lb threshold INFO", Logger::INFO, lb.getThreshold() );

		logger.info( "Nothing registered yet" );
		logger.addLogger( la );
		la.setExpected( "INFO\tRegistered A" );
		logger.info( "Registered A" );
		logger.addLogger( lb );
		la.setExpected( "INFO\tRegistered A+B" );
		lb.setExpected( "INFO\tRegistered A+B" );
		logger.info( "Registered A+B" );

		logger.setThreshold( Logger::WARNING );
		assert_eq( "Initial logger threshold WARNING", Logger::WARNING, logger.getThreshold() );
		assert_eq( "Initial la threshold WARNING", Logger::WARNING, la.getThreshold() );
		assert_eq( "Initial lb threshold WARNING", Logger::WARNING, lb.getThreshold() );

		la.setExpected( NULL );
		lb.setExpected( NULL );
		logger.info( "Registered A+B" );
	}

}
