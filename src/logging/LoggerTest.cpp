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

#include "TestLogger.hpp"

using namespace std;
using namespace unitpp;
using namespace logging;

namespace test {

	/** Tests the class {@link logging::Logger}. */
	class LoggerTest : public TestSuite {

		/** Test interface methods for reading. */
		void testLog ();

		public:

		/** Construct the test object. */
		LoggerTest ();

	} * loggerTest = new LoggerTest();	// automatically freed by unit++

	LoggerTest::LoggerTest () : TestSuite( "LoggerTest" ) {
		addTestMethod( "LoggerTest::testLog", this, &LoggerTest::testLog );
	}

	void LoggerTest::testLog () {
		const char
			*error = "ERROR\tWhorror",
			*warning = "WARNING\tWhoaning",
			*info = "INFO\tWhoafo",
			*debug = "DEBUG\tWobog";
		TestLogger logger;

		assert_eq( "Default threshold INFO", Logger::INFO, logger.getThreshold() );

		logger.setThreshold( Logger::ERROR );
		assert_eq( "New threshold ERROR", Logger::ERROR, logger.getThreshold() );
		logger.setExpected( NULL );
		logger.debug( "0" );
		logger.info( "1" );
		logger.warning( "2" );
		logger.setExpected( error );
		logger.error( "Whorror" );

		logger.setThreshold( Logger::WARNING );
		assert_eq( "New threshold WARNING", Logger::WARNING, logger.getThreshold() );
		logger.setExpected( NULL );
		logger.debug( "3" );
		logger.info( "4" );
		logger.setExpected( warning );
		logger.warning( "Whoaning" );
		logger.setExpected( error );
		logger.error( "Whorror" );

		logger.setThreshold( Logger::INFO );
		assert_eq( "New threshold INFO", Logger::INFO, logger.getThreshold() );
		logger.setExpected( info );
		logger.debug( "5" );
		logger.info( "Whoafo" );
		logger.setExpected( warning );
		logger.warning( "Whoaning" );
		logger.setExpected( error );
		logger.error( "Whorror" );

		logger.setThreshold( Logger::DEBUG );
		assert_eq( "New threshold DEBUG", Logger::DEBUG, logger.getThreshold() );
		logger.setExpected( debug );
		logger.debug( "Wobog" );
		logger.setExpected( info );
		logger.info( "Whoafo" );
		logger.setExpected( warning );
		logger.warning( "Whoaning" );
		logger.setExpected( error );
		logger.error( "Whorror" );
	}

}
