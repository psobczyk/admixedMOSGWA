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

#include "FileLogger.hpp"
#include "../TestSuite.hpp"
#include <string>
#include <fstream>
#include <cstdio>

using namespace std;
using namespace log;
using namespace unitpp;

namespace test {

	/** Tests the class {@link log::FileLogger}. */
	class FileLoggerTest : public TestSuite {

		/** Name template for temporary test directory. */
		static const char * const tmpDirnameTemplate;

		/** Name for temporary test file. */
		static const char * const tmpFilename;

		/** Prepare test setup.
		* @returns name of the test directory.
		*/
		string setUp ();

		/** Remove test setup. */
		void tearDown ( const string& tmpDirname, const string& tmpFilepath );

		/** Test logging. */
		void testLog ();

		public:

		/** Construct the test object. */
		FileLoggerTest ();

	} * fileLoggerTest = new FileLoggerTest();	// automatically freed by unit++

	const char * const FileLoggerTest::tmpDirnameTemplate = "filelogger-test.XXXXXX";
	const char * const FileLoggerTest::tmpFilename = "log.txt";

	FileLoggerTest::FileLoggerTest () : TestSuite( "FileLoggerTest" ) {
		addTestMethod( "FileLoggerTest::testLog", this, &FileLoggerTest::testLog );
	}

	string FileLoggerTest::setUp () {
		return createTmpDir( tmpDirnameTemplate );
	}

	void FileLoggerTest::tearDown ( const string& tmpDirname, const string& tmpFilepath ) {
		if ( 0 != remove( tmpFilepath.c_str() ) ) {
			perror( "Remove file failed" );
			assert_fail( "Remove tmp file" );
		}
		removeTmpDir( tmpDirname );
	}

	void FileLoggerTest::testLog () {
		const string
			tmpDirname( setUp() ),
			tmpFilepath( tmpDirname + '/' + tmpFilename );
		const char * const path = tmpFilepath.c_str();
		{
			FileLogger logger( path );
			logger.setLimit( Logger::WARNING );
			logger.info( "This goes nowhere: %u", 3u );
			logger.warning( "This goes to file: %u = %s", 12u, "12" );
			logger.error( "This goes to file, too: %d = %s", -1, "-1" );
		}
		{
			FileLogger logger( path );
			logger.setLimit( Logger::ERROR );
			logger.info( "This goes nowhere again: %u", 3u );
			logger.warning( "This now goes nowhere: %u = %s", 14u, "14" );
			logger.error( "This, finally, goes to file: %d = %s", -10, "-10" );
		}

		ifstream infile( path );
		string date, time, severity, message;

		// timestamps are tested in LoggerTest, so trusted here

		infile >> date >> time >> severity;
		getline( infile, message );
		assert_eq( "1st line severity", "WARNING", severity );
		assert_eq( "1st line message", "\tThis goes to file: 12 = 12", message );

		infile >> date >> time >> severity;
		getline( infile, message );
		assert_eq( "2nd line severity", "ERROR", severity );
		assert_eq( "2nd line message", "\tThis goes to file, too: -1 = -1", message );

		infile >> date >> time >> severity;
		getline( infile, message );
		assert_eq( "3rd line severity", "ERROR", severity );
		assert_eq( "3rd line message", "\tThis, finally, goes to file: -10 = -10", message );

		assert_false( "No 4th line", getline( infile, message ) ? 1 : 0 );

		tearDown( tmpDirname, tmpFilepath );
	}

}
