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

#include "TestSuite.hpp"
#include <cmath>
#include <cstdio>
#include <cstring>
#include <unistd.h>	// for rmdir(char[])
#include <cstdlib>	// for mkdtemp(char[])

using namespace std;
using namespace unitpp;

namespace test {

	TestSuite::TestSuite ( const string &title ) : suite ( title ) {
		suite::main().add( title, this );
	}

	string TestSuite::createTmpDir ( const char* tmpDirNameTemplate ) {
		const size_t bufferLength = strlen( tmpDirNameTemplate ) + 1;	// +1 for trailing 0
		vector<char> buffer( bufferLength );
		strncpy( buffer.data(), tmpDirNameTemplate, bufferLength );
		buffer.at( bufferLength - 1 ) = '\0';	// unnecessary unless data changes out of control
		assert_true( "Create test directory failed", NULL != mkdtemp( buffer.data() ) );
		const string tmpDirName( buffer.data() );
		return tmpDirName;
	}

	void TestSuite::removeTmpDir ( const string &tmpDirName ) {
		assert_eq( "Remove test directory failed", 0, rmdir( tmpDirName.c_str() ) );
	}

	void TestSuite::assert_close(
		const string& msg,
		const double expected,
		const double actual,
		const double tolerance
	) {
		char buffer[1024];
		snprintf(
			buffer,
			sizeof( buffer ),
			" expected close to %e and got %e, difference was %e, tolerance %e",
			expected,
			actual,
			actual-expected,
			tolerance
		);
		string xmsg = msg + buffer;
		assert_true( xmsg, fabs( expected - actual ) < tolerance );
	}

	void TestSuite::assert_nan (
		const string& msg,
		const double actual
	) {
		assert_true( msg, ::isnan( actual ) );
	}

	void TestSuite::assertData (
		const string message,
		const size_t countBytes,
		const char * const data,
		istream& input
	) {
		for ( size_t i = 0; i < countBytes; ++i ) {
			const int c = input.get();
			assert_eq( message.c_str(), (int) data[i], c );
		}
		const int c = input.get();
		assert_eq( message + " EOF expected", EOF, input.get() );
	}

	void TestSuite::assertText (
		const string message,
		const size_t countLines,
		const char * const lines[],
		istream& input
	) {
		vector<char> buffer;
		for ( size_t i = 0; i < countLines; ++i ) {
			const size_t bufferLength = strlen( lines[i] ) + 256;	// 256 for some unexpected extra characters and \0
			if ( bufferLength >= buffer.size() ) {
				buffer.resize( bufferLength, 0 );
			}
			input.getline( buffer.data(), bufferLength );
			assert_eq(
				message,
				string( lines[i] ),
				string( buffer.data(), strnlen( buffer.data(), bufferLength ) )
			);
		}
		assert_eq( message + " eof expected", EOF, input.get() );
	}

}
