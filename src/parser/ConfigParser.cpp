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

#include "../Exception.hpp"
#include "ConfigParser.hpp"
#include "SectionParslet.hpp"
#include "VariableParslet.hpp"
#include "BoolParslet.hpp"
#include "IntParslet.hpp"
#include "SizeTypeParslet.hpp"
#include "DoubleParslet.hpp"
#include "StringParslet.hpp"
#include "VectorStringParslet.hpp"
#include "ChoiceParslet.hpp"

namespace parser {

	ConfigParser::ConfigParser () {}

	void ConfigParser::skipWhitespace ( const char* &text ) {
		while ( isspace( *text ) ) ++text;
	}

	void ConfigParser::skipComment ( const char* &text ) {
		if ( '#' != *text ) return;
		while ( 0 != *text ) ++text;
	}

	void ConfigParser::parseLine (
		const char* &text,
		string &sectionName,
		const int lineNumber,
		const char* filename
	) {

		// Skip leading whitespace
		skipWhitespace( text );
		skipComment( text );

		// Skip blank or comment line
		if ( 0 == *text ) return;

		// Is it a [sectionname] line?
		SectionParslet sp( sectionName );
		const bool sectionHeading = sp.parse( text );
		if ( sectionHeading ) {
			skipWhitespace( text );
			skipComment( text );
		}

		// If the current sectionName is not known
		if ( 0 == sectionMaps.count( sectionName ) ) {
			throw Exception(
				"%s:%u:"
				" Unrecognised configuration section name \"%s\".",
				filename,
				lineNumber,
				sectionName.c_str()
			);
		}

		// Line information is either section heading or variable setting
		if ( ! sectionHeading ) {

			// Get the variable name
			string variableName;
			VariableParslet vp( variableName );
			if ( ! vp.parse( text ) ) {
				throw Exception(
					"%s:%u:"
					" Syntax error: Expected valid variable name, but got \"%s\".",
					filename,
					lineNumber,
					text
				);
			}
			skipWhitespace( text );

			// Expect an equals sign
			if ( '=' != *text ) {
				throw Exception(
					"%s:%u:"
					" Syntax error: Expected '=' after variable name \"%s\","
					" but got \"%s\".",
					filename,
					lineNumber,
					variableName.c_str(),
					text
				);
			} else {
				++text;
			}
			skipWhitespace( text );

			// Get the Parslet to parse the variable value
			map< const string, Parslet * > &sectionMap = *sectionMaps[ sectionName ];
			if ( 0 == sectionMap.count( variableName ) ) {
				throw Exception(
					"%s:%u:"
					" Unrecognized config variable name \"%s\""
					" in Section \"%s\".",
					filename,
					lineNumber,
					variableName.c_str(),
					sectionName.c_str()
				);
			}
			Parslet &parslet = *sectionMap[ variableName ];

			// Parse the variable value
			if ( ! parslet.parse( text ) ) {
				throw Exception(
					"%s:%u:"
					" Unrecognised value in Section \"%s\""
					" for %s variable name \"%s\": \"%s\".",
					filename,
					lineNumber,
					sectionName.c_str(),
					parslet.type(),
					variableName.c_str(),
					text
				);
			}
			skipWhitespace( text );
			skipComment( text );
		}

		// Anything left on the line?
		if ( 0 != *text ) {
			throw Exception(
				"%s:%u:"
				" Unrecognised characters: \"%s\".",
				filename,
				lineNumber,
				text
			);
		}
	}

	/** This function assumes as precondition
	* that the parslet has been allocated on the heap.
	* Otherwise, the destructor might run into trouble.
	*/
	void ConfigParser::declare ( const string &sectionName, const string &variableName, Parslet *parslet ) {
		if ( NULL == parslet ) {
			throw Exception(
				"Implementation error or out of memory:"
				" Null pointer to parslet for configuration section name \"%s\","
				" variable name \"%s\".",
				sectionName.c_str(),
				variableName.c_str()
			);
		}
		if ( 0 == sectionMaps.count( sectionName ) ) {
			sectionMaps[ sectionName ] = new map< const string, Parslet* >();
		}
		map< const string, Parslet* > &sectionMap = *sectionMaps[ sectionName ];
		if ( 0 == sectionMap.count( variableName ) ) {
			sectionMap[ variableName ] = parslet;
		} else {
			// This should not happen.
			// Free allocated parslet, because it would not be freed by ~ConfigParser
			delete parslet;
			throw Exception(
				"Implementation error:"
				" Duplicate configuration variable definition in section \"%s\","
				" variable name \"%s\".",
				sectionName.c_str(),
				variableName.c_str()
			);
		}
	}

	void ConfigParser::declare ( const string &sectionName, const string &variableName, bool &variable ) {
		declare( sectionName, variableName, new BoolParslet( variable ) );
	}

	void ConfigParser::declare ( const string &sectionName, const string &variableName, int &variable ) {
		declare( sectionName, variableName, new IntParslet( variable ) );
	}

	void ConfigParser::declare ( const string &sectionName, const string &variableName, size_t &variable ) {
		declare( sectionName, variableName, new SizeTypeParslet( variable ) );
	}

	void ConfigParser::declare ( const string &sectionName, const string &variableName, double &variable ) {
		declare( sectionName, variableName, new DoubleParslet( variable ) );
	}

	void ConfigParser::declare ( const string &sectionName, const string &variableName, string &variable ) {
		declare( sectionName, variableName, new StringParslet( variable ) );
	}
	void ConfigParser::declare ( const string &sectionName, const string &variableName, vector<string> &variable ) {
		declare( sectionName, variableName, new VectorStringParslet( variable ) );
	}


	void ConfigParser::declare ( const string &sectionName, const string &variableName, int &variable, const map< const string, int > &choices ) {
		declare( sectionName, variableName, new ChoiceParslet<int>( variable, choices ) );
	}

	void ConfigParser::parse( istream &input, const char* filename ) {
		string sectionName;
		string line;
		for ( int lineNumber = 1; getline( input, line ); ++lineNumber ) {
			const char* text = line.c_str();
			parseLine( text, sectionName, lineNumber, filename );
		}
	}

	ConfigParser::~ConfigParser () {
		for (
			map< const string, map< const string, Parslet* >* >::iterator iterator = sectionMaps.begin();
			iterator != sectionMaps.end();
			++iterator
		) {
			// iterator->first deallocate, too? TODO<BB>
			delete iterator->second; // Free the section map
		}
	}

}
