#include "ConfigParser.hpp"
#include "SectionParslet.hpp"
#include "VariableParslet.hpp"
#include "BoolParslet.hpp"
#include "IntParslet.hpp"
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

	bool ConfigParser::parseLine ( const char* &text, string &sectionName, const int lineNumber, const char* filename ) {

		// Skip leading whitespace
		skipWhitespace( text );
		skipComment( text );

		// Skip blank or comment line
		if ( 0 == *text ) return true;

		// Is it a [sectionname] line?
		SectionParslet sp( sectionName );
		const bool sectionHeading = sp.parse( text );
		if ( sectionHeading ) {
			skipWhitespace( text );
			skipComment( text );
		}

		// If the current sectionName is not known
		if ( 0 == sectionMaps.count( sectionName ) ) {
			// Warn only once: when section is opened
			if ( ( 1 == lineNumber ) || sectionHeading ) {
				cerr << filename << ":" << lineNumber << ": "
					<< "Warning: Config section name \""
					<< sectionName << "\" unknown: Section contents skipped."
					<< endl;
			}
			while ( 0 != *text ) ++text;
			return true;
		}

		// Line information is either section heading or variable setting
		if ( ! sectionHeading ) {

			// Get the variable name
			string variableName;
			VariableParslet vp( variableName );
			if ( ! vp.parse( text ) ) {
				cerr << filename << ":" << lineNumber << ": "
					<< "Syntax error: Expected valid variable name, but got \""
					<< text << "\"." << endl;
				return false;
			}
			skipWhitespace( text );

			// Expect an equals sign
			if ( '=' != *text ) {
				cerr << filename << "::" << lineNumber << ": "
					<< "Syntax error: Expected '=' after variable name \""
					<< variableName << "\", but got \""
					<< text << "\"." << endl;
				return false;
			} else {
				++text;
			}
			skipWhitespace( text );

			// Get the Parslet to parse the variable value
			map< const string, Parslet * > &sectionMap = *sectionMaps[ sectionName ];
			if ( 0 == sectionMap.count( variableName ) ) {
				cerr << filename << ":" << lineNumber << ": "
					<< "Warning: Config variable name \"" << variableName
					<< "\" unknown in Section \"" << sectionName
					<< "\": Variable setting skipped." << endl;
				while ( 0 != *text ) ++text;
				return true;
			}
			Parslet &parslet = *sectionMap[ variableName ];

			// Parse the variable value
			if ( ! parslet.parse( text ) ) {
				cerr << filename << ":" << lineNumber << ": "
					<< "Error: Unrecognised value for "
					<< parslet.type()
					<< " variable name \""
					<< variableName <<
					"\": \"" << text << "\"." << endl;
				return false;
			}
			skipWhitespace( text );
			skipComment( text );
		}

		// Anything left on the line?
		if ( 0 != *text ) {
			cerr << filename << ":" << lineNumber << ": "
				<< "Syntax error: Unrecognised characters: \""
				<< text << "\"." << endl;
			return false;
		}

		return true;
	}

	/** This function assumes as precondition
	* that the parslet has been allocated on the heap.
	* Otherwise, the destructor might run into trouble.
	*/
	void ConfigParser::declare ( const string &sectionName, const string &variableName, Parslet *parslet ) {
		if ( NULL == parslet ) {
			cerr << "Implementation error or out of memory:"
				<< " Null pointer to parslet for configuration section name \""
				<< sectionName << "\", variable name \""
				<< variableName << "\"." << endl;
			exit( 252 );
		}
		if ( 0 == sectionMaps.count( sectionName ) ) {
			sectionMaps[ sectionName ] = new map< const string, Parslet* >();
		}
		map< const string, Parslet* > &sectionMap = *sectionMaps[ sectionName ];
		if ( 0 == sectionMap.count( variableName ) ) {
			sectionMap[ variableName ] = parslet;
		} else {
			// This should not happen.
			cerr << "Implementation error:"
				<< " Duplicate configuration variable definition in section name \""
				<< sectionName << "\", variable name \""
				<< variableName << "\"." << endl;
			// Free allocated parslet, because it would not be freed by ~ConfigParser
			delete parslet;
			exit( 253 );	// Exit anyway
		}
	}

	void ConfigParser::declare ( const string &sectionName, const string &variableName, bool &variable ) {
		declare( sectionName, variableName, new BoolParslet( variable ) );
	}

	void ConfigParser::declare ( const string &sectionName, const string &variableName, int &variable ) {
		declare( sectionName, variableName, new IntParslet( variable ) );
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

	bool ConfigParser::parse( istream &input, const char* filename ) {
		bool ok = true;
		string sectionName;
		string line;
		for ( int lineNumber = 1; getline( input, line ); ++lineNumber ) {
			const char* text = line.c_str();
			ok &= parseLine( text, sectionName, lineNumber, filename );
		}
		return ok;
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
