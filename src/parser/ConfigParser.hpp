#ifndef _CONFIGPARSER_HPP_
#define _CONFIGPARSER_HPP_
#include <vector>
#include <stdlib.h>
#include <string>
#include <map>
#include <iostream>
#include "Parslet.hpp"

using namespace std;

/** Contains the configuration file parser.
* @author Bernhard Bodenstorfer
*/
namespace parser {

	/** Configures a set of variables from an ini-like configuration file */
	class ConfigParser {

		private:

		/** Known sections with their variable-accessors */
		map< const string, map< const string, class Parslet* >* > sectionMaps;

		/** Declare a configuration variable in the given section. */
		void declare ( const string &sectionName, const string &variableName, Parslet *parslet );

		/** Skip whitespace, if any */
		static void skipWhitespace ( const char* &text );

		/** Skip comment, if any */
		static void skipComment ( const char* &text );

		/** Parse one line of configuration.
		* @param filename used to formulate error messages
		* @param lineNumber used to formulate error messages
		* @returns whether an error has occurred
		*/
		bool parseLine ( const char* &text, string &sectionName, const int lineNumber, const char* filename );

		public:

		/** Construct a configuration parser with empty variable set */
		ConfigParser ();

		/** Declare a boolean configuration variable in the given section. */
		void declare ( const string &sectionName, const string &variableName, bool &variable );

		/** Declare an integer configuration variable in the given section. */
		void declare ( const string &sectionName, const string &variableName, int &variable );

		/** Declare a double precision numeric configuration variable in the given section. */
		void declare ( const string &sectionName, const string &variableName, double &variable );

		/** Declare a string configuration variable in the given section. */
		void declare ( const string &sectionName, const string &variableName, string &variable );

                /** Declare a vector<string> configuration variable in the given section. */
		void declare ( const string &sectionName, const string &variableName, vector<string> &variable );

		/** Declare an integer configuration variable represented by a choice of strings in the given section. */
		void declare ( const string &sectionName, const string &variableName, int &variable, const map< const string, int > &choices );


		/** Parse sections of name-value-pairs.
		* @param input to read the configuration text
		* @param filename used to formulate error messages
		* @returns whether an error has occurred
		*/
		bool parse ( istream &input, const char* filename );

		/** Destruct a configuration parser, free its internally allocated memory resources */
		~ConfigParser ();
	};

}

#endif /* _CONFIGPARSER_HPP_ */
