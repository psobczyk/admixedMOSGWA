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

#ifndef _CONFIGPARSER_HPP_
#define _CONFIGPARSER_HPP_
#include <vector>
#include <stdlib.h>
#include <string>
#include <map>
#include <iostream>
#include "Parslet.hpp"

/** Contains the configuration file parser.
* @author Bernhard Bodenstorfer
*/
namespace parser {

	/** Configures a set of variables from an ini-like configuration file */
	class ConfigParser {

		private:

		/** Known sections with their variable-accessors */
		std::map< const std::string, std::map< const std::string, class Parslet* >* > sectionMaps;

		/** Declare a configuration variable in the given section. */
		void declare (
			const std::string &sectionName,
			const std::string &variableName,
			Parslet *parslet
		);

		/** Skip whitespace, if any */
		static void skipWhitespace ( const char* &text );

		/** Skip comment, if any */
		static void skipComment ( const char* &text );

		/** Parse one line of configuration.
		* @param filename used to formulate error messages
		* @param lineNumber used to formulate error messages
		* @returns whether an error has occurred
		*/
		void parseLine (
			const char* &text,
			std::string &sectionName,
			const int lineNumber,
			const char* filename
		);

		public:

		/** Construct a configuration parser with empty variable set */
		ConfigParser ();

		/** Declare a boolean configuration variable in the given section. */
		void declare (
			const std::string &sectionName,
			const std::string &variableName,
			bool &variable
		);

		/** Declare an integer configuration variable in the given section. */
		void declare (
			const std::string &sectionName,
			const std::string &variableName,
			int &variable
		);

		/** Declare a double precision numeric configuration variable in the given section. */
		void declare (
			const std::string &sectionName,
			const std::string &variableName,
			double &variable
		);

		/** Declare a string configuration variable in the given section. */
		void declare (
			const std::string &sectionName,
			const std::string &variableName,
			std::string &variable
		);

                /** Declare a vector<string> configuration variable in the given section. */
		void declare (
			const std::string &sectionName,
			const std::string &variableName,
			std::vector<std::string> &variable
		);

		/** Declare an integer configuration variable represented by a choice of strings in the given section. */
		void declare (
			const std::string &sectionName,
			const std::string &variableName,
			int &variable,
			const std::map< const std::string, int > &choices
		);


		/** Parse sections of name-value-pairs.
		* @param input to read the configuration text
		* @param filename used to formulate error messages
		* @returns whether an error has occurred
		*/
		void parse ( std::istream &input, const char* filename );

		/** Destruct a configuration parser, free its internally allocated memory resources */
		~ConfigParser ();
	};

}

#endif /* _CONFIGPARSER_HPP_ */
