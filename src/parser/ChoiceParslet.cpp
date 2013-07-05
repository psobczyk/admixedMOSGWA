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

#ifndef _CHOICEPARSLET_CPP_
#define _CHOICEPARSLET_CPP_

#include "ChoiceParslet.hpp"
#include <string.h>
#include <iostream>

using namespace std;

namespace parser {

	template <class T> ChoiceParslet<T>::ChoiceParslet ( T &variable, const map < const string, T > &choices )
		: ParsletTemplate<T>( variable ), choices( choices ), typeString( "choice" )
	{
		if ( this->choices.empty() ) {
			cerr << "Implementation error: configuration choice with 0 possibilities." << endl;
		}
		typeString += "(";
		for (
			typename map< const string, T >::const_iterator iterator = this->choices.begin();
			iterator != this->choices.end();
			++iterator
		) {
			if ( iterator != this->choices.begin() ) {
				typeString += ",";
			}
			typeString += iterator->first;
		}
		typeString += ")";
	}

	template <class T> bool ChoiceParslet<T>::parse ( const char* &text ) {
		for (
			typename map< const string, T >::const_iterator iterator = choices.begin();
			iterator != choices.end();
			++iterator
		) {
			const char* choice = iterator->first.c_str();
			const int length = strlen( choice );
			if ( 0 == strncmp( choice, text, length ) ) {
				this->set( iterator->second );
				text += length;
				return true;
			}
		}
		return false;
	}

	// Possible space optimisation: factor dealing with ChoiceParslet::typeString out into a (second, non-template) superclass.
	template <class T> const char* ChoiceParslet<T>::type () {
		return typeString.c_str();
	}
}

#endif /* _CHOICEPARSLET_CPP_ */
