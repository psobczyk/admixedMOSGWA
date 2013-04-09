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
