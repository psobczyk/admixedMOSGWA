#include "VariableParslet.hpp"

namespace parser {

	VariableParslet::VariableParslet ( string &variable ) : NameParslet( variable ) {}

	bool VariableParslet::parse ( const char* &text ) {
		string name( parseName( text ) );
		if ( name.empty() ) return false;
		set( name );
		return true;
	}

	const char* VariableParslet::type () {
		return "boolean";
	}

}
