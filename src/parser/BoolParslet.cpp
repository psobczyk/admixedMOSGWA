#include "BoolParslet.hpp"
#include <string.h>

namespace parser {

	BoolParslet::BoolParslet ( bool &variable ) : ParsletTemplate<bool>( variable ) {}

	bool BoolParslet::parse ( const char* &text ) {
		if ( 0 == strncmp( text, "true", 4 ) ) {
			set( true );
			text += 4;
			return true;
		}
		if ( 0 == strncmp( text, "false", 5 ) ) {
			set( false );
			text += 5;
			return true;
		}
		return false;
	}

	const char* BoolParslet::type () {
		return "boolean";
	}

}
