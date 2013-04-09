#include "IntParslet.hpp"
#include <stdlib.h>

namespace parser {

	IntParslet::IntParslet ( int &variable ) : ParsletTemplate<int>( variable ) {}

	bool IntParslet::parse ( const char* &text ) {
		char* cursor = NULL;
		long l = strtol( text, &cursor, 10 );
		if ( NULL == cursor || text == cursor ) return false;
		if ( static_cast<int>( l ) != l ) return false;
		set( l );
		text = cursor;
		return true;
	}

	const char* IntParslet::type () {
		return "integer";
	}

}
