#include "DoubleParslet.hpp"
#include <stdlib.h>

namespace parser {

	DoubleParslet::DoubleParslet ( double &variable ) : ParsletTemplate<double>( variable ) {}

	bool DoubleParslet::parse ( const char* &text ) {
		char* cursor = NULL;
		double d = strtod( text, &cursor );
		if ( NULL == cursor || text == cursor ) return false;
		set( d );
		text = cursor;
		return true;
	}

	const char* DoubleParslet::type () {
		return "double precision numeric";
	}

}
