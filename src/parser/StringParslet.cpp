#include "StringParslet.hpp"

namespace parser {

	StringParslet::StringParslet ( string &variable ) : ParsletTemplate<string>( variable ) {}

	bool StringParslet::parse ( const char* &text ) {
		if ( '"' != *text ) return false;
		const char* endQuote = strrchr( text, '"' );
		if ( endQuote == text ) return false;
		set( string( text + 1, endQuote - text - 1 ) );
		text = endQuote + 1;
		return true;
	}

	const char* StringParslet::type () {
		return "string";
	}

}
