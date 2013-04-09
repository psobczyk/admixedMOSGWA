#include "NameParslet.hpp"

namespace parser {

	NameParslet::NameParslet ( string &variable ) : ParsletTemplate<string>( variable ) {}

	string NameParslet::parseName ( const char* &text ) {
		const char* cursor = text;
		if (
			( 'A' > *cursor || *cursor > 'Z' )
			&&
			( 'a' > *cursor || *cursor > 'z' )
		) return "";
		while (
			'0' <= *cursor && *cursor <= '9'
			||
			'A' <= *cursor && *cursor <= 'Z'
			||
			'a' <= *cursor && *cursor <= 'z'
			||
			'_' == *cursor
		) ++cursor;
		string name( text, cursor - text );
		text = cursor;
		return name;
	}

}
