#include "SectionParslet.hpp"

namespace parser {

	SectionParslet::SectionParslet ( string &variable ) : NameParslet( variable ) {}

	bool SectionParslet::parse ( const char* &text ) {
		if ( '[' != *text ) return false;
		const char* cursor = text + 1;
		string name( parseName( cursor ) );
		if ( name.empty() ) return false;
		if ( ']' != *cursor ) return false;
		set( name );
		text = cursor + 1;
		return true;
	}

	const char* SectionParslet::type () {
		return "section heading";
	}

}
