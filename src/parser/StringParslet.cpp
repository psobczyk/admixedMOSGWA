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
