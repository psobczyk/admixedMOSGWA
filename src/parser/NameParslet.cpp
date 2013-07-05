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
