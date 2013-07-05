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
