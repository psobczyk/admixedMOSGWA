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
