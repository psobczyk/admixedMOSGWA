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
