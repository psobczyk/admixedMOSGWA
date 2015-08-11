/********************************************************************************
 *	This file is part of the MOSGWA program code.				*
 *	Copyright ©2015, Bernhard Bodenstorfer.				*
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

#include "SizeTypeParslet.hpp"

namespace parser {

	SizeTypeParslet::SizeTypeParslet ( size_t &variable ) : ParsletTemplate<size_t>( variable ) {}

	bool SizeTypeParslet::parse ( const char* &text ) {
		char* cursor = NULL;
		const unsigned long l = strtoul( text, &cursor, 10 );
		for ( const char *c = text; c < cursor; ++c ) {
			if ( '-' == *c ) {
				return false;
			}
		}
		if ( NULL == cursor || text == cursor ) return false;
		if ( static_cast<size_t>( l ) != l ) return false;
		set( l );
		text = cursor;
		return true;
	}

	const char* SizeTypeParslet::type () {
		return "size_type";
	}

}
