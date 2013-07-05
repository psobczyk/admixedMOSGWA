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
