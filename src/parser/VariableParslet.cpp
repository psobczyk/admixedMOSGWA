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

#include "VariableParslet.hpp"

namespace parser {

	VariableParslet::VariableParslet ( string &variable ) : NameParslet( variable ) {}

	bool VariableParslet::parse ( const char* &text ) {
		string name( parseName( text ) );
		if ( name.empty() ) return false;
		set( name );
		return true;
	}

	const char* VariableParslet::type () {
		return "boolean";
	}

}
