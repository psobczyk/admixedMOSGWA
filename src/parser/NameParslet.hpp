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

#ifndef _NAMEPARSLET_HPP_
#define _NAMEPARSLET_HPP_

#include "ParsletTemplate.hpp"
#include <string>

using namespace std;

namespace parser {

	/** Abstract base class for parsing variable and section names */
	class NameParslet : public ParsletTemplate<string> {

		protected:

		/** Parse a variable name into a <code>string</code>. */
		string parseName ( const char* &text );

		public:

		/** Configure a <code>string</code> name of a variable.
		* @param variable to be configured.
		*/
		NameParslet ( string &variable );
	};

}

#endif /* _NAMEPARSLET_HPP_ */
