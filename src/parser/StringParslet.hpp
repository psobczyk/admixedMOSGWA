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

#ifndef _STRINGPARSLET_HPP_
#define _STRINGPARSLET_HPP_

#include "ParsletTemplate.hpp"

#include <string.h>
#include <string>

using namespace std;

namespace parser {

	/** Configures an <code>string</code> variable. */
	class StringParslet : public ParsletTemplate<string> {

		public:

		/** Construct an instance for a given <code>string</code> variable.
		* @param variable to be configured.
		*/
		StringParslet ( string &variable );

		virtual bool parse ( const char* &text );

		virtual const char* type ();
	};

}

#endif /* _STRINGPARSLET_HPP_ */
