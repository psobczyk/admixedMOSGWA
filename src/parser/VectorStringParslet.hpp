/********************************************************************************
 *	This file is part of the MOSGWA program code.				*
 *	Copyright ©2012–2013, Erich Dolejsi.					*
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

#ifndef _VECTORSTRINGPARSLET_HPP_
#define _VECTORSTRINGPARSLET_HPP_

#include "ParsletTemplate.hpp"
#include <vector>
#include <string.h>
#include <string>

using namespace std;

namespace parser {

	/** Configures an <code>string</code> variable. */
	class VectorStringParslet : public ParsletTemplate<vector<string> > {

		public:

		/** Construct an instance for a given <code>string</code> variable.
		* @param variable to be configured.
		*/
		VectorStringParslet (vector< string> &variable );

		virtual bool parse ( const char* &text );

		virtual const char* type ();

		/** Advance the text pointer to the next non-whitespace character. */
		void skipWhitespace ( const char* &text );
	};

}

#endif /* _VECTORSTRINGPARSLET_HPP_ */
