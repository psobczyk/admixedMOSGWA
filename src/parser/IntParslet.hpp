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

#ifndef _INTPARSLET_HPP_
#define _INTPARSLET_HPP_

#include "ParsletTemplate.hpp"

namespace parser {

	/** Configures an <code>int</code> variable. */
	class IntParslet : public ParsletTemplate<int> {

		public:

		/** Construct an instance for a given <code>int</code> variable.
		* @param variable to be configured.
		*/
		IntParslet ( int &variable );

		virtual bool parse ( const char* &text );

		virtual const char* type ();
	};

}

#endif /* _INTPARSLET_HPP_ */
