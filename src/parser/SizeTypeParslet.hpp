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

#ifndef _SIZETYPE_PARSLET_HPP_
#define _SIZETYPE_PARSLET_HPP_

#include "ParsletTemplate.hpp"
#include <cstdlib>

namespace parser {

	/** Configures an <code>int</code> variable. */
	class SizeTypeParslet : public ParsletTemplate<size_t> {

		public:

		/** Construct an instance for a given <code>int</code> variable.
		* @param variable to be configured.
		*/
		SizeTypeParslet ( size_t &variable );

		virtual bool parse ( const char* &text );

		virtual const char* type ();
	};

}

#endif /* _SIZETYPE_PARSLET_HPP_ */
