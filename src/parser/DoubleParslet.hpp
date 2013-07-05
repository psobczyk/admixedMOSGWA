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

#ifndef _DOUBLEPARSLET_HPP_
#define _DOUBLEPARSLET_HPP_

#include "ParsletTemplate.hpp"

namespace parser {

	/** Configures a <code>double</code> variable. */
	class DoubleParslet : public ParsletTemplate<double> {

		public:

		/** Construct an instance for a given <code>double</code> variable.
		* @param variable to be configured.
		*/
		DoubleParslet ( double &variable );

		virtual bool parse ( const char* &text );

		virtual const char* type ();
	};

}

#endif /* _DOUBLEPARSLET_HPP_ */
