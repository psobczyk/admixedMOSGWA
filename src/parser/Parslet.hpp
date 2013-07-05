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

#ifndef _PARSLET_HPP_
#define _PARSLET_HPP_

namespace parser {

	/** Abstract class to parse a value according to its expected type. */
	class Parslet {

		public:

		/** Parse from a string.
		* @param text is parsed
		* and upon success advanced to after the last character parsed
		* @returns whether anything has been be parsed,
		* i.e. whether <code>text</code> has been moved
		*/
		virtual bool parse ( const char* &text ) = 0;

		/** Informs about the expected variable type for log messages. */
		virtual const char* type () = 0;

		/** Empty placeholder for valuable subclass destructors. */
		virtual ~Parslet ();
	};

}

#endif /* _PARSLET_HPP_ */
