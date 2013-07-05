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

#ifndef _SECTIONPARSLET_HPP_
#define _SECTIONPARSLET_HPP_

#include "NameParslet.hpp"

namespace parser {

	/** Configures a <code>string</code> variable to hold a section name. */
	class SectionParslet : public NameParslet {

		public:

		/** Construct an instance for a given section name holder variable.
		* @param variable to be configured.
		*/
		SectionParslet ( string &variable );

		virtual bool parse ( const char* &text );

		virtual const char* type ();
	};

}

#endif /* _SECTIONPARSLET_HPP_ */
