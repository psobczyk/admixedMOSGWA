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

#ifndef _PARSLETTEMPLATE_HPP_
#define _PARSLETTEMPLATE_HPP_

#include "Parslet.hpp"

namespace parser {

	/** Abstract class template to help implementing {@link Parslet}. */
	template <class T> class ParsletTemplate : public Parslet {

		private:

		/** The variable to be configured. */
		T &variable;

		public:

		/** Construct an instance for a given variable.
		* @param variable to be configured.
		*/
		ParsletTemplate ( T &variable );

		/** Get the current variable value. */
		T get () const;

		/** Set the value of the configured variable. */
		void set ( const T value );
	};

}

// C++ compilers want template functions source code available when the template functions are used
#include "ParsletTemplate.cpp"

#endif /* _PARSLETTEMPLATE_HPP_ */
