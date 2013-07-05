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

#ifndef _CHOICEPARSLET_HPP_
#define _CHOICEPARSLET_HPP_

#include "ParsletTemplate.hpp"
#include <string>
#include <map>

using namespace std;

namespace parser {

	/** Configures a variable via syntactic choices taken from a {@link std::map}. */
	template <class T> class ChoiceParslet : public ParsletTemplate<T> {

		private:

		/** Maps possible syntactic choices to variable values */
		const map < const string, T > choices;

		/** Describes the kind of {@link Parslet} this is.
		* @see #type()
		*/
		string typeString;

		public:

		/** Construct an instance for a given variable and given map from choices to values.
		* @param variable to be configured.
		* @param choices maps syntactic keys to variable values
		*/
		ChoiceParslet ( T &variable, const map < const string, T > &choices );

		virtual bool parse ( const char* &text );

		virtual const char* type ();
	};

}

// C++ compilers want template functions source code available when the template functions are used
#include "ChoiceParslet.cpp"

#endif /* _CHOICEPARSLET_HPP_ */
