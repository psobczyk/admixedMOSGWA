/********************************************************************************
 *	This file is part of the MOSGWA program code.				*
 *	Copyright ©2013, Bernhard Bodenstorfer.					*
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

#ifndef UTIL_ARRAY_HPP
#define UTIL_ARRAY_HPP

#include <cstddef>	// for size_t
#include <cstdarg>

namespace util {

	/** Provides a way to pass an initialised array as constructor argument
	* in the member initialisation sequence of a constructor.
	* This class may become largely obsolete with C++ standard 2011,
	* where the package <code>initializer_list</code> does many tricks.
	* Let me know if somebody has an idea how to raise a compile-time error
	* when the constructor from ellipsis is given too few arguments,
	* i.e. how to make the template parameter <code>D</code>
	* determine the (minimum) number of arguments of the constructor.
	*/
	template<class E,size_t D> class Array {

		/** Holds the values. */
		E array[D];

		/** Helps construction. */
		void init ( const size_t begin, ::va_list args );

		public:

		/** Initialisation for <code>D=0</code>. */
		Array ();

		/** Initialise with given arguments for <code>D&gt;0</code>.
		* There must be at least <code>D</code> arguments.
		* The first of them is <code>initial</code>,
		* the following must be of the same type.
		*/
		Array ( E initial, ... );

		/** Initialise with arguments from a variable arguments list.
		* There must be at least <code>D</code> arguments
		* and they should be of the type <code>E</code>.
		*/
		Array ( ::va_list args );

		/** Provide the array for use.
		* @returns a pointer to the first element
		* of the <code>D</code>-dimensional array.
		*/
		operator E* ();
	};

}

// C++ compilers want template functions source code available when the template functions are used
#include "Array.tpl"

#endif	/* UTIL_ARRAY_HPP */
