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

#ifndef UTIL_REFERENCECOUNTED_HPP
#define UTIL_REFERENCECOUNTED_HPP

#include <stdlib.h>

namespace util {

	/** Furnish a given class with a reference counter.
	* As to the handle pattern, compare B.Stroustrup, "The C++ Programming Language", 3rd ed. p. 783.
	* @param T specifies the members of the array
	* @author Bernhard Bodenstorfer
	*/
	template <class T> class ReferenceCounted {

		protected:

		/** The stored <code>T</code> instance with the reference counter. */
		struct Instance {

			/** Reference counter.
			* 0 implies memory free.
			* 1 allows the handle to use optimising shortcuts because no respect to other handles is required.
			* 2 and may mean that copies must be allocated when non-constant methods are applied
			* (lazy-copy or copy-on-write).
			*/
			size_t referenceCount;

			private:

			/** The payload. */
			T object;

			public:

			/** Construct with given payload. */
			Instance ( const T& object );

			/** Get a reference to the object pointed to. */
			virtual operator T& ();

			/** Destruct and free all resources used. */
			virtual ~Instance ();

		} * instance;

		/** Construct from pointer to instance.
		* This is intended for use with subclasses of <code>Instance</code>.
		* An <code>Instance*</code> factory would not suffice,
		* because various input signatures would be needed in subclasses.
		*/
		ReferenceCounted ( Instance * instance );

		public:

		/** Construct a handle with given payload. */
		ReferenceCounted ( const T& object );

		/** Copy constructor with due respect to reference counting. */
		ReferenceCounted ( const ReferenceCounted& original );

		/** Assignment operator with due respect to reference counting. */
		virtual ReferenceCounted& operator= ( const ReferenceCounted &that );

		/** Get a reference to the object pointed to. */
		virtual operator T& ();

		/** Destruct a handle with due respect to reference counting. */
		virtual ~ReferenceCounted ();
	};

}

// Some C++ compilers want template functions source code available when the template functions are used:
// g++ does not implement the "export" keyword.
#include "ReferenceCounted.cpp"

#endif /* UTIL_REFERENCECOUNTED_HPP */
