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

#ifndef LINALG_AUTOVECTOR_HPP
#define LINALG_AUTOVECTOR_HPP

#include <ostream>
#include "Vector.hpp"

namespace linalg {

	class Vector;

	/** A cheap C++ wrapper for a gently resizing GSL vector.
	* As the logical and the physical sizes of the vector generally differ,
	* the logical vector is a view of a sub-vector of the physical storage array.
	*/
	class AutoVector : public Vector {

		/** The currently allocated memory block size. */
		size_t size;

		protected:

		/** Calculate memory requirement. */
		static size_t calculateSize ( const size_t dims );

		public:

		/** Construct a vector of given dimensions. */
		explicit AutoVector ( const size_t dims );

		/** Construct a vector of given dimensions with possibly some extra space for easy growth. */
		AutoVector ( const size_t dims, const size_t allocateDims );

		/** Copy constructor. */
		AutoVector ( const AutoVector& original );

		/** Assignment operator. */
		AutoVector& operator= ( const AutoVector& original );

		/** Resize to the given dimensions.
		* Values disappear when logical dimensions shrink.
		*/
		void exactSize ( const size_t dims );

		/** Resize to accommodate the given dimensions and possibly extra space for easy growth.
		* Values will disappear when dimensions shrink.
		*/
		void upSize ( const size_t dims );

		/** Delete the given coordinate.
		* Elements after that are shifted, if necessary.
		* If you only want to decrement the number of dimensions,
		* use {@link #upSize} or {@link #exactSize}.
		* @see #countDimensions
		*/
		void removeDimension ( const size_t dim );

		/** Destruct, particularly free the internal <code>double</code> array. */
		virtual ~AutoVector ();
	};

}

#endif	/* LINALG_AUTOVECTOR_HPP */
