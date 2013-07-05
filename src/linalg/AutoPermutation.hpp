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

#ifndef LINALG_AUTOPERMUTATION_HPP
#define LINALG_AUTOPERMUTATION_HPP

#include <ostream>
#include "Permutation.hpp"

namespace linalg {

	/** A cheap C++ wrapper for a gently resizing GSL permutation for pivoted QR decomposition.
	* As the logical and the physical sizes of the permutation generally differ,
	* the logical permutation uses a sub-array of the physical storage array.
	*
	* The interface is similar to that of {@link AutoVector}.
	*/
	class AutoPermutation : public Permutation {

		/** The currently allocated memory block size. */
		size_t size;

		protected:

		/** Calculate memory requirement. */
		static size_t calculateSize ( const size_t dims );

		public:

		/** Construct a vector of given dimensions. */
		explicit AutoPermutation ( const size_t dims );

		/** Construct a vector of given dimensions with possibly some extra space for easy growth. */
		AutoPermutation ( const size_t dims, const size_t allocateDims );

		/** Copy constructor. */
		AutoPermutation ( const AutoPermutation& original );

		/** Assignment operator. */
		AutoPermutation& operator= ( const AutoPermutation& original );

		/** Resize to the given dimensions.
		* Values disappear when logical dimensions shrink.
		*/
		void exactSize ( const size_t dims );

		/** Resize to accommodate the given dimensions and possibly extra space for easy growth.
		* Values will disappear when dimensions shrink.
		*/
		void upSize ( const size_t dims );

		/** Destruct, particularly free the internal <code>size_t</code> array. */
		virtual ~AutoPermutation ();
	};

}

#endif	/* LINALG_AUTOPERMUTATION_HPP */
