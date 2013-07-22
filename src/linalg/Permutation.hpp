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

#ifndef LINALG_PERMUTATION_HPP
#define LINALG_PERMUTATION_HPP

#include <string>
#include <ostream>
#include <gsl/gsl_permutation.h>

namespace linalg {

	/** A cheap C++ wrapper for a fixed-size GSL permutation,
	* which functions as a permutation view.
	* The class itself does not do any memory allocation or deallocation;
	* the subclass {@link AutoPermutation} takes care of that.
	* @see http://www.gnu.org/software/gsl/manual/html_node/Permutations.html
	*/
	class Permutation : protected gsl_permutation {

		/** Used to allow GSL-level access for coordinate-permutation. */
		friend class Vector;

		/** Used to allow GSL-level access for QRT-factorisation. */
		friend class Matrix;

		// Note: copy constructor and assignment operator are allowed
		// and equal default ones for a plain Permutation.

		protected:

		/** Set permutation from a given array.
		* The current pointer will be overwritten, so take care to avoid memory leaks!
		* If the Permutation is an {@link AutoPermutation},
		* this should be called only with base pointing to the current <code>permutation.base</code>
		* or when proper freeing of the current is guaranteed.
		* Similarly for the destructor to work,
		* an {@link AutoPermutation} should only get a <code>base</code> which has been malloc'ed.
		* It will be free'd by {@link AutoPermutation::~AutoPermutation()},
		* but not by {@link Permutation::~Permutation()}!
		* @returns whether any existing views might be invalidated by the resize
		*/
		bool gslInit ( const size_t dims, size_t *base );

		/** Construct a permutation of given dimensions from an array. */
		Permutation ( const size_t dims, size_t *base );

		public:

		/** Construct from a GSL permutation, but do not take ownership.
		* Data will be shared.
		* Hence, the gsl_permutation must continue to exist and eventually be properly freed.
		*/
		Permutation ( const gsl_permutation& permutation );

		/** Compare whether this equals that. */
		bool operator== ( const Permutation& that ) const;

		/** Initialise the permutation as identity. */
		void init ();

		/** Fill the permutation with values from an array.
		* The fill affects the logical portion of the permutation,
		* i.e. the array must be at least of size {@link countDimensions()}.
		*/
		void fill ( const size_t *array );

		/** Copy the values from another permutation.
		* Both permutations must have the same number of (logical) dimensions.
		*/
		void copy ( const Permutation& that );

		/** Get the number of logical dimensions. */
		size_t countDimensions () const;

		/** Get the element for a given dimension. */
		size_t get ( const size_t dim ) const;

		/** Set the element for a given dimension.
		* Note that this method allows to change a permutation into a more general integer mapping.
		*/
		void set ( const size_t dim, const size_t mappedDim );

		/** Swap the two elements for the given dimensions. */
		void swap ( const size_t dim1, const size_t dim2 );

		/** Create a string representation. Mainly for debugging. */
		std::string toString () const;

		/** Print to <code>cout</code>. Mainly for debugging. */
		void print () const;

		/** Placeholder destructor for {@link AutoPermutation::~AutoPermutation}. */
		virtual ~Permutation ();
	};

	/** Output a {@link Permutation}. */
	std::ostream& operator<< ( std::ostream& s, const Permutation& p );

}

#endif	/* LINALG_PERMUTATION_HPP */
