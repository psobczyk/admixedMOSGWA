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

#ifndef LINALG_AUTOMATRIX_HPP
#define LINALG_AUTOMATRIX_HPP

#include <ostream>
#include "Matrix.hpp"

namespace linalg {

	/** A cheap C++ wrapper for a gently resizing GSL matrix.
	* As the logical and the physical sizes of the matrix generally differ,
	* the logical matrix is a view of a sub-matrix of the physical storage array.
	*/
	class AutoMatrix : public Matrix {

		/** The currently allocated memory block size. */
		size_t size;

		protected:

		/** Calculate memory requirement. */
		static size_t calculateSize ( const size_t rows, const size_t cols );

		public:

		/** Construct a matrix of given dimensions. */
		AutoMatrix ( const size_t rows, const size_t cols );

		/** Construct a matrix of given dimensions with possibly some extra space for easy growth. */
		AutoMatrix (
			const size_t rows,
			const size_t cols,
			const size_t allocateRows,
			const size_t allocateCols
		);

		/** Copy constructor. */
		AutoMatrix ( const AutoMatrix& original );

		/** Assignment operator. */
		AutoMatrix& operator= ( const AutoMatrix& original );

		/** Resize to the given dimensions.
		* Values disappear when logical dimensions shrink.
		*/
		void exactSize ( const size_t rows, const size_t cols );

		/** Resize to accommodate the given dimensions and possibly extra space for easy growth.
		* Values will disappear when dimensions shrink.
		*/
		void upSize ( const size_t rows, const size_t cols );

		/** Add a row after the last and return it as modifiable vector.
		* The vector is only valid until the next modification of the matrix dimensions.
		*/
		Vector newRow ();

		/** Add a column after the last and return it as modifiable vector.
		* The vector is only valid until the next modification of the matrix dimensions.
		*/
		Vector newColumn ();

		/** Destruct, particularly free the internal <code>double</code> array. */
		virtual ~AutoMatrix ();
	};

}

#endif	/* LINALG_AUTOMATRIX_HPP */
