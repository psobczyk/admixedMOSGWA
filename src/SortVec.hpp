/********************************************************************************
 *	This file is part of the MOSGWA program code.				*
 *	Copyright ©2011–2013, Erich Dolejsi, Bernhard Bodenstorfer.		*
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

#ifndef SORT_VEC_HPP
#define SORT_VEC_HPP

#include <vector>
#include <algorithm>

using namespace std;

/** A sorted list of (size_t id, double value) pairs.
* Access by index, not by iterators.
*/
class SortVec {

	private:

	/** One pair to sort */
	struct SortItem {

		/** Used for SNP id */
		const size_t id;

		/** Used for P-value */
		const double value;

		/** Construct a pair */
		SortItem ( const size_t id, const double value );
	};

	friend bool order_function ( const SortVec::SortItem* i, const SortVec::SortItem* j );

friend bool order_function2 ( const SortVec::SortItem* i, const SortVec::SortItem* j );
	/** Vector of sortable pairs */
	vector <SortItem*> list_;

	/** Forbid use of the (default) copy constructor. */
	SortVec ( const SortVec& original );

	/** Forbid use of the (default) assignment operator. */
	SortVec& operator= ( const SortVec& original );

	protected:

	/** Clear all entries. */
	void clear ();

	public:

	/** Default constructor */
	SortVec ();

	/** Constructor with a given size */
	SortVec ( const size_t n );

	/** Constructor with arrays containing position and values, so that  ids[i], values[i] are a pair */
	SortVec ( const size_t n, const size_t ids[], const double values[], bool bigger = true );

	/** Destructor */
	~SortVec ();

	/** Clear and overwrite SortVec with arrays */
	void fillVec ( const size_t n, const size_t ids[], const double values[], bool bigger = true );

	/** Get position (or id_) for the k-th smallest value */
	size_t getId ( const size_t k ) const;

	/** Get the k-th smallest value */
	double getValue ( const size_t k ) const;
};

#endif	/* SORT_VEC_HPP */
