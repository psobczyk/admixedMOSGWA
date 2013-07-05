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

#include "AutoPermutation.hpp"
#include <string.h>
#include <assert.h>
#include "package.hpp"

using namespace std;

namespace linalg {

	size_t AutoPermutation::calculateSize ( const size_t dims ) {
		if ( 0 == dims ) return 0;	// guard against division by zero in assert below
		size_t required = dims * sizeof( size_t );
		assert( required / dims == sizeof( size_t ) );	// Guard against integer overflow
		return required;
	}

	AutoPermutation::AutoPermutation (
		const size_t dims
	) : Permutation(
		dims,
		// if not here, size would be initialised after parent class portion
		static_cast<size_t*>( malloc( size = calculateSize( dims ) ) )
	) {}

	AutoPermutation::AutoPermutation (
		const size_t dims,
		const size_t allocateDims
	) : Permutation(
		dims,
		static_cast<size_t*>( malloc( size = calculateSize( allocateDims ) ) )
	) {}

	AutoPermutation::AutoPermutation ( const AutoPermutation& original ) : Permutation(
		original.countDimensions(),
		static_cast<size_t*>( malloc( original.size ) )
	), size( original.size ) {
		const size_t used = calculateSize( countDimensions() );
		memcpy( data, original.data, used );
	}

	AutoPermutation& AutoPermutation::operator= ( const AutoPermutation& original ) {
		if ( this != &original ) {
			const size_t dims = original.countDimensions();
			const size_t used = calculateSize( dims );
			size_t* newArray;
			if ( used <= size ) {
				newArray = data;
			} else {
				newArray = static_cast<size_t*>( malloc( used ) );
				if ( NULL == newArray ) {
					throw bad_alloc();
				}
				free( data );
				data = newArray;
				size = used;
			}
			memcpy( data, original.data, used );
			gslInit( dims, newArray );
		}
		return *this;
	}

	void AutoPermutation::exactSize ( const size_t dims ) {
		const size_t
			required = calculateSize( dims ),
			used = calculateSize( countDimensions() );
		size_t *newArray;
		// Use realloc only if there is no risk to copy too many unused bytes
		// to get done in linear time with respect to min( used, required ).
		// I assume that realloc has constant overhead if size == required, i.e. no change.
		if ( required <= size || used > size >> 1 ) {
			newArray = static_cast<size_t*>( realloc( data, required ) );
			if ( 0 < dims && NULL == newArray ) {
				throw bad_alloc();
			}
		} else {
			newArray = static_cast<size_t*>( malloc( required ) );
			if ( NULL == newArray ) {
				throw bad_alloc();
			}
			memcpy( newArray, data, used );		// note: required > used
			free( data );
		}
		gslInit( dims, newArray );	// update GSL struct variables
		size = required;
	}

	/** The method guarantees that the physical dimensions will not shrink.
	* Moreover, if growth is necessary, it will be in geometrically increasing steps.
	*/
	void AutoPermutation::upSize ( const size_t dims ) {
		if ( size < calculateSize( dims ) ) {
			exactSize( upperPowerOf2( dims ) );
		}
		gslInit( dims, data );
	}

	AutoPermutation::~AutoPermutation () {
		free( data );
		data = NULL;
	}

}
