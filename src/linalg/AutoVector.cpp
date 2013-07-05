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

#include "AutoVector.hpp"
#include <string.h>
#include <assert.h>
#include "package.hpp"

using namespace std;

namespace linalg {

	size_t AutoVector::calculateSize ( const size_t dims ) {
		if ( 0 == dims ) return 0;	// guard against division by zero in assert below
		size_t required = dims * sizeof( double );
		assert( required / dims == sizeof( double ) );	// Guard against integer overflow
		return required;
	}

	AutoVector::AutoVector (
		const size_t dims
	) : Vector(
		dims,
		// if not here, size would be initialised after parent class portion
		static_cast<double*>( malloc( size = calculateSize( dims ) ) ),
		1u
	) {}

	AutoVector::AutoVector (
		const size_t dims,
		const size_t allocateDims
	) : Vector(
		dims,
		static_cast<double*>( malloc( size = calculateSize( allocateDims ) ) ),
		1u
	) {}

	AutoVector::AutoVector ( const AutoVector& original ) : Vector(
		original.countDimensions(),
		static_cast<double*>( malloc( original.size ) ),
		1u
	), size( original.size ) {
		const size_t used = calculateSize( countDimensions() );
		memcpy( vector.data, original.vector.data, used );
	}

	AutoVector& AutoVector::operator= ( const AutoVector& original ) {
		if ( this != &original ) {
			const size_t dims = original.countDimensions();
			const size_t used = calculateSize( dims );
			double* newArray;
			if ( used <= size ) {
				newArray = vector.data;
			} else {
				newArray = static_cast<double*>( malloc( used ) );
				if ( NULL == newArray ) {
					throw bad_alloc();
				}
				free( vector.data );
				vector.data = newArray;
				size = used;
			}
			memcpy( vector.data, original.vector.data, used );
			gslInit( dims, newArray );
		}
		return *this;
	}

	void AutoVector::exactSize ( const size_t dims ) {
		const size_t
			required = calculateSize( dims ),
			used = calculateSize( countDimensions() );
		double *newArray;
		// Use realloc only if there is no risk to copy too many unused bytes
		// to get done in linear time with respect to min( used, required ).
		// I assume that realloc has constant overhead if size == required, i.e. no change.
		if ( required <= size || used > size >> 1 ) {
			newArray = static_cast<double*>( realloc( vector.data, required ) );
			if ( 0 < dims && NULL == newArray ) {
				throw bad_alloc();
			}
		} else {
			newArray = static_cast<double*>( malloc( required ) );
			if ( NULL == newArray ) {
				throw bad_alloc();
			}
			memcpy( newArray, vector.data, used );	// note: required > used
			free( vector.data );
		}
		gslInit( dims, newArray );	// update GSL struct variables
		size = required;
	}

	/** The method guarantees that the physical dimensions will not shrink.
	* Moreover, if growth is necessary, it will be in geometrically increasing steps.
	*/
	void AutoVector::upSize ( const size_t dims ) {
		if ( size < calculateSize( dims ) ) {
			exactSize( upperPowerOf2( dims ) );
		}
		gslInit( dims, vector.data );
	}

	AutoVector::~AutoVector () {
		free( vector.data );
		vector.data = NULL;
	}

}
