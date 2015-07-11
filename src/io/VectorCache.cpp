/********************************************************************************
 *	This file is part of the MOSGWA program code.				*
 *	Copyright ©2015, Bernhard Bodenstorfer.				*
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

#include "VectorCache.hpp"

using namespace linalg;

namespace io {

	VectorCache::VectorCache ( const size_t dim, const size_t limit )
	: AutoMatrix( 0, dim ), limit( limit ), toCheck( lookup.begin() ) {
		// TODO: the vectors are currently stored as columns in the AutoMatrix.
		// When matrices will be column major, do not forget to un-transpose.
	}

	void VectorCache::store ( const size_t index, const Vector& vector ) {
		if ( 0 == limit ) return;
		const size_t size = lookup.size();
		size_t slot;	// free slot, to be determined
		bool hadToSearch = false;
		size_t whereToContinueSearch;
		LookupMap::iterator entry = lookup.find( index );
		if ( lookup.end() != entry ) {
			slot = entry->second.slot;	// case: overwrite cached vector
		} else if ( size < limit ) {
			slot = countRows();	// TODO<BB>: don't forget when transpose
			upSize( slot + 1, countColumns() );
		} else {
			hadToSearch = true;
			while ( true ) {
				if ( lookup.end() == toCheck ) {
					toCheck = lookup.begin();
				}
				// Divide the usage count by two and decide between continue or erase
				if ( toCheck->second.used >>= 1 ) {
					++toCheck;
				} else {
					whereToContinueSearch = toCheck->first;
					slot = toCheck->second.slot;	// this is the freed slot in the matrix
					lookup.erase( toCheck );
					break;
				}
			}
		}

		Vector target = rowVector( slot );
		target.copy( vector );

		Descriptor& descriptor = lookup[ index ];
		descriptor.slot = slot;
		descriptor.used = 1;

		if ( hadToSearch ) {
			toCheck = lookup.upper_bound( whereToContinueSearch );
		}
	}

	bool VectorCache::retrieve ( const size_t index, Vector& target ) {
		LookupMap::iterator entry = lookup.find( index );
		if ( lookup.end() != entry ) {
			++( entry->second.used );
			const Vector vector = rowVector( entry->second.slot );
			target.copy( vector );
			return true;
		}
		return false;
	}

}
