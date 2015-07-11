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

#ifndef IO_VECTORCACHE_HPP
#define IO_VECTORCACHE_HPP

#include "../linalg/AutoMatrix.hpp"
#include <map>

namespace io {

	/** Caches vectors in an {@link linalg::AutoMatrix}.
	* @author Bernhard Bodenstorfer
	*/
	class VectorCache : private linalg::AutoMatrix {

		/** Descriptor of cache entry. */
		struct Descriptor {
			/** Slot of the deposited vector in the storage matrix. */
			size_t slot;

			/** Indicates how often the vector has been used. More roughly means more often. */
			unsigned int used;
		};

		/** Can store by <code>index</code> all descriptors of cache entries. */
		typedef std::map<size_t,Descriptor> LookupMap;

		/** Stores by <code>index</code> all descriptors of cache entries. */
		LookupMap lookup;

		/** Limits the number of vectors to be cached. */
		const size_t limit;

		/** Index of the next vector to be checked for removal. */
		LookupMap::iterator toCheck;

		public:

		/** Allocate a cache for vectors of given equal dimension.
		* @param dim dimension of the vectors to be cached
		* @param limit how many vectors the cache may hold
		*/
		VectorCache ( const size_t dim, const size_t limit );

		/** Store a vector in the cache.
		* @param index of the vector in the totality of available vectors
		* @param vector of the vector in the totality of available vectors
		*/
		void store ( const size_t index, const linalg::Vector& vector );

		/** Try to load a vector from the cache into the given one.
		* @param index of the vector in the totality of available vectors
		* @param target where the cached values, if any, should be retrieved to
		* @returns whether the <code>index</code> has been found.
		* Otherwise, the <code>target</code> remains unchanged.
		*/
		bool retrieve ( const size_t index, linalg::Vector& target );
	};

}

#endif	/* IO_VECTORCACHE_HPP */
