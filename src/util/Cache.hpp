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

#ifndef UTIL_CACHE_HPP
#define UTIL_CACHE_HPP

#include <map>

/** Provides general utilities for use in more than one package.
* @author Bernhard Bodenstorfer
*/
namespace util {

	/** Provides cache logic. */
	template <class Key, class Value> class Cache {

		/** Holds the cached key-value pairs. */
		std::map<Key,Value> cache;

		public:

		/** Upon cache miss, provides values which should then be cached. */
		struct Retriever {

			/** Retrieve a value, when it is not in the cache. */
			virtual Value retrieve ( const Key& key ) = 0;
		};

		/** Give a reference to the object for a given key.
		* @param retriever to provide a value, when it is not in the cache.
		*/
		Value& get ( const Key& key, Retriever& retriever );
	};

}

// C++ compilers want template functions source code available when the template functions are used
#include "Cache.cpp"

#endif	/* UTIL_CACHE_HPP */
