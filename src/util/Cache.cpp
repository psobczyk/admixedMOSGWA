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

#ifndef UTIL_CACHE_CPP
#define UTIL_CACHE_CPP

#include "Cache.hpp"

using namespace std;

namespace util {

	/*
	* Check-and-insert idiom inspired by: Chris Yester-Young
	* @see http://stackoverflow.com/questions/3886593/how-to-check-if-stdmap-contains-a-key-without-doing-insert
	*/
	template <class Key, class Value> Value& Cache<Key,Value>::get ( const Key& key, Retriever& retriever ) {
		typename map<Key,Value>::iterator iterator( cache.lower_bound( key ) );
		if ( cache.end() == iterator || key < iterator->first ) {
			// not found
			Value retrievedValue( retriever.retrieve( key ) );
			// hinted insertion
			cache.insert( iterator, make_pair( key, retrievedValue ) );
			// point iterator to the added value
			--iterator;
		}
		Value& cachedValue( iterator->second );
		return cachedValue;
	}

}

#endif	/* UTIL_CACHE_CPP */
