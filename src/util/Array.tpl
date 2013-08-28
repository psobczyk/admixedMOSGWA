/********************************************************************************
 *	This file is part of the MOSGWA program code.				*
 *	Copyright ©2013, Bernhard Bodenstorfer.					*
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

#ifndef _UTIL_ARRAY_TPL_
#define _UTIL_ARRAY_TPL_

#include "Array.hpp"
#include <cassert>

namespace util {

	template<class E,size_t D> Array<E,D>::Array () {
		assert( 0 == D );
	}

	template<class E,size_t D> Array<E,D>::Array ( E initial, ... ) {
		assert( 0 < D );
		array[ 0 ] = initial;
		va_list args;
		va_start( args, initial );
		init( 1, args );
		va_end( args );
	}

	template<class E,size_t D> Array<E,D>::Array ( ::va_list args ) {
		init( 0, args );
	}

	template<class E,size_t D> void Array<E,D>::init ( const size_t begin, ::va_list args ) {
		for ( size_t i = begin; i < D; ++i ) {
			array[i] = va_arg( args, E );
		}
	}

	template<class E,size_t D> Array<E,D>::operator E* () {
		return array;
	}

}

#endif	/* _UTIL_ARRAY_TPL_ */
