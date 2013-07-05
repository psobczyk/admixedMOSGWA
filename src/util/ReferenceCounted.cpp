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

#ifndef UTIL_REFERENCECOUNTED_CPP
#define UTIL_REFERENCECOUNTED_CPP

#include "ReferenceCounted.hpp"

#include <assert.h>

namespace util {

	template <class T> ReferenceCounted<T>::Instance::Instance (
		const T& object
	) : referenceCount( 0 ), object( object ) {
	}

	template <class T> ReferenceCounted<T>::Instance::operator T& () {
		return object;
	}

	template <class T> ReferenceCounted<T>::Instance::~Instance () {
		assert( 0 == referenceCount );	// Otherwise destruction does not make sense
	}

	template <class T> ReferenceCounted<T>::ReferenceCounted (
		const T& object
	) : instance( new Instance( object ) ) {
		++instance->referenceCount;
	}

	template <class T> ReferenceCounted<T>::ReferenceCounted (
		Instance * instance
	) : instance( instance ) {
		++instance->referenceCount;
	}

	template <class T> ReferenceCounted<T>::ReferenceCounted (
		const ReferenceCounted &original
	) : instance( original.instance ) {
		++instance->referenceCount;
	}

	template <class T> ReferenceCounted<T>& ReferenceCounted<T>::operator= (
		const ReferenceCounted &that
	) {
		if ( instance != that.instance ) {
			if ( 0 == --instance->referenceCount ) {
				delete instance;
			}
			instance = that.instance;
			++instance->referenceCount;
		}
		return *this;
	}

	template <class T> ReferenceCounted<T>::operator T& () {
		return *instance;
	}

	template <class T> ReferenceCounted<T>::~ReferenceCounted () {
		if ( 0 == --instance->referenceCount ) {
			delete instance;
		}
		instance = NULL;
	}

}

#endif	/* UTIL_REFERENCECOUNTED_CPP */
