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
