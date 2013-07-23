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

#include "Permutation.hpp"
#include <iostream>
#include <sstream>
#include <string.h>
#include <assert.h>

using namespace std;

namespace linalg {

	bool Permutation::gslInit ( const size_t dims, size_t *base ) {
		bool invalidates
			= base != data
			|| dims < countDimensions();
		assert( 0 <= dims );
		data = base;
		size = dims;
		return invalidates;
	}

	Permutation::Permutation ( const size_t dims, size_t *base ) {
		gslInit( dims, base );
	}

	Permutation::Permutation ( const gsl_permutation& permutation ) : gsl_permutation( permutation ) {}

	bool Permutation::operator== ( const Permutation& that ) const {
		// later GSL:
		// return 0 != gsl_permutation_equal( this, &that );
		// earlier GSL: (not performance optimised, but anyway used only for testing)
		// TODO: make it more efficient
		const size_t dims = countDimensions();
		if ( dims != that.countDimensions() ) return false;
		for ( size_t dim = 0; dim < dims; ++dim ) {
			if ( get( dim ) !=  that.get( dim ) ) return false;
		}
		return true;
	}

	void Permutation::init () {
		gsl_permutation_init( this );
	}

	bool Permutation::isPermutation () const {
		const bool isP = gsl_permutation_valid( this );
		return isP;
	}

	void Permutation::fill ( const size_t *array ) {
		// A harmless snappy const_cast circumvents a complete class "ConstPermutation"
		Permutation that( countDimensions(), const_cast<size_t*>( array ) );
		copy( that );
	}

	void Permutation::copy ( const Permutation& that ) {
		gsl_permutation_memcpy( this, &that );
	}

	size_t Permutation::countDimensions () const {
		return size;
	}

	size_t Permutation::get ( const size_t dim ) const {
		return gsl_permutation_get( this, dim );
	}

	void Permutation::set ( const size_t dim, const size_t mappedDim ) {
		const size_t dimensions = countDimensions();
		assert( dim < dimensions );
		data[dim] = mappedDim;	// GSL has no set method for good reason.
	}

	void Permutation::swap ( const size_t dim1, const size_t dim2 ) {
		gsl_permutation_swap( this, dim1, dim2 );
	}

	/**
	* @see https://www.gnu.org/software/gsl/manual/html_node/Sorting-vectors.html
	* @see http://www.codecogs.com/library/computing/c/stdlib.h/qsort.php
	* @see http://stackoverflow.com/questions/4300896/how-portable-is-the-re-entrant-qsort-r-function-compared-to-qsort
	*/
	void Permutation::sort ( const Vector& vector ) {
		gsl_sort_vector_index( this, &vector.vector );
	}

	string Permutation::toString () const {
		ostringstream s;
		s << *this;
		return s.str();
	}

	void Permutation::print () const {
		cout << *this;
	}

	Permutation::~Permutation () {}

	ostream& operator<< ( ostream& s, const Permutation& p ) {
		for ( size_t dim = 0; dim < p.countDimensions(); ++dim ) {
			if ( 0 < dim ) s << '\t';
			s << p.get( dim );
		}
		s << endl;
		return s;
	}

}
