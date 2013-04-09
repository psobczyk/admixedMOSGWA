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

	string Permutation::toString () const {
		ostringstream s;
		s << *this;
		return s.str();
	}

	void Permutation::print () const {
		cout << *this;
	}

	Permutation::~Permutation () {}

	std::ostream& operator<< ( std::ostream& s, const Permutation& p ) {
		for ( size_t dim = 0; dim < p.countDimensions(); ++dim ) {
			if ( 0 < dim ) s << '\t';
			s << p.get( dim );
		}
		s << endl;
		return s;
	}

}
