#include "package.hpp"
#include <iostream>

using namespace std;

namespace linalg {

	/** If <code>required</code>==0, 0 is returned. */
	size_t upperPowerOf2 ( const size_t required ) {
		if ( 0 >= required ) return 0;
		size_t s = required - 1, p = 1;
		while ( s ) {
			s >>= 1;
			p <<= 1;
		}
		return p;
	}


	void printVector ( const Vector& v ) {
		cout << v;
	}

	void printMatrix ( const Matrix& m ) {
		cout << m;
	}

	void printPermutation ( const Permutation& p ) {
		cout << p;
	}

}
