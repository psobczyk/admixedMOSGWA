#ifndef _PARSLETTEMPLATE_CPP_
#define _PARSLETTEMPLATE_CPP_

#include "ParsletTemplate.hpp"

namespace parser {

	template <class T> ParsletTemplate<T>::ParsletTemplate ( T &variable ) : variable( variable ) {}

	template <class T> T ParsletTemplate<T>::get () const { return variable; }

	template <class T> void ParsletTemplate<T>::set ( const T value ) { variable = value; }

}

#endif /* _PARSLETTEMPLATE_CPP_ */
