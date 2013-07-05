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

#ifndef _PARSLETTEMPLATE_CPP_
#define _PARSLETTEMPLATE_CPP_

#include "ParsletTemplate.hpp"

namespace parser {

	template <class T> ParsletTemplate<T>::ParsletTemplate ( T &variable ) : variable( variable ) {}

	template <class T> T ParsletTemplate<T>::get () const { return variable; }

	template <class T> void ParsletTemplate<T>::set ( const T value ) { variable = value; }

}

#endif /* _PARSLETTEMPLATE_CPP_ */
