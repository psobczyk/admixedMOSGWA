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

#include "Id.hpp"

using namespace std;

namespace hdf5 {

	Id::Id ( const hid_t id, const string& name ) : id( id ), name( name ) {}

	hid_t Id::getId () const {
		return id;
	}

	const string& Id::getName () const {
		return name;
	}

	/** If the abstract destructor has no implementation, g++ yields:
	* "undefined reference to `hdf5::Id::~Id()'".
	*/
	Id::~Id () {}

}
