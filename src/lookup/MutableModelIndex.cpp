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

#include "MutableModelIndex.hpp"

using namespace lookup;

MutableModelIndex::MutableModelIndex ( const MutableModelIndex& original ) : ModelIndex( original ) {}

MutableModelIndex::MutableModelIndex ( const ModelIndex& original ) : ModelIndex( original ) {}

MutableModelIndex& MutableModelIndex::operator= ( const MutableModelIndex& original ) {
	return operator=( static_cast<const ModelIndex&>( original ) );
}

MutableModelIndex& MutableModelIndex::operator= ( const ModelIndex& original ) {
	ModelIndex::operator=( original );
	return *this;
}
