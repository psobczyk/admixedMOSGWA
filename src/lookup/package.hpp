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

#ifndef PACKAGE_LOOKUP_HPP
#define PACKAGE_LOOKUP_HPP

#include "ModelIndex.hpp"

/** Facilities to store already calculated {@link Model} selection criteria for later reference
* to avoid unnecessary recalculation.
* @author Bernhard Bodenstorfer
*/
namespace lookup {

	/** Output to <code>cout</code> for testing and debugging. */
	void printModelIndex ( const ModelIndex& mi );

}

#endif	/* PACKAGE_LOOKUP_HPP */
