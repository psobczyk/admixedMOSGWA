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

#ifndef BUILDINFO_HPP
#define BUILDINFO_HPP

/** Build information.
* The actual values are populated by the build system.
* @author Bernhard Bodenstorfer
*/
namespace buildinfo {

	/** The build-time as generated from the <code>cmake</code> makefile. */
	extern const char * timestamp;

}

#endif	/* BUILDINFO_HPP */
