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

#ifndef MOSGWA_TYPES_H
#define MOSGWA_TYPES_H

#include <stdint.h>

/** The index identifying an {@link SNP}. */
typedef uint32_t snp_index_t;

/** Maximum number of SNPs in a {@link Model}.
* Models will typically be much smaller.
*/
#define MAX_MODEL_SIZE	4294967295U

#endif	/* MOSGWA_TYPES_H */
