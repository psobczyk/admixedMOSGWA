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

#include "PlinkConstants.hpp"
#include <cmath>	// for nan(...)

namespace io {

	const char
		* const PlinkConstants::snpListExtension = ".bim",
		* const PlinkConstants::individualListExtension = ".fam",
		* const PlinkConstants::genotypeMatrixExtension = ".bed",
		* const PlinkConstants::covariateMatrixExtension = ".cov",
		* const PlinkConstants::phenotypeMatrixExtension = ".yvm";

	/** Mind http://pngu.mgh.harvard.edu/~purcell/plink/binary.shtml specifies in the long translation example that the bits are swapped.
	* E.g. 10 is stored with LSB 1! Therefore <code>genotypeTranslation[1]</code> yields <code>NaN</code>, not <code>genotypeTranslation[2]</code>!
	*/
	const double PlinkConstants::genotypeTranslation[4] = {
		-1.0,	// 00 homocygote
		::nan( "missing" ),	// 10 missing (stored as bit pattern 01 written with leading MSB)
		0.0,	// 01 heterocygote (stored as bit pattern 10 written with leading MSB)
		+1.0	// 11 homocygote
	};

	const char PlinkConstants::bedFileMagic[3] = {
		0x6c,
		0x1b,
		0x1	// SNP-majour mode (used for writing) is important for performance and so that write does not have to bit-fiddle with data already in the file
	};

}
