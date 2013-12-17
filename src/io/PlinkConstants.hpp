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

#ifndef IO_PLINKCONSTANTS_HPP
#define IO_PLINKCONSTANTS_HPP

namespace io {

	/** Shared constants for the plink file format. */
	namespace PlinkConstants {

		/** File extensions for the Plink files for SNPs, Individuals and genome. */
		extern const char
			* const snpListExtension,
			* const individualListExtension,
			* const genotypeMatrixExtension,
			* const covariateMatrixExtension,
			* const phenotypeMatrixExtension;

		/** Translation table from two genome bits to the number for the regression matrix entry. */
		extern const double genotypeTranslation[4];

		/** Binary genotype file magic numbers. */
		extern const char bedFileMagic[3];

	}

}

#endif	/* IO_PLINKCONSTANTS_HPP */
