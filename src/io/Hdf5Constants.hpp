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

#ifndef IO_HDF5CONSTANTS_HPP
#define IO_HDF5CONSTANTS_HPP

namespace io {

	/** Shared constants for the HDF5 file format. */
	namespace Hdf5Constants {

		/** Paths of the relevant objects in the HDF5 file. */
		extern const char
			* const snpListPath,
			* const individualListPath,
			* const genotypeMatrixPath,
			* const covariateListPath,
			* const covariateMatrixPath,
			* const phenotypeVectorPath;
	}

}

#endif	/* IO_HDF5CONSTANTS_HPP */
