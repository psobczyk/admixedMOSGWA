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

#include "Hdf5Constants.hpp"

namespace io {

	const char
		* const Hdf5Constants::snpListPath = "/single_nucleotide_polymorphisms",
		* const Hdf5Constants::individualListPath = "/individuals",
		* const Hdf5Constants::genotypeMatrixPath = "/genome_matrix",
		* const Hdf5Constants::covariateListPath = "/covariates",
		* const Hdf5Constants::covariateMatrixPath = "/covariate_matrix",
		* const Hdf5Constants::phenotypeVectorPath = "/phenotypes";

}
