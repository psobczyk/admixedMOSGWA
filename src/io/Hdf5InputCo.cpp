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

#include "Hdf5InputCo.hpp"
#include "../Exception.hpp"
#include <cassert>
#include <hdf5.h>

using namespace std;
using namespace linalg;

namespace io {

	const char
		* const Hdf5InputCo::covariateListPath = "/covariates",
		* const Hdf5InputCo::covariateMatrixPath = "/covariate_matrix";

	Hdf5InputCo::Hdf5InputCo ( const char * const filename )
		:
		Hdf5Input( filename ),
		covariateNames( fileId, covariateListPath ),
		covariatesTransposed( fileId, covariateMatrixPath )
	{
		// Tests for input data consistency
		if ( countIndividuals() != covariatesTransposed.countColumns() ) {
			throw Exception(
				"HDF5 input file \"%s\""
				" dataset has %u rows"
				" and transposed covariates dataset \"%s\" has minor %u dimensions"
				" but both numbers should be equal.",
				filename,
				countIndividuals(),
				covariateMatrixPath,
				covariatesTransposed.countColumns()
			);
		}
		if ( covariateNames.countDimensions() != covariatesTransposed.countRows() ) {
			throw Exception(
				"HDF5 input file \"%s\""
				" dataset \"%s\" has %u dimensions"
				" and dataset \"%s\" has %u dimensions"
				" but both numbers should be equal.",
				filename,
				covariateListPath,
				covariateNames.countDimensions(),
				covariateMatrixPath,
				covariatesTransposed.countRows()
			);
		}
	}

	size_t Hdf5InputCo::countCovariateVectors () {
		return covariateNames.countDimensions();
	}

	string Hdf5InputCo::getCovariateName ( const size_t covIndex ) {
		const string covariateId = covariateNames.readOne( covIndex );
		return covariateId;
	}

	void Hdf5InputCo::retrieveCovariatesIntoVector ( const size_t covIndex, linalg::Vector& v ) {
		const size_t
			rows = covariatesTransposed.countColumns(),	// equals rows
			covs = covariatesTransposed.countRows();	// equals covs
		assert( covIndex < covs );
		assert( rows == v.countDimensions() );
		vector<double> array( rows );
		covariatesTransposed.readRow( covIndex, array.data() );
		v.fill( array.data() );
	}

	Hdf5InputCo::~Hdf5InputCo () {}

}
