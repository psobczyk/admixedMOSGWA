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

#ifndef IO_HDF5INPUT_HPP
#define IO_HDF5INPUT_HPP

#include "InputAdapter.hpp"
#include "hdf5/StringList.hpp"
#include "hdf5/DoubleList.hpp"
#include "hdf5/DoubleTable.hpp"
#include <memory>

namespace io {

	/** Reads input data from HDF5 file format. */
	class Hdf5Input : public InputAdapter {

		protected:

		/** Holds the identifier of the opened HDF5-file. */
		hdf5::File file;

		/** Gets phenotype vector from HDF5. */
		hdf5::DoubleList phenotypes;

		/** Gets genotype matrix from HDF5. */
		hdf5::DoubleTable genotypesTransposed;

		/** Gets covariate matrix from HDF5. */
		std::auto_ptr<hdf5::DoubleTable> covariatesTransposedPtr;

		public:

		/** Set up the HDF5 file reading.
		* @param useCovariates indicates whether covariates should be read.
		*/
		Hdf5Input ( const char * const filename, const bool useCovariates = false );

		/** Copy the {@link countIndividuals} sized vector of genotype information for the given SNP into the given vector. */
		virtual void retrieveGenotypeVector ( const size_t snpIndex, linalg::Vector& v );

		/** Copy the {@link countIndividuals} sized vector of covariate information for the given covariate into the given vector. */
		virtual void retrieveCovariateVector ( const size_t covIndex, linalg::Vector& v );

		/** Copy a {@link countIndividuals} sized vector of phenotype information into the given vector. */
		virtual void retrievePhenotypeVector ( const size_t traitIndex, linalg::Vector& v );

		/** Declare access to be finished, release all resources. */
		virtual ~Hdf5Input ();
	};

}

#endif	/* IO_HDF5INPUT_HPP */
