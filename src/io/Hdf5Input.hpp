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

		/** Paths of the relevant objects in the HDF5 file. */
		static const char
			* const snpListPath,
			* const individualListPath,
			* const genotypeMatrixPath,
			* const covariateListPath,
			* const covariateMatrixPath,
			* const phenotypeVectorPath;

		/** Holds the identifier of the opened HDF5-file. */
		hdf5::File file;

		/** Gets phenotype data into MOSGWA. */
		hdf5::DoubleList phenotypes;

		/** Gets genotype data into MOSGWA. */
		hdf5::DoubleTable genotypesTransposed;

		/** Gets covariate matrix into MOSGWA. */
		std::auto_ptr<hdf5::DoubleTable> covariatesTransposedPtr;

		public:

		/** Set up the HDF5 file reading.
		* @param useCovariates indicates whether covariates should be read.
		*/
		Hdf5Input ( const char * const filename, const bool useCovariates = false );

		/** Copy the {@link countIndividuals} sized vector of phenotype information into the given vector. */
		virtual void retrievePhenotypeVector ( linalg::Vector& v );

		/** Copy the {@link countIndividuals} sized vector of genotype information for the given SNP into the given vector. */
		virtual void retrieveGenotypeVector ( const size_t snpIndex, linalg::Vector& v );

		/** Copy the {@link countIndividuals} sized vector of covariate information for the given covariate into the given vector. */
		virtual void retrieveCovariateVector ( const size_t covIndex, linalg::Vector& vector );

		/** Declare access to be finished, release all resources. */
		virtual ~Hdf5Input ();
	};

}

#endif	/* IO_HDF5INPUT_HPP */
