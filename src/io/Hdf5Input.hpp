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

#include "InputCo.hpp"
#include "Hdf5Object.hpp"

namespace io {

	/** Reads input data from HDF5 file format. */
	class Hdf5Input : virtual public Input {

		protected:

		/** Paths of the relevant objects in the HDF5 file. */
		static const char
			* const snpListPath,
			* const individualListPath,
			* const genotypeMatrixPath,
			* const phenotypeVectorPath;

		/** Holds the identifier of the opened HDF5-file. */
		Hdf5FileId fileId;

		/** Gets axis data names into MOSGWA. */
		Hdf5StringList
			snpNames,
			individualNames;

		/** Gets genotype data into MOSGWA. */
		Hdf5DoubleTable genotypesTransposed;

		/** Gets phenotype data into MOSGWA. */
		Hdf5DoubleList phenotypes;

		public:

		/** Set up the HDF5 file reading. */
		Hdf5Input ( const char * const filename );

		/** Return the number of SNPs in the data. */
		virtual size_t countSnps ();

		/** Return the number of individuals in the data. */
		virtual size_t countIndividuals ();

		/** Retrieve the data for the given SNP. */
		virtual SNP getSnp ( const size_t snpIndex );

		/** Retrieve the data for the given individual. */
		virtual Individual getIndividual ( const size_t individualIndex );

		/** Copy the {@link countIndividuals} sized vector of genotype information for the given SNP into the given vector. */
		virtual void retrieveGenotypesIntoVector ( const size_t snpIndex, linalg::Vector& v );

		/** Copy the {@link countIndividuals} sized vector of phenotype information into the given vector. */
		virtual void retrievePhenotypesIntoVector ( linalg::Vector& v );

		/** Declare access to be finished, release all resources. */
		virtual ~Hdf5Input ();
	};

}

#endif	/* IO_HDF5INPUT_HPP */
