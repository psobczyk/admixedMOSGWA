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

#ifndef IO_HDF5OUTPUT_HPP
#define IO_HDF5OUTPUT_HPP

#include "Output.hpp"
#include "hdf5/StringList.hpp"
#include "hdf5/DoubleList.hpp"
#include "hdf5/DoubleTable.hpp"

namespace io {

	/** Stores regression input data in HDF5 file format. */
	class Hdf5Output : public Output {

		protected:

		/** Relevant dimensions, i.e. maximum indices. */
		const size_t individualCount, snpCount, covariateCount;

		/** Holds the identifier of the opened HDF5-file. */
		hdf5::File file;

		/** Puts individual identifiers into HDF5. */
		hdf5::StringList individualList;

		/** Puts SNP identifiers into HDF5. */
		hdf5::StringList snpList;

		/** Puts covariate identifiers into HDF5. */
		hdf5::StringList covariateList;

		/** Puts phenotype vector into HDF5. */
		hdf5::DoubleList phenotypes;

		/** Puts genotype matrix into HDF5. */
		hdf5::DoubleTable genotypesTransposed;

		/** Puts covariate matrix into HDF5. */
		hdf5::DoubleTable covariatesTransposed;

		public:

		/** Set up the HDF5 file writing.
		* @param filename where to write to.
		*/
		Hdf5Output (
			const char * const filename,
			const size_t individualCount,
			const size_t snpCount,
			const size_t covariateCount
		);

		/** Store the data for the individuals.	 */
		virtual void setIndividuals ( const Individual * individuals );

		/** Store the <code>individualCount</code> sized vector of phenotype information. */
		virtual void storePhenotypeVector ( const linalg::Vector& v );

		/** Store the data for the SNPs.  */
		virtual void setSnps ( const SNP * snps );

		/** Store the <code>individualCount</code> sized vector of genotype information for the given SNP. */
		virtual void storeGenotypeVector ( const size_t snpIndex, const linalg::Vector& v );

		/** Set the names of the covariates. */
		virtual void setCovariates ( const std::string * covariates );

		/** Copy the {@link countIndividuals} sized vector of covariate information for the given covariate into the given vector. */
		virtual void storeCovariateVector ( const size_t covIndex, const linalg::Vector& v );

		/** Declare access to be finished, release all resources. */
		virtual ~Hdf5Output ();
	};

}

#endif	/* IO_HDF5OUTPUT_HPP */
