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

#ifndef IO_PLINKOUTPUT_HPP
#define IO_PLINKOUTPUT_HPP

#include "OutputAdapter.hpp"
#include "../linalg/AutoMatrix.hpp"
#include <fstream>

namespace io {

	/** Stores regression input data in <code>plink</code> file format. */
	class PlinkOutput : public OutputAdapter {

		/** Used to initialise genotype file with not-available. */
		static const char fillBytes[];

		linalg::AutoMatrix
			covariateMatrixTransposed,
			phenotypeMatrixTransposed;

		/** Common leading portion of all output files. */
		std::string filenameTrunk;

		/** Genotype data goes here. */
		std::ofstream genotypeStream;

		/** Actual genotype data is stored in bed-files from here onwards. */
		std::streampos genotypeOrigin;

		public:

		/** Construct with the necessary dimension information. */
		PlinkOutput (
			const char * const filenameTrunk,
			const size_t individualCount,
			const size_t snpCount,
			const size_t covariateCount,
			const size_t traitCount
		);

		/** Store the data for the SNPs. */
		virtual void setSnps ( const SNP * snps );

		/** Store the <code>individualCount</code> sized vector of genotype information for the given SNP. */
		virtual void storeGenotypeVector ( const size_t snpIndex, const linalg::Vector& v );

		/** Store the <code>individualCount</code> sized vector of covariate information for the given covariate. */
		virtual void storeCovariateVector ( const size_t covIndex, const linalg::Vector& v );

		/** Store the <code>individualCount</code> sized vector of phenotype information for the given trait. */
		virtual void storePhenotypeVector ( const size_t traitIndex, const linalg::Vector& v );

		/** Declare access to be finished,
		* write out all remaining data and release all allocated resources.
		*/
		virtual ~PlinkOutput ();
	};

}

#endif	/* IO_PLINKOUTPUT_HPP */
