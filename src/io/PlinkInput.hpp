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

#ifndef IO_PLINKINPUT_HPP
#define IO_PLINKINPUT_HPP

#include "InputCo.hpp"
#include <vector>
#include "../linalg/AutoVector.hpp"
#include "../linalg/AutoMatrix.hpp"

namespace io {

	/** Reads input data from PLink BIM, FAM and BED formatted files. */
	class PlinkInput : public InputCo {

		/** File extensions for the Plink files for SNPs, Individuals and genome. */
		static const char
			* const snpListExtension,
			* const individualListExtension,
			* const genotypeMatrixExtension,
			* const covariateMatrixExtension;

		/** Translation table from two genome bits to the number for the regression matrix entry. */
		static const double genotypeTranslation[];

		/** Stores all SNP characteristics in memory. */
		std::vector<SNP> snpList;

		/** Stores all Individual characteristics in memory. */
		std::vector<Individual> individualList;

		/** Stores all covariate names in memory. */
		std::vector<std::string> covariateList;

		/** Genome data matrix.
		* It is stored as transposed matrix in order to optimise memory access.
		* The vectors holding all individuals' data for one SNP should be in a contiguous piece of memory
		* so that the whole vector fits into cache.
		* It is these vectors which are added to the regression algorithm in the search phase.
		*/
		linalg::AutoMatrix genotypeMatrixTransposed;

		/** Phenotype data vector. */
		linalg::AutoVector phenotypeVector;

		/** Covariate matrix.
		* Similar to {@link #genomeMatrixTransposed},
		* it is stored as transposed matrix in order to optimise memory access.
		*/
		linalg::AutoMatrix covariateMatrixTransposed;

		public:

		/** Set up the Plink file reading. */
		PlinkInput ( const char* const filename );

		/** Return the number of SNPs in the data. */
		virtual size_t countSnps () const;

		/** Return the number of individuals in the data. */
		virtual size_t countIndividuals () const;

		/** Retrieve the data for the given SNP. */
		virtual SNP getSnp ( const size_t snpIndex );

		/** Retrieve the data for the given individual. */
		virtual Individual getIndividual ( const size_t individualIndex );

		/** Copy the {@link countIndividuals} sized vector of genotype information for the given SNP into the given vector. */
		virtual void retrieveGenotypesIntoVector ( const size_t snpIndex, linalg::Vector& vector );

		/** Copy the {@link countIndividuals} sized vector of phenotype information into the given vector. */
		virtual void retrievePhenotypesIntoVector ( linalg::Vector& vector );

		/** Return the number of covariate vectors in the data. */
		virtual size_t countCovariateVectors () const;

		/** Get the name of the given covariate. */
		virtual std::string getCovariateName ( const size_t covIndex ) const;

		/** Copy the {@link countIndividuals} sized vector of covariate information for the given covariate into the given vector. */
		virtual void retrieveCovariatesIntoVector ( const size_t covIndex, linalg::Vector& vector );

		/** Declare access to be finished, release all resources. */
		virtual ~PlinkInput ();
	};

}

#endif	/* IO_PLINKINPUT_HPP */
