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

#ifndef IO_OUTPUT_HPP
#define IO_OUTPUT_HPP

#include "../SNP.hpp"
#include "../Individual.hpp"
#include "../linalg/Vector.hpp"

namespace io {

	/** Interface to store regression data.
	* @see io::Input
	*/
	struct Output {

		/** Store the data for the individuals.
		* Note: The constructor should set the expected number of individuals.
		*/
		virtual void setIndividuals ( const Individual * individuals ) = 0;

		/** Store the <code>individualCount</code> sized vector of phenotype information. */
		virtual void storePhenotypeVector ( const linalg::Vector& v ) = 0;

		/** Store the data for the SNPs.
		* Note: The constructor should set the expected number of SNPs.
		*/
		virtual void setSnps ( const SNP * snps ) = 0;

		/** Store the <code>individualCount</code> sized vector of genotype information for the given SNP. */
		virtual void storeGenotypeVector ( const size_t snpIndex, const linalg::Vector& v ) = 0;

		/** Set the names of the covariates.
		* Note: The constructor should set the expected number of covariates.
		*/
		virtual void setCovariates ( const std::string * covariates ) = 0;

		/** Copy the <code>individualCount</code> sized vector of covariate information for the given covariate. */
		virtual void storeCovariateVector ( const size_t covIndex, const linalg::Vector& v ) = 0;

		/** Declare access to be finished, release all resources. */
		virtual ~Output ();
	};

}

#endif	/* IO_OUTPUT_HPP */
