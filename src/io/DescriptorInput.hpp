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

#ifndef IO_DESCRIPTORINPUT_HPP
#define IO_DESCRIPTORINPUT_HPP

#include "../SNP.hpp"
#include "../Individual.hpp"

/** Provides input/output exchange of regression data.
* @author Bernhard Bodenstorfer
*/
namespace io {

	/** Interface to provide the regression data input. */
	struct DescriptorInput {

		/** Return the number of individuals in the data. */
		virtual size_t countIndividuals () const = 0;

		/** Retrieve the descriptions of the individuals.
		* The pointer returned points into the <code>DescriptorInput</code> object
		* and must not be used after the lifetime of <code>this</code>.
		*/
		virtual const Individual* getIndividuals () const = 0;

		/** Return the number of SNPs in the data. */
		virtual size_t countSnps () const = 0;

		/** Retrieve the descriptions of the SNPs.
		* The pointer returned points into the <code>DescriptorInput</code> object
		* and must not be used after the lifetime of <code>this</code>.
		*/
		virtual const SNP* getSnps () const = 0;

		/** Return the number of covariate vectors in the data. */
		virtual size_t countCovariates () const = 0;

		/** Retrieve the names of the covariates.
		* The pointer returned points into the <code>DescriptorInput</code> object
		* and must not be used after the lifetime of <code>this</code>.
		*/
		virtual const std::string* getCovariates () const = 0;

		/** Return the number of traits i.e. phenotype vectors in the data. */
		virtual size_t countTraits () const = 0;

		/** Retrieve the names of the phenotype traits.
		* The pointer returned points into the <code>DescriptorInput</code> object
		* and must not be used after the lifetime of <code>this</code>.
		*/
		virtual const std::string* getTraits () const = 0;

		/** Declare access to be finished, release all resources. */
		virtual ~DescriptorInput ();
	};

}

#endif	/* IO_DESCRIPTORINPUT_HPP */
