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

#ifndef IO_INPUT_HPP
#define IO_INPUT_HPP

#include "../SNP.hpp"
#include "../Individual.hpp"
#include "../linalg/Vector.hpp"

/** Provides input/output exchange of regression data.
* @author Bernhard Bodenstorfer
*/
namespace io {

	/** Interface to provide the regression data input. */
	struct Input {

		/** Return the number of SNPs in the data. */
		virtual size_t countSnps () = 0;

		/** Return the number of individuals in the data. */
		virtual size_t countIndividuals () = 0;

		/** Retrieve the data for the given SNP. */
		virtual SNP getSnp ( const size_t snpIndex ) = 0;

		/** Retrieve the data for the given individual. */
		virtual Individual getIndividual ( const size_t individualIndex ) = 0;

		/** Copy the {@link countIndividuals} sized vector of genotype information for the given SNP into the given vector. */
		virtual void retrieveGenotypesIntoVector ( const size_t snpIndex, linalg::Vector& vector ) = 0;

		/** Copy the {@link countIndividuals} sized vector of phenotype information into the given vector. */
		virtual void retrievePhenotypesIntoVector ( linalg::Vector& vector ) = 0;

		/** Declare access to be finished, release all resources. */
		virtual ~Input ();
	};

}

#endif	/* IO_INPUT_HPP */
