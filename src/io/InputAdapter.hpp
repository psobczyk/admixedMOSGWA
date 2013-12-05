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

#ifndef IO_INPUTADAPTER_HPP
#define IO_INPUTADAPTER_HPP

#include "Input.hpp"
#include <string>
#include <vector>

namespace io {

	/** Provides some shared functionality to implement {@link Input} classes.
	* Use it if you want to program a new {@link Input} class.
	*/
	class InputAdapter : public Input {

		protected:

		/** Holds information about individuals. */
		std::vector<Individual> individuals;

		/** Holds SNP information. */
		std::vector<SNP> snps;

		/** Holds covariate information. */
		std::vector<std::string> covariates;

		/** Holds trait information. */
		std::vector<std::string> traits;

		/** Construct without descriptor information.
		* That information must be added later using the protected fields.
		*/
		InputAdapter ();

		public:

		/** Construct with the necessary descriptor information. */
		InputAdapter (
			const std::vector<Individual>& individuals,
			const std::vector<SNP>& snps,
			const std::vector<std::string>& covariates,
			const std::vector<std::string>& traits
		);

		/** Return the number of individuals in the data. */
		virtual size_t countIndividuals () const;

		/** Retrieve the descriptions of all individuals. */
		virtual const Individual * getIndividuals () const;

		/** Return the number of SNPs in the data. */
		virtual size_t countSnps () const;

		/** Retrieve the descriptions of all SNPs. */
		virtual const SNP * getSnps () const;

		/** Return the number of covariate vectors in the data. */
		virtual size_t countCovariates () const;

		/** Retrieve the names of all covariates. */
		virtual const std::string * getCovariates () const;

		/** Return the number of traits,  i.e. phenotype vectors in the data. */
		virtual size_t countTraits () const;

		/** Retrieve the names of all traits. */
		virtual const std::string * getTraits () const;

		/** Declare access to be finished, release all allocated resources. */
		virtual ~InputAdapter ();
	};

}

#endif	/* IO_INPUTADAPTER_HPP */
