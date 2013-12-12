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

#ifndef IO_IOHELPER_HPP
#define IO_IOHELPER_HPP

#include "DescriptorInput.hpp"
#include "DescriptorOutput.hpp"
#include "../linalg/Vector.hpp"
#include <string>
#include <vector>

namespace io {

	/** Holds data descriptor information to implement {@link Input} and {@kink Output} classes.
	* This class is extended to {@link InputAdapter} and {@link #OutputAdapter}.
	*/
	class IOHelper : public virtual DescriptorInput, public virtual DescriptorOutput {

		template<class T> static void set ( const T * sources, std::vector<T>& targets );

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
		IOHelper ();

		public:

		/** Construct with the necessary dimension information. */
		IOHelper (
			const size_t individualCount,
			const size_t snpCount,
			const size_t covariateCount,
			const size_t traitCount
		);

		/** Construct with the necessary descriptor information. */
		IOHelper (
			const std::vector<Individual>& individuals,
			const std::vector<SNP>& snps,
			const std::vector<std::string>& covariates,
			const std::vector<std::string>& traits
		);

		/** Return the number of individuals in the data. */
		virtual size_t countIndividuals () const;

		/** Retrieve the descriptions of all individuals. */
		virtual const Individual * getIndividuals () const;

		/** Store the data for the individuals.
		* Note: The constructor should set the expected number of individuals.
		*/
		virtual void setIndividuals ( const Individual * individuals );

		/** Return the number of SNPs in the data. */
		virtual size_t countSnps () const;

		/** Retrieve the descriptions of all SNPs. */
		virtual const SNP * getSnps () const;

		/** Store the data for the SNPs.
		* Note: The constructor should set the expected number of SNPs.
		*/
		virtual void setSnps ( const SNP * snps );

		/** Return the number of covariate vectors in the data. */
		virtual size_t countCovariates () const;

		/** Retrieve the names of all covariates. */
		virtual const std::string * getCovariates () const;

		/** Set the names of the covariates.
		* Note: The constructor should set the expected number of covariates.
		*/
		virtual void setCovariates ( const std::string * covariates );

		/** Return the number of traits,  i.e. phenotype vectors in the data. */
		virtual size_t countTraits () const;

		/** Retrieve the names of all traits. */
		virtual const std::string * getTraits () const;

		/** Set the names of the phenotype traits.
		* Note: The constructor should set the expected number of covariates.
		*/
		virtual void setTraits ( const std::string * traits );

		/** Declare access to be finished, release all allocated resources. */
		virtual ~IOHelper ();
	};

}

#endif	/* IO_IOHELPER_HPP */
