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

#ifndef IO_OUTPUTADAPTER_HPP
#define IO_OUTPUTADAPTER_HPP

#include "Output.hpp"
#include "IOHelper.hpp"

namespace io {

	/** Provides some shared functionality to implement {@link Output} classes.
	* Use it if you want to program a new {@link Output} class.
	*/
	class OutputAdapter : public Output, protected IOHelper {

		protected:

		/** Construct without descriptor information.
		* That information must be added later using the protected fields.
		*/
		OutputAdapter ();

		public:

		using IOHelper::setIndividuals;
		using IOHelper::setSnps;
		using IOHelper::setCovariates;
		using IOHelper::setTraits;

		/** Construct with the necessary dimension information. */
		OutputAdapter (
			const size_t individuals,
			const size_t snps,
			const size_t covariates,
			const size_t traits
		);

		/** Construct with the necessary descriptor information. */
		OutputAdapter (
			const std::vector<Individual>& individuals,
			const std::vector<SNP>& snps,
			const std::vector<std::string>& covariates,
			const std::vector<std::string>& traits
		);

		/** Declare access to be finished, release all allocated resources. */
		virtual ~OutputAdapter ();
	};

}

#endif	/* IO_OUTPUTADAPTER_HPP */
