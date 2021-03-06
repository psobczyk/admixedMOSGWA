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
#include "IOHelper.hpp"

namespace io {

	/** Provides some shared functionality to implement {@link Input} classes.
	* Use it if you want to program a new {@link Input} class.
	*/
	class InputAdapter : public Input, protected IOHelper {

		protected:

		/** Construct without descriptor information.
		* That information must be added later using the protected fields.
		*/
		InputAdapter ();

		public:

		using IOHelper::countIndividuals;
		using IOHelper::getIndividuals;
		using IOHelper::countSnps;
		using IOHelper::getSnps;
		using IOHelper::countCovariates;
		using IOHelper::getCovariates;
		using IOHelper::countTraits;
		using IOHelper::getTraits;

		/** Construct with the necessary descriptor information. */
		InputAdapter (
			const std::vector<Individual>& individuals,
			const std::vector<SNP>& snps,
			const std::vector<std::string>& covariates,
			const std::vector<std::string>& traits
		);

		/** Declare access to be finished, release all allocated resources. */
		virtual ~InputAdapter ();
	};

}

#endif	/* IO_INPUTADAPTER_HPP */
