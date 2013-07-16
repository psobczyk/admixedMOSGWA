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

#ifndef IO_INPUT_CO_HPP
#define IO_INPUT_CO_HPP

#include "Input.hpp"

namespace io {

	/** Interface adding covariables to regression data input. */
	struct InputCo : public Input {

		/** Return the number of covariate vectors in the data. */
		virtual size_t countCovariateVectors () const = 0;

		/** Copy the {@link countIndividuals} sized vector of covariate information for the given covariate into the given vector. */
		virtual void retrieveCovariatesIntoVector ( const size_t covIndex, linalg::Vector& vector ) = 0;

		/** Declare access to be finished, release all resources. */
		virtual ~InputCo ();
	};

}

#endif	/* IO_INPUT_CO_HPP */
