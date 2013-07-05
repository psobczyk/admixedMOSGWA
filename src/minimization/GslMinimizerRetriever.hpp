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

#ifndef MINIMIZATION_GSLMINIMIZERRETRIEVER_HPP
#define MINIMIZATION_GSLMINIMIZERRETRIEVER_HPP

#include "../util/Cache.hpp"
#include "GslMinimizerHolder.hpp"

namespace minimization {

	/** Provides <code>gsl_multimin_fdfminimizer</code>s for problems of given dimensionality as cache feed. */
	struct GslMinimizerRetriever : public util::Cache<const size_t, minimization::GslMinimizerHolder>::Retriever {

		/** Allocate the asked-for <code>gsl_multimin_fdfminimizer</code>. */
		virtual GslMinimizerHolder retrieve ( const size_t& dim );
	};

}

#endif	/* MINIMIZATION_GSLMINIMIZERRETRIEVER_HPP */
