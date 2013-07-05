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

#ifndef MINIMIZATION_GSLMINIMIZERHOLDER_HPP
#define MINIMIZATION_GSLMINIMIZERHOLDER_HPP

#include "../util/ReferenceCounted.hpp"
#include <gsl/gsl_multimin.h>

namespace minimization {

	/** Ties allocation/deallocation logic to a reference-counted <code>gsl_multimin_fdfminimizer*</code>. */
	class GslMinimizerHolder : public util::ReferenceCounted<gsl_multimin_fdfminimizer*> {

		protected:

		/** Provides <code>gsl_multimin_fdfminimizer</code> allocation and deallocation using GSL code. */
		struct Instance : public util::ReferenceCounted<gsl_multimin_fdfminimizer*>::Instance {

			/** Construct <code>gsl_multimin_fdfminimizer</code> for given problem dimension. */
			Instance ( const size_t dim );

			/** Destruct <code>gsl_multimin_fdfminimizer</code> using GSL code. */
			virtual ~Instance ();

		};

		public:

		/** Allocate a <code>gsl_multimin_fdfminimizer</code> to be referenced. */
		GslMinimizerHolder ( const size_t dim );
	};

}

#endif	/* MINIMIZATION_GSLMINIMIZERHOLDER_HPP */
