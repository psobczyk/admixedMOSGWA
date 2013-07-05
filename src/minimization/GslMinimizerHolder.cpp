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

#include "GslMinimizerHolder.hpp"

using namespace std;
using namespace util;

namespace minimization {

	GslMinimizerHolder::GslMinimizerHolder (
		const size_t dim
	) : ReferenceCounted<gsl_multimin_fdfminimizer*>(
		new GslMinimizerHolder::Instance( dim )
	) {}

	GslMinimizerHolder::Instance::Instance (
		const size_t dim
	) : ReferenceCounted<gsl_multimin_fdfminimizer*>::Instance(
		gsl_multimin_fdfminimizer_alloc(
			//gsl_multimin_fdfminimizer_steepest_descent,
			//gsl_multimin_fdfminimizer_conjugate_fr,
			//gsl_multimin_fdfminimizer_conjugate_pr,
			//gsl_multimin_fdfminimizer_vector_bfgs,
			gsl_multimin_fdfminimizer_vector_bfgs2,
			dim
		)
	) {}

	GslMinimizerHolder::Instance::~Instance () {
		gsl_multimin_fdfminimizer* &gslMinimizer( *this );
		gsl_multimin_fdfminimizer_free( gslMinimizer );
		gslMinimizer = NULL;
	}

}
