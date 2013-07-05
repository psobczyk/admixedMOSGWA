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

#include "Minimizer.hpp"
#include "GslMinimizerRetriever.hpp"
#include <assert.h>
#include <math.h>	// for NaN
#include <memory>	// for auto_ptr

using namespace std;
using namespace linalg;
using namespace util;

namespace minimization {

	Minimizer::Minimizer (
		const double initialStep,
		const double intermediateTolerance,
		const unsigned maxIterations
	) : initialStep( initialStep ), intermediateTolerance( intermediateTolerance ), maxIterations( maxIterations ) {
	}

	double Minimizer::minimize ( Minimizable& minimizable, Vector& x ) {
		const size_t dimensions = minimizable.countDimensions();
		assert( dimensions == x.countDimensions() );

		if ( 0 == dimensions ) {
			double result;
			minimizable.calculateFunction( x, result );
			return result;
		}

		// Cache-get the multidimensional minimizer
		GslMinimizerRetriever retriever;
		GslMinimizerHolder& gslMinimizerRef( gslMinimizers.get( dimensions, retriever ) );
		gsl_multimin_fdfminimizer* const gslMinimizer( gslMinimizerRef );

		// Prepare GSL minimization data struct
		gsl_multimin_function_fdf gslMinimizable;
		gslMinimizable.f = Minimizer::gslF;
		gslMinimizable.df = Minimizer::gslDF;
		gslMinimizable.fdf = Minimizer::gslFDF;
		gslMinimizable.n = dimensions;
		gslMinimizable.params = &minimizable;

		// Actual minimization
		gsl_multimin_fdfminimizer_set (
			gslMinimizer,
			&gslMinimizable,
			&x.vector,
			initialStep,
			intermediateTolerance
		);
		int status = GSL_CONTINUE;
		for (
			int iteration = 0;
			iteration < maxIterations && GSL_CONTINUE == status;
			++iteration
		) {
			const int gslError = gsl_multimin_fdfminimizer_iterate( gslMinimizer );
			if ( gslError ) break;
			status = gsl_multimin_test_gradient( gslMinimizer->gradient, 1e-3 );
		}
		if ( GSL_CONTINUE == status ) {
			return nan( "minimize failed" );
		} else {
			x.copy( Vector( *gslMinimizer->x ) );
			return gslMinimizer->f;
		}
	}

	Minimizer::~Minimizer () {
	}

	double Minimizer::gslF ( const gsl_vector* x, void* params ) {
		Minimizable& minimizable( *(Minimizable*) params );
		const Vector point( *x );
		double result;
		minimizable.calculateFunction( point, result );
		return result;
	}

	void Minimizer::gslDF ( const gsl_vector* x, void* params, gsl_vector* g ) {
		Minimizable& minimizable( *(Minimizable*) params );
		const Vector point( *x );
		Vector gradient( *g );
		minimizable.calculateDerivative( point, gradient );
	}

	void Minimizer::gslFDF ( const gsl_vector* x, void* params, double* f, gsl_vector* g ) {
		Minimizable& minimizable( *(Minimizable*) params );
		const Vector point( *x );
		Vector gradient( *g );
		minimizable.calculateFunctionAndDerivative( point, *f, gradient );
	}

}
