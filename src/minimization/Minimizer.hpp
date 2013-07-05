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

#ifndef MINIMIZATION_MINIMIZER_HPP
#define MINIMIZATION_MINIMIZER_HPP

#include "Minimizable.hpp"
#include "GslMinimizerHolder.hpp"
#include "../linalg/Vector.hpp"
#include "../util/Cache.hpp"

namespace minimization {

	/** Provides minimum search.
	* This class uses an internal pool of variables and is not thread-safe.
	* @see Minimizable
	* @see http://www.gnu.org/software/gsl/manual/html_node/Multidimensional-Minimization.html
	*/
	class Minimizer {

		private:

		static double gslF ( const gsl_vector* x, void* params );
		static void gslDF ( const gsl_vector* x, void* params, gsl_vector* g );
		static void gslFDF ( const gsl_vector* x, void* params, double* f, gsl_vector* g );

		protected:

		/** The initial step width for each call of {@link minimize}. */
		const double initialStep;

		/** The tolerance for intermediate partial minimisation results. */
		const double intermediateTolerance;

		/** How many iterations the minimum search loop should permit. */
		const unsigned maxIterations;

		/** Holds minimizers for different dimensions to be used and re-used. */
		util::Cache<const size_t, GslMinimizerHolder> gslMinimizers;

		public:

		/** Construct an instance.
		* @param initialStep specifies the initial search step distance
		* @param intermediateTolerance allows for imprecisions in intermediate loops
		* @param maxIterations bounds the number of iterations before giving up if no minimum is found
		*/
		Minimizer (
			const double initialStep,
			const double intermediateTolerance,
			const unsigned maxIterations
		);

		/** Minimize a given function, starting from the point x.
		* The dimensionalities of minimizable and x must match.
		* @param minimizable wraps the function to be minimised and its derivative.
		* @param x determines the starting point, and upon termination the point of the minimum.
		* @return the minimum, if one was found; <code>NaN</code> otherwise.
		*/
		double minimize ( Minimizable& minimizable, linalg::Vector& x );

		/** Frees all internally allocated resources. */
		~Minimizer ();
	};

}

#endif	/* MINIMIZATION_MINIMIZER_HPP */
