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

#ifndef MINIMIZATION_MINIMIZABLE_HPP
#define MINIMIZATION_MINIMIZABLE_HPP

#include "../linalg/Vector.hpp"

namespace minimization {

	/** Wraps a multi-parametric function and its derivative in a way suitable for minimum search.
	* @see Minimizer
	* @see http://www.gnu.org/software/gsl/manual/html_node/Multidimensional-Minimization.html
	*/
	class Minimizable {

		public:

		/** Get the number of logical dimensions required for the vector <code>x</code> in
		* {@link Minimizable::calculateFunction},
		* {@link Minimizable::calculateDerivative}
		* and
		* {@link Minimizable::calculateFunctionAndDerivative}.
		*/
		virtual size_t countDimensions () const = 0;

		/** Calculate the function. */
		virtual void calculateFunction ( const linalg::Vector& x, double& functionResult ) = 0;

		/** Calculate the gradient of the function. */
		virtual void calculateDerivative ( const linalg::Vector& x, linalg::Vector& derivativeResult ) = 0;

		/** Calculate both the function and its derivative.
		* This may be implemented more efficient
		* than separate calls to {@link Minimizable::calculateFunction} and {@link Minimizable::calculateDerivative},
		* which is the default implementation.
		*/
		virtual void calculateFunctionAndDerivative (
			const linalg::Vector& x,
			double& functionResult,
			linalg::Vector& derivativeResult
		);
	};

}

#endif	/* MINIMIZATION_MINIMIZABLE_HPP */
