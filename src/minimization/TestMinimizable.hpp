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

#ifndef MINIMIZATION_TESTMINIMIZABLE_HPP
#define MINIMIZATION_TESTMINIMIZABLE_HPP

#include "Minimizable.hpp"
#include "../linalg/AutoVector.hpp"

namespace test {

	/** A concretization for testing {@link minimization::Minimizable} and {@link minimization::Minimizer}. */
	class TestMinimizable : public minimization::Minimizable {

		private:

		/** Stores the minimum vector. */
		linalg::AutoVector minVec;

		public:

		/** Construct and initialize an instance. */
		TestMinimizable ( const linalg::Vector& minVec );

		/** @overrides Minimizable::countDimensions */
		virtual size_t countDimensions () const;

		/** Override the calculation of the function. */
		virtual void calculateFunction ( const linalg::Vector& x, double& functionResult );

		/** Override the calculation of the gradient. */
		virtual void calculateDerivative ( const linalg::Vector& x, linalg::Vector& derivativeResult );
	};

}

#endif	/* MINIMIZATION_TESTMINIMIZABLE_HPP */
