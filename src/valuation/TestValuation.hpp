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

#ifndef TEST_TESTVALUATION_HPP
#define TEST_TESTVALUATION_HPP

#include "Valuation.hpp"
#include <map>

namespace test {

	/** Test support implementation.
	* Useful to test client classes of {@link valuation::Valuation},
	* among which decorator <code>Valuation</code>s.
	*/
	class TestValuation : public valuation::Valuation {

		std::map<lookup::ModelIndex,double> mappings;

		public:

		/** Construct without mappings. */
		TestValuation ();

		/** Construct from given mappings. */
		TestValuation ( const std::map<lookup::ModelIndex,double>& mappings );

		/** Set the value for the gioven model. */
		void put ( const lookup::ModelIndex& modelIndex, const double value );

		/** Judge a given model. */
		virtual double valuate ( const lookup::ModelIndex& modelIndex );
	};

}

#endif	/* TEST_TESTVALUATION_HPP */
