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
