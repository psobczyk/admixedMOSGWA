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
