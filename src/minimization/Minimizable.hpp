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
