#include "Minimizable.hpp"

using namespace linalg;

namespace minimization {

	/** Default implementation is to just call
	* {@link Minimizable::calculateFunction} and {@link Minimizable::calculateDerivative} in sequence.
	*/
	void Minimizable::calculateFunctionAndDerivative (
		const Vector& x,
		double& functionResult,
		Vector& derivativeResult
	) {
		calculateFunction( x, functionResult );
		calculateDerivative( x, derivativeResult );
	}

}
