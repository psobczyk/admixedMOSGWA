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
