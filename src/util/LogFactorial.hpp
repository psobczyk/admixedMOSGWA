/********************************************************************************
 *	This file is part of the MOSGWA program code.				*
 *	Copyright ©2014, Bernhard Bodenstorfer.					*
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

#ifndef UTIL_LOGFACTORIAL_HPP
#define UTIL_LOGFACTORIAL_HPP

#include <cstdlib>
#include <vector>

namespace util {

	/** Calculates the logarithm of the factorial. */
	class LogFactorial {

		std::vector<double> cache;

		public:

		/** Construct and initialise. */
		LogFactorial ();

		/** Return \f$\log n!\f$. */
		double logFactorial ( const size_t n );

		/** Return \f$\log \left( n \choose k \right)\f$. */
		double logChoose ( const size_t n, const size_t k );
	};

}

#endif	/* UTIL_LOGFACTORIAL_HPP */
