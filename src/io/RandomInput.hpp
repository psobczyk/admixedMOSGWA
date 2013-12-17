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

#ifndef IO_RANDOMINPUT_HPP
#define IO_RANDOMINPUT_HPP

#include "InputAdapter.hpp"

namespace io {

	/** Generates random input data. This us useful as test data. */
	class RandomInput : public InputAdapter {

		/** Prefix plus index. */
		static std::string prefixedString ( const char prefix, const size_t i, const size_t upperBound );

		/** Generate a random vector of integers. */
		void retrieveVector (
			const size_t covIndex,
			linalg::Vector& vector,
			const int min,
			const int cases
		);

		public:

		/** Set up the random input dimensions. */
		RandomInput (
			const size_t individualCount,
			const size_t snpCount,
			const size_t covariateCount,
			const size_t traitCount
		);

		/** Generate a random vector for use as genotype data. */
		virtual void retrieveGenotypeVector ( const size_t snpIndex, linalg::Vector& vector );

		/** Generate a random vector for use as covariate data. */
		virtual void retrieveCovariateVector ( const size_t covIndex, linalg::Vector& vector );

		/** Generate a random vector for use as phenotype data. */
		virtual void retrievePhenotypeVector ( const size_t traitIndex, linalg::Vector& vector );
	};

}

#endif	/* IO_RANDOMINPUT_HPP */
