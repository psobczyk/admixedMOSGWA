/********************************************************************************
 *	This file is part of the MOSGWA program code.				*
 *	Copyright ©2015, Bernhard Bodenstorfer.				*
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

#ifndef IO_CACHEDINPUT_HPP
#define IO_CACHEDINPUT_HPP

#include "Input.hpp"
#include "VectorCache.hpp"

namespace io {

	/** Forwards requests to an <code>Input</code> and caches the vector type replies in memory. */
	class CachedInput : public Input {

		/** The wrapped input. */
		Input& input;

		/** Cache of genotype vectors.
		* Currently covariate matrix and phenotype are cached by <code>PlinkInput</code> itself.
		*/
		VectorCache genotypeCache;

		public:

		/** Wrap the source to be cached.
		* @param input source to be wrapped
		* @param limit how many vectors the cache may hold
		* TODO: Implement pre-fetch.
		*/
		CachedInput (
			Input& input,
			const size_t limit
		);

		/** Forwarded to the <code>input</code>. */
		virtual size_t countIndividuals () const;

		/** Forwarded to the <code>input</code>. */
		virtual const Individual* getIndividuals () const;

		/** Forwarded to the <code>input</code>. */
		virtual size_t countSnps () const;

		/** Forwarded to the <code>input</code>. */
		virtual const SNP* getSnps () const;

		/** Forwarded to the <code>input</code>. */
		virtual size_t countCovariates () const;

		/** Forwarded to the <code>input</code>. */
		virtual const std::string* getCovariates () const;

		/** Forwarded to the <code>input</code>. */
		virtual size_t countTraits () const;

		/** Forwarded to the <code>input</code>. */
		virtual const std::string* getTraits () const;

		/** Copy the {@link countIndividuals} sized vector of genotype information for the given SNP into the given vector.
		* Either from the cache or by forwarding the request to the wrapped <code>input</code>.
		*/
		virtual void retrieveGenotypeVector ( const size_t snpIndex, linalg::Vector& vector );

		/** Forwarded to the <code>input</code>. */
		virtual void retrievePhenotypeVector ( const size_t traitIndex, linalg::Vector& vector );

		/** Forwarded to the <code>input</code>. */
		virtual void retrieveCovariateVector ( const size_t covIndex, linalg::Vector& vector );
	};

}

#endif	/* IO_CACHEDINPUT_HPP */
