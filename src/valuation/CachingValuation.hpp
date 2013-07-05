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

#ifndef VALUATION_CACHINGVALUATION_HPP
#define VALUATION_CACHINGVALUATION_HPP

#include "Valuation.hpp"
#include "../util/Cache.hpp"

namespace valuation {

	/** Caching decorator for {@link valuation::Valuation}s.
	* Useful to avoid duplicate calculation.
	*/
	class CachingValuation : public Valuation, private util::Cache<lookup::ModelIndex,double> {

		/** Wraps a {@link valuation::Valuation} to access it as a {@link util::Cache::Retriever}. */
		class Retriever : public util::Cache<lookup::ModelIndex,double>::Retriever {

			/** From where to retrieve value upon cache miss. */
			Valuation& valuation;

			public:

			/** Reformulate a valuation as a retriever. */
			Retriever ( Valuation& valuation );

			/** Retrieve upon cache miss. */
			double retrieve ( const lookup::ModelIndex& modelIndex );

		} retriever;

		public:

		/** Decorate a given valuation. */
		CachingValuation ( Valuation& decorated );

		/** Judge a given model. */
		virtual double valuate ( const lookup::ModelIndex& modelIndex );
	};

}

#endif	/* VALUATION_CACHINGVALUATION_HPP */
