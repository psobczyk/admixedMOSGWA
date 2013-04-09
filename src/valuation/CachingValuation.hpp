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
