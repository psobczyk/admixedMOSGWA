#ifndef MINIMIZATION_GSLMINIMIZERRETRIEVER_HPP
#define MINIMIZATION_GSLMINIMIZERRETRIEVER_HPP

#include "../util/Cache.hpp"
#include "GslMinimizerHolder.hpp"

namespace minimization {

	/** Provides <code>gsl_multimin_fdfminimizer</code>s for problems of given dimensionality as cache feed. */
	struct GslMinimizerRetriever : public util::Cache<const size_t, minimization::GslMinimizerHolder>::Retriever {

		/** Allocate the asked-for <code>gsl_multimin_fdfminimizer</code>. */
		virtual GslMinimizerHolder retrieve ( const size_t& dim );
	};

}

#endif	/* MINIMIZATION_GSLMINIMIZERRETRIEVER_HPP */
