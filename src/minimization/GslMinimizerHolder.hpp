#ifndef MINIMIZATION_GSLMINIMIZERHOLDER_HPP
#define MINIMIZATION_GSLMINIMIZERHOLDER_HPP

#include "../util/ReferenceCounted.hpp"
#include <gsl/gsl_multimin.h>

namespace minimization {

	/** Ties allocation/deallocation logic to a reference-counted <code>gsl_multimin_fdfminimizer*</code>. */
	class GslMinimizerHolder : public util::ReferenceCounted<gsl_multimin_fdfminimizer*> {

		protected:

		/** Provides <code>gsl_multimin_fdfminimizer</code> allocation and deallocation using GSL code. */
		struct Instance : public util::ReferenceCounted<gsl_multimin_fdfminimizer*>::Instance {

			/** Construct <code>gsl_multimin_fdfminimizer</code> for given problem dimension. */
			Instance ( const size_t dim );

			/** Destruct <code>gsl_multimin_fdfminimizer</code> using GSL code. */
			virtual ~Instance ();

		};

		public:

		/** Allocate a <code>gsl_multimin_fdfminimizer</code> to be referenced. */
		GslMinimizerHolder ( const size_t dim );
	};

}

#endif	/* MINIMIZATION_GSLMINIMIZERHOLDER_HPP */
