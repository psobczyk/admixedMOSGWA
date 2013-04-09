#include "GslMinimizerHolder.hpp"

using namespace std;
using namespace util;

namespace minimization {

	GslMinimizerHolder::GslMinimizerHolder (
		const size_t dim
	) : ReferenceCounted<gsl_multimin_fdfminimizer*>(
		new GslMinimizerHolder::Instance( dim )
	) {}

	GslMinimizerHolder::Instance::Instance (
		const size_t dim
	) : ReferenceCounted<gsl_multimin_fdfminimizer*>::Instance(
		gsl_multimin_fdfminimizer_alloc(
			//gsl_multimin_fdfminimizer_steepest_descent,
			//gsl_multimin_fdfminimizer_conjugate_fr,
			//gsl_multimin_fdfminimizer_conjugate_pr,
			//gsl_multimin_fdfminimizer_vector_bfgs,
			gsl_multimin_fdfminimizer_vector_bfgs2,
			dim
		)
	) {}

	GslMinimizerHolder::Instance::~Instance () {
		gsl_multimin_fdfminimizer* &gslMinimizer( *this );
		gsl_multimin_fdfminimizer_free( gslMinimizer );
		gslMinimizer = NULL;
	}

}
