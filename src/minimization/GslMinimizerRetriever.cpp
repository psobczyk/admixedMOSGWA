#include "GslMinimizerRetriever.hpp"

namespace minimization {

	GslMinimizerHolder GslMinimizerRetriever::retrieve ( const size_t& dim ) {
		GslMinimizerHolder gslMinimizerHolder( dim );
		return gslMinimizerHolder;
	}

}
