#include "CachingValuation.hpp"

using namespace std;
using namespace lookup;
using namespace util;

namespace valuation {

	CachingValuation::Retriever::Retriever ( Valuation& valuation )
		: valuation( valuation ) {
	}

	double CachingValuation::Retriever::retrieve ( const ModelIndex& modelIndex ) {
		const double value = valuation.valuate( modelIndex );
		return value;
	}

	CachingValuation::CachingValuation ( Valuation& decorated )
		: retriever( decorated ) {
	}

	double CachingValuation::valuate ( const lookup::ModelIndex& modelIndex ) {
		const double value = get( modelIndex, retriever );
		return value;
	}

}
