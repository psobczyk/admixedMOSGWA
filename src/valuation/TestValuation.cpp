#include "TestValuation.hpp"
#include <cmath>
#include <algorithm>

using namespace std;
using namespace lookup;

namespace test {

	TestValuation::TestValuation () {}

	TestValuation::TestValuation ( const std::map<ModelIndex,double>& mappings )
		: mappings( mappings ) {
	}

	void TestValuation::put ( const lookup::ModelIndex& modelIndex, const double value ) {
		mappings[ modelIndex ] = value;
	}

	double TestValuation::valuate ( const lookup::ModelIndex& modelIndex ) {
		map<ModelIndex,double>::iterator iterator( mappings.lower_bound( modelIndex ) );
		if ( mappings.end() == iterator || modelIndex != iterator->first ) {
			return nan( "unset" );
		}
		return iterator->second;
	}

}
