#include "Search.hpp"

using namespace io;
using namespace lookup;

namespace search {

	Search::Search ( Input& input, lookup::ResultStore& resultStore ) : input( input ), resultStore( resultStore ) {}

	Search::~Search () {}

}
