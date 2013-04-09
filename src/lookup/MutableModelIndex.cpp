#include "MutableModelIndex.hpp"

using namespace lookup;

MutableModelIndex::MutableModelIndex ( const MutableModelIndex& original ) : ModelIndex( original ) {}

MutableModelIndex::MutableModelIndex ( const ModelIndex& original ) : ModelIndex( original ) {}

MutableModelIndex& MutableModelIndex::operator= ( const MutableModelIndex& original ) {
	return operator=( static_cast<const ModelIndex&>( original ) );
}

MutableModelIndex& MutableModelIndex::operator= ( const ModelIndex& original ) {
	ModelIndex::operator=( original );
	return *this;
}
