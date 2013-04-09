#include "ResultStore.hpp"

using namespace std;
using namespace lookup;

/** Note: In case of repeated entry of the same index,
* only one of the entries survives.
* It can then happen that the MSC for the index stored as optimal is out of sync with the optimal index.
*/
void ResultStore::setLinearMSC ( const ModelIndex& index, const double msc ) {
	linearMSC[ index ] = msc;
	if ( msc < optimalLinearMSC ) {
		optimalLinearModel = index;
		optimalLinearMSC = msc;
	}
}

double ResultStore::getLinearMSC ( const ModelIndex& index ) const {
	return const_cast<ResultStore*>( this )->linearMSC[ index ];
}

const ModelIndex ResultStore::getOptimalLinearModel () const {
	return optimalLinearModel;
}

/** Note: In case of repeated entry of the same index,
* only one of the entries survives.
* It can then happen that the MSC for the index stored as optimal is out of sync with the optimal index.
*/
void ResultStore::setLogisticMSC ( const ModelIndex& index, const double msc ) {
	logisticMSC[ index ] = msc;
	if ( msc < optimalLogisticMSC ) {
		optimalLogisticModel = index;
		optimalLogisticMSC = msc;
	}
}

double ResultStore::getLogisticMSC ( const ModelIndex& index ) const {
	return const_cast<ResultStore*>( this )->logisticMSC[ index ];
}

const ModelIndex ResultStore::getOptimalLogisticModel () const {
	return optimalLogisticModel;
}
