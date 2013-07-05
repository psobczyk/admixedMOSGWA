/********************************************************************************
 *	This file is part of the MOSGWA program code.				*
 *	Copyright ©2012–2013, Bernhard Bodenstorfer.				*
 *										*
 *	This program is free software; you can redistribute it and/or modify	*
 *	it under the terms of the GNU General Public License as published by	*
 *	the Free Software Foundation; either version 3 of the License, or	*
 *	(at your option) any later version.					*
 *										*
 *	This program is distributed in the hope that it will be useful,		*
 *	but WITHOUT ANY WARRANTY; without even the implied warranty of		*
 *	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.			*
 *	See the GNU General Public License for more details.			*
 ********************************************************************************/

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
