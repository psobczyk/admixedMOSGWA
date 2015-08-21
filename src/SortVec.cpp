/********************************************************************************
 *	This file is part of the MOSGWA program code.				*
 *	Copyright ©2011–2015, Erich Dolejsi, Bernhard Bodenstorfer.		*
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

#include "SortVec.hpp"
#include <algorithm>
#include <iostream>

using namespace std;

/** Determines the order between two SortItems:
* i < j iff the value of i is smaller than the value of j.
*/
bool order_function ( const SortVec::SortItem* i, const SortVec::SortItem* j ) {
	return ( i->value < j->value );
}
bool order_function2 ( const SortVec::SortItem* i, const SortVec::SortItem* j ){
	return ( i->value > j->value );
};
SortVec::SortItem::SortItem ( const size_t id, const double value ) : id( id ), value( value ) {}

SortVec::SortVec () {}

SortVec::SortVec ( const size_t n ) {
	list_.reserve( n );
}

SortVec::SortVec ( const size_t n, const size_t ids[], const double values[], bool bigger ) {
	fillVec( n, ids, values, bigger );
}

void SortVec::fillVec ( const size_t n, const size_t ids[], const double values[], bool bigger ) {
	clear();
	list_.reserve( n );
	for ( size_t i = 0; i < n; ++i ) {
		SortItem * entry = new SortItem( ids[i], values[i] );
		list_.push_back(entry);
	}
	if (bigger)
	sort( list_.begin(), list_.end(), order_function );
	else
	sort( list_.begin(), list_.end(), order_function2 );
}

size_t SortVec::getId ( const size_t k ) const {
	return list_.at(k)->id;
}

size_t SortVec::size () const {
	return list_.size();
}

double SortVec::getValue ( const size_t k ) const {
	return list_.at(k)->value;
}

void SortVec::clear () {
	for(
		vector<SortItem*>::iterator it = list_.begin();
		it < list_.end();
		++it
	) {
		delete *it;
	}
	list_.clear();
}

SortVec::~SortVec () {
	clear();
}
