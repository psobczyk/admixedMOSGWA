#include "SortVec.hpp"
#include <iostream>

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
