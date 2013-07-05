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

#include "ModelIndex.hpp"
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <algorithm>

using namespace std;
using namespace lookup;

ModelIndex::ReferenceCountedShared::ReferenceCountedShared ( const size_t size, const snp_index_t *snps )
	: referenceCount( 1 ), size( size ), snps( snps ) {}

ModelIndex::ReferenceCountedShared::~ReferenceCountedShared () {
	assert( 0 == referenceCount );	// Otherwise, why delete?
	free( const_cast<snp_index_t*>( snps ) );
	snps = NULL;
}

void ModelIndex::init ( const set<snp_index_t> &snpSet ) {
	const size_t size = snpSet.size();
	snp_index_t* const snps = static_cast<snp_index_t*>( malloc( size * sizeof( snp_index_t ) ) );
	copy( snpSet.begin(), snpSet.end(), snps );
	snpStruct = new ReferenceCountedShared( size, snps );
}

ModelIndex::ModelIndex () {
	const set<snp_index_t> emptySet;
	init( emptySet );
}

ModelIndex::ModelIndex ( const set<snp_index_t> &snpSet ) {
	init( snpSet );
}

ModelIndex::ModelIndex ( const vector<snp_index_t> &snpVec ) {
	// Use a set to sort and make unique
	set<snp_index_t> snpSet;
	copy( snpVec.begin(), snpVec.end(), inserter( snpSet, snpSet.end() ) );
	init( snpSet );
}

ModelIndex::ModelIndex ( const ModelIndex& original ) : snpStruct( original.snpStruct ) {
	++snpStruct->referenceCount;
}

ModelIndex& ModelIndex::operator= ( const ModelIndex& original ) {
	if ( this->snpStruct != original.snpStruct ) {
		if ( 0 >= --snpStruct->referenceCount ) {
			delete snpStruct;
		}
		++( snpStruct = original.snpStruct )->referenceCount;
	}
	return *this;
}

ModelIndex::ModelIndex ( const ModelIndex& original, const snp_index_t snp ) {
	const_iterator splitPoint = lower_bound( original.begin(), original.end(), snp );

	// whether snp is contained in original
	const bool alreadyContained = original.end() > splitPoint && snp == *splitPoint;
	const size_t size = alreadyContained
		? original.size() - 1
		: original.size() + 1;
	snp_index_t* const snps = static_cast<snp_index_t*>( malloc( size * sizeof( snp_index_t ) ) );

	// copy up to before the location for snp
	snp_index_t* insertionPoint = snps;
	copy( original.begin(), splitPoint, insertionPoint );
	insertionPoint += splitPoint - original.begin();

	// deal with the location of snp
	if ( alreadyContained ) {
		// skip it from original
		++splitPoint;
	} else {
		// insert it to copy
		*insertionPoint++ = snp;
	}

	// copy from after the location of snp to the end
	copy( splitPoint, original.end(), insertionPoint );

	snpStruct = new ReferenceCountedShared( size, snps );
}

ModelIndex::~ModelIndex () {
	if ( 0 >= --snpStruct->referenceCount ) {
		delete snpStruct;
	}
	snpStruct = NULL;
}

/** The internal object with higher reference count is shared.
* In case of a tie, <code>this</code> wins over <code>that</code>.
*/
void ModelIndex::shareInternals ( ModelIndex& that ) {
	if ( snpStruct->referenceCount < that.snpStruct->referenceCount ) {
		*this = that;
	} else {
		that = *this;
	}
}

/** Semantic const-ness of an immutable object means
* that it may be replaced by an equal object.
*/
int ModelIndex::compare ( const ModelIndex& that ) const {
	if ( snpStruct == that.snpStruct ) return 0;	// this was as fast as we like it
	if ( snpStruct->size < that.snpStruct->size ) return -1;
	if ( snpStruct->size > that.snpStruct->size ) return 1;

	// Remains the difficult case of equal size, but different structs
	int result = memcmp( snpStruct->snps, that.snpStruct->snps, snpStruct->size );
	if ( 0 == result ) {	// both should use the same snpStruct, so thus be it
		const_cast<ModelIndex*>( this )->shareInternals( const_cast<ModelIndex&>( that ) );
	}
	return result;
}

bool ModelIndex::operator== ( const ModelIndex& that ) const {
	return 0 == compare( that );
}

bool ModelIndex::operator!= ( const ModelIndex& that ) const {
	return 0 != compare( that );
}

bool ModelIndex::operator< ( const ModelIndex& that ) const {
	return 0 > compare( that );
}

bool ModelIndex::operator<= ( const ModelIndex& that ) const {
	return 0 >= compare( that );
}

bool ModelIndex::operator> ( const ModelIndex& that ) const {
	return 0 < compare( that );
}

bool ModelIndex::operator>= ( const ModelIndex& that ) const {
	return 0 <= compare( that );
}

bool ModelIndex::isSubsetOf ( const ModelIndex& that ) const {
	if ( snpStruct == that.snpStruct ) return true;
	if ( snpStruct->size > that.snpStruct->size ) return false;
	// Remains the case of less or equal size, but different structs
	return includes( that.begin(), that.end(), begin(), end() );
}

bool ModelIndex::isSupersetOf ( const ModelIndex& that ) const {
	if ( snpStruct == that.snpStruct ) return true;
	if ( snpStruct->size < that.snpStruct->size ) return false;
	// Remains the case of greater or equal size, but different structs
	return includes( begin(), end(), that.begin(), that.end() );
}

size_t ModelIndex::size () const { return snpStruct->size; }

ModelIndex::const_iterator ModelIndex::begin () const { return snpStruct->snps; }

ModelIndex::const_iterator ModelIndex::end () const { return snpStruct->snps + snpStruct->size; }

/** Comparison of two <code>snp_index_t</code> locations for {@link bsearch}. */
int snp_pos_cmp ( const void* a, const void* b ) {
	if ( * static_cast<const snp_index_t*>( a ) < * static_cast<const snp_index_t*>( b ) ) return -1;
	if ( * static_cast<const snp_index_t*>( a ) > * static_cast<const snp_index_t*>( b ) ) return 1;
	assert ( * static_cast<const snp_index_t*>( a ) == * static_cast<const snp_index_t*>( b ) );
	return 0;
}

bool ModelIndex::contains ( const snp_index_t snp ) const {
	return NULL != bsearch( &snp, snpStruct->snps, snpStruct->size, sizeof( snp_index_t ), &snp_pos_cmp );
}

ostream& lookup::operator<< ( ostream& s, const ModelIndex& m ) {
	s << '{';
	for ( ModelIndex::const_iterator iterator = m.begin(); iterator != m.end(); ++iterator ) {
		if ( m.begin() != iterator ) s << ',';
		s << *iterator;
	}
	s << '}';
	return s;
}
