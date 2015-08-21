/********************************************************************************
 *	This file is part of the MOSGWA program code.				*
 *	Copyright ©2012–2015, Bernhard Bodenstorfer.				*
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

#include "ScoreTestShortcut.hpp"
#include <assert.h>
#include <vector>

using namespace linalg;
using namespace lookup;
using namespace std;

/** Helps Firthimizer construction by temporarily allocating an <code>AutoVector</code>. */
static AutoVector getY( const MData& mData ) {
	AutoVector vector( mData.getIdvNo() );
	const_cast<MData&>( mData ).getY( vector );
	return vector;
}

ScoreTestShortcut::ScoreTestShortcut (
	const MData& mData
) : Firthimizer( getY( mData ) ), mData( mData ) {}

void ScoreTestShortcut::scoreTests (
	const Model& model,
	const size_t start,
	const size_t length,
	SortVec& sortVec
) {
	const ModelIndex index = model.getIndex();
	const size_t
		rows = mData.getIdvNo(),
		cols = 1 + index.size();	// 1 for leftmost column of ones
	AutoVector vec( rows );

	// track current number of columns
	size_t currentCols = countDimensions();

	// Empty xMat; from last to first is fastest
	while ( currentCols > 0 ) {
		removeColumn( --currentCols );
	}

	// X leftmost column is 1
	vec.fill( 1.0 );
	insertColumn( currentCols++, vec );

	// X columns 1 to modelsize
	for (
		ModelIndex::const_iterator iterator = index.begin();
		iterator < index.end();
		++iterator
	) {
		const size_t i = *iterator;
		// insertion at the end is fastest
		mData.getXcolumn( i, vec );
		insertColumn( currentCols++, vec );
	}

	// now all columns are there
	assert( cols == currentCols );

	// Set regression coefficients from what was calculated in Model
	// Precondition: Model already up-to-date
	// TODO<BB>: Add method to set whole vector
	AutoVector coefficients( cols );
	for ( size_t i = 0; i < cols; ++i ) {
		coefficients.set( i, model.getBeta( i ) );
	}
	setCoefficients( coefficients );

	// score tests for snps not yet in model
	const size_t stop = start + length;
	vector<size_t> snps( length );
	vector<double> scores( length );
	// Erich: Should this loop be parallelized?
	// BB: You can try, scoreTest() should be thread-safe.
	// But be careful with the conditional increment ++j.
	// I think it may not be worthwhile for the remaining experiments.
	size_t j = 0;
	for ( size_t i = start; i < stop; ++i ) {
		if ( ! index.contains( i ) ) {
			snps[j] = i;
			mData.getXcolumn( i, vec );
			scores[j] = scoreTest( vec );
			++j;
		}
	}
	sortVec.fillVec( j, &snps[0], &scores[0], false );
}
