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

void ScoreTestShortcut::scoreTests ( const Model& model, SortVec& sortVec ) {
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
		const snp_index_t i = *iterator;
		// insertion at the end is fastest
		mData.getXcolumn( i, vec );
		insertColumn( currentCols++, vec );
	}

	// now all columns are there
	assert( cols == currentCols );

	// Set regression coefficients from what was calculated in Model
	// Precondition: Model already up-to-date
	AutoVector coefficients( cols );
	for ( size_t i = 0; i < cols; ++i ) {
		coefficients.set( i, model.getBeta( i ) );
	}
	setCoefficients( coefficients );

	// score tests for snps not yet in model
	const int
		dataSize = mData.getSnpNo(),
		remainingSize = dataSize - index.size();
	vector<size_t> snps( remainingSize );
	vector<double> scores( remainingSize );
	// Erich: Should this loop be parallelized?
	// BB: You can try, scoreTest() should be thread-safe.
	// But be careful with the conditional increment ++j and use of vec.
	// I think it may not be worthwhile for the remaining experiments.
	for ( size_t i = 0, j = 0; i < dataSize; ++i ) {
		if ( ! index.contains( i ) ) {
			snps[j] = i;
			const_cast<MData&>( mData ).getXcolumn( i, vec );
			scores[j] = scoreTest( vec );
			++j;
		}
	}
	assert( remainingSize == snps.size() );
	assert( remainingSize == scores.size() );
	sortVec.fillVec( remainingSize, &snps[0], &scores[0], false );
}

/** TODO: unite this Erich-routine with the standard score test. */
size_t ScoreTestShortcut::scoreTests ( const Model& model, SortVec& sortVec, size_t start, size_t stop ) {
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
		const snp_index_t i = *iterator;
		// insertion at the end is fastest
		mData.getXcolumn( i, vec );
		insertColumn( currentCols++, vec );
	}

	// now all columns are there
	assert( cols == currentCols );

	// Set regression coefficients from what was calculated in Model
	// Precondition: Model already up-to-date
	AutoVector coefficients( cols );
	for ( size_t i = 0; i < cols; ++i ) {
		coefficients.set( i, model.getBeta( i ) );
	}
	setCoefficients( coefficients );

	// score tests for snps not yet in model
	const size_t remainingSize = stop - start + 1;
		// ERICH origninal: dataSize - index.size();
		// da liegt aber eigentlich immer ein SNP drinnen der schon im Model ist,
		// daher gibt es einen SNP zuviel fast immer!
	vector<size_t> snps( remainingSize );
	vector<double> scores( remainingSize );
	// Erich: Should this loop be parallelized?
	// BB: You can try, scoreTest() should be thread-safe.
	// But be careful with the conditional increment ++j.
	// I think it may not be worthwhile for the remaining experiments.
	size_t J = 0;
	for ( size_t i = start, j = 0; i < stop; ++i ) {
	//DEBUG	cerr<<"SCORE:SNP="<<i<<endl;
		if ( ! index.contains( i ) ) {
			snps[j] = i;
			mData.getXcolumn( i, vec );
			scores[j] = scoreTest( vec );
			J = j++;
		}
	}
	assert( remainingSize == snps.size() );
	assert( remainingSize == scores.size() );
	sortVec.fillVec( remainingSize, &snps[0], &scores[0], false );
	return J; //the last 
}
