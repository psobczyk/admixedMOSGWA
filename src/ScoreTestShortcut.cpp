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

ScoreTestShortcut::ScoreTestShortcut (
	const MData& mData
) : Firthimizer( const_cast<MData&>( mData ).getY() ), mData( mData ) {}

void ScoreTestShortcut::scoreTests ( const Model& model, SortVec& sortVec ) {
	const ModelIndex index = model.getIndex();
	const size_t
		rows = mData.getIdvNo(),
		cols = 1 + index.size();	// 1 for leftmost column of ones

	// track current number of columns
	size_t currentCols = countDimensions();

	// Empty xMat; from last to first is fastest
	while ( currentCols > 0 ) {
		removeColumn( --currentCols );
	}

	// X leftmost column is 1
	AutoVector ones( rows );
	ones.fill( 1.0 );
	insertColumn( currentCols++, ones );

	// X columns 1 to modelsize
	for (
		ModelIndex::const_iterator iterator = index.begin();
		iterator < index.end();
		++iterator
	) {
		const snp_index_t i = *iterator;
		// insertion at the end is fastest
		insertColumn( currentCols++, const_cast<MData&>( mData ).getXcolumn( i ) );
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
	vector<int> snps( remainingSize );
	vector<double> scores( remainingSize );
	// Erich: Should this loop be parallelized? - BB: You can try, scoreTest() should be thread-safe. But be careful with the conditional increment ++j. I think it may not be worhwhile at this point in time for the remaining experiments.
	for ( int i = 0, j = 0; i < dataSize; ++i ) {
		if ( ! index.contains( i ) ) {
			snps[j] = i;
			scores[j] = scoreTest( const_cast<MData&>( mData ).getXcolumn( i ) );
			++j;
		}
	}
	assert( remainingSize == snps.size() );
	assert( remainingSize == scores.size() );
	sortVec.fillVec( remainingSize, &snps[0], &scores[0], false );
}
int ScoreTestShortcut::scoreTests ( const Model& model, SortVec& sortVec,int  start, int stop ) {
	const ModelIndex index = model.getIndex();
	const size_t
		rows = mData.getIdvNo(),
		cols = 1 + index.size();	// 1 for leftmost column of ones

	// track current number of columns
	size_t currentCols = countDimensions();

	// Empty xMat; from last to first is fastest
	while ( currentCols > 0 ) {
		removeColumn( --currentCols );
	}

	// X leftmost column is 1
	AutoVector ones( rows );
	ones.fill( 1.0 );
	insertColumn( currentCols++, ones );

	// X columns 1 to modelsize
	for (
		ModelIndex::const_iterator iterator = index.begin();
		iterator < index.end();
		++iterator
	) {
		const snp_index_t i = *iterator;
		// insertion at the end is fastest
		insertColumn( currentCols++, const_cast<MData&>( mData ).getXcolumn( i ) );
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
		dataSize =stop,//ERICH original: mData.getSnpNo(),
		remainingSize = stop-start+1; //ERICH origninal dataSize - index.size();
	                                     //da liegt aber eigentlich immer ein SNP drinnen der schon im Model ist, daher gbit es einen SNP zuviel fast immer!
	vector<int> snps( remainingSize );
	vector<double> scores( remainingSize );
	// Erich: Should this loop be parallelized? - BB: You can try, scoreTest() should be thread-safe. But be careful with the conditional increment ++j. I think it may not be worhwhile at this point in time for the remaining experiments.
	int J=0;
	for ( int i = start, j = 0  ; i < stop+1; ++i ) {
	//DEBUG	cerr<<"SCORE:SNP="<<i<<endl;
		if ( ! index.contains( i ) ) {
			snps[j] = i;
			scores[j] = scoreTest( const_cast<MData&>( mData ).getXcolumn( i ) );
			J=j;++j;
		}
	}
	assert( remainingSize == snps.size() );
	assert( remainingSize == scores.size() );
	sortVec.fillVec( remainingSize, &snps[0], &scores[0], false );
	return J; //the last 
}
