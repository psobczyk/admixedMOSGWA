/********************************************************************************
 *	This file is part of the MOSGWA program code.				*
 *	Copyright ©2012–2013, Bernhard Bodenstorfer, Erich Dolejsi		*
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

#ifndef SCORETESTSHORTCUT_HPP
#define SCORETESTSHORTCUT_HPP

#include "Firthimizer.hpp"
#include "Model.hpp"
#include "MData.hpp"
#include "SortVec.hpp"

/** A quick and dirty implementation of score tests based on the guts of the yet unfinished Firthimizer. */
class ScoreTestShortcut : private Firthimizer {

	private:

	/** The common data to which all score tests refer. */
	const MData& mData;

	public:

	/** Prepares an object with given data. */
	ScoreTestShortcut ( const MData& mData );

	/** Calculates the absolute values of the score tests for all SNPs not yet in the model.
	* @param model with regard to which score tests of the other SNPs are performed.
	* @param sortVec will be updated to contain the SNP indices sorted by test result.
	*/
	void scoreTests ( const Model& model, SortVec& sortVec );
	size_t scoreTests ( const Model& model, SortVec& sortVec, size_t, size_t );


};

#endif	/* SCORETESTSHORTCUT_HPP */
