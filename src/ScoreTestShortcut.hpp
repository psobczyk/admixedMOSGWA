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
        int scoreTests ( const Model& model, SortVec& sortVec,int, int );


};

#endif	/* SCORETESTSHORTCUT_HPP */
