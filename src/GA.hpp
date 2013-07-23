/********************************************************************************
 *	This file is part of the MOSGWA program code.				*
 *	Copyright ©2011–2013, Artur Gola, Bernhard Bodenstorfer.		*
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

#ifndef GA_HPP
#define GA_HPP

#include <iostream>
#include <vector>
#include <iostream>
#include <set>
#include <string>
#include <sstream>
#include "Log.hpp"
#include "Pool.hpp"
#include "MData.hpp"
#include "Model.hpp"

using namespace std;

ostream &operator<< ( ostream &out, vector<snp_index_t> &v );

class GA {

private:

  /** set of exclusived SNPs */
  set<snp_index_t> exclusivedSNP;
  
  
	/** "table" of correration values for snps in all models */
	static vector< vector<snp_index_t> > correlations;

	/** data for models */
	MData data;

	/** models */
	Model **models;

	/** number of models */
	unsigned int modelsNo;

	/** pool of models.
	* PoolItem has got information about snps - the position-no of the SNPs in the Model and msc of the model
	*/
	set<PoolItem> pool;

	/** number of iteration with no progress */
	unsigned int maxNoProgressIter;

	/** probabilitiy of crossover */
	double pCross;

	/** probabilitiy of mutation */
	double pMutation;

	/** tournament size */
	unsigned int tournamentSize;

	/** correlation threshold - used in recombination and local improvement functions */
	double correlationThreshold;

	/** snps for correlation are from range [snp - correlationRange, snp + correlationRange] */
	int correlationRange;

	/** mutates model */
	void mutation ( Model *aModel, double pMutation, double threshold );

	/** makes one child from two parenst */
	Model* recombination ( const Model & s1, const Model & s2 );

	/** finds and returns index to the worst model in the population */
	unsigned int findTheWorst () const;

	/** tournament selection */
	Model* tournamentSelection ( unsigned int tournamentSize ) const;

	/** find the best, the worst and average models */
	void statistics ( double &theBest, double &average, double &theWorst );

	/** correlation for snp on snpPosition position - relative position of the snp at modelsnps_ vector */
	vector<snp_index_t> stronglyCorrelatedSnpsAt ( Model *model, const int& snpPosition, const double& threshold, int correlationRange );

	/** correlation for snp */
	vector<snp_index_t> stronglyCorrelatedSnps ( Model *model, const int& snp, const double& threshold, int correlationRange );

	/** local inprovement of a given model */
	void localImprovement ( Model *model, double threshold, int correlationRange );
	
	/** computes and writes to file (*_pProb.txt) posterior probalibities of models */
	void computePosteriorProbabilitys();

public:

	GA (
		unsigned int modelsNo_,
		unsigned int maxNoProgressIter_,
		double pCross_,
		double pMutation_,
		unsigned int tournamentSize_,
		double correlationThreshold_ = 0.7,
		int correlationRange = 500
	);

	/** writes pool on the screen. WARNING Be carefull, pool may be very large */
	void writePool () const;

	/** writes pool to log file */
	void writePoolToLog () const;

	/** runs genethic algorithm */
	void run ();

	~GA ();

  void selectModel ( Model& selectedModel );
  
  void calculateIndividualTests();
  
	/** TESTING for testing only */
	void tests ();

	/** TESTING for testing only */
	void testRecombination ();

	/** TESTING for testing only */
	void testTournament ();

	/** TESTING for testing only */
	void testMutation ();

	/** TESTING for testing only */
	void testsCorrelation ();
  
  void printExclusivedSNP();
};

#endif
