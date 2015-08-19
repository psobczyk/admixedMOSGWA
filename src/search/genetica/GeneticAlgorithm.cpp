/********************************************************************************
 *	This file is part of the MOSGWA program code.				*
 *	Copyright ©2015, Bernhard Bodenstorfer.					*
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

#include "GeneticAlgorithm.hpp"
#include "GA.hpp"
#include <sstream>
#include <fstream>
#include <string>
#include <cmath>

using namespace std;
using namespace logging;

namespace genetica {

	GeneticAlgorithm::GeneticAlgorithm () {}

	void GeneticAlgorithm::run () {
		logger->info( "Start genetic algorithm search" );

		const unsigned int modelsNo_ = parameter->modelsNo;
		const unsigned int maxNoProgressIter_ = parameter->maxNoProgressIter;
		const double pCross_ = parameter->pCross;
		const double pMutation_ = parameter->pMutation;
		const unsigned int tournamentSize_ = parameter->tournamentSize;
		const double correlationThreshold_ = parameter->correlationThreshold;
		const int correlationRange_ = parameter->correlationRange;
		const double regionMinCorrelation_ = parameter->regionMinCorrelation;
		const int B_ = parameter->B;

		logger->info( "modelsNo_ = %u", parameter->modelsNo );
		logger->info( "maxNoProgressIter_ = %u", parameter->maxNoProgressIter );
		logger->info( "pCross_ = %f", parameter->pCross );
		logger->info( "pMutation_ = %f", parameter->pMutation );
		logger->info( "tournamentSize_ = %u", parameter->tournamentSize );
		logger->info( "correlationThreshold_ = %f", parameter->correlationThreshold );
		logger->info( "correlationRange_ = %d", parameter->correlationRange );
		logger->info( "regionMinCorrelation_ = %f", parameter->regionMinCorrelation );
		logger->info( "B_ = %d", parameter->B );
		logger->info( "causalModelFilename = \"%s\"", parameter->causalModelFilename.c_str() );

		string old_out_file_name;
		if ( parameter->outNo >= 0 ) {
			old_out_file_name = parameter->out_file_name;
		}

		const time_t time_start = time(NULL);
		timespec ts;
		clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &ts);
		time_t now;
		time(&now);
		srand(now);

		map<snp_index_t, int> mapSNPCausal_ind;
		vector< multiset<long double> > tabCausalPost;
		vector< multiset<long double> > tabCausalPost_b;

		GA ga(
			modelsNo_,
			maxNoProgressIter_,
			pCross_,
			pMutation_,
			tournamentSize_,
			B_,
			parameter->modelsFilename,
			correlationThreshold_,
			correlationRange_,
			regionMinCorrelation_
		);
		ga.run();

		stringstream ss( "GA time: " );
		ga.writePoolToFile(ss);

		ga.initCausalPost( mapSNPCausal_ind );
		tabCausalPost.resize( mapSNPCausal_ind.size() );
		tabCausalPost_b.resize( mapSNPCausal_ind.size() );

		stringstream ssModels;
		ssModels << "";
		ga.computePosteriorProbability(
			ssModels,
			mapSNPCausal_ind,
			tabCausalPost,
			tabCausalPost_b
		);
		writePosterior(
			(parameter->out_file_name + "_post_12.txt").c_str(),
			mapSNPCausal_ind,
			tabCausalPost,
			1
		);

		const time_t time_end = time(NULL);
		ss << (
			time_end > time_start
				? sec2time( time_end - time_start )
				: sec2time( 24 * 3600 + time_end - time_start )
			);
		logger->info( "%s", ss.str().c_str() );
		ga.writePoolToFile(ss);

		ofstream Pmi_YsortFile;
		Pmi_YsortFile.exceptions ( ofstream::eofbit | ofstream::failbit | ofstream::badbit );
		try {
			Pmi_YsortFile.open(
				( parameter->out_file_name + "_PjMi_YsortFile.txt" ).c_str(),
				fstream::out | fstream::trunc
			);
			Pmi_YsortFile << ss.str() << endl << ssModels.str() << endl;
			Pmi_YsortFile.flush();
			Pmi_YsortFile.close();
		} catch ( ofstream::failure e ) {
			logger->error( "Could not write PMi_Y-sorted File: %s", e.what() );
		}
		logger->info(
			"Posterior probalibities of SNPs are in the file: \"%s\"",
			( parameter->out_file_name + "_PjMi_YsortFile.txt" ).c_str()
		);
		logger->info(
			"A report is in the file: \"%s\"",
			( parameter->out_file_name + "_recognized_SNPs.txt" ).c_str()
		);
		logger->info( "%s", ss.str().c_str() );
		if ( parameter->outNo >= 0) {
			parameter->out_file_name = old_out_file_name;
		}
		logger->info( "Finish elaborated greedy search" );
	}

}
