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
#include <time.h> 

#include "Log.hpp"
#include "Pool.hpp"
#include "MData.hpp"
#include "Model.hpp"
#include "search/Search.hpp"

ostream &operator<< ( ostream &out, vector<snp_index_t> &v );



struct TPOWER_FDR 
{ 
  long double POWER;
  long double FDR;
  unsigned int FDcount;
  snp_index_t badSNP;       // the best SNP of posteriori <= 0.5
  long double posteriorBad;
  long double posteriorBadCluster;
};

/** calculates the difference between two times */
void time_diff(int &hour, int & min, int &sec, int &nano, timespec &ts, timespec &te);

/** outputs the map on the screen */
template <typename T1, typename T2>
ostream &operator<< (ostream &out, map<T1, T2> &v)
{
  for (typename map<T1, T2>::iterator it = v.begin(); it != v.end(); it++)
    out << "[" << it->first << " -> " << it->second << "]" << endl;
  return out;
}

/** outputs the vector on the screen */
template <typename T>
ostream &operator<< (ostream &out, vector<T> &v)
{
  out << "[";
  typename vector<T>::iterator it = v.begin();
  if (v.size() > 0)
  {
    out << *it;
    it++;
  }
  for (; it != v.end(); it++)
    out << ", " << *it;
  return out << "]";
}

/** outputs the set on the screen */
template <typename T>
ostream &operator<< (ostream &out, set<T> &v)
{
  out << "{";
  typename set<T>::iterator it = v.begin();
  if (v.size() > 0)
  {
    out << *it;
    it++;
  }
  for (; it != v.end(); it++)
    out << ", " << *it;
  return out << "}";
}


class THeri_PMiY
{
public:
  long double h2M;         // heriability
  long double PMi_Y;  
  THeri_PMiY(long double h2M = 0.0L, long double PMi_Y = 0.0L): h2M(h2M), PMi_Y(PMi_Y) {}
  bool operator < (const THeri_PMiY & t) const {return h2M < t.h2M;}
};


class GA //: public search::Search 
{
private:
  map<int, int> cs;
  const double correlationTh;
  
  TBitset exclusivedSNP;   // set of exclusived SNPs - for initial population
  TBitset goodSNPs;        // 
  int goodSNPsNo;          // number of SNPs wchich p-val are > 0.1
  
  ofstream poolFilePart;
  
  /** "table" of correration values for snps in all models */
  static vector< vector<snp_index_t> > correlations;

  /** Data for models */
  MData data;

  /** Models */
  Model **models;

  /** The number of models */
  unsigned int modelsNo;

  /** The pool of models */
  set<PoolItem> pool;

  /** The number of iteration with no progress */
  unsigned int maxNoProgressIter;

  /** The probabilitiy of crossover */
  double pCross;

  /** The probabilitiy of mutation */
  double pMutation;

  /** The tournament size */
  unsigned int tournamentSize;

  /** The correlation threshold - it's used in recombination and local improvement function */
  double correlationThreshold;

  /** SNPs for the correlation are from range [snp - correlationRange, snp + correlationRange] */
  int correlationRange;

  /** Mutates the model */
  void mutation ( Model *aModel, double pMutation, double threshold );

  /** Makes one child from two parents */
  Model* recombination ( const Model & s1, const Model & s2 );

  /** Finds and returns index to the worst model in the population */
  unsigned int findTheWorst () const;

  /** The tournament selection */
  Model* tournamentSelection ( unsigned int tournamentSize ) const;

  /** The stop criterion.*/
  int B;
  
  /** Finds the best, the worst and average models */
  void statistics ( double &theBest, double &average, double &theWorst );

  /** correlation for snp on snpPosition position - relative position of the snp at modelsnps_ vector */
  vector<snp_index_t> stronglyCorrelatedSnpsAt ( Model *model, const int& snpPosition, const double& threshold, int correlationRange );

  /** correlation for snp */
  vector<snp_index_t> stronglyCorrelatedSnps ( Model *model, const int& snp, const double& threshold, int correlationRange );

  /** local inprovement of a given model */
  void old_localImprovement ( Model *model, double threshold, int correlationRange );
  
  /** new version - works not too good */
  void localImprovement ( Model &model, double threshold, int correlationRange );

  set<snp_index_t> testSet;  
  
  long double RSSo;
  
  /** Pointer to real model */
  Model *realModel;  
  
  map<snp_index_t, int> mapSNPid_label;    // Map of SNP_id -> cluster label - to calculate POWER, FDR, ...
  
  map<snp_index_t, int> mapSNPid_labelLO;    // Map of SNP_id -> cluster label - for the local improvement
  
  map<int, long double> mapLabel_PmiY;     // Maps of label -> posterior of whole cluster (the sum) - to calculate POWER, FDR, ...
 
  vector< set<snp_index_t> >vectClust;      // An element of this vector is a claster - the set of snp_id
  
  map<int, int> mapLabel_ind;               // Maps of label -> index. The index of the vector of set
  
  vector<int> clusterLabel;                 // clusterLabel[i] is a label of i-th cluster
                                            // clusterLabel[i] - the label of i-th cluster, vectClust[i] - the set of snp of the i-th cluster
  static map<int, int> mapLabel_count;       // zlicza etykiety wykrytych klastrów dla wielu symulacji - dane rzeczywiste
                                            
  
  /** Returns number of models which have got better msc value than the parameter mscVal */
  int isInNBestModels(double mscValue);
  
  static map<snp_index_t, int> recognizedSNPs;
  
  /**
   * @brief Prepares vectClust, mapLabel_ind, clusterLabel
   */
  void makeVectCluster(map<snp_index_t, int>& mapSNPid_label);
  
  void readClusterLabel(map<snp_index_t, int>& mapSNPid_label, const string& fileName, map<string, int>& SNP_Names_Id, map<int, long double>* mapLabel_PmiY = 0);
  
  void toPool(const Model *model, char c);
  
public:

  GA (
    unsigned int modelsNo_,
    unsigned int maxNoProgressIter_,
    double pCross_,
    double pMutation_,
    unsigned int tournamentSize_,
    int B_,
    string fileName,
    double correlationThreshold_,
    int correlationRange,
    bool statisticsOnly = false
  );

  /** writes pool on the screen. WARNING Be carefull, pool may be very large */
  void writePool () const;

  /** writes pool to log file */
  void writePoolToFile (stringstream &ssModels) const;

  /** runs genethic algorithm */
  void run ();

  ~GA ();

  void selectModel ( Model& selectedModel );
  
  //void calculateIndividualTests();
  
  /** computes and writes to file (*_pProb.txt) posterior probalibities of models */

  void computePosteriorProbability(stringstream &ssModels, map<snp_index_t, int> &mapSNPCausal_ind, vector< multiset<long double> > &tabCausalPost,
                                    vector< multiset<long double> > &tabCausalPost_b, long double minPosterior);
  
  void poolReader(string fileName, stringstream& sGATime);//, map<snp_index_t, int> *recognizedSNPs = 0);
  
  /** TESTING for testing only */
//  void tests ();

  

  /** TESTING for testing only */
  //void testsCorrelation ();
  
  //void printExclusivedSNP();
  
  
  void piMassExtract(const string &fileName, string &outFileName);//, map<snp_index_t), int> &recognizedSNPs) ;
//  void piMassExtract_old_before_2013_05_01(int from, int to);
  
  int readInitialPop(string fileName, set<PoolItem> &population, int popSize);
  
  void calculateClusterPosterior(string fileName, long double minPosterior = 0.001);
  
  void calculateClusterPosterior(vector<snp_index_t> &clusterSNP, vector<long double> &clusterSNPposterior);
  
  /** the old first method, no cluster, only correlation between SNPs */
  void calculatePOWER_FDR(set<snp_index_t> &mySnps, vector<snp_index_t> &realSNPs, long double &POWER, long double &FDR, unsigned int & FDcount, set<snp_index_t> &trueSNPs);
  
  void calculatePOWER_FDR_clustGA(set<snp_index_t> &mySnps, vector<snp_index_t> &realSNPs, long double &POWER, long double &FDR, unsigned int & FDcount, set<snp_index_t> &trueSNPs);
  
  void calculatePOWER_FDR_clust(set<snp_index_t> &mySnps, vector<snp_index_t> &causalSNPs, TPOWER_FDR &powerFDR, map<snp_index_t, long double> &mapSNPid_Pmi_Y, 
                                snp_index_t badSNP, map<snp_index_t, int> &recognizedSNPs, map<snp_index_t, int> &mapSNPCausal_ind,
                                vector< multiset<long double> > &tabCausalPost);
  
       
  
  void calculatePOWER_FDR_clust(set<snp_index_t> &mySnps, vector<snp_index_t> &causalSNPs, TPOWER_FDR &powerFDR, map<snp_index_t, long double> &mapSNPid_Pmi_Y,
                                map<snp_index_t, int> &mapSNPCausal_ind, vector< multiset<long double> > &tabCausalPost_b);
  
  
  void calculatePOWER_FDR_clust_sum(vector<snp_index_t> &causalSNPs, TPOWER_FDR &powerFDR, map<snp_index_t, long double> &mapSNPid_Pmi_Y, 
                                map<snp_index_t, int> &recognizedSNPs, map<snp_index_t, int> &mapSNPCausal_ind, vector< multiset<long double> > &tabCausalPost);  

  void calculatePOWER_FDR_clust_max(vector<snp_index_t> &causalSNPs, TPOWER_FDR &powerFDR, map<snp_index_t, long double> &mapSNPid_Pmi_Y, 
                                map<snp_index_t, int> &recognizedSNPs, map<snp_index_t, int> &mapSNPCausal_ind, vector< multiset<long double> > &tabCausalPost);
  
  
  void checkCorrelation(set<snp_index_t> &mySnps, vector<snp_index_t> &realSNPs);
  
  /** Read the real model from file */
  void modelReader(string fileName, vector<snp_index_t> &v);

  void setRecognisedSNPs(); //  
  
  static map<snp_index_t, int> getRecognisedSNPs() {return recognizedSNPs;}
  
  double aBestCorrelatedSNP(snp_index_t aSNP, set<snp_index_t> & snps, snp_index_t &bestCorrelatedSNP) const;
  
  // tworzy klastry do obliczenia POWER dla piMassa
  void makeClusers(vector<set<snp_index_t> > &tab);
  
  void piMassCluserPOWER(const string &fileName, const string &outFileName, map<snp_index_t, int> &mapSNPCausal_ind, vector< multiset<long double> > &tabCausalPost);
  
//  void calculatePOWER_FDR_clust(const set<snp_index_t> &mySnps, const vector<snp_index_t> &causalSNPs, TPOWER_FDR &powerFDR, const map<snp_index_t, long double> &Pmi_Y);
  
  void coumputeHeritability(stringstream &sp_sort, const vector<long double> &diff, vector<THeri_PMiY> &tab);
  
  void initCausalPost( map<snp_index_t, int> &mapSNPCausal_ind );
  
  void makeCausalClasters(string fileName, vector<snp_index_t> causalSNPs);  
  
  void saveLabelCount(const string &fileName);
  
  void setCausalCount(map<snp_index_t, int> &map_SNP2count);
  
};

void writePosterior(string fileName, map<snp_index_t, int> &mapSNPCausal_ind, vector< multiset<long double> > &tabCausalPost, int size);

#endif

