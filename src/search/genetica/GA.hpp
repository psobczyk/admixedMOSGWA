/********************************************************************************
 *	This file is part of the MOSGWA program code.				*
 *	Copyright ©2011–2015, Artur Gola, Bernhard Bodenstorfer.		*
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

#include "../../Pool.hpp"
#include "../../MData.hpp"
#include "../../Model.hpp"
#include "../../search/Search.hpp"

ostream &operator<< ( std::ostream &out, std::vector<snp_index_t> &v );

namespace genetica {

struct TSNP_Info 
{
  snp_index_t SNP;
  unsigned int Chr;
  int pos;
  bool operator < (const TSNP_Info & s) const
  {
    return Chr < s.Chr ? true : (Chr > s.Chr ? false : pos < s.pos);
  }
  friend std::ostream& operator << (std::ostream &out, const TSNP_Info & t)
  {
    return out << "SNP: " << t.SNP << ", Chr: " << t.Chr << ", pos: " << t.pos;
  }
};

struct TRegionSet_info
{
  set<TSNP_Info>* s;
  
  bool operator < (const TRegionSet_info &t) const;
};


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
std::string stemp_diff(const tm* t_start, const tm* t_end);
std::string sec2time(const time_t &t);


/** outputs the map on the screen */
template <typename T1, typename T2>
std::ostream &operator<< (std::ostream &out, std::map<T1, T2> &v)
{
  for (typename std::map<T1, T2>::iterator it = v.begin(); it != v.end(); it++)
    out << "[" << it->first << " -> " << it->second << "]" << endl;
  return out;
}

/** outputs the vector on the screen */
template <typename T>
std::ostream &operator<< (std::ostream &out, std::vector<T> &v)
{
  out << "[";
  typename std::vector<T>::iterator it = v.begin();
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
std::ostream &operator<< (std::ostream &out, std::set<T> &v)
{
  out << "{";
  typename std::set<T>::iterator it = v.begin();
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
  std::map<int, int> cs;
  const double correlationTh;

  TBitset exclusivedSNP;   // a set of exclusived SNPs - for the initial population
  TBitset goodSNPs;        // a set of SNPs which p-val are > 0.1
  int goodSNPsNo;          // a number of SNPs wchich p-val are > 0.1
  
  std::ofstream poolFilePart;
  
  /** "table" of correration values for snps in all models */
  static std::vector< std::vector<snp_index_t> > correlations;

  /** Data for models */
  MData data;

  /** Models */
  Model **models;

  /** The number of models */
  unsigned int modelsNo;

  /** The pool of models */
  std::set<PoolItem> pool;

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
  
  double regionMinCorrelation;  // threshold for a region method

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
  
  /** SNPs calculatee by MOSGWA. poolReader() sets up this variable */
  std::vector<snp_index_t> v_mosgwa;

  /**  allows to save a pool when it will have 30K (100K) models */
  bool p_30K,
       p_100K;
  time_t time_start; // the start time for a pool of size 30K or 100K

  bool isCluster;
       
  
  /** Finds the best, the worst and average models */
  void statistics ( double &theBest, double &average, double &theWorst );

  /** correlation for snp on snpPosition position - relative position of the snp at modelsnps_ vector */
  std::vector<snp_index_t> stronglyCorrelatedSnpsAt ( Model *model, const int& snpPosition, const double& threshold, int correlationRange );

  /** correlation for snp */
  std::vector<snp_index_t> stronglyCorrelatedSnps ( Model *model, const int& snp, const double& threshold, int correlationRange );

  /** local inprovement of a given model */
  void old_localImprovement ( Model *model, double threshold, int correlationRange );
  
  /** new version - works not too good */
  void localImprovement ( Model &model, double threshold, int correlationRange );

  std::set<snp_index_t> testSet;  
  
  long double RSSo;
  
  /** Pointer to real model */
  Model *realModel;  
  
  std::map<snp_index_t, int> mapSNPid_label;     // Map of SNP_id -> cluster label - to calculate POWER, FDR, ...
  
  std::map<snp_index_t, int> mapSNPid_labelLO;   // Map of SNP_id -> cluster label - for the local improvement
  
  std::map<int, long double> mapLabel_PmiY;      // Maps of label -> posterior of whole cluster (the sum) - to calculate POWER, FDR, ...
 
  std::vector< set<snp_index_t> >vectClust;      // An element of this vector is a claster - the set of snp_id
  
  std::map<int, int> mapLabel_ind;               // Maps of label -> index. The index of the vector of set
  
  std::vector<int> clusterLabel;                 // clusterLabel[i] is a label of i-th cluster
                                            // clusterLabel[i] - the label of i-th cluster, vectClust[i] - the set of snp of the i-th cluster
                                            
  static std::map<int, int> mapLabel_count;      // It is needed to count the labels of clusters for many runs of GA.
                                            
  
  /** Returns number of models which have got better msc value than the parameter mscVal */
  int isInNBestModels(double mscValue);
  
  static std::map<snp_index_t, int> recognizedSNPs_Region;
  static std::map<snp_index_t, int> recognizedSNPs_bestGA;
  static std::map<snp_index_t, int> recognizedSNPs_mosgwa;
  static std::map<snp_index_t, int> recognizedSNPs_posterioriModel;
  static std::map<snp_index_t, int> recognizedSNPs_clusterMax;  
  static std::map<snp_index_t, int> recognizedSNPs_clusterSum;  
  static std::map<snp_index_t, int> recognizedSNPs_piMass;  
  
  /**
   * @brief Prepares vectClust, mapLabel_ind, clusterLabel
   */
  void makeVectCluster(std::map<snp_index_t, int>& mapSNPid_label);
  
  void readClusterLabel(std::map<snp_index_t, int>& mapSNPid_label, const std::string& fileName, std::map<std::string, int>& SNP_Names_Id, std::map<int, long double>* mapLabel_PmiY = 0);
  
  void toPool(const Model *model, char c);
  
public:

  GA (
    unsigned int modelsNo_,
    unsigned int maxNoProgressIter_,
    double pCross_,
    double pMutation_,
    unsigned int tournamentSize_,
    int B_,
    std::string fileName,
    double correlationThreshold_,
    int correlationRange,
    double regionMinCorrelation,
    bool statisticsOnly = false
  );

  /** writes pool on the screen. WARNING Be carefull, pool may be very large */
  void writePool () const;

  /** writes pool to log file */
  void writePoolToFile (std::stringstream &ssModels, std::string postfix = "") const;

  /** runs genethic algorithm */
  void run ();

  ~GA ();

  void selectModel ( Model& selectedModel );
  
  /** computes and writes to file (*_pProb.txt) posterior probalibities of models */
  void computePosteriorProbability(std::stringstream &ssModels, std::map<snp_index_t, int> &mapSNPCausal_ind, std::vector< std::multiset<long double> > &tabCausalPost,
                                    std::vector< std::multiset<long double> > &tabCausalPost_b);//, long double minPosterior);
  
  void poolReader(std::string fileName, std::stringstream& sGATime, int real_modelsNo);//, std::map<snp_index_t, int> *recognizedSNPs = 0);
  
  /** TESTING for testing only */
//  void tests ();

  

  /** TESTING for testing only */
  //void testsCorrelation ();
  
  //void printExclusivedSNP();
  
  
  void piMassExtract(const std::string &fileName, std::string &outFileName);//, std::map<snp_index_t), int> &recognizedSNPs) ;
//  void piMassExtract_old_before_2013_05_01(int from, int to);
  
  int readInitialPop(std::string fileName, std::set<PoolItem> &population, int popSize);
  
  void calculateClusterPosterior(std::string fileName, long double minPosterior = 0.001);
  
  void calculateClusterPosterior(std::vector<snp_index_t> &clusterSNP, std::vector<long double> &clusterSNPposterior);
  
  
//  void calculatePOWER_FDR_clustGA(std::set<snp_index_t> &mySnps, std::vector<snp_index_t> &causalSNPs, long double &POWER, long double &FDR, unsigned int & FDcount, std::set<snp_index_t> &TP_SNPs);
  void calculatePOWER_FDR_clustGA(std::set<snp_index_t> &mySnps, std::vector<snp_index_t> &causalSNPs, TPOWER_FDR &powerFDR, std::set<snp_index_t> &TP_SNPs, std::map<snp_index_t, int> &recognizedSNPs);


  // Na potrzeby 2-go artykułu
  void calculatePOWER_FDR_clustGA_2ndArticle(std::set<snp_index_t> &mySnps, std::vector<snp_index_t> &realSNPs, long double &POWER, long double &FDR, 
                                    unsigned int & FDcount, std::set<snp_index_t> &trueSNPs);
  
  
  /** the old first method, no cluster, only correlation between SNPs */
  void calculatePOWER_FDR(std::set<snp_index_t> &mySnps, std::vector<snp_index_t> &realSNPs, long double &POWER, long double &FDR, unsigned int & FDcount, std::set<snp_index_t> &trueSNPs);
  
/* 
  void calculatePOWER_FDR_clust(std::set<snp_index_t> &mySnps, std::vector<snp_index_t> &causalSNPs, TPOWER_FDR &powerFDR, std::map<snp_index_t, long double> &mapSNPid_Pmi_Y, 
                                snp_index_t badSNP, std::map<snp_index_t, int> &recognizedSNPs, std::map<snp_index_t, int> &mapSNPCausal_ind,
                                std::vector< std::multiset<long double> > &tabCausalPost);
*/  
       
  
//  void calculatePOWER_FDR_clust(std::set<snp_index_t> &mySnps, std::vector<snp_index_t> &causalSNPs, TPOWER_FDR &powerFDR, std::map<snp_index_t, long double> &mapSNPid_Pmi_Y,
//                                std::map<snp_index_t, int> &mapSNPCausal_ind, std::vector< std::multiset<long double> > &tabCausalPost_b);
  
  
  void calculatePOWER_FDR_clust_sum(std::vector<snp_index_t> &causalSNPs, TPOWER_FDR &powerFDR, std::map<snp_index_t, long double> &mapSNPid_Pmi_Y, 
                                std::map<snp_index_t, int> &recognizedSNPs, std::map<snp_index_t, int> &mapSNPCausal_ind, std::vector< std::multiset<long double> > &tabCausalPost);  

  void calculatePOWER_FDR_clust_max(std::vector<snp_index_t> &causalSNPs, TPOWER_FDR &powerFDR, std::map<snp_index_t, long double> &mapSNPid_Pmi_Y, 
                                std::map<snp_index_t, int> &recognizedSNPs, std::map<snp_index_t, int> &mapSNPCausal_ind, std::vector< std::multiset<long double> > &tabCausalPost);
  
/** test */  
//void calculate_clusters(std::vector<snp_index_t> &causalSNPs, TPOWER_FDR &powerFDR, std::map<snp_index_t, long double> &mapSNPid_Pmi_Y, 
//                                std::map<snp_index_t, int> &recognizedSNPs, std::map<snp_index_t, int> &mapSNPCausal_ind, std::vector< std::multiset<long double> > &tabCausalPost, bool clusterSum);
void calculate_clusters(std::map<snp_index_t, long double> &mapSNPid_Pmi_Y, std::vector<snp_index_t> &modelSNPs, bool clusterSum, std::string method_Name = "");

  
  void checkCorrelation(std::set<snp_index_t> &mySnps, std::vector<snp_index_t> &realSNPs);
  
  /** Read the real model from file */
  void modelReader(std::string fileName, std::vector<snp_index_t> &v);

  void setRecognisedSNPs(); //  
  
  static std::map<snp_index_t, int> getRecognizedSNPs_Region() {return recognizedSNPs_Region;}
  static std::map<snp_index_t, int> getRecognizedSNPs_bestGA() {return recognizedSNPs_bestGA;}  
  static std::map<snp_index_t, int> getRecognisedSNPs_mosgwa() {return recognizedSNPs_mosgwa;}
  static std::map<snp_index_t, int> getRecognisedSNPs_posterioriModel() {return recognizedSNPs_posterioriModel;}
  static std::map<snp_index_t, int> getRecognisedSNPs_clusterMax() {return recognizedSNPs_clusterMax;}
  static std::map<snp_index_t, int> getRecognisedSNPs_clusterSum() {return recognizedSNPs_clusterSum;}
  static std::map<snp_index_t, int> getRecognisedSNPs_piMass() {return recognizedSNPs_piMass;}
  
//  double theBestCorrelatedSNP(snp_index_t aSNP, std::set<snp_index_t> & snps, snp_index_t &bestCorrelatedSNP) const;
  
  // tworzy klastry do obliczenia POWER dla piMassa
  void makeClusers(std::vector<std::set<snp_index_t> > &tab);
  
  void piMassCluserPOWER(const std::string &fileName, const std::string &outFileName, std::map<snp_index_t, int> &mapSNPCausal_ind, std::vector< std::multiset<long double> > &tabCausalPost);
  
//  void calculatePOWER_FDR_clust(const std::set<snp_index_t> &mySnps, const std::vector<snp_index_t> &causalSNPs, TPOWER_FDR &powerFDR, const std::map<snp_index_t, long double> &Pmi_Y);
  
  void coumputeHeritability(std::stringstream &sp_sort, const std::vector<long double> &diff, std::vector<THeri_PMiY> &tab);
  
  void initCausalPost( std::map<snp_index_t, int> &mapSNPCausal_ind );
  
  void makeCausalClasters(std::string fileName, std::vector<snp_index_t> causalSNPs);  
  
  void saveLabelCount(const std::string &fileName);
  
  void setCausalCount(std::map<snp_index_t, int> &map_SNP2count);

//  void stronglyCorrelatedSnpsCluster(const double& threshold = 0.5);  
  void stronglyCorrelatedSnpsCluster(const std::vector<snp_index_t> & tabSNPs, const double& threshold );
  
  void saveRecognizedSNPinfo(const std::vector<snp_index_t> &v, const std::string &method_Name, std::map<snp_index_t, long double> &mapSNPid_Pmi_Y);
  
  
  void writePoolofSize(bool &flag, int maxSize);
  void regionStrategy(const std::string method_Name, Model &model, std::multimap<long double, snp_index_t>& Pmi_Ysort, std::map<snp_index_t, long double> &mapSNPid_Pmi_Y, double th, std::set<TRegionSet_info> &setNewRegion);
  void calculatePOWER_FDR_clustRegion(std::vector<snp_index_t> &realSNPs, TPOWER_FDR &powerFDR, std::set<TRegionSet_info> &setNewRegion, std::map<snp_index_t, int>& recognizedSNPs);
  void test_region();
};
//  void calculatePOWER_FDR_clustGA(std::set<snp_index_t> &mySnps, std::vector<snp_index_t> &causalSNPs, TPOWER_FDR &powerFDR, std::set<snp_index_t> &TP_SNPs);
void writePosterior(std::string fileName, std::map<snp_index_t, int> &mapSNPCausal_ind, std::vector< std::multiset<long double> > &tabCausalPost, int size);

}

#endif

