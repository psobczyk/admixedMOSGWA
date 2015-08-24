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

#ifndef SEARCH_MEMETICA_MA_HPP
#define SEARCH_MEMETICA_MA_HPP

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

ostream &operator<< ( std::ostream &out, std::vector<size_t> &v );

namespace memetica {
  
struct TSNP_Info 
{
  size_t SNP;
  unsigned int Chr;
  unsigned int pos;
  long double posterior;
  
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

/** Converts the seconds to the time format */
std::string sec2time(const time_t &t);


/** outputs a map on the screen */
template <typename T1, typename T2>
std::ostream &operator<< (std::ostream &out, std::map<T1, T2> &v)
{
  for (typename std::map<T1, T2>::iterator it = v.begin(); it != v.end(); it++)
    out << "[" << it->first << " -> " << it->second << "]" << endl;
  return out;
}

/** outputs a vector on the screen */
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
  {
    out << ", " << *it;
  }
  return out << "]";
}

/** outputs a set on the screen */
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
  {
    out << ", " << *it;
  }
  return out << "}";
}


struct THeri_PMiY
{
public:
  long double h2M;         // a heriability
  long double PMi_Y;       // a posterior
  THeri_PMiY(long double h2M = 0.0L, long double PMi_Y = 0.0L): h2M(h2M), PMi_Y(PMi_Y) {}
  bool operator < (const THeri_PMiY & t) const {return h2M < t.h2M;}
};


class MA  : public search::Search
{
private:
  
  /** a set of exclusived SNPs - for the initial population */
  TBitset exclusivedSNP;   
  
  /** a set of SNPs which p-val are > 0.1 */
  TBitset goodSNPs;        
  
  /** a number of SNPs wchich p-val are > 0.1 */
  int goodSNPsNo;          
  
  /** the "table" of correration values for snps in all models */
  static std::vector< std::vector<size_t> > correlations; 

  /** Data for models */
  MData data;

  /** Models */
  Model **models;

  /** The number of models */
  unsigned int modelsNo;

  /** The pool of models */
  std::set<PoolItem> pool;

  /** The maximum number of iterations without the progress */
  unsigned int maxNoProgressIter;

  /* The maximum size of a pool */
  unsigned int maxPoolSize;  
  
  /** The probabilitiy of crossover */
  double pCross;

  /** The probabilitiy of mutation */
  double pMutation;

  /** The tournament size */
  unsigned int tournamentSize;

  /** The correlation threshold - it's used in a recombination and a local improvement function */
  double correlationThreshold;

  /** SNPs for the correlation are from range [snp - correlationRange, snp + correlationRange] */
  int correlationRange;
  
  double regionMinCorrelation;  // threshold for a region method

  /** For the stop criterion.*/
  int B;
  
  /** The result model */
  Model resultModel;
  
  /** SNPs calculatee by MOSGWA (a stepwise selection). poolReader() function also sets up this variable */
  std::vector<size_t> v_mosgwa;

  std::set<size_t> testSet;  

  long double RSSo;

  
  bool isResult;
  
  /** Mutates the model */
  void mutation ( Model *aModel, double pMutation, double threshold );

  /** Makes one child from two parents */
  Model* recombination ( const Model & s1, const Model & s2 );

  /** Finds and returns the index of the worst model in the population */
  unsigned int findTheWorst () const;

  /** The tournament selection */
  Model* tournamentSelection ( unsigned int tournamentSize ) const;

  /** Finds the best, the worst and an average model */
  void statistics ( double &theBest, double &average, double &theWorst );

  /** calculates a correlation for a snp on a snpPosition position - relative position of the snp at modelsnps_ vector */
  std::vector<size_t> stronglyCorrelatedSnpsAt ( Model *model, const int& snpPosition, const double& threshold, int correlationRange );

  /** calculates correlation for a snp */
  std::vector<size_t> stronglyCorrelatedSnps ( Model *model, const int& snp, const double& threshold, int correlationRange );

  /** a local inprovement of a given model */
  void localImprovement ( Model *model, double threshold, int correlationRange );
  
  /** Returns number of models which have got better msc value than the parameter mscVal */
  int isInNBestModels(double mscValue);

  /** writes a model to the pool */
  void toPool(const Model *model, char c);

  /** computes the heritabilities */
  void coumputeHeritability(std::stringstream &sp_sort, const std::vector<long double> &diff, std::vector<THeri_PMiY> &tab);
  
  /**  selects a model by the region strategy method */
  void regionStrategy(std::stringstream & ss_Region, Model &model, std::multimap<long double, size_t>& Pmi_Ysort, std::map<size_t, long double> &mapSNPid_Pmi_Y, double th, std::set<TRegionSet_info> &setNewRegion);

  /** selects model for the initial population) */
  void selectModel ( Model& selectedModel );
  
  /** computes the results and writes them to the output files */
  void computeAndWriteResults(); 
  
    
public:

  /* construct the MA object */
  MA ();
  /** writes a pool to a file */
  void writePoolToFile (std::stringstream &ssModels, std::string postfix = "") const;

  /** runs the memetic algorithm */
  virtual void run ();

  /** destructs the MA object */
  virtual ~MA ();
  
  /** */
  virtual const Model * result();
};

}
#endif
