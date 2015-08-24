#include <algorithm>
#include <iterator>
#include <iomanip>
#include <cstdlib>                     
#include <set>
#include <list>
#include <cassert>
#include <fstream>

#include "MA.hpp"  
#include "../../SNP.hpp"
#include "../../logging/Logger.hpp"  
#include "../../Parameter.hpp"  
#include "../../Exception.hpp"

using namespace logging;


namespace memetica {
  
  
const int minModelSize = 3;                    // the minimal model size for the initial population
const double p_value_threshold = 0.1;          // MA uses only SNPs which p_value is smalle than p_value_threshold

vector<vector<size_t> > MA::correlations; // a correlations vector

/**
 * @brief Prints a given vector of SNPs on the screen
 * @param v - a vector to print
 * @returns a reference to ostream object
 */ 
ostream &operator<< (ostream &out, vector<size_t> &v)
{
  out << "[";
  vector<size_t>::iterator it = v.begin();
  if (v.size() > 0)
  {
    out << *it;
    it++;
  }
  for (; it != v.end(); it++)
    out << ", " << *it;
  return out << "]";
}

bool TRegionSet_info::operator < (const TRegionSet_info &t) const
{
  set<TSNP_Info>::iterator it_2 = t.s->begin();
  set<TSNP_Info>::iterator it = s->begin();
  for (; it != s->end() && it_2 != t.s->end(); ++it, ++it_2)
  {
    TSNP_Info s1 = *it;
    TSNP_Info s2 = *it_2;
    if (s1.Chr == s2.Chr)
    {
      if (s1.pos < s2.pos)
        return true;
      else
        if (s1.pos > s2.pos)
          return false;
    }
    else    
      return s1.Chr < s2.Chr;   
  }
  if (it == s->end() && it_2 != t.s->end())
    return true;  // this->s is a subset of t.s
  else
    return false;
}

/**
 * @brief Returns a time staple as a string
 */
string timeStemple()
{
  const time_t ltime=time(NULL);          //get current calendar time
  const tm* now = localtime(&ltime);      // struct for the day, year....
  return int2strPadWith( now->tm_hour, 2, '0' ) + ":" + int2strPadWith( now->tm_min, 2, '0' ) 
         + ":" + int2strPadWith( now->tm_sec, 2, '0' );  
}

/**
 * @brief Converts the seconds to the time format
 * @param t - time to convert
 * @result a converted time as a string
 */
string sec2time(const time_t &t)
{
  int hour = t / 3600;
  int sec = t - hour * 3600;
  int min = sec / 60;
  sec = sec % 60;
  return int2strPadWith( hour, 2, '0' ) + ":" + int2strPadWith( min, 2, '0' ) 
         + ":" + int2strPadWith( sec, 2, '0' );    
}


/**
 * @brief Creates a progress bar. It is helpfull to choose the values for the stop criterion
 * @param current - the current value of progress
 * @param max - the maximum value of progress
 * @param size - the width of the progress bar
 * @param first - the char which indicates a progress
 * @param second - the char which indicates no progress
 * @return a progress bar as a string
 */
string progressBar(unsigned int current, unsigned int max, unsigned int size, char first = '>', char second = '-')
{
  unsigned int percent = current * 100.0 / max + 0.5;
  unsigned int pos = size * current / max;
  string bar;
  bar = "";
  for (unsigned int i = 0; i < size; ++i)  
    if (i < pos)
      bar += first;
    else
      bar += second;
  if (current > max)
  {
    current = max;
  }
  stringstream returnBar;
  returnBar << " complete [" << bar << "] " << setw(3) << percent << "%\t";
  return returnBar.str();
}

/**
 * @brief normCDF function 
 */
double normCDF(double z)
{
  if (z > 6.0) return 1.0;
  if (z < -6.0) return 0.0;
  
  double b1 = 0.31938153,
         b2 = -0.356563782,
         b3 = 1.781477937,
         b4 = -1.821255978,
         b5 = 1.330274429,
         p  = 0.2316419,
         c2 = 0.3989423;
         
  double a = fabs(z);
  double t = 1.0/(1.0 + a * p);
  double b = c2 * exp((-z) * (z/2.0));
  double n = ((((b5 * t + b4) * t + b3) * t + b2) * t + b1) * t;
  n = 1.0 - b * n;
  if (z < 0.0)
    n = 1.0 - n;
  return n;
}


/**
 * @brief Constructs the Memetic Algoritm object.
 * All parameters come from config file or take default values
*/
MA::MA()
 : data(), 
   models(0),
   modelsNo(parameter->modelsNo),
   maxNoProgressIter(parameter->maxNoProgressIter),
   maxPoolSize(parameter->maxPoolSize),
   pCross(parameter->pCross),
   pMutation(parameter->pMutation),
   tournamentSize(parameter->tournamentSize),
   correlationThreshold(parameter->correlationThreshold),
   correlationRange(parameter->correlationRange),
   regionMinCorrelation(parameter->regionMinCorrelation),
   B(parameter->B),
   resultModel(data),
   isResult(false)
{
  logger->info( "Start memetic algorithm search" );
  logger->info( "modelsNo = %u", modelsNo );
  logger->info( "B = %d", B );  
  logger->info( "maxNoProgressIter = %u", maxNoProgressIter );
  logger->info( "pCross = %f", pCross );
  logger->info( "pMutation = %f", pMutation );
  logger->info( "tournamentSize = %u", tournamentSize );
  logger->info( "correlationThreshold = %f", correlationThreshold );
  logger->info( "correlationRange = %d", correlationRange );
  logger->info( "regionMinCorrelation = %f", regionMinCorrelation );
  logger->info( "maxPoolSize = %d", maxPoolSize );

  time_t now;
  time(&now);
  srand(now);                                  // runs the random generator

  std::ofstream file;  
  file.open((parameter->out_file_name + "_Report.txt").c_str(), ios::trunc);
  file.close();        

  if (modelsNo <= 0)
  {
    throw Exception("A number of models (%d) is too small. Set the number of models greater than 0", modelsNo);
  }
  if ((unsigned int) (B) > modelsNo)
  {
    throw Exception ("A parameter B (B = %d) must be greater or equal to the number of models (%d).\r\nSet up a valid value of B in your config file.", B, modelsNo);
  }

  data.calculateIndividualTests();

  Model m0(data);
  m0.computeRegression();
  data.setLL0M( m0.getMJC() );
  
  if ((unsigned int) (parameter->ms_MaximalSNPsMultiForwardStep) > data.getIdvNo()/4)
  {
    parameter->ms_MaximalSNPsMultiForwardStep = data.getIdvNo()/4;
  }
  
  if (parameter->silent == false)
  {
    logger->info("maxInitModelSize: %u, individuals: %u", parameter->ms_MaximalSNPsMultiForwardStep, data.getIdvNo());
  }
  
  correlations.resize(data.getSnpNo());
  exclusivedSNP.reset();
  goodSNPs.reset();
  goodSNPsNo = 0;

  for (unsigned int i = 0; i < data.getSnpNo(); ++i)
  {
    if ( data.getSingleMarkerTestAt( data.getOrderedSNP(goodSNPsNo) ) < p_value_threshold )
    {
      goodSNPs[data.getOrderedSNP(goodSNPsNo)] = true;
      ++goodSNPsNo;
    }
    else
      break;
  }
  
  if (parameter->silent == false)
  {
    logger->info("The number of SNPs which the p-value is greater than %d: %d", p_value_threshold, goodSNPsNo);
  }
  
  vector<size_t> v;
  m0.computeMSC();
  RSSo = m0.getMJC();
  double h2_M = (RSSo - m0.getMJC()) / RSSo;
  unsigned int pSize = pool.size(); 
  v = m0.getModelSnps();
  PoolItem p(v, m0.getMSC(), h2_M, pSize + 1, 'I');
  pool.insert(p);
  
  try
  {
    models = new Model*[modelsNo];
  }
  catch (bad_alloc &ba)
  {
    throw Exception ("MA(...): Not enought memory to create %u models.\r\nTry to use a smaller number of models", modelsNo);
  }
  
  vector<size_t>::iterator it;
  vector<size_t> snps;
  if (parameter->silent == false)
  {
    logger->debug("The maximum model size is: %u", parameter->maximalModelSize);
  }
  // creates the first model by the stepwise method 
  try
  {
    models[0] = new Model(data);                       
    models[0]->computeRegression();
    data.setLL0M(models[0]->getMJC());
    data.selectModel(
      models[0],
      parameter->PValueBorder,                         
      parameter->maximalModelSize,                     
      Parameter::selectionCriterium_mBIC_firstRound   
    );
    data.selectModel(models[0], 5000, parameter->maximalModelSize, Parameter::selectionCriterium_mBIC2); 
    vector<size_t> v = models[0]->getModelSnps();
    v_mosgwa = v;
    for (unsigned int i = 0; i < v.size(); ++i)
    {
      exclusivedSNP[ v[i] ] = true;
    }
  }
  catch (bad_alloc &ba)
  {
    throw Exception ("MA(...): Not enought memory to create a %uth model: ", 0);
  }
  // creates the first half of populations
  for (unsigned int i = 1; i < modelsNo/2; i++)
  {
    try
    {
      models[i] = new Model(data);
      selectModel(*models[i]);
    }
    catch (bad_alloc &ba)
    {
      throw Exception ("MA(...): Not enought memory to create a %uth model.", i);      
    }
  }
  // creates the second half of population
  for (unsigned int i = modelsNo/2; i < modelsNo; i++)
  {
    try
    {
      models[i] = new Model(data);
    }
    catch (bad_alloc &ba)
    {
      throw Exception ("MA(...): Not enought memory to create a %uth model.", i);
    }
    set<size_t> initSNP;
    int randSNP;
    unsigned int randSNPno = 5;  // the maximum model size
    while ( initSNP.size() < randSNPno )
    { 
      randSNP = rand() % goodSNPsNo;
      if (exclusivedSNP[randSNP] == true)
      {
        continue;
      }
      initSNP.insert(randSNP);
      exclusivedSNP[randSNP] = true;
    }
    models[i]->createFromSNPs(initSNP);
  }
  // writes the models to the pool file
  stringstream ss;  
  for (unsigned int i = 0; i < modelsNo; i++)
  {
    if (models[i]->computeRegression() == false)
    {
      models[i]->computeMSCfalseRegression();
    }
    else
    {
      models[i]->computeMSC();
    }
    double h2_M = (RSSo - models[i]->getMJC()) / RSSo;
    unsigned int pSize = pool.size(); 
    v = models[i]->getModelSnps();
    PoolItem p(v, models[i]->getMSC(), h2_M, pSize + 1, 'I');
    if (i == 0)
    {
      toPool(models[i], 'S');
    }
    else  
    {
      toPool(models[i], 'I');
    }
    ss << endl << (i + 1) << "). " << p;
  }
  logger->debug("MEMETIC ALGORITHM - THE INITIAL POPULATION:");
  logger->debug("%s", ss.str().c_str());
}

/**
 * @brief Writes a model to the pool
 * @param model - a model
 * @param c - an information where the model was created: 
 *                           'I' - a model was created during creating an initial population
 *                           'S' - the model was created by a greedy search of MOSGWA
 *                           'B' - a model was created in a backward step of a reproduction
 *                           'F' - a model was created in a forward step of a reproduction
 *                           'L' - a model was created in a local improvement step
 *                           'A' - a model was created by adding a SNP in a mutation step
 *                           'E' - a model was created by erasing a SNP in a mutation step
 *                           'C' - a model was createt by changing a SNP in a mutation step
 */
void MA::toPool(const Model *model, char c)
{
  if (pool.size() < maxPoolSize)
  {
    vector<size_t> v = model->getModelSnps();
    double h2_M = (RSSo - model->getMJC()) / RSSo;
    unsigned int pSize = pool.size(); 
    PoolItem p(v, model->getMSC(), h2_M, pSize + 1, c);
    pool.insert(p);
  }
}

/**
 * @brief Mutates a given model. The model is modified.
 * @param aModel pointer to a model
 * @param pMutation a probability of mutation
 * @param threshold - a threshold for the correlated snps. If the threshold is less or equal 0 this method does not take into account correlated snps
 */
void MA::mutation(Model *model, double pMutation, double threshold)
{
  int snpsNo = data.getSnpNo();
  vector<size_t> snpsVector = model->getModelSnps();
  vector<size_t>::iterator it;
  set<size_t> snpsSet(snpsVector.begin(), snpsVector.end());
  set<size_t>::iterator itSet;
  vector<size_t> v;
  int oneSnp;
  int maxTrial = 1000;
  if (random() % 2)
  {  // adds a gen to the model
    if (threshold > 0)
    {
      for (it = snpsVector.begin(); it != snpsVector.end(); ++it)
      {
        if (correlations[*it].size() == 0)
          correlations[*it] = stronglyCorrelatedSnps(model, *it, threshold, correlationRange);
        snpsSet.insert(correlations[*it].begin(), correlations[*it].end());
      } // snpsSet - set of SNPs from the model and all SNPs which are correlatad with them 
    }
    int noTrial = 0;
    do
    {
      oneSnp = random() % snpsNo;
      itSet = snpsSet.find(oneSnp);
      ++noTrial;
    }
    while (itSet != snpsSet.end() && noTrial < maxTrial);
    if (noTrial < maxTrial)
    {  
      model->addSNPtoModel(oneSnp);
      if (model->computeRegression())
        model->computeMSC();
      else
        model->computeMSCfalseRegression();
    }  
    toPool(model, 'A');  
  }
  else
  {  
    if (model->getModelSize() > 1)
    { // erases a gen from the model
      oneSnp = random() % model->getModelSize();
      model->removeSNPfromModel(oneSnp);
      if (model->computeRegression() == false)
      {
        model->computeMSCfalseRegression();
      }
      else 
      {
        model->computeMSC();
      }
      toPool(model, 'E');  
    }
    else
    { // changes a gen to anoter (if the model has got only one SNP, then this SNP is changed to another)
      if (model->getModelSize() != 0)
      {
        model->removeSNPfromModel(0);
      }
      oneSnp = rand() % goodSNPsNo;
      model->addSNPtoModel(data.getOrderedSNP(oneSnp));
      if (model->computeRegression() == false)
      {
        model->computeMSCfalseRegression();
      }
      else 
      {
        model->computeMSC();
      }
      toPool(model, 'C');  
    }
  }
}

/**
 * @brief Makes a child model from two parents models. MSC is calculated
 * @param s1 - the first parent model
 * @param s2 - the second parent model
 * @returns a pointer to a child model
 * WARNING Allocates memory for new Model. 
 */                                                               
Model* MA::recombination(const Model & s1, const Model & s2) 
{
  vector<size_t> v1 = s1.getModelSnps();       
  vector<size_t> v2 = s2.getModelSnps();
    
  set<size_t> s_1(v1.begin(), v1.end());       
  set<size_t> s_2(v2.begin(), v2.end());
  
  set<size_t> SI;                              // an intersection of s1 and s2
  set_intersection(s_1.begin(), s_1.end(), s_2.begin(), s_2.end(), insert_iterator<set<size_t> >(SI, SI.begin()));
  set<size_t> SU;                              // an union of s1 and s2
  set_union(s_1.begin(), s_1.end(), s_2.begin(), s_2.end(), insert_iterator<set<size_t> >(SU, SU.begin()));
  set<size_t> SD;                              // a symmetric difference of s1 and s2
  set_difference(SU.begin(), SU.end(), SI.begin(), SI.end(), insert_iterator<set<size_t> >(SD, SD.begin()));

  set<size_t>::iterator it;
  int addedSNP;
  set<size_t> tempSD = SD;
  vector<size_t> v;

  // intersection and forward selection with symmeric difference
  Model* modelForward = new Model(data);
  modelForward->createFromSNPs(SI);
  if (modelForward->computeRegression())
  {
    modelForward->computeMSC();
  }
  else
  {
    modelForward->computeMSCfalseRegression();
  }
  double bestMSC = modelForward->getMSC();
  do
  {
    addedSNP = -1;
    for (it = tempSD.begin(); it != tempSD.end(); ++it)
    {
      modelForward->addSNPtoModel(*it);
      if (modelForward->computeRegression())
      {
        modelForward->computeMSC();
      }
      else
      {
        modelForward->computeMSCfalseRegression();
      }
      toPool(modelForward, 'F');  
      if (modelForward->getMSC() < bestMSC)
      {
        addedSNP = *it;
        bestMSC = modelForward->getMSC();
      }
      modelForward->removeSNPValFromModel(*it);
    }
    if (addedSNP >= 0)
    {
      tempSD.erase(addedSNP);
      modelForward->addSNPtoModel(addedSNP);
      if (modelForward->computeRegression())
      {
        modelForward->computeMSC();
      }
      else
      {
        modelForward->computeMSCfalseRegression();
      }
    }
  }
  while (addedSNP >= 0);
  if (modelForward->computeRegression())
  {
    modelForward->computeMSC();
  }
  else
  {
    modelForward->computeMSCfalseRegression();
  }
  if (SD.size() > data.getIdvNo() / 2)           
  {
    return modelForward;
  }
  // an union and the backward selection form the symetric difference     
  Model* modelBackward = new Model(data); 
  modelBackward->createFromSNPs(SU);       
  if (modelBackward->computeRegression())
  {
    modelBackward->computeMSC();
  }
  else
  {
    modelBackward->computeMSCfalseRegression();
  }
  set<size_t> bestBackward = SU;
  set<size_t> currentSNPs = SU;
  bestMSC = modelBackward->getMSC();
  int erasedSNP;
  tempSD = SD;
  do
  {
    erasedSNP = -1;
    for (it = tempSD.begin(); it != tempSD.end(); ++it)
    {
      modelBackward->removeSNPValFromModel(*it);
      if (modelBackward->computeRegression())
      {
        modelBackward->computeMSC();
      }
      else
      {
        modelBackward->computeMSCfalseRegression();
      }
      toPool(modelBackward, 'B');  
      if (modelBackward->getMSC() < bestMSC)
      {
        bestMSC = modelBackward->getMSC();
        erasedSNP = *it;
      }
      modelBackward->addSNPtoModel(*it);
    }
    if (erasedSNP >= 0)
    {
      tempSD.erase(erasedSNP);
      modelBackward->removeSNPValFromModel(erasedSNP);
      if (modelBackward->computeRegression())
      {
        modelBackward->computeMSC();
      }
      else
      {
        modelBackward->computeMSCfalseRegression();
      }
    }
  }
  while (erasedSNP >= 0);  
  
  if (modelBackward->computeRegression())
  {
    modelBackward->computeMSC();
  }
  else
  {
    modelBackward->computeMSCfalseRegression();
  }
  if (modelBackward->getMSC() < modelForward->getMSC())
  {
    delete modelForward;
    return modelBackward;
  }
  else
  {
    delete modelBackward;
    return modelForward;
  }
}

/**
 * @brief Finds an index to the worst model in the population
 * @returns an index to the worst model in the population
 */
unsigned int MA::findTheWorst() const
{
  unsigned int theWorst = 0;
  for (unsigned int i = 1; i < modelsNo; i++)
  {
    if (models[i]->getMSC() > models[theWorst]->getMSC())
    {
      theWorst  = i;
    }
  }
  return theWorst;  
}

/**
 * @brief Makes tournament selection.
 * @param tournamentSize - size of tournament
 * @returns pointer to a winner model.
 */ 
Model* MA::tournamentSelection(unsigned int tournamentSize) const
{
   if (tournamentSize <= 1)
   {
     throw Exception ("The tournament size (%u) must be greater than 0.\n\rSet up the valid value in your config file.", tournamentSize);
   }
   set<int> tournament;
   set<int>::iterator it;
   int randModel;
   int theBest = 0;
   unsigned int i = 0;
   
   while (i < tournamentSize)
   {
     randModel = rand() % modelsNo;
     it = tournament.find(randModel);
     if (it == tournament.end())
     {
       tournament.insert(randModel);
       if (i == 0 || models[randModel]->getMSC() < models[theBest]->getMSC())
       {
          theBest = randModel;
       }  
       ++i;
     }
   }
   return models[theBest];
}

/**
 * @brief Sorting function for the Pool 
 * @returns boolean value of comparison a < b
 * @note Takes into acount: first  - size of strings
 *                          second - alphabetical order
 */
bool sort_string(string a, string b)
{
  if (a.length() == b.length())
  {
     return a < b;
  }
  return a.length() < b.length();
}


/**
 * @brief Runs the Memetic Algoritm
 */
void MA::run()
{
  unsigned int MAcount = 0;
  unsigned int iterCount = 0;
  Model *firstParent = 0,                  // parents for recombination 
        *secondParent = 0;
  Model *childModel = 0;
  vector<size_t> v;
  double best, adv, worst;
  statistics(best, adv, worst);
  stringstream ssBegin;
  //cout << endl;
  ssBegin << "Initial population, the worst: " << worst << ", adv: " << adv << ", the best: " << best << endl;
  logger->info("%s", ssBegin.str().c_str());

  while (MAcount < maxNoProgressIter &&  pool.size() < maxPoolSize)
  {
    firstParent = tournamentSelection(tournamentSize);
    do  
    {
      secondParent = tournamentSelection(tournamentSize);
    } while (firstParent == secondParent);
    
    if (random() < pCross * RAND_MAX) 
    { // a recombination
      childModel = recombination(*firstParent, *secondParent); // allocates memory   
      if (random() < pMutation * RAND_MAX)
      {
        mutation(childModel, pMutation, correlationThreshold);   
      }
    }
    else
    { // only mutation
      childModel = new Model(data);
      if (firstParent->getMSC() > secondParent->getMSC())
      {
        *childModel = *firstParent;
      }
      else
      {
        *childModel = *secondParent;
      }
      mutation(childModel, pMutation, correlationThreshold);
    }
    if (childModel->computeRegression())
    {
      childModel->computeMSC();
    }
    else
    {
      childModel->computeMSCfalseRegression();
    }
    localImprovement(childModel, correlationThreshold, correlationRange);
    unsigned int theWorst = findTheWorst();
    int rankB = isInNBestModels(childModel->getMSC());
    if (childModel->getMSC() < models[theWorst]->getMSC())
    { // we have a better model than the worst model in the population
      unsigned int toChange = 0;
      // We don't want to have two or more the same models in the population. This code prevents it.
      while (toChange < modelsNo && *childModel != *models[toChange])
      {
         ++toChange;
      }
      if (toChange == modelsNo && childModel->getModelSize() > 0)
      {  
/*        if (parameter->silent == false)
        {
          cout << endl << "We have got a better model " << endl;
        }
*/        delete models[theWorst];
        models[theWorst] = childModel;
        childModel = 0;      
        if (parameter->silent == false)
        {
          stringstream ss;
          ss << endl << "iter count: " << iterCount << ":" << endl;
          for (unsigned int m = 0; m < modelsNo; ++m)
          {
            ss << m << ") " << *models[m] << endl;
          }
          logger->info( "%s", ss.str().c_str() );
        }  
        if (rankB < B)
        {
          MAcount = 0;
        }
      }
      else
      {
        delete childModel;
        childModel = 0;
        ++MAcount;
      }
    } 
    else
    {
      ++MAcount;
      if (childModel != 0)
      {
        delete childModel;
      }
      childModel = 0;
    }  
    ++iterCount;
    if (iterCount % 100 == 0)
    {
      double best, adv, worst;
      statistics(best, adv, worst);
      stringstream ss;
/*      cout << timeStemple();
      if ((MAcount + 0.0)/maxNoProgressIter >=  (pool.size() + 0.0)/maxPoolSize)
      {
        cout << progressBar(MAcount, parameter->maxNoProgressIter, 30);
      }
      else
      {
        cout << progressBar(pool.size(), maxPoolSize, 30);
      }
      cout << "iterations: " << setw(5) << iterCount << ", MSC of the worst model: " << setprecision(10) 
           << worst << ", ave: " << adv << ", the best: " << best << ", the pool size: " << setw(6) << pool.size() << endl;
      if (parameter->silent == false) 
      {
        for (unsigned int m = 0; m < modelsNo; ++m)
        {
          cout << m << ") " << *models[m] << endl;
        }     
      }
*/      
      if ((MAcount + 0.0)/maxNoProgressIter >=  (pool.size() + 0.0)/maxPoolSize)
      {
        ss << progressBar(MAcount, parameter->maxNoProgressIter, 30) << " iteration: " << setw(5) << iterCount << ", the worst: " << setprecision(10) 
        << worst << ", ave: " << adv << ", the best: " << best << ", poolSize: " << setw(6) << pool.size();// << endl;
      }
      else
      {
        ss << progressBar(pool.size(), maxPoolSize, 30) << " iteration: " << setw(5) << iterCount << ", the worst: " << setprecision(10) 
        << worst << ", ave: " << adv << ", the best: " << best << ", pool size: " << setw(6) << pool.size();// << endl;
      }
      logger->info("%s", ss.str().c_str());
    }
  }
  statistics(best, adv, worst);
  logger->info("At the end of MA, after %u ienerations, msc of the worst model: %d, adv: %d, the best: %u, the pool size: %u", iterCount, worst, adv, best, pool.size() );

  stringstream ssFinal;
  for (unsigned int i = 0; i < modelsNo; i++)
  {
    if (models[i]->computeRegression() == false)
    {
      models[i]->computeMSCfalseRegression();
    }
    else
    {
      models[i]->computeMSC();
    }
    ssFinal << endl << (i + 1) << "] " << *models[i] << "\tsize: " << models[i]->getModelSize();
  }
  logger->info("The end of MA");
  computeAndWriteResults();
}

/**
 * @brief Finds the best, the worst and the average value of mBIc2 
 * @param theBest - a reference to a variable of the best msc
 * @param average -  a reference to a variable of the average msc
 * @param theWorst - a reference to a variable of the worst msc
 */
void MA::statistics(double &theBest, double &average, double &theWorst) 
{
  double sum = 0.0;
  theBest = theWorst = sum = models[0]->getMSC();
  double msc;
  for (unsigned int i = 1; i < modelsNo; ++i)
  {
    msc = models[i]->getMSC();
    if (msc < theBest)
    {
      theBest = msc;
    }
    else
    {
      if (msc > theWorst)
      {
        theWorst = msc;
      }
    }
    sum += msc;  
  }
  average = sum / modelsNo;
}

/**
 * @brief The destructor
 */
MA::~MA()
{
  if (models != 0)
  {
//    cout << "END" << endl;
    for (unsigned int i = 0; i < modelsNo; ++i)
    {
      if (models[i])
      {
        delete models[i];
      }
    }
    delete []models;    
    models = 0;
  }
}

/**
 * @brief Computes a correlation for a given SNP (referred by an index 'snp') and a given threshold value
 * @param model - a model
 * @param snp - an index of SNP in the 'model' parameter
 * @param threshold - a threshold value of a correlation. 
 * @return a vector of SNPs which are correlated with 'snp' above a threshold value.
 */
vector<size_t> MA::stronglyCorrelatedSnps(Model *model, const int& snp, const double& threshold, int correlationRange )
{
  vector<size_t> snps = model->getModelSnps();
  vector<size_t>::iterator it;
  it = find(snps.begin(), snps.end(), snp);
  if ((unsigned int) (it - snps.begin()) < snps.size())  
  {
    return stronglyCorrelatedSnpsAt(model, it - snps.begin(), threshold, correlationRange);
  }
  else
  {
    snps.clear();  // a snp is not in the model, no correlated snps, so this method returns the empty vector
    return snps;
  }
}  

/**
 * @brief Computes a correlation for a snp on a given position with a given threshold value
 * @param model - a model to compute a correlation
 * @param snpPosition - it is a relative position of a SNP at the model
 * @param threshold - a threshold value of correlation. 
 * @return a vector of snps with a correlation above a threshold
 */
vector<size_t> MA::stronglyCorrelatedSnpsAt(Model *model, const int& snpPosition, const double& threshold, int correlationRange )
{
  if (snpPosition < 0 || snpPosition >= model->getModelSize())
  {
    throw Exception ("stronglyCorrelatedSnps(): A SNP at the position: %u is out of range [0, %u]", snpPosition, model->getModelSize());
  }
  stringstream ss;                 // to save the output
  multimap<double, int> StrongCor; // to sort the correlated SNPs

  double abscor;                   // abs|Correlation| of two SNPs
  unsigned int j;                           
  
  unsigned int rangeFrom;
  if (model->getSNPat(snpPosition) - correlationRange < 0)
  {
    rangeFrom = 0;
  }
  else
  {
    rangeFrom = model->getSNPat(snpPosition) - correlationRange;
  }
  if (rangeFrom < 0)
  {
    rangeFrom = 0;
  }
  unsigned int rangeTo = model->getSNPat(snpPosition) + correlationRange;
  if (rangeTo > data.getSnpNo())
  {
    rangeTo = data.getSnpNo();
  }
  // searches for the strongly correlated SNPs
  for (j = rangeFrom; j < rangeTo; j++)
  {
    abscor = fabs( data.computeCorrelation( model->getSNPat(snpPosition), j ) ); // computes a correlation between a SNP which is from the model and SNP_j
    if (abscor  >= threshold)                                                    // inserts if a correlation is big enough
    {
      StrongCor.insert(pair<double, int>(abscor, j));
    }
  }
  vector<size_t> correratedSnps(StrongCor.size());
  int i = 0;
  for (multimap<double, int>::reverse_iterator it = StrongCor.rbegin(); it != StrongCor.rend(); it++, ++i)
  {
    correratedSnps[i] = (*it).second;
  }
  return correratedSnps;
}  

/**
 * @brief Makes a local improvement of a given model
 * @param model - a model to improve
 * @param threshold - a minimal correlation value
 * @param correlationRange - a correlation range. For a given SNP, the corelations are computed in the range [SNP - correlationRange, SNP + correlationRange]
 */
void MA::localImprovement(Model *model, double threshold, int correlationRange) 
{
  vector <size_t> snps = model->getModelSnps();
  vector <size_t>::iterator it;                        // snp to improvement
  vector<size_t> v;
  v = model->getModelSnps();
  vector<size_t> v_pool;
  if (model->computeRegression())
  {
    model->computeMSC();
  }
  else
  {
    model->computeMSCfalseRegression();
  }
  double bestMSC = model->getMSC();
  int changedSnp;     
  vector<size_t>::iterator itCorrelation;
  for (it = snps.begin(); it != snps.end(); ++it)  // for every snp in the model
  {                                                // tests snp in it
    if (correlations[*it].size() == 0)
    {
      correlations[*it] = stronglyCorrelatedSnps(model, *it, threshold, correlationRange);
    }
    changedSnp = -1;
    for (itCorrelation = correlations[*it].begin(); itCorrelation != correlations[*it].end(); ++itCorrelation)
    {
      model->removeSNPValFromModel(*it);            // removes snp from model           
      model->addSNPtoModel(*itCorrelation);         // adds correlated snp to model    
      
      v = model->getModelSnps();
      testSet.clear(); 
      testSet.insert(v.begin(), v.end());
      if (v.size() != testSet.size())
      {  
        model->createFromSNPs(testSet);
      }  
      if (model->computeRegression())
      {
        model->computeMSC();
      }
      else
      {
        model->computeMSCfalseRegression();                          // calculates MSC
      }
      toPool(model, 'L');  
      if (model->getMSC() < bestMSC)
      {
        bestMSC = model->getMSC();
        changedSnp = *itCorrelation;
      }
      model->removeSNPValFromModel(*itCorrelation);
      model->addSNPtoModel(*it);
    }
    if (changedSnp >= 0)
    {
      model->removeSNPValFromModel(*it);
      model->addSNPtoModel(changedSnp); 
      v = model->getModelSnps();      
      testSet.clear(); 
      testSet.insert(v.begin(), v.end());
      if (v.size() != testSet.size())
      {  
        model->createFromSNPs(testSet);
      }  
      if (model->computeRegression())
      {
        model->computeMSC();
      }
      else
      {
        model->computeMSCfalseRegression();            
      }
    }
  }
  if (model->computeRegression())
  {
    model->computeMSC();
  }
  else
  {
    model->computeMSCfalseRegression();
  }
}

/**
 * @brief Writes the pool of models to a file
 * @param ssTime - the runtime of MA
 * @param postfix - a text which is added to the end of the file name.
 */
void MA::writePoolToFile(stringstream &ssTime, string postfix) const
{
  set<PoolItem>::iterator poolIt;
  int i = 1;
  ofstream poolFile;
  poolFile.exceptions ( ofstream::eofbit | ofstream::failbit | ofstream::badbit ); // checks if poolfile can be written
  try
  {
    poolFile.open( ( parameter->out_file_name + "_pool" + postfix + ".txt" ).c_str(),  fstream::out | fstream::trunc );
    for (poolIt = pool.begin(); poolIt != pool.end(); poolIt++, ++i)
    {
      PoolItem p1 = *poolIt;
      poolFile << p1 << "\n";
    }  
    poolFile.flush();
    poolFile << "END\t" << endl << endl
       << "LEGEND: " << endl
       << "I - a model was created during creating an initial population" << endl
       << "S - the model was created by a greedy search of MOSGWA" << endl
       << "B - a model was created in a backward step of a reproduction" << endl
       << "F - a model was created in a forward step of a reproduction" << endl
       << "L - a model was created in a local improvement step" << endl
       << "A - a model was created by adding a SNP in a mutation step" << endl
       << "E - a model was created by erasing a SNP in a mutation step" << endl
       << "C - a model was createt by changing a SNP in a mutation step" << endl << endl
       << "poolSize\t" << pool.size()  << endl << endl
       << "modelsNo\t" << parameter->modelsNo << endl
       << "maxNoProgressIter\t" << parameter->maxNoProgressIter << endl
       << "maxPoolSize\t" << parameter->maxPoolSize << endl
       << "pCross\t" << parameter->pCross << endl
       << "pMutation\t " << parameter->pMutation << endl
       << "tournamentSize\t" << parameter->tournamentSize << endl
       << "correlationThreshold\t" << parameter->correlationThreshold << endl
       << "correlationRange\t" << parameter->correlationRange << endl
       << "regionMinCorrelation\t" <<  parameter->regionMinCorrelation << endl
       << "B\t" << parameter->B << endl << endl
       << "plink_files\t" << parameter->in_files_plink << endl
       << "use_extra_yvalues\t" << parameter->y_value_extra_file << endl  
       << "trait_index\t" << parameter->in_values_int << endl             
       << "multi_forward_step_max\t" << parameter->ms_MaximalSNPsMultiForwardStep << endl << endl;
    poolFile << ssTime.str() << endl;   
    poolFile.flush();
    poolFile.close();
  }
  catch (ofstream::failure e)
  {
    cerr << "Could not write Pool-File" << endl; 
  }
}

/**
 * @brief Computes a heritabilities
 * @param sp_sort - to return an text information about h2_Low and h2_Hight
 * @param diff - a vector of diff values of each model in the pool
 * @param tab - a vector of heritability and posterior of each model in the pool 
 */
void MA::coumputeHeritability(stringstream &sp_sort, const vector<long double> &diff, vector<THeri_PMiY> &tab)
{
  long double h2M_advL = 0.0L, 
              h2M_advH = 0.0L;  
  long double h1 = 0.0,
              h  = 0.5,
              h2 = 1.0;
  double Phi_h1 = 0.0,
         Phi = 0.0,
         Phi_h2 = 0.0;
  for (unsigned int i = 0; i < pool.size(); ++i)
  {
    if ( diff[i] == 0.0 || tab[i].h2M == 0.0)
    {
      Phi_h1 += 1.0 * tab[i].PMi_Y;
    }
    else  
    {
      Phi_h1 += normCDF((h1 - tab[i].h2M) / diff[i]) * tab[i].PMi_Y;
    }
    if ( diff[i] == 0.0 || tab[i].h2M == 0.0)
    {
      Phi += 1.0 * tab[i].PMi_Y;
    }
    else  
    {
      Phi += normCDF((h  - tab[i].h2M) / diff[i]) * tab[i].PMi_Y;
    }
    if ( diff[i] == 0.0 || tab[i].h2M == 0.0)
    {
      Phi_h2 += 1.0  * tab[i].PMi_Y;
    }
    else
    {
      Phi_h2 += normCDF((h2 - tab[i].h2M) / diff[i]) * tab[i].PMi_Y;
    }
  }
  double Phi_05 = Phi;
  while (fabs(Phi - 0.05) > 0.0001 && h1 != h2) // fabs(h1 - h2) > 1.0e-50)
  {
    if (Phi > 0.05)
    {
      h2 = h;
    }
    else
    {
      h1 = h;
    }
    h = (h2 + h1)/2.0;
    Phi = 0.0;
    for (unsigned int i = 0; i < pool.size(); ++i)
    {
      if ( diff[i] == 0.0 || tab[i].h2M == 0.0)
      {
        Phi += 1.0  * tab[i].PMi_Y;
      }
      else
      {
        Phi += normCDF((h - tab[i].h2M) / diff[i]) * tab[i].PMi_Y;
      }
    }  
  }
  h2M_advL = h;
  Phi = Phi_05;
  h = 0.5;
  h1 = 0.0;
  h2 = 1.0;
  while (fabs(Phi - 0.95) > 0.0001  && h1 != h2) // fabs(h1 - h2) > 1.0e-50)
  {
    if (Phi < 0.95)
    {
      h1 = h;
    }
    else
    {
      h2 = h;
    }
    h = (h2 + h1)/2.0;
    Phi = 0.0;
    for (unsigned int i = 0; i < pool.size(); ++i)
    {
      if ( diff[i] == 0.0 || tab[i].h2M == 0.0)
      {
        Phi += 1.0  * tab[i].PMi_Y;
      }
      else
      {
        Phi += normCDF((h - tab[i].h2M) / diff[i]) * tab[i].PMi_Y;
      }
    }  
  }
  h2M_advH = h;
  sp_sort << ", (" << setprecision(3) << h2M_advL << ", " << setprecision(3) << h2M_advH << ")";
}

/**
 * @brief Computes the results and writes them to the output files. 
 */
void MA::computeAndWriteResults() 
{
  long double minMsc, msc;
  vector<long double> sortPool(pool.size());
  int nPool = 0;
  vector<size_t> snps;
  for (set<PoolItem>::iterator it = pool.begin(); it != pool.end(); ++it, ++nPool)
  {   
    snps = it->getPoolItem();                                  // takes snps from a model
    msc = it->getMsc();
    sortPool[nPool] = msc;
  }
  sort(sortPool.begin(), sortPool.end());
  minMsc = sortPool[0];

  // calculate $ \sum _{M \in Models} { \exp {-mBIC(M)} $ Latex
  long double sum = 0.0L;
  for (unsigned  int sP = 0; sP < sortPool.size(); ++sP)
  {
    msc = sortPool[sP];
    sum += exp((minMsc - msc)/2.0L);  
  }

  // calculates $ P(M_i | Y)$
  vector<THeri_PMiY> tab(pool.size());
  unsigned int i;
  long double sum_PMi_Y = 0.0L;
  long double h2M_adv = 0.0L;
  unsigned int n = data.getIdvNo ();
  unsigned int k;
  vector<long double> diff(pool.size());
  
  i = 0;
  for (set<PoolItem>::iterator  it = pool.begin(); it != pool.end(); ++it, ++i)
  {
    msc = (*it).getMsc();
    tab[i].PMi_Y = exp((minMsc - msc)/2.0L) / sum;
    tab[i].h2M = (*it).getHeritability();
    sum_PMi_Y += tab[i].PMi_Y;
    k = (*it).getModelSize();
    diff[i] = sqrt( (4.0 * tab[i].h2M * (1 - tab[i].h2M) * (1 - tab[i].h2M) * (n - k - 1.0) * (n - k - 1.0))  / ((n * n - 1) * (3.0 + n)) );
    h2M_adv += (tab[i].h2M * tab[i].PMi_Y);
  }
  // calculates heritability
  stringstream sp_heritability;
  sp_heritability << "Estimated heritability plus Credible interval: " << setprecision(3) << h2M_adv;
  coumputeHeritability(sp_heritability, diff, tab); 
  
  stringstream ss_Report;
  ss_Report << "MA_MOSGWA used the selection criterion ";
  if (parameter->selectionCriterium == parameter->selectionCriterium_BIC)
  {
    ss_Report << "BIC" << endl;
  }
  if (parameter->selectionCriterium == parameter->selectionCriterium_EBIC)
  {
    ss_Report << "EBIC" << endl;
  }
  if (parameter->selectionCriterium == parameter->selectionCriterium_mBIC)
  {
    ss_Report << "mBIC" << endl;
  }
  if (parameter->selectionCriterium == parameter->selectionCriterium_mBIC2)
  {
    ss_Report << "mBIC2" << endl;
  }
  ss_Report << endl;

  // calculates a posterior of all SNPs $ P(m_i | Y) $
  sort(tab.begin(), tab.end());  
  set<size_t> calculated_j;         // a set of calculated SNPs - to don't calcuate it again
  map<size_t, long double> Pmi_Y;
  multimap<long double, size_t> Pmi_Ysort;
  size_t aSNP;
  
  set<TSNP_Info> set_All_SNPs;
  TSNP_Info aSNP_info;
  i = 0;
  for (set<PoolItem>::iterator it = pool.begin(); it != pool.end(); ++it, ++i)   // for every model in the pool
  {
    snps = (*it).getPoolItem();                                                  // takes snps from a model
    msc = (*it).getMsc();
    for (unsigned int j = 0; j < snps.size(); ++j)                               // for every snp from the model
    {
      aSNP = snps[j];
      if (Pmi_Y.find(aSNP) == Pmi_Y.end())
      {  
        Pmi_Y.insert(pair<size_t, long double>(aSNP, 0.0L)); 
      }  
      Pmi_Y[aSNP] = Pmi_Y[aSNP] + exp((minMsc - msc)/2.0L);
      aSNP_info.SNP = aSNP;
      aSNP_info.Chr = (data.getSNP(aSNP)).getChromosome();
      aSNP_info.pos = (data.getSNP(aSNP)).getBasePairPosition();
      set_All_SNPs.insert(aSNP_info);
    }  
  }
  
  for (map<size_t, long double>::iterator itMap = Pmi_Y.begin(); itMap != Pmi_Y.end(); ++itMap)
  {
    (*itMap).second /= sum;
    Pmi_Ysort.insert(pair<long double, size_t>( (*itMap).second, (*itMap).first));
  }  

  set<size_t> mySnps;    // a recognized SNPs with a posteriori greater than 0.5

  multimap<long double, size_t>::reverse_iterator itM;
  for (itM = Pmi_Ysort.rbegin(); itM != Pmi_Ysort.rend() && (*itM).first > 0.5; itM++)
  {
    mySnps.insert((*itM).second);
  }
 
  set<size_t> sPosterior;
  for (multimap<long double, size_t>::reverse_iterator itM = Pmi_Ysort.rbegin(); itM != Pmi_Ysort.rend() && (*itM).first > 0.5; itM++)
  {
    sPosterior.insert((*itM).second);
  }
  vector<size_t>  vPosterior(sPosterior.begin(), sPosterior.end());  
  Model modelRegion(data);                          // a model for a region strategy
  set<TRegionSet_info> setNewRegion;

  stringstream ss_Detail;
  regionStrategy(ss_Detail, modelRegion, Pmi_Ysort, Pmi_Y, regionMinCorrelation, setNewRegion);
  ss_Report << "The best model found by MA includes the following ";
  if (modelRegion.getModelSize() == 0 || modelRegion.getModelSize() == 1)
  {
    ss_Report << modelRegion.getModelSize() << " SNP:" << endl << endl;
  }
  else
  {
    ss_Report << modelRegion.getModelSize() << " SNPs:" << endl << endl;
  }
  for (int i = 0; i < modelRegion.getModelSize(); ++i)
  {
    ss_Report << modelRegion.getSNPId(i) << endl;
  }
  resultModel = modelRegion;  // set up the result model
  isResult = true;
  ss_Report << endl;
  ss_Report << "Model selection criterion: " << modelRegion.getMSC() << endl << endl;
  ss_Report << sp_heritability.str() << endl << endl;
  ss_Report << "Details on detected SNPs and SNPs in LD: (compact report)" << endl << endl;
  ss_Report << ss_Detail.str() << endl;

  for (set<TRegionSet_info>::iterator itS = setNewRegion.begin(); itS != setNewRegion.end(); ++itS)
  {
    TRegionSet_info s = *itS;
    delete s.s;
  }

  ofstream  reportFile;
  reportFile.exceptions ( ofstream::eofbit | ofstream::failbit | ofstream::badbit ); 
  try
  {
    reportFile.open( ( parameter->out_file_name + "_Report.txt" ).c_str(),  fstream::out | fstream::trunc );
    reportFile << ss_Report.str() << endl;
    reportFile.flush();
    reportFile.close();
  }
  catch (ofstream::failure e)
  {
    cerr << "Could not write '...Report.txt' File" <<endl;
  }  

  stringstream ss_Best_SNPs;
  stringstream ss_All_SNPs;
  ss_All_SNPs << setw(17) << "SNP Id" << "\t" << setw(4) << "Chr" << "\t" << setw(12) << "Pos" << "\t" << setw(12) << "posterior" << "\t" << endl;  
  ss_Best_SNPs << setw(17) << "SNP Id" << "\t" << setw(4) << "Chr" << "\t" << setw(12) << "Pos" << "\t" << setw(12) << "posterior" << "\t" << endl;  
  
  multimap<long double, size_t> mmap_Pmi_SNP;
  for (set<TSNP_Info>::iterator it = set_All_SNPs.begin(); it != set_All_SNPs.end(); ++it)
  {
    ss_All_SNPs << setw(17) << data.getSNP( it->SNP ).getSnpId() << "\t" << setw(4) << (data.getSNP(it->SNP )).getChromosome() << "\t" << setw(12) 
                << (data.getSNP( it->SNP)).getBasePairPosition() << "\t" << setw(12) << Pmi_Y[ it->SNP ] << endl;      
    if (Pmi_Y[ it->SNP ] > regionMinCorrelation)
    {
       mmap_Pmi_SNP.insert( pair<long double, size_t> (Pmi_Y[ it->SNP ], it->SNP) );
    }
  }

  for (multimap<long double, size_t>::iterator it = mmap_Pmi_SNP.begin(); it != mmap_Pmi_SNP.end(); ++it)
  {
    ss_Best_SNPs << setw(17) << data.getSNP( it->second ).getSnpId() << "\t" << setw(4) << (data.getSNP(it->second )).getChromosome() << "\t" << setw(12) 
                 << (data.getSNP( it->second)).getBasePairPosition() << "\t" << setw(12) << Pmi_Y[ it->second ] << endl;      
  }
  
  reportFile.exceptions ( ofstream::eofbit | ofstream::failbit | ofstream::badbit ); 
  try
  {
    reportFile.open( ( parameter->out_file_name + "_Best_SNPs.txt" ).c_str(),  fstream::out | fstream::trunc );
    reportFile << ss_Best_SNPs.str() << endl;
    reportFile.flush();
    reportFile.close();
  }
  catch (ofstream::failure e)
  {
    cerr << "Could not write '...Best_SNPs.txt' File" <<endl;
  }  
  
  reportFile.exceptions ( ofstream::eofbit | ofstream::failbit | ofstream::badbit ); 
  try
  {
    reportFile.open( ( parameter->out_file_name + "_All_SNPs.txt" ).c_str(),  fstream::out | fstream::trunc );
    reportFile << ss_All_SNPs.str() << endl;
    reportFile.flush();
    reportFile.close();
  }
  catch (ofstream::failure e)
  {
    cerr << "Could not write '...All_SNPs.txt' File" <<endl;
  }  
}


/**-------------------------------------------------------------------------
 * @brief Model selection for MA. It works like selectModel() but makes only a stepwise selection.
 * @param model - a reference to a returned model 
 * -------------------------------------------------------------------------
 */
void MA::selectModel ( Model& model ) {
  logger->debug( "MA: Model Selection started: " );
  if (model.computeRegression())
    model.computeMSC();
  else
    model.computeMSCfalseRegression();
  model.makeMultiForwardStep(0, Parameter::selectionCriterium_BIC, NULL, &exclusivedSNP, &goodSNPs );
} 


/** 
 * @brief Returns number of models which have got better msc value than the mscVal parameter
 * @param mscValue - a msc value of the checked model 
 */
int MA::isInNBestModels(double mscValue)
{
  int n = 0;
  for (unsigned int i = 0; i < modelsNo; i++)
    if (models[i]->getMSC() < mscValue)
      ++n; 
  return n;
}

/**
 * @brief selects SNPs for a model based on a posterior of SNP and a correlation
 * @param ss_Region - to return an information 
 * @param model - to return a selected model
 * @param Pmi_Ysort - a multimap containing of sorted posterior values and SNPs
 * @param mapSNPid_Pmi_Y - a multimap containing SNPs and their posterior values
 * @param th - a threshold value of a posterior
 * @param setNewRegion - a returned set containing information about the selected SNPs.
*/
void MA::regionStrategy(stringstream & ss_Region, Model &model, std::multimap<long double, size_t>& Pmi_Ysort, 
                        std::map<size_t, long double> &mapSNPid_Pmi_Y, double th, std::set<TRegionSet_info> &setNewRegion)
{
  map<size_t, unsigned int> mapSNP2ind;   
  vector<vector<double> > tabCorr;
  vector<size_t> SNPs;                    
  vector< set<size_t> > tabRegionSNPs;    
  vector<double> SNPs_posteriori;              
  size_t aSNP;
  double absCorr;
  unsigned int i = 0;
  multimap<long double, size_t>::reverse_iterator itM;
  vector<double> vv;
  vector<size_t> v_all;
  set<size_t>s;
  ofstream  outFile;
  for (itM = Pmi_Ysort.rbegin(); itM != Pmi_Ysort.rend() && (*itM).first > th; ++itM)
  {
    aSNP = itM->second;
    mapSNP2ind.insert(pair<size_t, unsigned int> (aSNP, i));
    tabCorr.push_back(vv);
    SNPs.push_back(aSNP);
    SNPs_posteriori.push_back(mapSNPid_Pmi_Y[aSNP]);
    tabRegionSNPs.push_back(s);
    if (i != 0)
    {
      tabCorr[i].resize(i);
    }
    for (unsigned int k = 0; k < i; k++)
    {
      absCorr = fabs( data.computeCorrelation( SNPs[i], SNPs[k] ) ); // compute correlation between model SNP and SNP j     
      tabCorr[ mapSNP2ind[aSNP] ][k] = absCorr;
    }
    ++i;
  }
  
  // for each SNP which is correlated each other we sum up a posteriori probability
  for (i = 0; i < tabCorr.size(); ++i)
  {
    aSNP = SNPs[i];
    for (unsigned int k = 0; k < i; ++k)
    {
      if (tabCorr[i][k] > 0.5  && (data.getSNP( SNPs[i] )).getChromosome() == (data.getSNP( SNPs[k] )).getChromosome() )  
      {                                                // 3 conditions: 1. correlations > 0.5, 2. the same getChromosome, a distance is smaller than 10000000
        if ( fabs((data.getSNP(SNPs[i])).getBasePairPosition() - (data.getSNP(SNPs[k])).getBasePairPosition()) < 10000000) 
        {
          SNPs_posteriori[i] += mapSNPid_Pmi_Y[ SNPs[k] ];
          SNPs_posteriori[k] += mapSNPid_Pmi_Y[ SNPs[i] ];
          tabRegionSNPs[i].insert(k);
          tabRegionSNPs[k].insert(i);
        }
      }
    }
  }
  vector<double> SNPs_posteriori_orig = SNPs_posteriori;
  
  vector<size_t> v;
  v_all.clear();
  TSNP_Info aSNP_info;
  set<TSNP_Info> aRegion;
  set<TSNP_Info> *aRegionPtr;
  stringstream ssInfo;

  for (i = 0; i < SNPs.size(); ++i)
  {
    if (SNPs_posteriori[i] > 0.5)  // If a posteriori probability of a region is greater than 0.5 we report this region.
    { 
      v.push_back(SNPs[i]);
      for (set<size_t>::iterator it = tabRegionSNPs[i].begin(); it != tabRegionSNPs[i].end(); ++it)
      {
          if ( mapSNPid_Pmi_Y[ SNPs[*it] ] > mapSNPid_Pmi_Y[ v[v.size() -1] ] )
          {
            v[v.size() -1] = SNPs[*it];    // wybór najlepszego SNPa z danego rejonu, żeby reprezentował ten rejon w modelu
          }
      }
    }
  }    
  ///// a compact raport
  SNPs_posteriori = SNPs_posteriori_orig;

  setNewRegion.clear();
  for (i = 0; i < SNPs.size(); ++i)
  {
    if (SNPs_posteriori[i] > 0.5)  // If a posteriori probability of a region is greater than 0.5 we report this region.
    { 
      aRegionPtr = new set<TSNP_Info>;
      aRegion.clear();
      aSNP_info.SNP = SNPs[i];
      aSNP_info.Chr = (data.getSNP(SNPs[i])).getChromosome();
      aSNP_info.pos = (data.getSNP(SNPs[i])).getBasePairPosition();
      aRegion.insert(aSNP_info);
      aRegionPtr->insert(aSNP_info);
      for (set<size_t>::iterator it = tabRegionSNPs[i].begin(); it != tabRegionSNPs[i].end(); ++it)
      {
        aSNP_info.SNP = SNPs[*it];
        aSNP_info.Chr = (data.getSNP(SNPs[*it])).getChromosome();
        aSNP_info.pos = (data.getSNP(SNPs[*it])).getBasePairPosition();
        aRegion.insert(aSNP_info);
        aRegionPtr->insert(aSNP_info);
      }
      TRegionSet_info anItem;
      anItem.s = aRegionPtr;
      setNewRegion.insert(anItem);
      for (set<TSNP_Info>::iterator it = aRegion.begin(); it != aRegion.end(); ++it)
      {
        SNPs_posteriori[ mapSNP2ind[it->SNP] ] = -1;
      }
    }
  }

  ss_Region << setw(15) << "SNP Id" << "\t" << setw(4) << "Chr" << "\t" << setw(12) << "Position" << "\t" << setw(8) << "Corr" << "\t" <<  setw(12) << "posteriori" << "\t" 
                 << setw(12) << "posteriori (region)" << endl << endl;         
  double myCorr;       
  for (set<TRegionSet_info>::iterator itS = setNewRegion.begin(); itS != setNewRegion.end(); ++itS)
  {
    TRegionSet_info s = *itS;
    set<TSNP_Info> *aSNP = s.s;
    
    // looking for the best SNP to set a star (*)
    set<TSNP_Info>::iterator itBest = aSNP->begin();
    size_t bestSNP = itBest->SNP;
    for (set<TSNP_Info>::iterator it = aSNP->begin(); it != aSNP->end(); ++it)
    {
      if (mapSNPid_Pmi_Y[it->SNP] >  mapSNPid_Pmi_Y[bestSNP])
      {
        bestSNP = it->SNP;
      }
    }
    for (set<TSNP_Info>::iterator it = aSNP->begin(); it != aSNP->end(); ++it)
    {
      if (bestSNP == it->SNP)
      {
        myCorr = 1.0;
      }
      else
      {
        if (mapSNP2ind[bestSNP] < mapSNP2ind[ it->SNP])
        {
          myCorr = tabCorr[ mapSNP2ind[it->SNP] ][ mapSNP2ind[ bestSNP] ];
        }
        else
        {
          myCorr = tabCorr[ mapSNP2ind[bestSNP] ][ mapSNP2ind[ it->SNP] ];
        }
      }
      ss_Region << setw(15) << (bestSNP == it->SNP ? (string("*") + data.getSNP( it->SNP ).getSnpId()).c_str() : data.getSNP( it->SNP ).getSnpId()) << "\t" << setw(4) << it->Chr << "\t" 
                << setw(12) << it->pos << "\t" << setw(7) << fixed << setprecision(3) << myCorr << "\t" << setw(12) << fixed << setprecision(3) << mapSNPid_Pmi_Y[ it->SNP ] << "\t" 
                << setw(12) << fixed << setprecision(3) << SNPs_posteriori_orig[ mapSNP2ind[it->SNP] ] << endl;        
    }  
    ss_Region << endl;
  }

  if (v.size() > 0)
  {  
    model.addManySNP(v);
  }
  if (model.computeRegression())
  {
    model.computeMSC();
  }
  else
  {
    model.computeMSCfalseRegression();
  }
}

/**
 * @brief Returns the result model
 */
const Model* MA::result() {
  if (isResult == true)
  {
    return &resultModel;
  }
  else 
  {
    throw Exception("Call result() only aftern run()!");
  }
}

}
