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

#include "GA.hpp"  
#include <algorithm>
#include <iterator>
#include <iomanip>
#include <cstdlib>                     // for random number
#include "SNP.hpp"
#include <set>

const int minModelSize = 3;            // minimal model size for initial population

vector<vector<snp_index_t> > GA::correlations; // correlations vector

/**
 * @brief For printing a given integer vector on the screen
 * @param v a vector to print
 * @returns reference to ostream object
 */ 
ostream &operator<< (ostream &out, vector<snp_index_t>&v)
{
  out << "[";
  vector<snp_index_t>::iterator it = v.begin();
  if (v.size() > 0)
  {
    out << *it;
    it++;
  }
  for (; it != v.end(); it++)
    out << ", " << *it;
  return out << "]";
}

/**
 * @brief For printing a given integer set on the screen
 * @param v a set to print
 * @returns reference to ostream object
 */ 
ostream &operator<< (ostream &out, set<snp_index_t>&v)
{
  out << "{";
  set<snp_index_t>::iterator it = v.begin();
  if (v.size() > 0)
  {
    out << *it;
    it++;
  }
  for (; it != v.end(); it++)
    out << ", " << *it;
  return out << "}";
}



/**
 * @brief Constructor of Genethic Algoritm object
 * @param modelsNo number of models
 * @param maxNoProgressIter number of iteration without progress - stop criterion 
 * @param pCross probabilities of crossover
 * @param pMutation probabilities of mutation
 * @param tournamentSize size of each tournament
 * @param correlationRange snps for correlation are from range [snp - correlationRange, snp + correlationRange]
 * @param correlationThreshold correlation threshold used in local improvement and mutation function
*/
GA::GA(unsigned int modelsNo_, unsigned int maxNoProgressIter_, double pCross_, double pMutation_, unsigned int tournamentSize_, 
       double correlationThreshold_, int correlationRange_)
:modelsNo(modelsNo_), maxNoProgressIter(maxNoProgressIter_), pCross(pCross_), pMutation(pMutation_), 
 tournamentSize(tournamentSize_), correlationThreshold(correlationThreshold_), correlationRange(correlationRange_) 
{
  int firstModelSize;                 // Size of the first model. This model is created with all data (no snps are excluded)
  vector<snp_index_t>::iterator it;
  vector<snp_index_t> snps;

  vector<snp_index_t> v;
  stringstream ss;

  correlations.resize(data.getSnpNo());
  srand(0);
  if ( !parameter.imp_is_imputated)
  {
    cout << "imp_in_imputated" << endl;
    data.imputateMissingValues();  
    //data.writeBEDfile();
  }
  data.calculateIndividualTests();                                                              
  try
  {
    models = new Model*[modelsNo];
  }
  catch (bad_alloc &ba)
  {
    cerr << "Not enought memory for GA(), models number: " << modelsNo << endl;
    exit(-1);
  }
  
  int n;
  int snpsSize = data.getSnpNo();;
  int randSNP;
   
  for (unsigned int i = 0; i < modelsNo; i++)
  {
    try
    {
      models[i] = new Model(data);
    }
    catch (bad_alloc &ba)
    {
      cerr << "Not enought memory for GA(), for model number: " << i << endl;
      exit(-1);   
    }
    selectModel(*models[i]);
    if (i == 0)
    {
      firstModelSize = models[i]->getModelSize();
      if (firstModelSize < minModelSize)
        firstModelSize = minModelSize;  
    }
    if (models[i]->getModelSize() == 0)    // adds to model random snps
    {                                      // TESTING this code haven't been tested
       n = 0; 
       while (n < minModelSize)
       {
         randSNP = random() % snpsSize;
         snps = models[i]->getModelSnps();
         it = find(snps.begin(), snps.end(), randSNP);
         if (it == snps.end())            // random snp does not exsist in the model
         {
           models[i]->addSNPtoModel(randSNP);  
           ++n;
         }
       }
         
    }
    if (models[i]->computeRegression() == false)
      models[i]->computeMSCfalseRegression();
    else
      models[i]->computeMSC();
    v = models[i]->getModelSnps();
    PoolItem p(v, models[i]->getMSC());
    pool.insert(p);
    ss << endl << (i + 1) << "). " << p;
   
  }
  v.clear();
  printLOG("GENETICS ALGORITHM - INITIAL POPULATION:");
  printLOG(ss.str());
}

/**
 * @brief Mutates given model. The model is modified.
 * @param aModel pointer to a model
 * @param pMutation a probability of mutation
 * @param threshold - threshold for correlated snps. If threshold <=0 then this method does not take into account correlated snps
 * WARNING Does not calculete MSC
 */
void GA::mutation(Model *model, double pMutation, double threshold)
{
  if (random() < pMutation * RAND_MAX)
  {  
    int snpsNo = data.getSnpNo();
    vector<snp_index_t> snpsVector = model->getModelSnps();
    vector<snp_index_t>::iterator it;
    set<snp_index_t> snpsSet(snpsVector.begin(), snpsVector.end());
    set<snp_index_t>::iterator itSet;
    vector<snp_index_t> v;
    int oneSnp;
    int maxBadRegression = 1000;
    if (random() % 2)
    {  // add gen
      if (threshold > 0)
      {
        for (it = snpsVector.begin(); it != snpsVector.end(); ++it)
        {
           if (correlations[*it].size() == 0)
             correlations[*it] = stronglyCorrelatedSnps(model, *it, threshold, correlationRange);
           snpsSet.insert(correlations[*it].begin(), correlations[*it].end());
        }  
      }
      
      int noBadRegression = 0;
      do
      {
        do
        {
          oneSnp = random() % snpsNo;
          itSet = snpsSet.find(oneSnp);
        }
        while (itSet != snpsSet.end());
        model->addSNPtoModel(oneSnp);
/*        if (model->computeRegression())
          noBadRegression = -1;
        else
          ++noBadRegression;
*/        
      }  
      while (noBadRegression < 0 || noBadRegression > maxBadRegression);  
      if (noBadRegression > maxBadRegression)
      {
        model->removeSNPValFromModel(oneSnp);
//        model->computeRegression();
//        cerr << "Can not make mutation (add gen) for " << *model << endl;
//        char c; cout << "Press a key.. "; cin >> c;
      }
      if (model->computeRegression())
        model->computeMSC();
      else
        model->computeMSCfalseRegression();
      v = model->getModelSnps();
      PoolItem p(v, model->getMSC());
      pool.insert(p);
    }
    else
    {  
      if (model->getModelSize() > 1)
      { // erase gen
        oneSnp = random() % model->getModelSize();
        model->removeSNPfromModel(oneSnp);
        if (model->computeRegression() == false)
          model->computeMSCfalseRegression();
        else 
          model->computeMSC();
        v = model->getModelSnps();
        PoolItem p(v, model->getMSC());
        pool.insert(p);
      }
      else
      { // change gen to anoter (if model has got only one snp, then this snp is changed to another.
        model->removeSNPfromModel(0);
        oneSnp = random() % snpsNo;
        model->addSNPtoModel(oneSnp);
        if (model->computeRegression() == false)
          model->computeMSCfalseRegression();
        else 
          model->computeMSC();
        v = model->getModelSnps();
        PoolItem p(v, model->getMSC());
        pool.insert(p);
      }
    }
  }
}

/**
 * @brief Makes a child model from two parents models. MSC is calculated
 * @param s1 - parent model
 * @param s2 - parent model
 * @returns pointer to child Model object
 * WARNING Allocates memory for new Model. 
 */                                                               
Model* GA::recombination(const Model & s1, const Model & s2) 
{
  vector<snp_index_t> v1 = s1.getModelSnps();       // get vector of snp's
  vector<snp_index_t> v2 = s2.getModelSnps();
    
  set<snp_index_t> s_1(v1.begin(), v1.end());       // set of snp's
  set<snp_index_t> s_2(v2.begin(), v2.end());
  
  set<snp_index_t> SI;                              // intersection of s1 and s2
  set_intersection(s_1.begin(), s_1.end(), s_2.begin(), s_2.end(), insert_iterator<set<snp_index_t> >(SI, SI.begin()));
  set<snp_index_t> SU;                              // union of s1 and s2
  set_union(s_1.begin(), s_1.end(), s_2.begin(), s_2.end(), insert_iterator<set<snp_index_t> >(SU, SU.begin()));
  set<snp_index_t> SD = SU;                         // symmetric difference of s1 and s2
  set_difference(SU.begin(), SU.end(), SI.begin(), SI.end(), insert_iterator<set<snp_index_t> >(SD, SD.begin()));
  //cout << "SI = " << SI << endl;
  //cout << "SU = " << SU << endl;
   // intersection and forward selection with symmeric difference
  Model* modelForward = new Model(data);

  modelForward->createFromSNPs(SI);
  if (modelForward->computeRegression())
    modelForward->computeMSC();
  else
    modelForward->computeMSCfalseRegression();
  double bestMSC = modelForward->getMSC();
  set<snp_index_t>::iterator it;
  int addedSNP;
  set<snp_index_t> tempSD = SD;
  vector<snp_index_t> v;
  //cout << endl << endl << "forward step: " << *modelForward << endl;
  //cout << "bestMSC = " << bestMSC << endl;
  do
  {
    addedSNP = -1;
    for (it = tempSD.begin(); it != tempSD.end(); ++it)
    {
      modelForward->addSNPtoModel(*it);
      if (modelForward->computeRegression())
        modelForward->computeMSC();
      else
      {
        modelForward->computeMSCfalseRegression();
      }  
      v = modelForward->getModelSnps();
      //cout << "addedSNP: " << *it << ", model = "  << *modelForward <<  endl;
      //char c; cout << "press a key..."; cin >> c;
      PoolItem p(v, modelForward->getMSC());
      pool.insert(p);
      if (modelForward->getMSC() < bestMSC)
      {
        addedSNP = *it;
        bestMSC = modelForward->getMSC();
        //cout << "new best, addedSNP: " << addedSNP << ", model: " << *modelForward  << endl;
        //char c; cout << "press a key..."; cin >> c;
      }
      modelForward->removeSNPValFromModel(*it);
    }
    if (addedSNP >= 0)
    {
      tempSD.erase(addedSNP);
      modelForward->addSNPtoModel(addedSNP);
      
      if (modelForward->computeRegression())
        modelForward->computeMSC();
      else
        modelForward->computeMSCfalseRegression();
      //cout << " new forward model: " << *modelForward << endl;
      //char c; cout << "press a key..."; cin >> c;
    }
  }
  while (addedSNP >= 0);

   if (SD.size() > data.getIdvNo() / 2)           ////!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! only Individual/2.0
    return modelForward;
  
  
  // union and backward selection form symetric difference     
  Model* modelBackward = new Model(data); 
  modelBackward->createFromSNPs(SU);       
  if (modelBackward->computeRegression())
    modelBackward->computeMSC();
  else
  {
    modelBackward->computeMSCfalseRegression();
    //cout << "computeMSCfalseRegression : initial modelBackward" << endl;
      //  char c; cout << "press a key..."; cin >> c;
  }  
 
  bestMSC = modelBackward->getMSC();
  int erasedSNP;
  tempSD = SD;
  //cout << endl << endl << "backward step: " << endl;
  //cout << "bestMSC = " << bestMSC << ", model: " << *modelBackward << endl;
  do
  {
    erasedSNP = -1;
    for (it = tempSD.begin(); it != tempSD.end(); ++it)
    {
      modelBackward->removeSNPValFromModel(*it);
      if (modelBackward->computeRegression())
        modelBackward->computeMSC();
      else
      {
        modelBackward->computeMSCfalseRegression();
//        cout << "computeMSCfalseRegression : initial modelBackward" << endl;
  //      char c; cout << "press a key..."; cin >> c;
      } 
      v = modelBackward->getModelSnps();
      PoolItem p(v, modelBackward->getMSC());
      pool.insert(p);
      //cout << "erasedSNP: " << *it << ", model = " << *modelBackward << endl;
      //char c; cout << "press a key..."; cin >> c;
      if (modelBackward->getMSC() < bestMSC)
      {
        bestMSC = modelBackward->getMSC();
        erasedSNP = *it;
//        cout << "newBest, erasedSNP: " << erasedSNP  <<  ", model = " << *modelBackward << endl;
  //      char c; cout << "press a key..."; cin >> c;
      }
      modelBackward->addSNPtoModel(*it);
    }
    if (erasedSNP >= 0)
    {
      tempSD.erase(erasedSNP);
      modelBackward->removeSNPValFromModel(erasedSNP);
      if (modelBackward->computeRegression())
        modelBackward->computeMSC();
      else
        modelBackward->computeMSCfalseRegression();

    //  cout << "new Bakward model: erasedSNP: " << erasedSNP  << *modelBackward << endl;
      //char c; cout << "press a key..."; cin >> c;
    }
  }
  while (erasedSNP >= 0);  

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
 * @brief Finds index to the worst model in the population
 * @returns index to the worst model in the population
 */
unsigned int GA::findTheWorst() const
{
  unsigned int theWorst = 0;
  for (unsigned int i = 1; i < modelsNo; i++)
    if (models[i]->getMSC() > models[theWorst]->getMSC())
      theWorst  = i;
  return theWorst;  
}

/**
 * @brief Makes tournament selection.
 * @param tournamentSize - size of tournament
 * @returns pointer to winner model.
 */ 
Model* GA::tournamentSelection(unsigned int tournamentSize) const
{
   if (tournamentSize <= 0)
   {
     cerr << "tournament size must be greater than 0" << endl;
     exit(-1);
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
 * @brief Sorting function for Pool 
 * @returns boolean value of comparison a < b
 * @note Takes into acount: first  - size of strings
 *                          second - alphabetical order
 */
bool sort_string(string a, string b)
{
  if (a.length() == b.length())
     return a < b;
  if (a.length() < b.length())
    return true;
  else 
    return false;    
}

/**
 * @brief Writes pool on the screen
 */
void GA::writePool() const
{
  cout << "-------------------------------------------------------------------------" << endl;
  cout << "Pool:" << endl;;
  set<PoolItem>::iterator poolIt;
  int i = 1;
  for (poolIt = pool.begin(); poolIt != pool.end(); poolIt++, ++i)
  {
    PoolItem p1 = *poolIt;
    cout << i << "]" << p1 << ", size: " << p1.getModelSize() << endl << endl;
  }
}

/**
 * @brief Runs the genetics algoritm
 */
void GA::run()
{
  unsigned int GAcount = 0;
  unsigned int iterCount = 0;
  Model *firstParent = 0,                  // parents for recombination 
        *secondParent = 0;
  double randVal;
  Model *childModel = 0;
  vector<snp_index_t> v;
  ofstream	poolFile;
  poolFile.exceptions ( ofstream::eofbit | ofstream::failbit | ofstream::badbit ); // checks if poolfile can be written
  try
  {
    poolFile.open( ( parameter.out_file_name + "_pool.txt" ).c_str(),  fstream::out );
    for (unsigned int i = 0; i < modelsNo; ++i)
      poolFile << (i + 1) << "] " << *models[i] << endl;
    poolFile.flush();
  }
  catch (ofstream::failure e)
  {
    cerr << "Could not write Pool-File" <<endl;
    return;
  }

  while (GAcount < maxNoProgressIter)
  {
    cout << "\riter count = " << setw(5) << iterCount 
         << ", GA count = " << setw(4) << GAcount << endl;
    if (iterCount == 19)
    {
      char c; cout << "wciśnij 19"; cin >> c;
    }

    firstParent = tournamentSelection(tournamentSize);
    randVal = random();    
    if (randVal < pCross * RAND_MAX) 
    { // recombination
      do 
      {
        secondParent = tournamentSelection(tournamentSize);
      } while (firstParent == secondParent);
      childModel = recombination(*firstParent, *secondParent); // allocates memory   
      //cout << "after recomb.: " << childModel << " <> " << *childModel << endl;
    }
    else
    { // only mutation
      childModel = new Model(data);
      *childModel = *firstParent;
    }
    if (iterCount == 19)
      cout << "before mutation: " << childModel << " <> " << *childModel << endl;
    mutation(childModel, pMutation, correlationThreshold);
    if (childModel->computeRegression())
      childModel->computeMSC();
    else
      childModel->computeMSCfalseRegression();
    if (iterCount == 19)
      cout << "after mutation: " << childModel << " <> " << *childModel << endl;
    localImprovement(childModel, correlationThreshold, correlationRange);
    if (iterCount == 19)
      cout << "after local im: " << childModel << " <> " << *childModel << endl;
    //char ch; cout << " press a key ";    cin >> ch;
    unsigned int theWorst = findTheWorst();
    if (childModel->getMSC() < models[theWorst]->getMSC())
    { // we have better model then the worst model in the population
      unsigned int toChange = 0;
      // We don't want to have two or more the same models in population. Those code prevent it.
      while (toChange < modelsNo && childModel->getMSC() != models[toChange]->getMSC())
         ++toChange;
      if (toChange == modelsNo && childModel->getModelSize() > 0)
      { 
        cout << endl << "We have got a better model " << endl;
        delete models[theWorst];
        models[theWorst] = childModel;
        v = childModel->getModelSnps();
        PoolItem p(v, childModel->getMSC());
        if (pool.find(p) == pool.end())
        {
          p = PoolItem(v, childModel->getMSC());
          unsigned int poolSize = pool.size();
          pool.insert(p);    
          if (pool.size() > poolSize)
          { // the child model was added to the pool
            
          }
          poolFile << (pool.size() + 1) << "] " << *childModel << endl;
          poolFile.flush();
        }
        childModel = 0;      
        stringstream ss;
        ss << "iter count: " << iterCount << ":" << endl;
        for (unsigned int m = 0; m < modelsNo; ++m)
          ss << m << ") " << *models[m] << endl;
        printLOG(ss.str());
        GAcount = 0;
      }
      else
      {
        delete childModel;
        ++GAcount;
        //cout << "childModel is in population!!!!!" << endl;
      }
    } 
    else
    {
      ++GAcount;
      delete childModel;
    }  
    ++iterCount;
    if (iterCount % 100 == 0)
    {
      double best, adv, worst;
      statistics(best, adv, worst);
      stringstream ss;
      cout << endl;
      ss << "After " << iterCount << " generations, " << "the worst: " << worst << ", adv: " << adv << ", the best: " << best << endl;
      printLOG(ss.str());
    }
  }
  double best, adv, worst;
  statistics(best, adv, worst);
  stringstream ss;
  cout << endl;
  ss << "At the end of GA, after " << iterCount << " generations, " << "the worst: " << worst << ", adv: " << adv << ", the best: " << best << endl;
  printLOG(ss.str());
  writePoolToLog();
  try 
  {
    poolFile.close();
  }
  catch (ofstream::failure e)
  {
    cerr << "Could not close Pool-File" <<endl;
  }
  
  computePosteriorProbabilitys();
  
}

/**
 * @brief Finds the best, the worst and average value of mBIc2 
 * @param theBest reference to variable of the best msc
 * @param average reference to variable of the average msc
 * @param theWorst reference to variable of the worst msc
 * Statistics are returned with reference variables
 */
void GA::statistics(double &theBest, double &average, double &theWorst) 
{
  double sum = 0.0;
  theBest = theWorst = sum = models[0]->getMSC();
  double msc;
  for (unsigned int i = 1; i < modelsNo; ++i)
  {
    msc = models[i]->getMSC();
    if (msc < theBest)
      theBest = msc;
    else
      if (msc > theWorst)
        theWorst = msc;
    sum += msc;  
  }
  average = sum / modelsNo;
}

/**
 * @brief Destructor
 */
GA::~GA()
{
   if (models != 0)
   {
      cerr << "in destructor" << endl;
      for (unsigned int i = 0; i < modelsNo; ++i)
        if (models[i])
          delete models[i];
      delete []models;    
      models = 0;
   }
}

/**
 *  @brief For testing
 */
void GA::testTournament()
{
   Model *m = new Model(data);
   m = tournamentSelection(tournamentSize);
   //cout << "after tournamenSelection" << endl;
   vector<snp_index_t> v = m->getModelSnps();
   cout << "model from tournament" << v << endl;
   delete m;
}

/**
 *  @brief For testing recombination
 */
void GA::testRecombination()
{
  /*cout << "----------------------------------------------------------------" << endl;
  cout << "Pool:" << endl;;
  set<PoolItem>::iterator poolIt;
  int i = 1;
  for (poolIt = pool.begin(); poolIt != pool.end(); poolIt++, ++i)
  {
    PoolItem p1 = *poolIt;
    cout << i << "]\t" << p1 << ",\t size: " << p1.getModelSize() << endl << endl;
  }
  */
  /*
  //set<snp_index_t>snps;
  int ind_1 = 0,
      ind_2 = 1;
  if (models[ind_1]->computeRegression())
    models[ind_1]->computeMSC();
  else 
    models[ind_1]->computeMSCfalseRegression();
  cout << endl;
  cout << "model[" << ind_1 << "]: " << *models[ind_1] << endl << endl;
  
  
  if (models[ind_2]->computeRegression())
    models[ind_2]->computeMSC();
  else 
    models[ind_2]->computeMSCfalseRegression();
  cout << "model[" << ind_2 << "]: " << *models[ind_2] << endl << endl;
  
  Model *child = recombination(*models[ind_1], *models[ind_2]) ;
  
  if (child->computeRegression())
    child->computeMSC();
  else 
    child->computeMSCfalseRegression();
  cout << "after recombination: " << *child << endl;
  
  return;
  */
  snp_index_t temp_1[] = //{1291, 555, 483, 813, 891, 675, 42, 501, 814, 735, 1160, 1009, 625, 807, 71, 597, 486, 277, 1004, 957, 89, 228, 225, 541, 1197, 411, 351, 172, 151, 192, 1149, 363, 1152, 644, 1059, 713, 630, 305, 765, 317, 1273, 883};
  {1291, 42, 71, 89, 151, 172, 192, 225, 228, 277, 305, 317, 351, 363, 411, 483, 486, 501, 541, 555, 597, 625, 630, 644, 675, 713, 735, 765, 807, 813, 814, 883, 891, 957, 1004, 1009, 1059, 1149, 1152, 1160, 1197, 1273};
  //{555, 228, 483, 813, 891, 675, 42, 501, 814, 735, 1160, 1009, 625, 807, 71, 597, 486, 277, 1004, 957, 89, 225, 541, 1197, 411, 351, 172, 151, 192, 1149, 363, 1152, 644, 1059, 713, 630, 305, 765, 317, 1273, 1291, 883}; // błąd w localImprovement
  
  int sizeTemp_1 = sizeof(temp_1)/sizeof(temp_1[0]);
  vector<snp_index_t> snps;
  snps.assign(temp_1, temp_1 + sizeTemp_1);
  Model m(data);
  m.addManySNP ( snps );
  if (  m.computeRegression() == true )
    m.computeMSC();
  else
  {
    m.computeMSCfalseRegression();
    cout << "false:" << endl;
  }  
  
  cout << endl << "before localImprovement M: " << m << endl;
  char c; cout << "Press key ..."; cin >> c;
  localImprovement(&m, correlationThreshold, correlationRange);
  cout << endl << "a localImprovement M: " << m << endl;
  
  /*
  cout << "Models for recombinations: " << endl;
  cout << "model" << model0 << ": " << models[model0]->getMSC() << endl;
  cout << m1 << endl;
  cout << "model" << model1 << ": " << models[model1]->getMSC() << endl;
  cout << m2 << endl;
  Model *child = recombination(*models[model0], *models[model1]);  
  m1 = child->getModelSnps();
  cout << "child model" << m1 << endl
       << "msc: " << child->getMSC() << endl;
   */    
}

/**
 *  @brief For testing mutation
 */
void GA::testMutation()
{
  vector<snp_index_t> m1;
  for (unsigned i = 0; i < modelsNo; i++)
  {
    m1= models[i]->getModelSnps();
    cout << "Models before mutation: " << endl
         << m1 << ", msc: " << models[i]->getMSC() << endl;
    mutation(models[i], 1.0, correlationThreshold);
    m1 = models[i]->getModelSnps();
    if (models[i]->computeRegression())
      models[i]->computeMSC();
    else
      models[i]->computeMSCfalseRegression();
    cout << "Models after mutation: " << endl
         << m1 << ", msc: " << models[i]->getMSC() << endl; 
  }
}

/**
 * @brief Computes correlation for a given snp with given threshold value
 * @param model - model to compute correlation
 * @param snp - snp at vector modelSnps_
 * @param threshold - threshold value of correlation. 
 * @return vector of snps with a correlation above threshold
 */
vector<snp_index_t> GA::stronglyCorrelatedSnps(Model *model, const int& snp, const double& threshold, int correlationRange )
{
  vector<snp_index_t> snps = model->getModelSnps();
  vector<snp_index_t>::iterator it;
  it = find(snps.begin(), snps.end(), snp);
  //cout << "FIND for snp: " << snp << " ---" << it - snps.begin() << endl;
  if ((unsigned int) (it - snps.begin()) < snps.size())  // snp belongs to the model - computes correlation
    return stronglyCorrelatedSnpsAt(model, it - snps.begin(), threshold, correlationRange);
  else
  {
    snps.clear();  // snp is not in model, no correlated snps, so method returns empty vector
    return snps;
  }
}  

/**
 * @brief Computes correlation for a snp on given position with given threshold value
 * @param model - model to compute correlation
 * @param snpPosition - snpPostition is relative position of the snp at vector modelSnps_
 * @param threshold - threshold value of correlation. 
 * @return vector of snps with a correlation above threshold
 */
vector<snp_index_t> GA::stronglyCorrelatedSnpsAt(Model *model, const int& snpPosition, const double& threshold, int correlationRange )
{
  if (snpPosition < 0 || snpPosition >= model->getModelSize())
  {
    cerr << "ERROR: In stronglyCorrelatedSnps()" << endl
         << "snp: " << snpPosition << " is out of range [0, " << model->getModelSize() << "]" << endl;
    exit(-1); 
  }
  stringstream ss;                 // to save output
  multimap<double, int> StrongCor; // to sort correlated SNPs
  int count;                       // count strongly correlated SNPs for a given SNP
  double abscor;                   // |Correlation| for two SNPs
  count = 0;
  int rangeFrom = model->getSNPat(snpPosition) - correlationRange;
  if (rangeFrom < 0)
    rangeFrom = 0;
  int rangeTo = model->getSNPat(snpPosition) + correlationRange;
  if (rangeTo > data.getSnpNo())
    rangeTo = data.getSnpNo();
  // search for strongly correlated SNPs
  for ( size_t j = rangeFrom; j < rangeTo; j++ ) {
    abscor = fabs( data.computeCorrelation( model->getSNPat(snpPosition), j ) ); // compute correlation between model SNP and SNP j
    if (abscor  >= threshold) // add if correlation is big enough
    {
      StrongCor.insert(pair<double, int>(abscor, j));
      count++;
    }
  }
  vector<snp_index_t> correratedSnps(StrongCor.size());
  // save strongly correlated SNPs for output 
  //ss << model->getSNPId(snpPosition) << " has an absolute correlation of at least "<< threshold  <<" with the following "<< count << " SNPs:\n";
  int i = 0;
  for (multimap<double, int>::reverse_iterator it = StrongCor.rbegin(); it != StrongCor.rend(); it++, ++i)
  {
     /*ss << "\t" << (data.getSNP((*it).second))->getSnpId() 
        << "\t" << (data.getSNP((*it).second))->getChromosome() 
        << "\t" << (data.getSNP((*it).second))->getBasePairPosition() 
        << "\t" << (*it).first << "\n";
     */
    correratedSnps[i] = (*it).second;
  }
  // output to screen
  //cout << ss.str();
  return correratedSnps;
}  

/**
 *  @brief For testing correlation
 */
void GA::testsCorrelation()
{
  int no = 0;
  int oneSnp = 1;
  vector<snp_index_t>corSNP = stronglyCorrelatedSnpsAt(models[no], oneSnp, correlationThreshold, correlationRange);
  stringstream ss;                 // to save output
  ss << "correlatad snps: " << endl
       << corSNP << endl;
  ss << "in test(): " << endl;     
  ss << "snp: " << models[no]->getSNPId(oneSnp) << ", ==, " << (data.getSNP(models[no]->getSNPat(oneSnp)))->getSnpId() << endl;
  ofstream OUT;                    // output model to file
  for (vector<snp_index_t>::iterator it = corSNP.begin(); it != corSNP.end(); it++)     
  {
    ss << "\t" << (data.getSNP(*it))->getSnpId() 
       << "\t" << (data.getSNP(*it))->getChromosome() 
       << "\t" << (data.getSNP(*it))->getBasePairPosition() 
       << endl;
  }
  cout << ss.str() << endl;
  try
  {
    OUT.open( "test_Corr.txt", ios::out);
    OUT << ss.str();
    OUT.close();
  }
  catch( ofstream::failure ) 
  {
    printLOG ("Could not write the correlation file test_Corr.txt");
  }
  //models[no]->printStronglyCorrelatedSnps(0.7);
}



/**
 *  @brief For testing GA
 */
void GA::tests()
{
  set<snp_index_t>snps;
  
  snp_index_t temp_1[] = {42, 71, 89, 151, 172, 192, 225, 228, 277, 305, 317, 351, 363, 411, 483, 486, 501, 541, 555, 597, 625, 630, 644, 675, 713, 735, 765, 807, 813, 814, 883, 891, 957, 1004, 1009, 1059, 1149, 1152, 1160, 1197, 1273, 1291};
  /*{200, 361, 610, 1362, 2476, 3572, 4298, 4371, 4464, 4744, 4758, 5058, 6124, 6390, 6567, 6703, 7769, 9260, 9766, 10604, 11693, 12878, 12929, 13111, 13504, 13976, 14526, 14617, 14821, 14937, 17314, 19204, 19386, 19698, 19889, 20063, 20081, 21348, 22409, 23152, 23605, 24171, 25094, 26599, 26798, 28595, 29741, 30462, 30849, 32951, 35611, 35653, 36170, 39239, 39316, 39874, 40245, 40739, 41305, 41372, 42149, 42203, 42848, 42955, 43702, 44506, 45572, 46045, 47992, 48247, 48260, 49033, 49077, 50853, 50926, 51529, 51606, 52822, 53875, 54557, 56242, 56662, 56831, 57633, 57969, 59173, 60625, 61587, 62179, 62576, 64536, 68174, 69341, 72086, 72310, 73198, 76018, 76744, 76993, 77040, 77438, 78659, 79166, 79782, 80980, 81158, 82323, 82942, 82982, 84187, 86199, 86455, 86768, 88162, 90042, 90088, 90561, 90796, 91501, 91511, 91517, 91636, 92060, 93456, 93835, 94081, 95063, 96671, 96686, 96798, 97212, 97741, 97971, 99127, 100125, 100772, 101117, 102312, 102355, 103316, 104629, 104950, 105676, 107146, 108287, 108866, 109415, 110026, 110039, 111924, 115474, 118036};
   * */
  
  /*{977, 2436, 2798, 3901, 3915, 3917, 3922, 7206, 11635, 12008, 12010, 13852, 15819, 16360, 17181, 17756, 18860, 19610, 19612, 20883, 20886, 22403, 22498, 23523, 23856, 24524, 30868, 34602, 34726, 35385, 37135, 37728, 39424, 40554, 41875, 43053, 48540, 49719, 51170, 51283, 61525, 62532, 64585, 65590, 65596, 66747, 67832, 73235, 73362, 74708, 74714, 74853, 75726, 77894, 79927, 80862, 83062, 83728, 86337, 89785, 91779, 91933, 96524, 99712, 99718, 101886, 103016, 108369, 108905, 109238, 109602, 110968, 113678, 115921, 115924, 116724, 118754, 118758, 118768, 119519, 123274, 124201, 125095, 125466, 128897, 136759, 139537, 141558, 141559, 143212, 144021, 145311, 147902, 150609, 150612, 150621, 150622, 150915, 150991, 151171, 153793, 153795, 154987, 157828, 159554, 160448, 160717, 162130, 163137, 173997, 174006, 174420, 177309, 179211, 183721, 186369, 186372, 187784, 187792, 188629, 188652, 188936, 190321, 191083, 191086, 191225, 193308, 194265, 195142, 197521, 198470, 198616, 199707, 200059, 200103, 204291, 205222, 205501, 207003, 207130, 210207, 213889, 213965, 214200, 216724, 217145, 217149, 218078, 220355, 222561, 222574, 224177};
   */
  snp_index_t temp_2[] = 
  //{1, 10, 20};
  {71, 89, 151, 172, 192, 225, 228, 277, 305, 317, 351, 363, 411, 483, 486, 501, 541, 555, 597, 625, 630, 644, 675, 713, 735, 765, 807, 813, 814, 883, 891, 957, 1004, 1009, 1059, 1149, 1152, 1160, 1197, 1273, 1291, 42};
  /*{3650, 3903, 3926, 18859, 19318, 20886, 23845, 30872, 30992, 31613, 35382, 35517, 37131, 37897, 40560, 41828, 45348, 48525, 51279, 51596, 66738, 70054, 70411, 74710, 74720, 77748, 80863, 83435, 84483, 86341, 91783, 98477, 98781, 105857, 113670, 118752, 140047, 143215, 150611, 150622, 151455, 153243, 153790, 157833, 160131, 162131, 165104, 170179, 174539, 177308, 179091, 187360, 187776, 188632, 188937, 191092, 191100, 193441, 194267, 194424, 195143, 195147, 197518, 204021, 204290, 208737, 209620, 215305, 217151, 217558, 220354, 222553, 223502, 225279, 226673, 227110, 228121};
   */
  int sizeTemp_1 = sizeof(temp_1)/sizeof(temp_1[0]);
  int sizeTemp_2 = sizeof(temp_2)/sizeof(temp_2[0]);
  
  vector<snp_index_t> snpsV;
  snpsV.assign(temp_1, temp_1 + sizeTemp_1);
  vector<snp_index_t> snpsV2;
  snpsV2.assign(temp_2, temp_2 + sizeTemp_2);
  
  Model m1(data);
  Model m2(data);
  
  cout << endl << endl;
  m1.addManySNP(snpsV);
  if (m1.computeRegression())
    m1.computeMSC();
  else
    m1.computeMSCfalseRegression();
  //cout << "m1: " << m1 << endl;
     
  snps.insert(temp_2, temp_2 + sizeTemp_2);
  m2.createFromSNPs(snps);
  if (m2.computeRegression())
    m2.computeMSC();
  else
    m2.computeMSCfalseRegression();
  cout << "m2: " << m2 << endl;
  
  m2.removeSNPfromModel(0);
  if (m2.computeRegression())
  {
    m2.computeMSC();
    cout << "regresion = true" << endl;
  }  
  else
    m2.computeMSCfalseRegression();
  cout << "after remove 1st snp: " << m2 << endl;
  m2.addSNPtoModel(42);
  if (m2.computeRegression())
  {
    m2.computeMSC();
    cout << "regresion = true" << endl;
  }  
  else
    m2.computeMSCfalseRegression();
  cout << "after add snp " << 42 << ": " << m2 << endl;
  
  /*
  set<snp_index_t> linDependent;
  for (int i = 0; i < sizeTemp_1; i++)
  {
    m.addSNPtoModel(temp_1[i]);
    if (m.computeRegression() == false)
    { 
      linDependent.insert(temp_1[i]);
      cerr << "can't coumpute regresion after add snp " << temp_1[i] << endl;
      cerr << "model: " << endl << m << endl;
      m.removeSNPValFromModel(temp_1[i]);
      //exit(-1);
    }
    else
    {
      
      cout << i << "] " << m << endl;
      if (i >= 87 && i <= 89)
      {  
        m.computeMSC(0, true);  
        char c; cout << i << "!!! " << "press a key... "; cin >> c;
      }  
      else
        m.computeMSC();  
      
    }
  }
  if (m.computeRegression() == false)
  {
      cerr << "can't coumpute regresion after add all snps " << endl;
      cerr << "model: " << endl << m << endl;
      exit(-1);
  }
  else
  {
    if (m.computeRegression())
      m.computeMSC();
    else
      m.computeMSCfalseRegression();
    cout << "first model: " << m << endl;
    char c; cout << "press a key.. "; cin >> c;
  }
  if (linDependent.size() != 0)
  {
    set<snp_index_t>tempSet;
    cout << "linear dependent snps: " << linDependent << endl;
    while (linDependent.size() != 0)
    {
      tempSet = linDependent;
      m.clearModel();
      for (set<snp_index_t>::iterator it = tempSet.begin(); it != tempSet.end(); ++it)
      {
        m.addSNPtoModel(*it);
        if (m.computeRegression() == false)
        { 
          cerr << "can't coumpute regresion after add snp " << *it << endl;
          cerr << "model: " << endl << m << endl;
          m.removeSNPValFromModel(*it);
        }
        else
        {
          linDependent.erase(*it);    
        }
      }
      if (m.computeRegression())
        m.computeMSC();
      else
        m.computeMSCfalseRegression();
      cout << "next model: " << m << endl;
      char c; cout << "press a key.. "; cin >> c;
    }
  }
  else
  {
    return;
  }
  char c; cout << "press a key.. "; cin >> c;
  if (m.computeRegression())
    m.computeMSC();
  else
    m.computeMSCfalseRegression();
  cout << "complete model: " << m << endl;
  models[1]->clearModel();
  snps.insert(temp_1, temp_1 + sizeTemp_1);
  //models[1]->createFromSNPs(snps);
  vector<snp_index_t>snpVect;
  snpVect.assign(temp_1, temp_1 + sizeTemp_1);
  models[1]->addManySNP(snpVect);
  models[1]->printModelNew();
  if (models[1]->computeRegression())
    models[1]->computeMSC();
  else
    models[1]->computeMSCfalseRegression();
  cout << "model 1: " << *models[1] << endl;
  
  Model* modelBackward = new Model(data); 
  modelBackward->clearModel();
  modelBackward->createFromSNPs(snps);       
  if (modelBackward->computeRegression() == false)
  {
    cerr << "(modelBackward) can't computeRegression for " << *modelBackward << endl;
    exit(-1);
  }
  modelBackward->computeMSC();
  cout << "MODEL backward: " << *modelBackward << endl;
  models[2]->clearModel();
  snps.clear();
  snps.insert(temp_2, temp_2 + sizeTemp_2);
  models[2]->createFromSNPs(snps);
  if (models[2]->computeRegression())  
  {
    models[2]->computeMSC();
    cout << "model 2: " << *models[2] << endl;
  }  
  else 
  {
    cout << "Could not compute regresion for model:\n\r" << *models[2] << endl;
  }
 
  //testsCorrelation();
  
  
  //v = m.getModelSnps();
  //cout << "empty model: " << snps << endl << "msc: " << msc_ << endl;
  //char c; cout << "wcisnij "; cin >> c;
  
    //localImprovement(models[1], correlationThreshold, correlationRange);
//  cout << " after  impr: " << *models[1] << endl;
 // testMutation();
  */
}

/**
 * @brief Local improvement of a given model
 * @param model - improvement model
 * @param threshold - minimal correlation value
 * @param correlationRange - correlation for snp is compute for snps from range [snp - correlationRange, snp + correlationRange]
 * @note Test all snps which are correlated with tested snp with correlation above threshold.
 *       If any correlated snp has got msc better then tested snp then tested snp is changed to correlated snp.
 *       Test all snps in the model
 */
void GA::localImprovement(Model *model, double threshold, int correlationRange) 
{
  vector <snp_index_t> snps = model->getModelSnps();
  vector <snp_index_t>::iterator it;                        // snp to improvement
  
  vector<snp_index_t> v;
  v = model->getModelSnps();
  vector<snp_index_t> v_pool;
  if (model->computeRegression())
    model->computeMSC();
  else
    model->computeMSCfalseRegression();
  double bestMSC = model->getMSC();
  cout << "model to improvement: " << endl << *model << endl;
  int changedSnp;     
  vector<snp_index_t>::iterator itCorrelation;
  for (it = snps.begin(); it != snps.end(); ++it)  // for every snp in the model
  {                                                // tests snp in it
    if (correlations[*it].size() == 0)
      correlations[*it] = stronglyCorrelatedSnps(model, *it, threshold, correlationRange);
    cout << "correlation for snp: " << *it << endl;
    cout << correlations[*it] << ": (" << correlations[*it].size() << ")" <<endl;
    changedSnp = -1;
    for (itCorrelation = correlations[*it].begin(); itCorrelation != correlations[*it].end(); ++itCorrelation)
    {
      model->removeSNPValFromModel(*it);            // removes snp from model           
      model->addSNPtoModel(*itCorrelation);         // adds correlated snp to model    
      if (model->computeRegression())
        model->computeMSC();
      else
        model->computeMSCfalseRegression();                          // calculates MSC
      cout << "removed: " << *it << ", added: " << *itCorrelation << " model: " << *model << endl;
      v_pool = model->getModelSnps();
      PoolItem p(v_pool, model->getMSC());
      pool.insert(p);
      cout << "snp: " << *itCorrelation << ", msc: " << model->getMSC() << endl;
      if (model->getMSC() < bestMSC)
      {
        bestMSC = model->getMSC();
        changedSnp = *itCorrelation;
        cout << "=> changed snp: " << *it << " with snp: " << *itCorrelation << endl;
        cout << "=> best model: " << *model << endl;
        
      }
      model->removeSNPValFromModel(*itCorrelation);
      model->addSNPtoModel(*it);
      /*
      if ((itCorrelation - correlations[*it].begin()) % 40 == 0)
      {
        char a;
        cout << "press key and Enter ";
        cin >> a;
      }
      */
    }
    if (changedSnp >= 0)
    {
      model->removeSNPValFromModel(*it);
      cout << "remove snp: " << *it << endl;
      cout << "add snp: " << changedSnp << endl;
      model->addSNPtoModel(changedSnp); 
      if (model->computeRegression())
        model->computeMSC();
      else
        model->computeMSCfalseRegression();
      
      cout << "<==> best model: " << *model << endl;
    }
    /*
    char c;
    cout << "for next snp press any key..";
    cin >> c;
    */
  }
  if (model->computeRegression())
    model->computeMSC();
  else
    model->computeMSCfalseRegression();
  cout << "at the end of localImprovement, m = " << *model << endl;
}

/**
 * @brief Writes the pool of models to log file
 */
void GA::writePoolToLog() const
{
  stringstream ss; 
  ss << "-------------------------------------------------------------------------" << endl;
  ss << "                     --                                    Pool:                            --" << endl;
  ss << "                     -------------------------------------------------------------------------" << endl;
  set<PoolItem>::iterator poolIt;
  int i = 1;
  ofstream	poolFile;
  poolFile.exceptions ( ofstream::eofbit | ofstream::failbit | ofstream::badbit ); // checks if poolfile can be written
  try
  {
    poolFile.open( ( parameter.out_file_name + "_pool.txt" ).c_str(),  fstream::out | fstream::trunc );
    for (poolIt = pool.begin(); poolIt != pool.end(); poolIt++, ++i)
    {
      PoolItem p1 = *poolIt;
      ss << i << ") " << p1 << ", size: " << p1.getModelSize() << endl;
      //poolFile << i << ") " << p1 << ", size: " << p1.getModelSize() << endl;
    }  
    poolFile.flush();
  }
  catch (ofstream::failure e)
  {
    cerr << "Could not write Pool-File" <<endl;
  }
  printLOG(ss.str());
  
  poolFile << ss.str() << endl;
  poolFile.close();
}

/**
 * @brief Computes and writes to file (*_pProb.txt) posterior probalibities of models
 */
void GA::computePosteriorProbabilitys()
{
  set<PoolItem>::iterator it = pool.begin();
  long double minMsc, msc, maxMsc; 
  // find minimal msc
  if (it != pool.end())
  {
    minMsc = (*it).getMsc();
    maxMsc = minMsc;
    ++it;
  }
  else
  {
    cerr << "Pool is empty" << endl;
    exit(-1);
    
  }
  for (  ; it != pool.end(); ++it)
  {
    msc = (*it).getMsc();
    if (msc < minMsc)           
      minMsc = msc;
    if (msc > maxMsc)
      maxMsc = msc;
  }
  char c; cout << "minMSC = " << minMsc << ", maxMsc = " << maxMsc << endl; cin >> c;
  //minMsc /= 2.0;
  //maxMsc /= 2.0;
  // calculate $ \sum _{M \in Models} { \exp {-mBIC(M)} $ Latex
  long double sum = 0.0;
  for (it = pool.begin(); it != pool.end(); ++it)
  {
    msc = (*it).getMsc();
    sum += exp((-msc + minMsc)/2.0);  
  }
  cout << " sumM = " << sum << endl; cin >> c;
  // calculate $ P(M_i | Y)$
  vector<long double> PMi_Y(pool.size());
  unsigned int i = 0;
  long double sum_Pmi_Y = 0.0;
  
  stringstream ss; 
  
  ofstream  PMi_YFile;
  PMi_YFile.exceptions ( ofstream::eofbit | ofstream::failbit | ofstream::badbit ); // checks if poolfile can be written
  try
  {
    PMi_YFile.open( ( parameter.out_file_name + "_PMi_YFile.txt" ).c_str(),  fstream::out | fstream::trunc );
    for (it = pool.begin(); it != pool.end(); ++i, ++it)
    {
      msc = (*it).getMsc();
      PMi_Y[i] = exp((-msc + minMsc)/2.0) / sum;
      ss << "PMi_Y[" << i << "] = " << PMi_Y[i] << endl; 
      sum_Pmi_Y += PMi_Y[i];
    }
    ss << "sum Pmi_Y: " <<  sum_Pmi_Y << endl; cout << "wciśnij "; cin >> c;
    PMi_YFile << ss;
    PMi_YFile.flush();
  }
  catch (ofstream::failure e)
  {
    cerr << "Could not write PMi_Y-File" <<endl;
  }
  printLOG(ss.str());
  
  PMi_YFile << ss.str() << endl;
  PMi_YFile.close();
  
  
  // calculate $ P(m_i | Y) $
  set<snp_index_t> calculated_j;         // set of calculated snps - to don't calcuate it again
  vector<snp_index_t> snps;
  multimap<snp_index_t, long double> Pmi_Y;
  multimap<long double, snp_index_t> Pmi_Ysort;
  long double sum_PM_Y;
  //cout << "pool size: " << pool.size() << endl; 
  //cout << "wciśnij "; cin >> c;
  
  
    
  
  set<PoolItem>::iterator itTemp;
  for (it = pool.begin(); it != pool.end(); ++it)                // for every model in the pool
  {
    snps = (*it).getPoolItem();                                  // takes snps from a model
    sum_PM_Y = 0.0;
    for (unsigned int j = 0; j < snps.size(); ++j)                            // for every snp from the model
    {
      if (calculated_j.find(snps[j]) == calculated_j.end())      // if there is no calculation for the snp
      {
        i = 0;
        for (itTemp = pool.begin(); itTemp != pool.end(); ++itTemp, ++i)
        {
          PoolItem spi = *itTemp;
          if (spi.findSns(snps[j]))
          {
            sum_PM_Y += PMi_Y[i];
          }
        }
        Pmi_Y.insert(pair<int, long double>(snps[j], sum_PM_Y));
        Pmi_Ysort.insert(pair<long double, int>(sum_PM_Y, snps[j]));
        cout << snps[j] << " -> " << sum_PM_Y << endl;
        calculated_j.insert(snps[j]);
      }
    }
  }
  
  stringstream sp;   
  ss << "Porteriori probabilities: " << endl; cout << "Porteriori probabilities: " << endl; cout << "wciśnij "; cin >> c;
  
  // sort Pmi_Y
  
  
  for (multimap<snp_index_t, long double>::iterator itM = Pmi_Y.begin(); itM != Pmi_Y.end(); itM++)
  {
    sp << (*itM).first << ": " << (*itM).second << endl;
    cout << (*itM).first << ": " << (*itM).second << endl;
  }
  
  
  printLOG(sp.str());  
  
  ofstream  Pmi_YFile;
  Pmi_YFile.exceptions ( ofstream::eofbit | ofstream::failbit | ofstream::badbit ); // checks if poolfile can be written
  try
  {
    Pmi_YFile.open( ( parameter.out_file_name + "_PjMi_YFile.txt" ).c_str(),  fstream::out | fstream::trunc );
    Pmi_YFile << sp.str() << endl;; 
    PMi_YFile.flush();
    PMi_YFile.close();
  }
  catch (ofstream::failure e)
  {
    cerr << "Could not write PMi_Y-File" <<endl;
  }
  
  
  stringstream sp_sort;   
  for (multimap<long double, snp_index_t>::iterator itM = Pmi_Ysort.begin(); itM != Pmi_Ysort.end(); itM++)
  {
    sp_sort << (*itM).first << ": " << (*itM).second << endl;
    //cout << (*itM).first << ": " << (*itM).second << endl;
  }
  
  ofstream  Pmi_YsortFile;
  Pmi_YsortFile.exceptions ( ofstream::eofbit | ofstream::failbit | ofstream::badbit ); // checks if poolfile can be written
  try
  {
    Pmi_YsortFile.open( ( parameter.out_file_name + "_PjMi_YsortFile.txt" ).c_str(),  fstream::out | fstream::trunc );
    Pmi_YsortFile << sp_sort.str() << endl;; 
    Pmi_YsortFile.flush();
    Pmi_YsortFile.close();
  }
  catch (ofstream::failure e)
  {
    cerr << "Could not write PMi_Y-sorted File" <<endl;
  }
  
}


/**-------------------------------------------------------------------------
 * @brief Model selection for GA. It works like selectModel() but makes only stepwise selection.
 * @param SelectedModel model for selection
 * -------------------------------------------------------------------------
 */
void GA::selectModel ( Model& model ) {
  printLOG( "GA: Model Selection started: " );
  int PValueBorder = data.getSnpNo() - 1; 
  // determines the SNPs with p-Value < Parameter::ms_MaximalPValueForwardStep
  // this are tested in the Forward Step
  while ( PValueBorder >= 0 && data.getSnp_order_Value(PValueBorder) > parameter.ms_MaximalPValueForwardStep ) {
    --PValueBorder;
  }

  // compute the log-likelihood of the 0 Model
  
  if (model.computeRegression())
    model.computeMSC();
  else
    model.computeMSCfalseRegression();
  int dummy=0;
  int *startIndex=&dummy;
  model.makeMultiForwardStep(0, 1, startIndex, &exclusivedSNP ); // the 0 is the default value for PValueBorder, 1 for BIC


//    model.computeRegression();
  model.printModel("After Multi-Forward AG:");
  model.printModelNew();
  model.printModelInR();
} 

void GA::calculateIndividualTests()
{
  
  printLOG("GA: Start Individual Tests");
  
  Model SingleSNP( data ); // create Model with current MData
  size_t*  SNPList = new size_t[data.getSnpNo()]; // to store information for SortVec snp_order_
  double* TestStat = new double[data.getSnpNo()]; // to store information for SortVec snp_order_
  
  for ( size_t i = 0; i < data.getSnpNo(); ++i )
  {
    SNPList[i] = i; // for sorting, store positon 
    TestStat[i]= SingleSNP.computeSingleRegressorTest(i); // compute p-value of single marker test, and store for sorting
    data.setSingleMarkerTestAt(i, TestStat[i]); // store singel marker test in SNP
  }
  data.fillSnp_order_Vec(data.getSnpNo(), SNPList, TestStat); // sort the SNPs w.r.t ascending p-values in snp_order_
  data.printSnpOrder();  
  delete[] SNPList;
  delete[] TestStat;

}

void GA::printExclusivedSNP()
{
  cout << "{";
  for (set<snp_index_t>::iterator it = exclusivedSNP.begin(); it != exclusivedSNP.end(); ++it)
    if (it == exclusivedSNP.begin())
      cout << *it;
    else  
      cout << ", " << *it;
  cout << "}" << endl;  
}
