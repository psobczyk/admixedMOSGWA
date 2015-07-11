#include <algorithm>
#include <iterator>
#include <iomanip>
#include <cstdlib>                     // for random number
#include "SNP.hpp"
#include <set>
#include <list>
#include <cassert>
#include "GA.hpp"  

const bool printInfo = false;

const int minModelSize = 3;            // minimal model size for initial population

vector<vector<snp_index_t> > GA::correlations; // correlations vector


map<snp_index_t, int> GA::recognizedSNPs_Region;
map<snp_index_t, int> GA::recognizedSNPs_bestGA;
map<snp_index_t, int> GA::recognizedSNPs_mosgwa;
map<snp_index_t, int> GA::recognizedSNPs_posterioriModel;
map<snp_index_t, int> GA::recognizedSNPs_clusterMax;  
map<snp_index_t, int> GA::recognizedSNPs_clusterSum;  
map<snp_index_t, int> GA::recognizedSNPs_piMass;  



map<int, int> GA::mapLabel_count;


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

bool compareModels(const Model *m1, const Model *m2)
{
  return m1->getMSC() < m2->getMSC();
}

/**
 * @brief Returns a tame staple as a string
 */
string timeStemple()
{
  const time_t ltime=time(NULL);          //get current calendar time
  const tm* now = localtime(&ltime);      // struct for the day, year....
  return int2strPadWith( now->tm_hour, 2, '0' ) + ":" + int2strPadWith( now->tm_min, 2, '0' ) 
         + ":" + int2strPadWith( now->tm_sec, 2, '0' );  
}

string sec2time(const time_t &t)
{
  int hour = t / 3600;
  int sec = t - hour * 3600;
  int min = sec / 60;
  sec = sec % 60;
  return int2strPadWith( hour, 2, '0' ) + ":" + int2strPadWith( min, 2, '0' ) 
         + ":" + int2strPadWith( sec, 2, '0' );    
}


string stemp_diff(const tm* t_start, const tm* t_end)
{
  unsigned int start_sec = t_start->tm_sec + t_start->tm_min * 60 + t_start->tm_hour * 3600;
  unsigned int end_sec = t_end->tm_sec + t_end->tm_min * 60 + t_end->tm_hour * 3600;
  unsigned int diff_sec;
  if (start_sec > end_sec)
    diff_sec = 24 * 3600 + end_sec - start_sec;
  else
    diff_sec = end_sec - start_sec;
  int hour = diff_sec / 3600;
  //sec -= hour * 3600;
  int sec = diff_sec % 3600;
  int min = sec / 60;
  sec = sec % 60;
  
  return int2strPadWith( hour, 2, '0' ) + ":" + int2strPadWith( min, 2, '0' ) + ":" + int2strPadWith( sec, 2, '0' );    
}

/**
 * @brief Creates a progress bar. It is helpfull to set the stop criterion
 */
string progressBar(int current, int max, int size, char first = '>', char second = '-')
{
  int percent = current * 100.0 / max + 0.5;
  int pos = size * current / max;
  string bar;
  bar = "";
  for (int i = 0; i < size; ++i)  
    if (i < pos)
      bar += first;
    else
      bar += second;
  if (current > max)
  {
    cerr << "progressBar(), current > max!!" << endl;
    exit(-1);
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

/*
string stemp_diff(const tm* t_end, const tm* t_start)
{
  unsigned int start_sec = t_start->tm_sec + t_start->tm_min * 60 + t_start->tm_hour * 3600;
  unsigned int end_sec = t_end->tm_sec + t_end->tm_min * 60 + t_end->tm_hour * 3600;
  unsigned int diff_sec;
  if (start_sec > end_sec)
    diff_sec = 24 * 3600 + end_sec - start_sec;
  else
    diff_sec = end_sec - start_sec;
  
  return int2strPadWith( now->tm_hour, 2, '0' ) + ":" + int2strPadWith( now->tm_min, 2, '0' ) + ":" + int2strPadWith( now->tm_sec, 2, '0' );    
}
*/
/** @brief For calculatation a difference between two times 
 * @param hour, min, sec, nano - calculated time
 * @param tp - start time
 * @param te - end time
 */
void time_diff(int &hour, int & min, int &sec, int &nano, timespec &ts, timespec &te)
{
  if (ts.tv_nsec > te.tv_nsec)
  {
     nano = 1000000000 + te.tv_nsec - ts.tv_nsec;
     if (te.tv_sec == 0)
       te.tv_sec = 24 * 3600;
     --te.tv_sec;
  }
  else
    nano = te.tv_nsec - ts.tv_nsec;
  sec = te.tv_sec - ts.tv_sec;
  hour = sec / 3600;
  sec = sec % 3600;
  min = sec / 60;
  sec = sec % 60;
  nano = nano / 1000000;
  if (nano % 10 >= 5)
    nano += 10;
  nano = nano / 10;
}

/**
 * @brief Prints a given vector of SNPs on the screen
 * @param v - a vector to print
 * @returns reference to ostream object
 */ 

ostream &operator<< (ostream &out, vector<snp_index_t> &v)
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
 * @brief Constructs the Genethic Algoritm object
 * @param modelsNo - the number of models
 * @param maxNoProgressIter - the number of iterations without progress - stop criterion 
 * @param pCross - probabilities of crossover
 * @param pMutation - probabilities of mutation
 * @param tournamentSize - size of each tournament
 * @param correlationRange - snps for correlation are from range [snp - correlationRange, snp + correlationRange]
 * @param correlationThreshold - correlation threshold which is used in local improvement and mutation function
*/
GA::GA(unsigned int modelsNo_, unsigned int maxNoProgressIter_, double pCross_, double pMutation_, unsigned int tournamentSize_, 
       int B_, string fileName, double correlationThreshold_, int correlationRange_, double regionMinCorrelation, bool statisticsOnly)
:correlationTh(0.5), models(0), modelsNo(modelsNo_), maxNoProgressIter(maxNoProgressIter_), pCross(pCross_), pMutation(pMutation_), 
 tournamentSize(tournamentSize_), correlationThreshold(correlationThreshold_), correlationRange(correlationRange_), regionMinCorrelation(regionMinCorrelation), B(B_), realModel(0)
{
//  char c; cout << "statisticsOnly: " << statisticsOnly << ", press a key.."; cin >> c;
  time_start = time(NULL);
  p_30K = false;
  p_100K = false;
  isCluster = false;
  ofstream file;  
  file.open((parameter.out_file_name + "_Report.txt").c_str(), ios::trunc);
  file.close();

  assert(modelsNo > 0);
  assert((unsigned int) (B) <= modelsNo);

  data.calculateIndividualTests();

//  test_region();
//  exit(0);
/*  vector<snp_index_t> tab;
  tab.push_back(697067); sprawdzenie
  tab.push_back(599750 );
  stronglyCorrelatedSnpsCluster(tab, .5 );
  exit(0);  
*/  
//double absCorr = fabs( data.computeCorrelation( SNPs[i], SNPs[k] ) ); // compute correlation between model SNP and SNP j  
  
/*  data.writeBEDfile();
  cout << "Cacluates only IndividualTests()" << endl; 
  exit(0);
*/  
  
  Model m0(data);
  m0.computeRegression();
  data.setLL0M( m0.getMJC() );
  
  if (printInfo) printLOG("Initial population start");
  
  /*
  int myTab[] = {242029, 302459};
  int myTabSize = sizeof(myTab)/sizeof(myTab[0]);
  vector<snp_index_t> tabSNPs(myTab, myTab + myTabSize);
  stronglyCorrelatedSnpsCluster(tabSNPs, 0.5 );
  cout << "end of Correlated" << endl;
//  exit(0);
  */
  
  map<string, int> SNP_Names_Id;
  for (unsigned int i = 0; i < data.getSnpNo(); ++i ) 
  {
    SNP_Names_Id.insert( pair<string, snp_index_t>(data.getSNP( i ).getSnpId(), i) );
  }

  string clusterFile_POWER( (parameter.in_files_plink + ".clu").c_str() );
  readClusterLabel(mapSNPid_label, clusterFile_POWER, SNP_Names_Id, &mapLabel_PmiY);
  //stronglyCorrelatedSnpsCluster();
//  exit(0);
  if ((unsigned int) (parameter.ms_MaximalSNPsMultiForwardStep) > data.getIdvNo()/4)
    parameter.ms_MaximalSNPsMultiForwardStep = data.getIdvNo()/4;
  cout << "maxInitModelSize: " << parameter.ms_MaximalSNPsMultiForwardStep << ", individuals: " << data.getIdvNo() << endl;
  
  poolFilePart.exceptions ( ofstream::eofbit | ofstream::failbit | ofstream::badbit ); // checks if poolfile can be written
  if (parameter.silent == false)
  {
    try
    {
      if (parameter.silent == false)
        poolFilePart.open( ( parameter.out_file_name + "_" + int2str(parameter.in_values_int) + "_pool_part.txt" ).c_str(),  fstream::out );
    }
    catch (...)
    {
      cerr << "Can't open file: " << (parameter.out_file_name + "_" + int2str(parameter.in_values_int) + "_pool_part.txt").c_str() << endl;
      exit(-1);
    }
  }
  
  correlations.resize(data.getSnpNo());
/* 07.07  if ( !parameter.imp_is_imputated)
  {
    cout << "imp_in_imputated" << endl;
    data.imputateMissingValues();  
    data.writeBEDfile();
  }
*/
  exclusivedSNP.reset();
  goodSNPs.reset();
  goodSNPsNo = 0;
  double p_value_threshold = 0.1;
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
  
/*
  vector<snp_index_t> causalSNPs;             // snps from the real model
  modelReader(parameter.causalModelFilename, causalSNPs);  
  cout << "causal model: " << causalSNPs << endl;
  Model m(data);
  //vector<snp_index_t> vM((realSNPs_), (realSNPs_ + realSNPsize - 1));
  m.addManySNP(causalSNPs);
  if (m.computeRegression())
    m.computeMSC();
  else
  {
    m.computeMSCfalseRegression();
    cerr << "FALSE_regression" << endl;
  }
  
  cout << m << endl;
  m.printModelInMatlab ();
  m.printModelInR ();
//  data.printSelectedSNPsInMatlab(causalSNPs);
  exit(0);
*/
  
  
  
  cout << "goodSNPsNo: " << goodSNPsNo << endl;

  
 /*
  {  
    string tabNameCSDA[] = 
//    {"rs9525262", "rs10048748", "rs10937559", "rs8028606"};                // hmm25278 CSDA
//      { "rs9525181", "rs10490450", "rs6441934", "rs2044109", "rs17455546"}; // hmm26651 CSDA
    {"rs9525262", "rs17386102", "rs1370718", "rs2044109", "rs2819755", "rs13021147"}; //hmm32074 CSDA
//    {"rs9525262", "rs10490450", "rs6441934", "rs2044109", "rs17455546", "rs3761945", "rs17326215"}; // hmm34610
    const int tabNameCSDA_size = (sizeof(tabNameCSDA)/sizeof(tabNameCSDA[0]));
    vector<string> output(tabNameCSDA, tabNameCSDA + tabNameCSDA_size);
    for (unsigned int i = 0; i < output.size(); ++i ) {
      output[i] = tabNameCSDA[i]; //getSNPId(i);
    }
    data.printSelectedSNPsInMatlab(output, "CSDA");
    vector<snp_index_t> csdaSNPs;
    data.findSNPIndex(output, csdaSNPs);
    Model m(data);
    if (csdaSNPs.size() > 0)
      m.addManySNP(csdaSNPs);
    if (m.computeRegression())
      m.computeMSC();
    else
    {
      m.computeMSCfalseRegression();
      cerr << "FALSE_regression" << endl;
    }
    cout << m << endl;
    m.printModel( "csda_hmm25278-S", Parameter::selectionCriterium_mBIC2 );
    
    
//    int tabMA[] = {91647, 534917, 575306, 690497};                         // hmm25278-S
//    int tabMA[] = {91196, 163722, 222895, 328440, 329108, 331097, 534260}; // hmm26651-S
      int tabMA[] = {48142, 65058, 69328, 330235, 389936, 473715, 534523, 690589}; // hmm32074
//    int tabMA[] = {22407, 122680, 138994, 199909, 276045, 330332, 535351, 576665}; // hmm34610-S
    const int tabMA_size = (sizeof(tabMA)/sizeof(tabMA[0]));
    vector<snp_index_t> maSNPs(tabMA, tabMA + tabMA_size);
    Model ma(data);
    if (maSNPs.size() > 0)
      ma.addManySNP(maSNPs);
    if (ma.computeRegression())
      ma.computeMSC();
    else
    {
      ma.computeMSCfalseRegression();
      cerr << "FALSE_regression" << endl;
    }
    cout << "MA: " << ma << endl;
    ma.printModel( "the_best_hmm32074-S", Parameter::selectionCriterium_mBIC2 );
    
    exit(0);
  }
*/
    // zapisanie błędnie obliczonej korelacji dla klastrów  
  /*
  {
    int snps[] = {535149, 91264, 155080, 575381};
    const int snps_size = sizeof(snps)/sizeof(snps[0]);
    
    vector<snp_index_t> v(snps, snps + snps_size);
    Model m(data);
    m.addManySNP(v);
    if (m.computeRegression())
      m.computeMSC();
    else
    {
      m.computeMSCfalseRegression();
      cerr << "FALSE_regression" << endl;
    }
    
    cout << m << endl;
    m.printModelInMatlab ();
   // saveRecognizedSNPinfo(v, "CSDA");
  }
//  exit(0);
*/
  
  
/*  
//  data.writeBEDfile( 0, 0 );
//  data.writeBEDfilePlink();
  cout << "printSNPs ();: " << endl;
  data.printSNPs ();
  cout << " void printGenoData () const;" << endl;
  ofstream genoFile((parameter.out_file_name + "_geno.txt").c_str());
  data.printGenoDataAG (genoFile);
  genoFile.close();
  cout << "OK" << endl;
  exit(0);
*/  
  
  vector<snp_index_t> v;
  m0.computeMSC();
  RSSo = m0.getMJC();
  double h2_M = (RSSo - m0.getMJC()) / RSSo;
  unsigned int pSize = pool.size(); 
  v = m0.getModelSnps();
  PoolItem p(v, m0.getMSC(), h2_M, pSize + 1, 'I');
  pool.insert(p);
  if (pSize < pool.size()  && parameter.silent == false)
  {
    poolFilePart << p << "\n";
  }
  
  try
  {
    models = new Model*[modelsNo];
  }
  catch (bad_alloc &ba)
  {
    cerr << "Not enought memory for GA(), models number: " << modelsNo << endl;
    exit(-1);
  }

  if (statisticsOnly == true)
  { //<TODO>  sprawdzić, czy potrzeba przydzielać pamięć dla models;
    for (unsigned int i = 0; i < modelsNo; ++i)
      models[i] = 0;
    return;
  }  
  // for reading an initial population from a file
  set<PoolItem> population;   
  if (fileName != "")  // reads an initial population to a set of models
    readInitialPop(fileName, population, modelsNo);
  set<PoolItem> ::iterator itPop = population.begin();
  
  vector<snp_index_t>::iterator it;
  vector<snp_index_t> snps;
  cout << "parameter.maximalModelSize: " << parameter.maximalModelSize << endl;
//  {char c; cout << "Pres a key.."; cin >> c;}
  stringstream ss;
  // the first model is created by stepwise method
  
  try
  {
    models[0] = new Model(data);
    if (fileName != "")  // reads models from a file. 
    {
      if (itPop != population.end())
      {
        v.clear();
        v = itPop->getPoolItem();                             
        v_mosgwa = v;   
        if (v.size() > 0)
          models[0]->addManySNP(v);
        ++itPop;
      }
    }
    else                 // stepwise selection, mBIC2
    {
      models[0]->computeRegression();
      data.setLL0M(models[0]->getMJC());
//      cout << "Stepwise setLLOM: " << models[0]->getMJC() << endl;;
      data.selectModel(
        models[0],
        parameter.PValueBorder,                         // PValueBorder=350;
        parameter.maximalModelSize,                     // maximalModelSize
        Parameter::selectionCriterium_mBIC_firstRound   // uses different expected_causal_SNPs
      );
      data.selectModel(
        models[0], 
        5000,
//        350,
//        parameter.PValueBorder, 
        parameter.maximalModelSize, 
        Parameter::selectionCriterium_mBIC2
      ); 
//      cout << "First model (Stepwise): " << *models[0] << endl;
//      char chr; cout << "Press a key..."; cin >> chr;
      vector<snp_index_t> v = models[0]->getModelSnps();
      v_mosgwa = v;
      cout << "MOSGWA: " << v << endl;
      for (unsigned int i = 0; i < v.size(); ++i)
      {
        exclusivedSNP[ v[i] ] = true;
      }
//      {cout << "Press a key.."; char c; cin >> c;}
    }  
  }
  catch (bad_alloc &ba)
  {
    cerr << "Not enought memory to create a model " << 0 << " (in GA())" <<  endl;
    exit(-1);   
  }
  // for the first half of populations
  for (unsigned int i = 1; i < modelsNo/2; i++)
  {
    try
    {
      models[i] = new Model(data);
      if (fileName != "")  // reads models from a file. 
      {
        if (itPop != population.end())
        {
          v.clear();
          v = itPop->getPoolItem(); 
          if (v.size() > 0)
            models[i]->addManySNP(v);
          ++itPop;
        }
      }
      else                 // multi forward selection with BIC     
      {
        selectModel(*models[i]);
//        {cout << i << "> Press a key.."; char c; cin >> c;}  
      }  
    }
    catch (bad_alloc &ba)
    {
      cerr << "Not enought memory to create a model " << i << " (in GA())" <<  endl;
      exit(-1);   
    }
  }
//  cout << "parameter.ms_MaximalSNPsMultiForwardStep: " << parameter.ms_MaximalSNPsMultiForwardStep << endl;
//  char c; cout << "Press a key.."; cin >> c;
  // for the second half of population
  for (unsigned int i = modelsNo/2; i < modelsNo; i++)
  {
    try
    {
      models[i] = new Model(data);
    }
    catch (bad_alloc &ba)
    {
      cerr << "Not enought memory to create a model " << i << " (in GA())" <<  endl;
      exit(-1);   
    }
    if (fileName != "")  // reads model from a file. 
    {
      if (itPop != population.end())
      {
        v.clear();
        v = itPop->getPoolItem();  
        if (v.size() > 0)
          models[i]->addManySNP(v);
        ++itPop;
      }
    }
    else                 // creates a random model
    {
      set<snp_index_t> initSNP;
      int randSNP;
      unsigned int randSNPno = 5; //models[0]->getModelSize();
/*      if (randSNPno > 2)
        randSNPno = models[0]->getModelSize();
      if (randSNPno < 2)
        randSNPno = 50;   // Pozostają drobne niuanse do rozwiązania, np. warunek, żeby liczba SNPów w populacji była mniejsza od liczby goodSNPsNo
*/        
      while ( initSNP.size() < randSNPno )
      { 
        randSNP = rand() % goodSNPsNo;
        if (exclusivedSNP[randSNP] == true)
          continue;
        initSNP.insert(randSNP);
        exclusivedSNP[randSNP] = true;
      }
      models[i]->createFromSNPs(initSNP);
    }
    
  }
  // writes the models to the pool file
  cout << endl;
  for (unsigned int i = 0; i < modelsNo; i++)
  {
    if (models[i]->computeRegression() == false)
      models[i]->computeMSCfalseRegression();
    else
      models[i]->computeMSC();
    
    double h2_M = (RSSo - models[i]->getMJC()) / RSSo;
    unsigned int pSize = pool.size(); 
    v = models[i]->getModelSnps();
    PoolItem p(v, models[i]->getMSC(), h2_M, pSize + 1, 'I');
    if (i == 0)
      toPool(models[i], 'S');
    else  
      toPool(models[i], 'I');
    ss << endl << (i + 1) << "). " << p;
    cout << i << "] " << *models[i] << "\tsize: " << models[i]->getModelSize() << endl;
  }
  printLOG("GENETICS ALGORITHM - INITIAL POPULATION:");
  printLOG(ss.str());
  if (printInfo) printLOG("Initial population END");
  cout << endl;
}

/**
 * @brief Writes a model to the pool
 * @param model - a model to write
 * @param c - an information where the model was created: I - initial population, R - recombination, L - local improvement, A - add-mutation etc.
 */
void GA::toPool(const Model *model, char c)
{
  vector<snp_index_t> v = model->getModelSnps();
  double h2_M = (RSSo - model->getMJC()) / RSSo;
  unsigned int pSize = pool.size(); 
  PoolItem p(v, model->getMSC(), h2_M, pSize + 1, c);
  pool.insert(p);
  if (pSize < pool.size()  && parameter.silent == false)
  {
    poolFilePart << p << "\n";
  }
}

/**
 * @brief Mutates a given model. The model is modified.
 * @param aModel pointer to a model
 * @param pMutation a probability of mutation
 * @param threshold - a threshold for the correlated snps. If the threshold is less or equal 0 this method does not take into account correlated snps
 */
void GA::mutation(Model *model, double pMutation, double threshold)
{
  int snpsNo = data.getSnpNo();
  vector<snp_index_t> snpsVector = model->getModelSnps();
  vector<snp_index_t>::iterator it;
  set<snp_index_t> snpsSet(snpsVector.begin(), snpsVector.end());
  set<snp_index_t>::iterator itSet;
  vector<snp_index_t> v;
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
    writePoolofSize(p_30K, 30000);
    writePoolofSize(p_100K, 100000);
  }
  else
  {  
    if (model->getModelSize() > 1)
    { // erases a gen from the model
      oneSnp = random() % model->getModelSize();
      model->removeSNPfromModel(oneSnp);
      if (model->computeRegression() == false)
        model->computeMSCfalseRegression();
      else 
        model->computeMSC();
      toPool(model, 'E');
      writePoolofSize(p_30K, 30000);
      writePoolofSize(p_100K, 100000);
    }
    else
    { // changes a gen to anoter (if the model has got only one SNP, then this SNP is changed to another)
      if (model->getModelSize() != 0)
        model->removeSNPfromModel(0);
      oneSnp = rand() % goodSNPsNo;
      model->addSNPtoModel(data.getOrderedSNP(oneSnp));
      if (model->computeRegression() == false)
        model->computeMSCfalseRegression();
      else 
        model->computeMSC();
      toPool(model, 'C');
      writePoolofSize(p_30K, 30000);
      writePoolofSize(p_100K, 100000);
    }
  }
}

/**
 * @brief Makes a child model from two parents models. MSC is calculated
 * @param s1 - parent model
 * @param s2 - parent model
 * @returns pointer to a child model
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
  set<snp_index_t> SD;//                            // symmetric difference of s1 and s2
  set_difference(SU.begin(), SU.end(), SI.begin(), SI.end(), insert_iterator<set<snp_index_t> >(SD, SD.begin()));

  set<snp_index_t>::iterator it;
  int addedSNP;
  set<snp_index_t> tempSD = SD;
  vector<snp_index_t> v;

  // intersection and forward selection with symmeric difference
  Model* modelForward = new Model(data);
  modelForward->createFromSNPs(SI);
  if (modelForward->computeRegression())
    modelForward->computeMSC();
  else
    modelForward->computeMSCfalseRegression();
  double bestMSC = modelForward->getMSC();
  do
  {
    addedSNP = -1;
    for (it = tempSD.begin(); it != tempSD.end(); ++it)
    {
      modelForward->addSNPtoModel(*it);
      if (modelForward->computeRegression())
        modelForward->computeMSC();
      else
        modelForward->computeMSCfalseRegression();
      toPool(modelForward, 'F');
      writePoolofSize(p_30K, 30000);
      writePoolofSize(p_100K, 100000);
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
        modelForward->computeMSC();
      else
        modelForward->computeMSCfalseRegression();
    }
  }
  while (addedSNP >= 0);
  if (modelForward->computeRegression())
    modelForward->computeMSC();
  else
    modelForward->computeMSCfalseRegression();
  
  if (SD.size() > data.getIdvNo() / 2)           // !!! 
    return modelForward;
  
  // union and backward selection form symetric difference     
  Model* modelBackward = new Model(data); 
  modelBackward->createFromSNPs(SU);       
  if (modelBackward->computeRegression())
    modelBackward->computeMSC();
  else
    modelBackward->computeMSCfalseRegression();
  set<snp_index_t> bestBackward = SU;
  set<snp_index_t> currentSNPs = SU;
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
        modelBackward->computeMSC();
      else
        modelBackward->computeMSCfalseRegression();
      toPool(modelBackward, 'B'); 
      writePoolofSize(p_30K, 30000);
      writePoolofSize(p_100K, 100000);
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
        modelBackward->computeMSC();
      else
        modelBackward->computeMSCfalseRegression();
    }
  }
  while (erasedSNP >= 0);  
  
  if (modelBackward->computeRegression())
    modelBackward->computeMSC();
  else
    modelBackward->computeMSCfalseRegression();
  
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
 * @returns pointer to a winner model.
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
 * @brief Sorting function for the Pool 
 * @returns boolean value of comparison a < b
 * @note Takes into acount: first  - size of strings
 *                          second - alphabetical order
 */
bool sort_string(string a, string b)
{
  if (a.length() == b.length())
     return a < b;
  return a.length() < b.length();
}

/**
 * @brief Writes pool on the screen
 */
void GA::writePool() const
{
  cout << "Pool:" << endl;;
  set<PoolItem>::iterator poolIt;
  int i = 1;
  for (poolIt = pool.begin(); poolIt != pool.end(); poolIt++, ++i)
  {
    cout << i << ") " << *poolIt << endl;
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
  Model *childModel = 0;
  vector<snp_index_t> v;
  double best, adv, worst;
  statistics(best, adv, worst);
  stringstream ssBegin;
  cout << endl;
  ssBegin << "Initial population, the worst: " << worst << ", adv: " << adv << ", the best: " << best << endl;
  printLOG(ssBegin.str());
  p_30K = true;
  p_100K = true;


  while (GAcount < maxNoProgressIter)
  {
    if (parameter.silent == false)
    {
      cout << "\riter count = " << setw(5) << iterCount 
           << ", GA count = " << setw(4) << GAcount << ", -> " << setw(5) << GAcount * 100.0 / parameter.maxNoProgressIter << "%... done            ";
      cout.flush();     
    }
    

    firstParent = tournamentSelection(tournamentSize);
    do  
    {
      secondParent = tournamentSelection(tournamentSize);
    } while (firstParent == secondParent);
    if (random() < pCross * RAND_MAX) 
    { // recombination
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
        *childModel = *firstParent;
      else
        *childModel = *secondParent;
      mutation(childModel, pMutation, correlationThreshold);
    }
    if (childModel->computeRegression())
      childModel->computeMSC();
    else
      childModel->computeMSCfalseRegression();
    //localImprovement(*childModel, correlationThreshold, correlationRange);
    old_localImprovement(childModel, correlationThreshold, correlationRange);
    //old_localImprovement(childModel, 0.3, correlationRange);
    unsigned int theWorst = findTheWorst();
    int rankB = isInNBestModels(childModel->getMSC());
    if (childModel->getMSC() < models[theWorst]->getMSC())
    { // we have a better model than the worst model in the population
      unsigned int toChange = 0;
      // We don't want to have two or more the same models in the population. This code prevents it.
      //while (toChange < modelsNo && childModel->getMSC() != models[toChange]->getMSC())
      while (toChange < modelsNo && *childModel != *models[toChange])
         ++toChange;
      if (toChange == modelsNo && childModel->getModelSize() > 0)
      {  
        
/*        if (childModel->getMSC() < 34128.82)
        {
          for (unsigned int m = 0; m < modelsNo; ++m)
            cout << m << ") " << setprecision(15) << *models[m] << "<>" << setprecision(20) << models[m]->getMSC() << endl;
          cout << "child: " << setprecision(15) << *childModel  << "<>" << setprecision(20) << childModel->getMSC() << endl;
          char c; cout << "Press a key..."; cin >> c;
        }
*/        
        
        
        if (parameter.silent == false)
          cout << endl << "We have got a better model " << endl;
        delete models[theWorst];
        models[theWorst] = childModel;
        childModel = 0;      
        if (parameter.silent == false)
        {
          stringstream ss;
          ss << endl << "iter count: " << iterCount << ":" << endl;
          for (unsigned int m = 0; m < modelsNo; ++m)
            ss << m << ") " << *models[m] << endl;
          printLOG(ss.str());
        }  
        if (rankB < B)
          GAcount = 0;
      }
      else
      {
        delete childModel;
        childModel = 0;
        ++GAcount;
      }
    } 
    else
    {
      ++GAcount;
      if (childModel != 0)
        delete childModel;
      childModel = 0;
    }  
    ++iterCount;
    if (iterCount % 100 == 0)
    {
      double best, adv, worst;
      statistics(best, adv, worst);
      stringstream ss;
      if (parameter.silent == false)
        cout << endl;
      else
      {
        cout << timeStemple() << progressBar(GAcount, parameter.maxNoProgressIter, 30) << "iterations: " << setw(5) << iterCount << ", the worst: " << setprecision(10) 
             << worst << ", ave: " << adv << ", the best: " << best << ", poolSize: " << setw(6) << pool.size() << endl;
        for (unsigned int m = 0; m < modelsNo; ++m)
           cout << m << ") " << *models[m] << endl;
             
      }       
      ss << progressBar(GAcount, parameter.maxNoProgressIter, 30) << " iteration: " << setw(5) << iterCount << ", the worst: " << setprecision(10) 
         << worst << ", ave: " << adv << ", the best: " << best << ", poolSize: " << setw(6) << pool.size();// << endl;
      printLOG(ss.str());
    }
  }
  statistics(best, adv, worst);
  stringstream ss;
  if (parameter.silent == false)
    cout << endl;
  ss << "At the end of GA, after " << iterCount << " generations, " << "the worst: " << worst << ", adv: " << adv << ", the best: " << best;// << endl;
  cout << endl << "At the end of GA, after " << iterCount << " generations, " << "the worst: " << worst << ", adv: " << adv << ", the best: " << best << endl;
  printLOG(ss.str());
  if (parameter.silent == false)
     poolFilePart << "END" << endl;
  
  stringstream ssFinal;
//  ssFinal << endl << endl;
  
  for (unsigned int i = 0; i < modelsNo; i++)
  {
    if (models[i]->computeRegression() == false)
      models[i]->computeMSCfalseRegression();
    else
      models[i]->computeMSC();
    ssFinal << endl << (i + 1) << "] " << *models[i] << "\tsize: " << models[i]->getModelSize();
  }
  printLOG("GENETICS ALGORITHM - FINAL POPULATION:");
  printLOG(ssFinal.str());

  
}

/**
 * @brief Finds the best, the worst and the average value of mBIc2 
 * @param theBest reference to variable of the best msc
 * @param average reference to variable of the average msc
 * @param theWorst reference to variable of the worst msc
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
   if (realModel != 0)
   {
     delete realModel;
     realModel = 0;
   }
   if (parameter.silent == false)
     poolFilePart.close();
}

/**
 * @brief Computes correlation for a given SNP and a given threshold value
 * @param model - model to compute correlation
 * @param snp - SNP at vector modelSnps_
 * @param threshold - threshold value of correlation. 
 * @return vector of SNPs with a correlation above threshold
 */
vector<snp_index_t> GA::stronglyCorrelatedSnps(Model *model, const int& snp, const double& threshold, int correlationRange )
{
  vector<snp_index_t> snps = model->getModelSnps();
  vector<snp_index_t>::iterator it;
  it = find(snps.begin(), snps.end(), snp);
  if ((unsigned int) (it - snps.begin()) < snps.size())  // snp belongs to the model - computes correlation
    return stronglyCorrelatedSnpsAt(model, it - snps.begin(), threshold, correlationRange);
  else
  {
    snps.clear();  // a snp is not in the model, no correlated snps, so this method returns the empty vector
    return snps;
  }
}  

/**
 * @brief Computes correlation for a snp on given position with given threshold value
 * @param model - model to compute correlation
 * @param snpPosition - snpPostition is a relative position of a SNP at the model
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
  stringstream ss;                 // to save the output
  multimap<double, int> StrongCor; // to sort the correlated SNPs

  double abscor;                   // abs|Correlation| of two SNPs
  unsigned int j;                           
  
  unsigned int rangeFrom;
  if (model->getSNPat(snpPosition) - correlationRange < 0)
    rangeFrom = 0;
  else
    rangeFrom = model->getSNPat(snpPosition) - correlationRange;
  if (rangeFrom < 0)
    rangeFrom = 0;
  unsigned int rangeTo = model->getSNPat(snpPosition) + correlationRange;
  if (rangeTo > data.getSnpNo())
    rangeTo = data.getSnpNo();
  // search for strongly correlated SNPs
  for (j = rangeFrom; j < rangeTo; j++)
  {
    abscor = fabs( data.computeCorrelation( model->getSNPat(snpPosition), j ) ); // compute correlation between model SNP and SNP j
    if (abscor  >= threshold) // add if correlation is big enough
    {
      StrongCor.insert(pair<double, int>(abscor, j));
    }
  }
  vector<snp_index_t> correratedSnps(StrongCor.size());
  int i = 0;
  for (multimap<double, int>::reverse_iterator it = StrongCor.rbegin(); it != StrongCor.rend(); it++, ++i)
  {
    correratedSnps[i] = (*it).second;
  }
  return correratedSnps;
}  

/**
 * @brief Makes local improvement of a given model
 * @param model - a model to improve
 * @param threshold - minimal correlation value
 * @param correlationRange - a correlation range. For a given SNP, corelations are computed in range [SNP - correlationRange, SNP + correlationRange]
 */
void GA::old_localImprovement(Model *model, double threshold, int correlationRange) 
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
  int changedSnp;     
  vector<snp_index_t>::iterator itCorrelation;
  for (it = snps.begin(); it != snps.end(); ++it)  // for every snp in the model
  {                                                // tests snp in it
    if (correlations[*it].size() == 0)
      correlations[*it] = stronglyCorrelatedSnps(model, *it, threshold, correlationRange);
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
        model->computeMSC();
      else
        model->computeMSCfalseRegression();                          // calculates MSC
      toPool(model, 'L');  
      writePoolofSize(p_30K, 30000);
      writePoolofSize(p_100K, 100000);
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
        model->computeMSC();
      else
        model->computeMSCfalseRegression();            
    }
  }
  if (model->computeRegression())
    model->computeMSC();
  else
    model->computeMSCfalseRegression();
}

/**
 * @brief Writes the pool of models to a file
 * @param ssTime - the runtime of GA
 */
void GA::writePoolToFile(stringstream &ssTime, string postfix) const
{
  
  
  set<PoolItem>::iterator poolIt;
  int i = 1;
  ofstream poolFile;
  poolFile.exceptions ( ofstream::eofbit | ofstream::failbit | ofstream::badbit ); // checks if poolfile can be written
  try
  {
    poolFile.open( ( parameter.out_file_name + "_pool" + postfix + /*+ int2str(parameter.in_values_int)*/ + ".txt" ).c_str(),  fstream::out | fstream::trunc );
    poolFile << ssTime.str() << endl;
//    cout << "Pool file postfix: " << postfix << endl;
    for (poolIt = pool.begin(); poolIt != pool.end(); poolIt++, ++i)
    {
      PoolItem p1 = *poolIt;
      poolFile << p1 << "\n";
//      ss << p1 << "\n";
    }  
    poolFile.flush();
    poolFile << "END\t" << endl
       << "poolSize\t" << pool.size()  << endl
       << "modelsNo\t" << parameter.modelsNo << endl
       << "maxNoProgressIter\t" << parameter.maxNoProgressIter << endl
       << "pCross\t" << parameter.pCross << endl
       << "pMutation\t " << parameter.pMutation << endl
       << "tournamentSize\t" << parameter.tournamentSize << endl
       << "correlationThreshold\t" << parameter.correlationThreshold << endl
       << "correlationRange\t" << parameter.correlationRange << endl
       << "regionMinCorrelation_ =" <<  parameter.regionMinCorrelation << endl
       << "causalModelFilename\t" << parameter.causalModelFilename << endl
       << "B\t" << parameter.B << endl << endl
       << "plink_files\t" << parameter.in_files_plink << endl
       << "use_extra_yvalues\t" << parameter.y_value_extra_file << endl  // use_extra_yvalues", y_value_extra_file
       << "trait_index\t" << parameter.in_values_int << endl                  // trait_index", in_values_int
       << "multi_forward_step_max\t" << parameter.ms_MaximalSNPsMultiForwardStep << endl;
       
      
       
/*       
    cout << "END\t" << endl
       << "poolSize\t" << pool.size()  << endl
       << "modelsNo\t" << parameter.modelsNo << endl
       << "maxNoProgressIter\t" << parameter.maxNoProgressIter << endl
       << "pCross\t" << parameter.pCross << endl
       << "pMutation\t " << parameter.pMutation << endl
       << "tournamentSize\t" << parameter.tournamentSize << endl
       << "correlationThreshold\t" << parameter.correlationThreshold << endl
       << "correlationRange\t" << parameter.correlationRange << endl
       << "causalModelFilename\t" << parameter.causalModelFilename << endl
       << "B\t" << parameter.B << endl << endl
       << "plink_files\t" << parameter.in_files_plink << endl
       << "use_extra_yvalues\t" << parameter.y_value_extra_file << endl  // use_extra_yvalues", y_value_extra_file
       << "trait_index\t" << parameter.in_values_int << endl                  // trait_index", in_values_int
       << "multi_forward_step_max\t" << parameter.ms_MaximalSNPsMultiForwardStep << endl;
*/       
    poolFile.flush();
  }
  catch (ofstream::failure e)
  {
    cerr << "Could not write Pool-File" <<endl; 
  }
//  poolFile << ss.str() << endl;
  
  poolFile.close();
}

/**
 * @brief Computes a heritability
 * @param sp_sort - to return an text information about h2_Low, h2_Hight and a pool size
 * @param diff - a vector of diff values of each model in the pool
 * @param tab - a vector of heritability and P(Mi/Y) of each model in the pool 
 */
void GA::coumputeHeritability(stringstream &sp_sort, const vector<long double> &diff, vector<THeri_PMiY> &tab)
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
      Phi_h1 += 1.0 * tab[i].PMi_Y;
    else  
      Phi_h1 += normCDF((h1 - tab[i].h2M) / diff[i]) * tab[i].PMi_Y;
    if ( diff[i] == 0.0 || tab[i].h2M == 0.0)
      Phi += 1.0 * tab[i].PMi_Y;
    else  
      Phi += normCDF((h  - tab[i].h2M) / diff[i]) * tab[i].PMi_Y;
    if ( diff[i] == 0.0 || tab[i].h2M == 0.0)
      Phi_h2 += 1.0  * tab[i].PMi_Y;
    else
      Phi_h2 += normCDF((h2 - tab[i].h2M) / diff[i]) * tab[i].PMi_Y;
  }
  double Phi_05 = Phi;
  while (fabs(Phi - 0.05) > 0.0001 && h1 != h2) // fabs(h1 - h2) > 1.0e-50)
  {
    if (Phi > 0.05)
      h2 = h;
    else
      h1 = h;
    h = (h2 + h1)/2.0;
    Phi = 0.0;
    for (unsigned int i = 0; i < pool.size(); ++i)
    {
      if ( diff[i] == 0.0 || tab[i].h2M == 0.0)
        Phi += 1.0  * tab[i].PMi_Y;
      else
        Phi += normCDF((h - tab[i].h2M) / diff[i]) * tab[i].PMi_Y;
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
      h1 = h;
    else
      h2 = h;
    h = (h2 + h1)/2.0;
    Phi = 0.0;
    for (unsigned int i = 0; i < pool.size(); ++i)
    {
      if ( diff[i] == 0.0 || tab[i].h2M == 0.0)
        Phi += 1.0  * tab[i].PMi_Y;
      else
        Phi += normCDF((h - tab[i].h2M) / diff[i]) * tab[i].PMi_Y;
    }  
  }
  h2M_advH = h;
  sp_sort << "h2_Low: " << h2M_advL << "\th2_Hight: " << h2M_advH << "\tpoolSize: " << pool.size() << endl;
}

/**
 * @brief Computes and writes to file (*_pProb.txt) posterior probalibities of models
 */
void GA::computePosteriorProbability(stringstream &ssModels, map<snp_index_t, int> &mapSNPCausal_ind, vector< multiset<long double> > &tabCausalPost,
                                    vector< multiset<long double> > &tabCausalPost_b)//, long double minPosterior)
{
  if (isCluster == true)
    makeVectCluster(mapSNPid_label);
  long double minMsc, msc;//, maxMsc; 
  vector<long double> sortPool(pool.size());
  int nPool = 0;
  vector<snp_index_t> snps;
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
  for (unsigned         int sP = 0; sP < sortPool.size(); ++sP)
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
  stringstream sp_sort;
  sp_sort << "h2M_adv: " << h2M_adv << "\t";
  coumputeHeritability(sp_sort, diff, tab); 
  sort(tab.begin(), tab.end());
  // calculates $ P(m_i | Y) $
  set<snp_index_t> calculated_j;         // set of calculated snps - to don't calcuate it again
  map<snp_index_t, long double> Pmi_Y;
  multimap<long double, snp_index_t> Pmi_Ysort;
  snp_index_t aSNP;
  i = 0;
  for (set<PoolItem>::iterator it = pool.begin(); it != pool.end(); ++it, ++i)                // for every model in the pool
  {
    snps = (*it).getPoolItem();                                       // takes snps from a model
    msc = (*it).getMsc();
    for (unsigned int j = 0; j < snps.size(); ++j)                    // for every snp from the model
    {
      aSNP = snps[j];
      if (Pmi_Y.find(aSNP) == Pmi_Y.end())
      {  
        Pmi_Y.insert(pair<snp_index_t, long double>(aSNP, 0.0L)); 
      }  
      Pmi_Y[aSNP] = Pmi_Y[aSNP] + exp((minMsc - msc)/2.0L);
    }  
  }
  
  for (map<snp_index_t, long double>::iterator itMap = Pmi_Y.begin(); itMap != Pmi_Y.end(); ++itMap)
  {
    (*itMap).second /= sum;
    Pmi_Ysort.insert(pair<long double, snp_index_t>( (*itMap).second, (*itMap).first));
  }  

  // Power, FDR, Clusters, etc
  set<snp_index_t> mySnps;    
  vector<long double> clusterSNPposterior; // prawd. posteriori dla wynikowego modelu posteriori zbudowanego na postawie klastrów
                                           //
  multimap<long double, snp_index_t>::reverse_iterator itM;
  for (itM = Pmi_Ysort.rbegin(); itM != Pmi_Ysort.rend() && (*itM).first > 0.5; itM++)
  {
    mySnps.insert((*itM).second);
  }
 
  // POWER, FDR, FD count
  TPOWER_FDR powerFDR;
    powerFDR.POWER = 0.0;
    powerFDR.FDR = 0;
    powerFDR.FDcount = 0;
    powerFDR.badSNP = 0;       // the best SNP of posteriori <= 0.5
    powerFDR.posteriorBad = 0.0;
    powerFDR.posteriorBadCluster = 0.0;  
  set<snp_index_t> sPosterior;
  for (multimap<long double, snp_index_t>::reverse_iterator itM = Pmi_Ysort.rbegin(); itM != Pmi_Ysort.rend() && (*itM).first > 0.5; itM++)
    sPosterior.insert((*itM).second);
  vector<snp_index_t>  vPosterior(sPosterior.begin(), sPosterior.end());  
  int bestGA = 0;
  for (unsigned int i = 1; i < modelsNo; ++i)
    if (models[i]->getMSC() < models[bestGA]->getMSC())
      bestGA = i;
  set<snp_index_t> trueSNPs;                        // the set of only true positive snps
  set<snp_index_t> trueSNPsBest;                    // the set of only true positive snps - for the best model in the pool
  Model modelRegion(data);                          // model for a region strategy
  set<TRegionSet_info> setNewRegion;
  if (parameter.causalModelFilename.length() > 0 && isCluster)   // The simulation. We have the causalModel
  {
    vector<snp_index_t> realSNPs;                   // snps from the casual model
    modelReader(parameter.causalModelFilename, realSNPs);
    
    Model *mx = new Model(data);
    if (realSNPs.size() > 0)
      mx->addManySNP(realSNPs);
    if (mx->computeRegression())
      mx->computeMSC();
    else
      mx->computeMSCfalseRegression();
    ssModels << "causalModel: " << *mx << endl;
    delete mx; mx = 0;
//    vector<snp_index_t> clust_max;
    calculatePOWER_FDR_clust_max(realSNPs, powerFDR, Pmi_Y, recognizedSNPs_clusterMax, mapSNPCausal_ind, tabCausalPost);

    sp_sort << "POWER_max: " << powerFDR.POWER;
    sp_sort << "\tFDR_max: " << powerFDR.FDR;
    sp_sort << "\tFDcount_max: " << powerFDR.FDcount << endl;
    sp_sort << "badSNP_max: " << powerFDR.badSNP;
    sp_sort << "\tposteriori_bedSNP_max: " << powerFDR.posteriorBad;
    sp_sort << "\tposteriori_cluster_bedSNP_max: " << powerFDR.posteriorBadCluster << endl;
   
//    vector<snp_index_t> clust_sum;
    calculatePOWER_FDR_clust_sum(realSNPs, powerFDR, Pmi_Y, recognizedSNPs_clusterSum, mapSNPCausal_ind, tabCausalPost);

    sp_sort << "POWER_sum: " << powerFDR.POWER;
    sp_sort << "\tFDR_sum: " << powerFDR.FDR;
    sp_sort << "\tFDcount_sum: " << powerFDR.FDcount << endl;
    sp_sort << "badSNP_sum: " << powerFDR.badSNP;
    sp_sort << "\tposteriori_bedSNP_sum: " << powerFDR.posteriorBad;
    sp_sort << "\tposteriori_cluster_bedSNP_sum: " << powerFDR.posteriorBadCluster << endl;
    
    vector<snp_index_t> v_GA = models[bestGA]->getModelSnps();  // for the best model of the pool
    set<snp_index_t> GA_snps(v_GA.begin(), v_GA.end());         // set of snp's
    //      calculatePOWER_FDR_clustGA(GA_snps, realSNPs, powerFDR, Pmi_Y, mapSNPCausal_ind, tabCausalPost_b);
    //      calculatePOWER_FDR_clustGA(GA_snps, , set<snp_index_t> &TP_S NPs)
    set<snp_index_t> TP_SNPs;
    calculatePOWER_FDR_clustGA(GA_snps, realSNPs, powerFDR, TP_SNPs, recognizedSNPs_bestGA);
    sp_sort << "POWER_bestGA: " << powerFDR.POWER;
    sp_sort << "\tFDR_bestGA: " << powerFDR.FDR;
    sp_sort << "\tFDcount_bestGA: " << powerFDR.FDcount << endl;
    
    set<snp_index_t> mosgwa_snps(v_mosgwa.begin(), v_mosgwa.end());         // set of snp's
    //      calculatePOWER_FDR_clust(mosgwa_snps, realSNPs, powerFDR, Pmi_Y, mapSNPCausal_ind, tabCausalPost_b);
    calculatePOWER_FDR_clustGA(mosgwa_snps, realSNPs, powerFDR, TP_SNPs, recognizedSNPs_mosgwa);
    sp_sort << "POWER_MOSGWA: " << powerFDR.POWER;
    sp_sort << "\tFDR_MOSGWA: " << powerFDR.FDR;
    sp_sort << "\tFDcount_MOSGWA: " << powerFDR.FDcount << endl;
 
    calculatePOWER_FDR_clustGA(sPosterior, realSNPs, powerFDR, TP_SNPs, recognizedSNPs_posterioriModel);
    sp_sort << "POWER_POSTERIORI: " << powerFDR.POWER;
    sp_sort << "\tFDR_POSTERIORI: " << powerFDR.FDR;
    sp_sort << "\tFDcount_POSTERIORI: " << powerFDR.FDcount << endl;
    
    regionStrategy("Region:", modelRegion, Pmi_Ysort, Pmi_Y, regionMinCorrelation, setNewRegion);
    cout << "Region: " << modelRegion << endl;
    ssModels << "Region: " << modelRegion << endl;
    calculatePOWER_FDR_clustRegion(realSNPs, powerFDR, setNewRegion, recognizedSNPs_Region);
    sp_sort << "POWER_REGION: " << powerFDR.POWER;
    sp_sort << "\tFDR_REGION: " << powerFDR.FDR;
    sp_sort << "\tFDcount_REGION: " << powerFDR.FDcount << endl;
    
  }
  else
  {
    cout << "casual snps: N/A" << endl;

    sp_sort << "POWER_max: " << powerFDR.POWER;
    sp_sort << "\tFDR_max: " << powerFDR.FDR;
    sp_sort << "\tFDcount_max: " << powerFDR.FDcount << endl;
    sp_sort << "badSNP_max: " << powerFDR.badSNP;
    sp_sort << "\tposteriori_bedSNP_max: " << powerFDR.posteriorBad;
    sp_sort << "\tposteriori_cluster_bedSNP_max: " << powerFDR.posteriorBadCluster << endl;

    sp_sort << "POWER_sum: " << powerFDR.POWER;
    sp_sort << "\tFDR_sum: " << powerFDR.FDR;
    sp_sort << "\tFDcount_sum: " << powerFDR.FDcount << endl;
    sp_sort << "badSNP_sum: " << powerFDR.badSNP;
    sp_sort << "\tposteriori_bedSNP_sum: " << powerFDR.posteriorBad;
    sp_sort << "\tposteriori_cluster_bedSNP_sum: " << powerFDR.posteriorBadCluster << endl;

    sp_sort << "POWER_bestGA: " << powerFDR.POWER;
    sp_sort << "\tFDR_bestGA: " << powerFDR.FDR;
    sp_sort << "\tFDcount_bestGA: " << powerFDR.FDcount << endl;

    sp_sort << "POWER_REGION: " << powerFDR.POWER;
    sp_sort << "\tFDR_REGION: " << powerFDR.FDR;
    sp_sort << "\tFDcount_REGION: " << powerFDR.FDcount << endl;

    Model *mx = new Model(data);
    if (mx->computeRegression())
      mx->computeMSC();
    else
      mx->computeMSCfalseRegression();
    ssModels << "causalModel: " << *mx << endl;
    delete mx;
//    ssModels << "causalModel: [], msc: 0.0" << endl;
    regionStrategy("Region:", modelRegion, Pmi_Ysort, Pmi_Y, regionMinCorrelation, setNewRegion);
    cout << "Region: " << modelRegion << endl;
    ssModels << "Region: " << modelRegion << endl;
  }

  for (set<TRegionSet_info>::iterator itS = setNewRegion.begin(); itS != setNewRegion.end(); ++itS)
  {
    TRegionSet_info s = *itS;
    delete s.s;
  }
  
/*  // CSDA models  
  {
    string tabNameCSDA[] = 
//     {"rs9525262", "rs10048748", "rs10937559", "rs8028606"};  // hmm25278 CSDA
//     {"rs9525181", "rs10490450", "rs6441934", "rs2044109", "rs17455546"}; //  hmm26651-S CSDA
//      {"rs9525262", "rs17386102", "rs1370718", "rs2044109", "rs2819755", "rs13021147"}; //  hmm32074-S CSDA
      {"rs9525262", "rs10490450", "rs6441934", "rs2044109", "rs17455546", "rs3761945", "rs17326215"};  // hmm34610-S CSDA

    vector<string> SNPs_list(tabNameCSDA, tabNameCSDA + (sizeof(tabNameCSDA)/sizeof(tabNameCSDA[0])));
    data.printSelectedSNPsInMatlab(SNPs_list, "_#CSDA#");

    vector<snp_index_t> v;
    data.findSNPIndex(SNPs_list, v);
//    int snps[] = {535149, 91264, 155080, 575381};
 //   const int snps_size = sizeof(snps)/sizeof(snps[0]);
    
//    vector<snp_index_t> v(snps, snps + snps_size);
    Model m(data);
    m.addManySNP(v);
    if (m.computeRegression())
      m.computeMSC();
    else
    {
      m.computeMSCfalseRegression();
      cerr << "FALSE_regression" << endl;
    }
    
    cout << m << endl;
    m.printModelInMatlab ("#CSDA#");
    saveRecognizedSNPinfo(v, "CSDA");
  }
  
*/  
  
  ssModels << "bestGA: " << *models[bestGA] << endl;

  vector<snp_index_t> vI = models[bestGA]->getModelSnps();
  if (isCluster == true)
    saveRecognizedSNPinfo(vI, "bestGA", Pmi_Y);
//  models[bestGA]->printModelInMatlab ("#bestGA#");
  if (isCluster == true)
    saveRecognizedSNPinfo(v_mosgwa, "mosgwa", Pmi_Y);
  Model  *mX = new Model(data);
  if (v_mosgwa.size() > 0)
  {
    mX->addManySNP(v_mosgwa);
  }
  if (mX->computeRegression())
    mX->computeMSC();
  else
    mX->computeMSCfalseRegression();    
  ssModels << "mosgwa: " << *mX << endl;
  cout << "mosgwa: " << *mX << endl;
//  mX->printModelInMatlab ("#mosgwa#");
  delete mX;
  
  cout << "bestGA: " << *models[bestGA] << ", models[0]: " << *(models[0]) << endl;  
  
//  set<snp_index_t> sPosterior;
  sp_sort << "---------------------------------------" << endl;
  for (multimap<long double, snp_index_t>::reverse_iterator itM = Pmi_Ysort.rbegin(); itM != Pmi_Ysort.rend() && (*itM).first > 0.5; itM++)
  {
    sp_sort << setw(20) << setprecision(16) << (*itM).first << " " << setw(7) << (*itM).second << " " << setw(10) << data.getSNP( (*itM).second ).getSnpId() << endl;
//    sPosterior.insert((*itM).second);
  }
  
  Model *mx = new Model(data);
  if (vPosterior.size() > 0)
     mx->addManySNP(vPosterior);
  if (mx->computeRegression())
     mx->computeMSC();
  else
     mx->computeMSCfalseRegression();
  ssModels << "posterioriModel " << *mx << endl;
  if (isCluster == true)
    saveRecognizedSNPinfo(vPosterior, "posterioriModel", Pmi_Y);
  if (parameter.silent == false)
    cout << "posterioriMode: " << *mx << endl;    
  delete mx;  
  mx = 0;
  if (isCluster == true)   
  {
    vector<snp_index_t> modelSNPs_max;
    vector<snp_index_t> modelSNPs_sum;
    bool clusterSum = true;
    calculate_clusters(Pmi_Y, modelSNPs_max, !clusterSum, "cluster_MAX");
    calculate_clusters(Pmi_Y, modelSNPs_sum, clusterSum, "cluster_SUM");
    Model *mx = new Model(data);
    if (modelSNPs_max.size() > 0)
      mx->addManySNP(modelSNPs_max);
    if (mx->computeRegression())
      mx->computeMSC();
    else
      mx->computeMSCfalseRegression();
    cout << "cluster_max: " << *mx << endl;
    ssModels << "cluster_max: " << *mx << endl;
    delete mx;    
    mx = new Model(data);
    if (modelSNPs_sum.size() > 0)
      mx->addManySNP(modelSNPs_sum);
    if (mx->computeRegression())
      mx->computeMSC();
    else
      mx->computeMSCfalseRegression();
    cout << "cluster_sum: " << *mx << endl;
    ssModels << "cluster_sum: " << *mx << endl;
    delete mx;
  }
  else
  {
    
  }
  
/*  
    powerFDR.POWER = 0.0;
    powerFDR.FDR = 0;
    powerFDR.FDcount = 0;
    powerFDR.badSNP = 0;       // the best SNP of posteriori <= 0.5
    powerFDR.posteriorBad = 0.0;
    powerFDR.posteriorBadCluster = 0.0;    
*/  
  sp_sort << "---------------------------------------" << endl;
  for (multimap<long double, snp_index_t>::reverse_iterator itM = Pmi_Ysort.rbegin(); itM != Pmi_Ysort.rend(); itM++)
  {
    sp_sort << setw(20) << setprecision(16) << (*itM).first << ": " << setw(7) << (*itM).second << " " << setw(10) << data.getSNP( (*itM).second ).getSnpId() << endl;
  }
  ssModels << sp_sort.str();  
}

/**
 * @brief Calculates the statistics (e.g POWER, FDR) of the given SNPs 
 * @param mySnp       - SNPs to test, (the result model)
 * @param realSNPs    - the causal SNPs
 * @param POWER       - returns the POWER
 * @param FDR         - returns the false discovery rate
 * @param FDcount     - returns the number of the false discovery SNPs. WARNING if FDcount is equal 0, FDposterior is unsetted
 * @param trueSNPs    - returnes the set of true discovery SNPs
 **/

void GA::calculatePOWER_FDR(set<snp_index_t> &mySnps, vector<snp_index_t> &realSNPs, long double &POWER, long double &FDR, unsigned int & FDcount, set<snp_index_t> &trueSNPs)
{
  int sumTP = 0;
  int sumFP = 0;
  double abscor;
  map<snp_index_t, int> snpTP;
  set<snp_index_t> causualSNP(realSNPs.begin(), realSNPs.end());
  set<snp_index_t> snp2check = mySnps;
  snp_index_t aSNP;
  for (vector<snp_index_t>::iterator itS = realSNPs.begin(); itS != realSNPs.end(); ++itS)
  {
    aSNP = *itS;
    if ( snp2check.find(aSNP) != snp2check.end() )
    {
      ++sumTP;
      trueSNPs.insert(aSNP);
      if (recognizedSNPs_piMass.size() > 0)
        recognizedSNPs_piMass[aSNP]++;
      snp2check.erase(aSNP);
      causualSNP.erase(aSNP);
    }
  }
  bool TP;
  for (set<snp_index_t>::iterator it = snp2check.begin(); it != snp2check.end(); ++it)
  { // for each snp 
    aSNP = *it;   
    TP = false;
    for (set<snp_index_t>::iterator itC = causualSNP.begin(); itC != causualSNP.end(); ++itC)
    {
      abscor = fabs( data.computeCorrelation( *itC, aSNP) ); // computes the correlation between the two snps
      if (abscor > correlationThreshold)
      {
        ++sumTP;
        trueSNPs.insert(*itC);
        causualSNP.erase(*itC);
        if (recognizedSNPs_piMass.size() > 0)
          recognizedSNPs_piMass[*itC]++;
        TP = true;
        break;
      } 
    }
    if (TP == false)
    {
      ++sumFP;
    }
  }  
  POWER = (sumTP + 0.0) / realSNPs.size();
  if (sumFP + sumTP > 0)
  {
    FDR = (sumFP + 0.0) / (sumFP + sumTP);
    FDcount = sumFP;
  }
  else
  {
    FDR = 0.0;
    FDcount = 0;
//    cerr << "Something wrong: sumFP + sumTP = 0!" << endl;
//    char ch; cout << "Press any key..."; cin >> ch;
  }  
}



/** 
 *  Jeżeli SNP z modelu i przyczynowy SNP należy do tego samego klastra, to mamy TP
 *  w przeciwnym razie FP
 *  Jeżeli dwa lub więcej SNPów należy do tego samego klastra, który nie jest przyczynowy, to mamy 
 *  jeden False Positive.
 */
void GA::calculatePOWER_FDR_clustGA(set<snp_index_t> &mySnps, vector<snp_index_t> &causalSNPs, TPOWER_FDR &powerFDR, set<snp_index_t> &TP_SNPs, map<snp_index_t, int> &recognizedSNPs)
{
  const bool TEST = false;
  TP_SNPs.clear();
  assert(causalSNPs.size() > 0);
  
//    assert(mySnps.size() > 0);
  if (mySnps.size() <= 0)
  {
    powerFDR.POWER = 0;
    powerFDR.FDR = 0;
    powerFDR.FDcount = 0;
    return;
  }
  
  set<int> causalLabels;          // etykiety SNPów przyczynowych
  set<int> FP_LabelSet;          // etykiety klastrów zaklasyfikowanych jako False Positive
  
  map<int, snp_index_t> mapCausalLabel2SNP;
  for (vector<snp_index_t>::iterator it = causalSNPs.begin(); it != causalSNPs.end(); ++it)
  {
    causalLabels.insert(mapSNPid_label[*it]);
    mapCausalLabel2SNP.insert(pair<int, snp_index_t> (mapSNPid_label[*it], *it));
  }
  if (TEST) cout << "causal labels: " << causalLabels << endl;
  if (TEST) cout << "GA labesl    : ";
  
  snp_index_t aSNP_Label;
  for (set<snp_index_t>::iterator it = mySnps.begin(); it != mySnps.end(); ++it)
  {
    aSNP_Label = mapSNPid_label[*it];
    if (causalLabels.find(aSNP_Label) != causalLabels.end())
    {
      TP_SNPs.insert( mapCausalLabel2SNP[aSNP_Label] );
      if (TEST) cout << "+" << aSNP_Label << ", ";
    }
    else
    {
      FP_LabelSet.insert(aSNP_Label);
      if (TEST) cout << "-" << aSNP_Label << ", ";
    }
  }
  for (set<snp_index_t>::iterator it = TP_SNPs.begin(); it != TP_SNPs.end(); ++it)
    ++recognizedSNPs[*it];
  int sumTP = TP_SNPs.size();
  int sumFP = FP_LabelSet.size();
  
  powerFDR.POWER = (sumTP + 0.0) / causalLabels.size();
  if (sumFP + sumTP > 0)
  {
    powerFDR.FDR = (sumFP + 0.0) / (sumFP + sumTP);
  }
  else
  {
//    cerr << "Something wrong: sumFP + sumTP = 0!" << endl;
    powerFDR.FDR = 0;
  }  
  powerFDR.FDcount = sumFP;
  if (TEST) cout << endl 
                 << "POWER: " << powerFDR.POWER << endl
                 << "FDR:   " << powerFDR.FDR << endl
                 << "FD:    " << powerFDR.FDcount << endl;
/*  
  for (set<int>::iterator it = trueLabels.begin(); it != trueLabels.end(); ++it)
  {
    for (unsigned int i = 0; i < causalSNPs.size(); ++i)
      if (mapSNPid_label[ causalSNPs[i] ] == *it)
      {
        trueSNPs.insert( causalSNPs[i] );
        break;
      }
  }
  */
}

// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
// Na potrzeby 2-go artykułu
void GA::calculatePOWER_FDR_clustGA_2ndArticle(set<snp_index_t> &mySnps, vector<snp_index_t> &realSNPs, long double &POWER, long double &FDR, 
                                    unsigned int & FDcount, set<snp_index_t> &trueSNPs)
{
  assert(realSNPs.size() > 0);
  assert(mySnps.size() > 0);

  set<int> causalLabels;  // etykiety SNPów przyczynowych
  set<int> trueLabels;    // etykiety SNPów zaklasyfikowanych jako True Positive
  set<int> mySNPsLabels;  // etykiety SNPów wskazanych przez GA
  
  for (set<snp_index_t>::iterator it = mySnps.begin(); it != mySnps.end(); ++it)
    mySNPsLabels.insert(mapSNPid_label[*it]);
  
  for (vector<snp_index_t>::iterator it = realSNPs.begin(); it != realSNPs.end(); ++it)
    causalLabels.insert(mapSNPid_label[*it]);
  
  int sumTP = 0;
  int sumFP = 0;
  
  cout << "causal labels: " << causalLabels << endl;
  cout << "GA labesl    : " << mySNPsLabels << endl;
  for (set<int>::iterator it = mySNPsLabels.begin(); it != mySNPsLabels.end(); ++it)
  {
    if (causalLabels.find(*it) != causalLabels.end())
    {
      ++sumTP;
      trueLabels.insert(*it);
    }
    else
      ++sumFP;
  }
  POWER = (sumTP + 0.0) / causalLabels.size();
  if (sumFP + sumTP > 0)
  {
    FDR = (sumFP + 0.0) / (sumFP + sumTP);
  }
  else
  {
//    cerr << "Something wrong: sumFP + sumTP = 0!" << endl;
    FDR = 0;
  }  
  FDcount = sumFP;
  
  for (set<int>::iterator it = trueLabels.begin(); it != trueLabels.end(); ++it)
  {
    for (unsigned int i = 0; i < realSNPs.size(); ++i)
      if (mapSNPid_label[ realSNPs[i] ] == *it)
      {
        trueSNPs.insert( realSNPs[i] );
        break;
      }
  }
}

/**-------------------------------------------------------------------------
 * @brief Model selection for GA. It works like selectModel() but makes only a stepwise selection.
 * @param model a model of selection
 * -------------------------------------------------------------------------
 */
void GA::selectModel ( Model& model ) {
  printLOG( "GA: Model Selection started: " );
  if (model.computeRegression())
    model.computeMSC();
  else
    model.computeMSCfalseRegression();
  model.makeMultiForwardStep(0, Parameter::selectionCriterium_BIC, NULL, &exclusivedSNP, &goodSNPs );
} 

/** 
 * @brief Reads the pool of models from a given file 
 * @param fileName the name of the data file
 */
void GA::poolReader(string fileName, stringstream& sGATime, int real_modelsNo)
{
  pool.clear();
  vector<snp_index_t> vBest;   // to remember the best model in the pool
  long double bestMSC = -1.0;
  pool.clear();
  vector<snp_index_t> v;
  int snp;
  string s;
  long double h2;
  int size;
  char chr;
  string line;
  ifstream  poolFile;
  ofstream  poolFileOut;
  unsigned int id;
  double msc;
  stringstream ss;
  poolFile.open( fileName.c_str(), ios::in );
  cout << "Try open the file: " << fileName.c_str() << endl;
  bool isTime = true;
  int initial_population_Count = 0;
  bool isStepwise = false;
  vector <snp_index_t> temp_mosgwa;
  if ( poolFile.is_open() )
  {
    cout << "pool file is opened" << endl;
    while ( ! poolFile.eof() ) 
    {
      if (isTime == true)
      {
        chr = poolFile.peek();
        if (chr == 'G')
        {
          poolFile >> s;
          sGATime << s << " ";
          poolFile >> s;
          sGATime << s << " ";
          poolFile >> s;
          sGATime << s << endl;
          cout << sGATime.str() << endl;
        }   
        isTime = false;
      }  
      char char_id;
      poolFile >> char_id;  // P
      if (char_id == 'I')
        ++initial_population_Count;
      if (char_id == 'E')  // after the last data line in the file is a line 'END'
      {
        chr = poolFile.peek();
        if (chr == 'N')
          break;
      }
      if (poolFile.eof())
        break;
      poolFile >> id;   // id
      poolFile >> line; // "msc:"
      poolFile >> msc;  // msc
      if (poolFile.fail() == true)
      {
        cerr << "Something gone wrong. Incorect file format. Wrong msc value: " << msc << endl
             << "id: " << id << endl;
        exit(-1);
      }
      poolFile >> line;  // "size:"
      poolFile >> size;  // size
      if (poolFile.fail() == true)
      {
        cerr << "Something gone wrong. Incorect file format. Wrong the model size: " << size << endl
             << "id: " << id << endl;
        exit(-1);
      }      
      poolFile >> line; // "h2_M:"
      poolFile >> h2;   // h2_M
      if (poolFile.fail() == true)
      {
        cerr << "Something gone wrong. Incorect file format. Wrong the h2_M value: " << h2 << endl
             << "id: " << id << endl;
        exit(-1);
      }            
      poolFile >> chr;  // [
      chr = poolFile.peek();   // check for zero model "[]"
      while (chr != ']')
      {
         poolFile >> snp;
         if (poolFile.fail() == true)
         {
           cerr << "Something gone wrong. Incorect file format, snp:" << snp << endl
                << "id: " << id << endl;
           exit(-1);
         }
         v.push_back(snp);
         chr = poolFile.peek();
         if (chr != ']')
           poolFile >> s;  //  ", "
      }
      poolFile >> s; // '],' for zero model
      PoolItem p(v, msc, h2, id, char_id);
      if (id == 2)
      {
        v_mosgwa = v;
      }
      if (char_id == 'S')
      {
        isStepwise = true;
        temp_mosgwa = v;
      }
      if (msc < bestMSC || bestMSC < 0.0)
      {
        bestMSC = msc;
        vBest = v;
      }
      pool.insert(p); 
      
      testSet.clear();
      testSet.insert(v.begin(), v.end());
      if (v.size() != testSet.size())
      {
        cerr << "ERROR: the same snps in poolReader(), info: " << endl << v << ", msc: " << msc << endl;
        cerr << testSet << endl;
        char c; cerr << "Press a key... "; cin >> c;
      }
      v.clear();
    }
    poolFile.close();
  }  
  else
  {
    cerr << "file: " << fileName << " is not open" << endl;
    exit (-1);
  }
  cout << "bestMSC: " << bestMSC << ", vBest = " << vBest << endl;
  cout << "initial_population_Count: " << initial_population_Count << endl
       << "modelsNo: " << real_modelsNo << endl;
  if (isStepwise == true)
  {
    v_mosgwa = temp_mosgwa;
  }
  else
    if (initial_population_Count <= real_modelsNo)
    {
      v_mosgwa.clear();
      cout << "Zero model of MOSGWA: " << v_mosgwa << endl;
    }
  if (models == 0)
  {
    models = new Model*[modelsNo];
    for (unsigned int i = 0; i < modelsNo; ++i)
      models[i] = 0;
  }  
  for (unsigned int i = 0; i < modelsNo; ++i)
  {
    if (models[i])
      delete models[i];
    models[i] = new Model(data);
    if (vBest.size() > 0)
      models[i]->addManySNP(vBest);
    if (models[i]->computeRegression())
      models[i]->computeMSC();
    else
    {
      cerr << "poolReader(). FAIL computeRegression for models[i]" << endl;
      exit(-1);
    }    
  }
  cout << "poolReader best GA: " << *models[0] << endl; 
  cout << "end poolReader" << endl;
}

/** 
 * @brief Reads the causal model from a file (the file name is in *.conf)
 * @param fileName - the file name where the model data is stored
 */
void GA::modelReader(string fileName, vector<snp_index_t> &v)
{
  int snp;
  int modelSize;
  string line, s;
  ifstream poolFile;
  stringstream ss;
  v.clear();
  poolFile.open( fileName.c_str(), ios::in ); 
  poolFile.clear();
  if (poolFile.is_open())
  {
    getline(poolFile, line);
    getline(poolFile, line);
    ss.str(line);
    ss.clear();
    ss >> line >> modelSize;
    v.resize(modelSize);
    getline(poolFile, line);
    int i = 0;
    while ( ! poolFile.eof() ) 
    {
      getline(poolFile, line);
      ss.str(line);
      if (line == "")
        continue;
      ss.clear();
      ss >> snp;
      if (ss.fail())
      {
        unsigned int loc = line.find("Intercept", 0);
        if( loc != string::npos ) 
          ss.clear();
        else
        {
          cerr << "Wrong file format, Intercept" << endl;
          exit(-1);
        }
      }
      else
        v[i++] = snp;
    }
  }
  else
  {
    cerr << "GA::modelReader() function: Can not open file: " << fileName << endl;
    exit(-1);
  }
  if (realModel != 0)
    delete realModel;
  realModel = new Model(data);
  if (v.size() > 0)
    realModel->addManySNP(v);
  if (realModel->computeRegression())
    realModel->computeMSC();
  else
  {
    cerr << "something wrong with real model" << endl;
    exit (-1);
  }
  delete realModel;
  realModel = 0;
}

/**
 * @brief Extracts some information from piMass files
 * @note not used for the model selection
 */
void GA::piMassExtract(const string &fileName, string &outFileName)
{
  int intVal;
  double dVal;
  string line, s;  
  ofstream outFile;
  ofstream modelsFile;
  outFile.open((outFileName + "piMass_PowerFRD.txt").c_str(), ios::app); 
  modelsFile.open((outFileName + "piMassModels.txt").c_str(), ios::app);
  
  stringstream sp_sort;
  map<int, unsigned int> snp;  // liczebność poszczególnych snp'ów
  ifstream  aFile;
  stringstream ss;
  map<snp_index_t, long double> Pmi_Y;
  multimap<long double, snp_index_t> Pmi_Ysort;
  aFile.open( (fileName + "/pref.mcmc.txt").c_str(), ios::in ); 
  aFile.clear();    
  cout << "a file: " << (fileName + "/pref.mcmc.txt").c_str() << " is opened" << endl;
  getline(aFile, line);
  stringstream sInfo;
  while (!aFile.eof())
  {
    aFile >> s;                     // rs
    aFile >> intVal;                // chr
    aFile >> intVal;                // position
    aFile >> dVal;                  // postc
    --intVal;
    Pmi_Y.insert(pair<snp_index_t, long double> (intVal, dVal));
    aFile >> dVal;                  // postb
    aFile >> dVal;                  // beta
    aFile >> dVal;                  // betarb
  }  // reads id and posterior. (odczyt id i poster z jednego pliku)
  aFile.close();
  cout << "Pmi_Y.size = " << Pmi_Y.size() << endl;
  for (map<snp_index_t, long double>::iterator itMap = Pmi_Y.begin(); itMap != Pmi_Y.end(); ++itMap)
  {
    Pmi_Ysort.insert(pair<long double, snp_index_t>( (*itMap).second, (*itMap).first));
  }
  set<snp_index_t> mySnps;    
  for (multimap<long double, snp_index_t>::reverse_iterator itM = Pmi_Ysort.rbegin(); itM != Pmi_Ysort.rend() && (*itM).first > 0.5; itM++)
  {
    cout << setw(15) << setprecision(11) << (*itM).first << ": " << setw(7) << (*itM).second << " " << setw(10) << data.getSNP( (*itM).second ).getSnpId() << endl;
    mySnps.insert((*itM).second);
  }
  
  Model piMassModel(data);
  vector<snp_index_t> v_piMass(mySnps.begin(), mySnps.end());
  if (v_piMass.size() > 0)
    piMassModel.addManySNP(v_piMass);
  if (piMassModel.computeRegression())
    piMassModel.computeMSC();
  else
  {
    cerr << "something wrong with piMass model - piMassExtract" << endl;
    exit (-1);
  }
  cout << "piMass model: " << piMassModel << endl;
  modelsFile << piMassModel << endl;
  
  Model *mx = 0;
  if (parameter.causalModelFilename.length() > 0)
  {
    vector<snp_index_t> realSNPs;             // snps from real model
    modelReader(parameter.causalModelFilename, realSNPs);  
    
    long double POWER, FDR;
    unsigned int FDcount;
    if (realSNPs.size() > 0)
    {
      mx = new Model(data);
      if (realSNPs.size() > 0)
        mx->addManySNP(realSNPs);
      if (mx->computeRegression())
        mx->computeMSC();
      else
        mx->computeMSCfalseRegression();
      cout << "causalModel: " << *mx << endl;    
      delete mx;  
      mx = 0;
      set<snp_index_t> trueSNPs;
      calculatePOWER_FDR(mySnps, realSNPs, POWER, FDR, FDcount, trueSNPs);//, Pmi_Y, tab_Prob);
      cout << "POWER: " <<  POWER <<  ", FDR: " << FDR << ", FDcount: " << FDcount << endl;
      sp_sort << POWER << "\t" << FDR << "\t" << FDcount << endl;
    }
    else
    {
      sp_sort << "\tPOWER: " << "NA";
      sp_sort << "\tFPD: " << "NA" << endl;
      cerr << "The real model size is zero!" << endl;
    }
  }
  else
  {
    cout << "casual snps: N/A" << endl;
    sp_sort << "\tPOWER: " << "NA";
    sp_sort << "\tFPD: " << "NA" << endl;    
  }
  // number of SNPs (liczebność snpów)
  for (set<snp_index_t>::iterator it = mySnps.begin(); it != mySnps.end(); it++)
  {
    if (snp.find(*it) != snp.end())
    {
      ++snp[*it];
    }
    else
    {
      snp.insert(pair<int, unsigned int>(*it, 1));
    }
  }
  outFile << sp_sort.str() << endl;    
  outFile.flush();
  outFile.close();
  
  ofstream file;      
  file.open((outFileName + "piMassModels_result_snp.txt").c_str(), ios::out);  
  for (map<int, unsigned int>::iterator it = snp.begin(); it != snp.end(); ++it)
    file << (*it).first << "\t" << (*it).second << endl;
  file.close();  
}

/**
 * @brief Reads the initial population of GA form a fila
 * <TODO check the input data format>
 */
int GA::readInitialPop(string fileName, set<PoolItem> &population, int popSize)
{
  population.clear();
  vector<snp_index_t> v;
  int snp;
  string s;
  int id;
  int size;
  char chr;
  string line;
  double msc;
  ifstream poolFile;
  stringstream ss;
  poolFile.open( fileName.c_str(), ios::in );
  poolFile.clear();
  if ( poolFile.is_open())
  {
    while ( ! poolFile.eof() && popSize) 
    {
      poolFile >> chr;     // id
      poolFile >> id;      // number
      poolFile >> s;   // ", msc:"
      poolFile >> msc;   
      //cout << "msc " << msc << endl;
      if (poolFile.fail() == true)
      {
        cerr << "Something gone wrong. Incorect file format. Wrong msc value: " << msc << endl;
        exit(-1);
      }
      poolFile >> s;  //  ", size:"
      poolFile >> size;
      //cout << "size " << size << endl;
      if (poolFile.fail() == true)
      {
        cerr << "Something gone wrong. Incorect file format. Wrong model size value: " << size << endl;
        exit(-1);
      }
      long double h;
      poolFile >> s;  // , h2_M:
      poolFile >> h;
      if (poolFile.fail() == true)
      {
        cerr << "Something gone wrong. Incorect file format. Wrong model heraitability value: " << size << endl;
        exit(-1);
      }
      
      if (poolFile.fail() == true)
        break;
      poolFile >> chr;      // [
      chr = poolFile.peek();   // check for zero model "[]"
      if (chr != ']')
      {
        poolFile >> snp;
        if (poolFile.fail() == true)
        {
          cerr << "Something gone wrong. Incorect file format" << endl;
          exit(-1);
        }
        while (poolFile.fail() != true)
        {
          v.push_back(snp);
          poolFile >> s;  //  ", "
          poolFile >> snp;
        }
      }
      else
      {
        poolFile >> s; // '],' for zero model
      }
      poolFile.clear();
      
      PoolItem p(v, msc, h, id, 'I');
      population.insert(p); 
      v.clear();
      --popSize;      
    }
    poolFile.close();
  } 
  else 
  {
    cerr << "something gone wrong" << endl;
    exit(1);
  }
  return 0;
}

/** 
 * @brief Returns number of models which have got better msc value than the mscVal parameter
 * @param mscValue - msc value of the checked model 
 */
int GA::isInNBestModels(double mscValue)
{
  int n = 0;
  for (unsigned int i = 0; i < modelsNo; i++)
    if (models[i]->getMSC() < mscValue)
      ++n; 
  return n;
}


/**
 * @brief Reads data from a file fileName, prepares data and makes calculations
 *        (Wyznacza model na podstawie klastrów z plików *sort)
 * @param fileName - the file name of posterior data  // (nazwa pliku z danymi - sortPosteriori)
 * @param minPosterior - the minimal posterior value, default 0.01
 * @note not used for model selection
 */
void GA::calculateClusterPosterior(string fileName, long double minPosterior)
{
  vector<snp_index_t> clusterSNP; // wynikowy model posteriori zbudowany na postawie klastrów
  vector<long double> clusterSNPposterior; // prawd. posteriori dla wynikowego modelu posteriori zbudowanego na postawie klastrów
  ifstream file;
  file.open( fileName.c_str(), ios::in );
  file.clear();
  vector<long double> posteriorVect;
  long double val,  
              valPost; // value of postrerior
  string line;
  if ( file.is_open() )
  {
    while (line.compare("---------------------------------------") != 0)
    {
      getline(file, line);
      //cout << line << endl;
    }
    file >> valPost;                       //  posteriori
    while ( valPost >= minPosterior)
    {
      clusterSNPposterior.push_back(valPost);
      file >> line;                    //  :
      file >> val;                     //  
      if (valPost >= 0.5)
      { // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! SPRAWDZIĆ !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
//        originalSNP.push_back(val);
 //       cout << "post: " << valPost << "\tsnp: " << val << endl;
      }
      clusterSNP.push_back(val);       //  snp_id
      file >> line;                    //  snp name
      file >> valPost;                     //  posteriori
    } 
  }  
  else
  {
    cerr << "calcuateClusterPosterior(), file not open: " << fileName << endl;
    exit(-1);
  }  
  file.close();
  //clusterSNPposterior = originalSNPposterior;
  cout << "clusterSNP: " << clusterSNP << endl;; // wynikowy model posteriori zbudowany na postawie klastrów
  cout << "clusterSNPposterior:: " << clusterSNPposterior << endl; // prawd. posteriori dla wynikowego modelu posteriori zbudowanego na postawie klastrów
//  cout << "originalSNP: " << originalSNP << endl; // oryginalne model
//  cout << "originalSNPposterior: " << originalSNPposterior << endl;
  /*
   * clasterSNP zawiera snpy o prawd. a posteriori > minPosterior
   * originalSNP zawiera snpy o prawd. a posteriori > 0.5
   * clusterSNPposterior zawiera wartości a posteriori > minPosterior
   * originalSNPposterior zawiera wartości a posteriori > minPosterior
   **/
  calculateClusterPosterior(clusterSNP, clusterSNPposterior);//, originalSNP, originalSNPposterior);
}


/*
 * clasterSNP zawiera snpy o prawd. a posteriori > minPosterior
 * originalSNP zawiera snpy o prawd. a posteriori > 0.5
 * clusterSNPposterior zawiera wartości a posteriori > minPosterior
 * originalSNPposterior zawiera wartości a posteriori > minPosterior
 * @brief Calculates the cluster model
 * There is the cluster model in the clusterSNP - get back
 */
void GA::calculateClusterPosterior(vector<snp_index_t> &clusterSNP, vector<long double> &clusterSNPposterior)
{
  vector<long double> originalSNPposterior;
  originalSNPposterior = clusterSNPposterior;
//  cout << "inside <<<<" << endl;
//  cout << "clusterSNP: " << clusterSNP << endl;; // wynikowy model posteriori zbudowany na postawie klastrów
//  cout << "clusterSNPposterior:: " << clusterSNPposterior << endl; // prawd. posteriori dla wynikowego modelu posteriori zbudowanego na postawie klastrów
  //cout << "originalSNP: " << originalSNP << endl; // oryginalne model
  //cout << "originalSNPposterior: " << originalSNPposterior << endl << endl;
  
  int *tabCorrelatedWith = new int[clusterSNP.size()]; // wskazuje, z kŧórym snp'em (indeks w clusterSNP) jest skorelowany
  long double correlation;
  int N = clusterSNP.size();
  double *table = new double[N * N];     // tablica z wartościami korelacji pomiędzy poszczególnymi snp'ami
                                         // pod główną przekątną
  for (int i = 0; i < N * N; ++i)
    table[i] = 0.0;
  
  for (int r = clusterSNP.size() - 1; r >= 0; --r)  // kolumna
  {
    tabCorrelatedWith[r] = -1;
    for (int c = r - 1; c >= 0; --c) // wiersz
    {
      correlation = fabs( data.computeCorrelation( clusterSNP[r], clusterSNP[c]) );
      table[r * N + c] = correlation;
      if (correlation > correlationTh)
      {
         clusterSNPposterior[c] += clusterSNPposterior[r];
      }
    }
  }
  ofstream outFile;
  outFile.open((parameter.out_file_name + "_cluster_" + int2str(parameter.in_values_int) + ".txt").c_str(), ios::app);   
  outFile << "snp_id\t org. poster\t new poster.\tcorrelated with\tid";
  for (int i = 0; i < N; ++i)
  {
    outFile << "\t" << i;
  }
  outFile << endl;
  outFile.flush();
  //cout << endl << "to file <<< " << endl << endl;
  for (int i = 0; i < N; ++i)
  {
    outFile << clusterSNP[i] << "\t" << originalSNPposterior[i] << "\t" << clusterSNPposterior[i];
  //  cout << clusterSNP[i] << "\t" << originalSNPposterior[i] << "\t" << clusterSNPposterior[i];
    if (clusterSNPposterior[i] - originalSNPposterior[i] > 0)
    {
      outFile << "\t" << clusterSNPposterior[i] - originalSNPposterior[i] << "\t" << i;
      //cout << "\t" << clusterSNPposterior[i] - originalSNPposterior[i] << "\t" << i;
    }
    else
    {
      outFile << "\t  " << "\t" << i;
//      cout  << "\t  " << "\t" << i;
    }
    for (int j = 0; j < N; ++j)
    {
      outFile << "\t" << table[i*N + j];
//      cout  << "\t" << table[i*N + j];
      outFile.flush();
    }
    outFile << endl;
  }  
//  cout << endl << endl;
  outFile.close();
  outFile.open((parameter.out_file_name + "_cluster_msc_" + int2str(parameter.in_values_int) + ".txt").c_str(), ios::app);   
  vector<snp_index_t> snps;
  Model *mx = new Model(data);
  for (unsigned int i = 0; i < originalSNPposterior.size(); ++i)
    if (originalSNPposterior[i] > 0.5)
     snps.push_back(clusterSNP[i]);
  if (snps.size() > 0)  
    mx->addManySNP(snps);
  if (mx->computeRegression())
  {
    mx->computeMSC();
  }  
  else
  {
    mx->computeMSCfalseRegression();
    cerr << "MSCfalseRegression,calculateClusterPosterior(4p)" << *mx << endl;    
    exit(-1);
  }
  
  vector<snp_index_t> v = mx->getModelSnps();
  outFile << v << "\t" << mx->getMSC() << "\t";
  delete mx; 
  snps.clear();
  mx = new Model(data);
  for (unsigned int i = 0; i < clusterSNPposterior.size(); ++i)
  {
    if (clusterSNPposterior[i] > 0.5)
    {
      snps.push_back(clusterSNP[i]);
    }
  }
  if (snps.size() > 0)
    mx->addManySNP(snps);
  if (mx->computeRegression())
  {
    mx->computeMSC();
  }  
  else
  {
    mx->computeMSCfalseRegression();
    cerr << "MSCfalseRegression, calculateClusterPosterior(4p)" << *mx << endl;    
  }
  clusterSNP = snps;  // zwrócenie modelu obliczonego przy pomocy klastrów
  v = mx->getModelSnps();
  outFile << mx->getMSC() << "\t" << v << endl;
  outFile.close();

  /*
  if (parameter.causalModelFilename.length() > 0)
  {
    vector<snp_index_t> realSNPs;             // snps from real model
    modelReader(parameter.causalModelFilename, realSNPs);  
    
    double POWER, FDR;
    if (realSNPs.size() > 0)
    {
      outFile.open((parameter.out_file_name + "_cluster_msc_POWE_FDR" + int2str(parameter.in_values_int) + ".txt").c_str(), ios::app);   
      set<snp_index_t> mySnps(snps.begin(), snps.end());
      //mySnps.assigne(snps.begin(), snps.end());
      calculatePOWER_FDR(mySnps, realSNPs, POWER, FDR);
      outFile << setprecision(14) <<  mx->getMSC() << "\t" << POWER << "\t" << FDR << endl;
      outFile.close();
    }
    else
    {
      cerr << "The real model size is zero!" << endl;
    }
  }
  else
  {
    outFile.open((parameter.out_file_name + "_cluster_msc_POWE_FDR" + int2str(parameter.in_values_int) + ".txt").c_str(), ios::app);   
    outFile << "casual snps: N/A" << endl;
    outFile.close();
  }  
  */
  delete [] table;
  delete [] tabCorrelatedWith;
}

/**
 * Oblicza korelacje i zapisuje wyniki w tabeli:
 * snp  casualSNP  maxCorr
 */
void GA::checkCorrelation(set<snp_index_t> &mySnps, vector<snp_index_t> &realSNPs)
{
  int size = mySnps.size();
  double *snp = new double[size];
  double *causal = new double[size];
  double *correlation = new double[size];
  sort(realSNPs.begin(), realSNPs.end());
  int aSNP;
  double maxCorrelation;
  snp_index_t causalSNPmax = 0;
  double absCorr;
  int counter = 0;
  ofstream outFile;  
  outFile.open((parameter.out_file_name + "_All_correlation_" + int2str(parameter.in_values_int) + ".txt").c_str(), ios::app);
  outFile << " \t";
  for (unsigned int i = 0; i < realSNPs.size(); ++i)
    outFile << realSNPs[i] << "\t";
  outFile << endl;
  for (set<snp_index_t>::iterator it = mySnps.begin(); it != mySnps.end(); ++it)
  {
    aSNP = *it; 
    maxCorrelation = -1.0;
    outFile << aSNP << "\t";
    for (unsigned int i = 0; i < realSNPs.size(); ++i)
    {
      absCorr = fabs( data.computeCorrelation( realSNPs[i], aSNP) ); // compute correlation between snps
      outFile << absCorr << "\t";
      if (absCorr > maxCorrelation)
      {
        maxCorrelation = absCorr;
        causalSNPmax = realSNPs[i];
      }
    }
    outFile << maxCorrelation << endl;
    snp[counter] = aSNP;
    causal[counter] = causalSNPmax;
    correlation[counter] = maxCorrelation;
    cout << snp[counter] << "\t" << causal[counter] << "\t" << correlation[counter] << endl;
    counter++;
  }
  outFile.close();
  //ofstream outFile;  
  outFile.open((parameter.out_file_name + "_correlation_" + int2str(parameter.in_values_int) + ".txt").c_str(), ios::app);
  
  if (outFile.is_open() == true)
  {
    for (int i = 0; i < size; ++i)
      outFile << snp[i] << "\t" << causal[i] << "\t" << correlation[i] << endl;
  }
  else
  {
    cerr << "Can not open the file: " << endl;
    exit(-1);
  }  
  cout << "file: " << (parameter.out_file_name + "_correlation_" + int2str(parameter.in_values_int) + ".txt").c_str() << endl;
}




void GA::setRecognisedSNPs()
{

  recognizedSNPs_Region.clear();
  recognizedSNPs_bestGA.clear();
  recognizedSNPs_mosgwa.clear();
  recognizedSNPs_posterioriModel.clear();
  recognizedSNPs_clusterMax.clear();
  recognizedSNPs_clusterSum.clear();

  if (parameter.causalModelFilename.length() > 0)
  {
    vector<snp_index_t> realSNPs;             // snps from real model
    modelReader(parameter.causalModelFilename, realSNPs);   
    for (unsigned int i = 0; i < realSNPs.size(); ++i)
    {
      recognizedSNPs_Region.insert(pair<snp_index_t, int>(realSNPs[i], 0));
      recognizedSNPs_bestGA.insert(pair<snp_index_t, int>(realSNPs[i], 0));
      recognizedSNPs_mosgwa.insert(pair<snp_index_t, int>(realSNPs[i], 0));
      recognizedSNPs_posterioriModel.insert(pair<snp_index_t, int>(realSNPs[i], 0));
      recognizedSNPs_clusterMax.insert(pair<snp_index_t, int>(realSNPs[i], 0));
      recognizedSNPs_clusterSum.insert(pair<snp_index_t, int>(realSNPs[i], 0));
    } 
  }
}
/*
double GA::theBestCorrelatedSNP(snp_index_t aSNP, set<snp_index_t> & snps, snp_index_t &bestCorrelatedSNP) const
{
  double max_correlation = 0.0;
  double abscor;
  for (set<snp_index_t>::iterator itC = snps.begin(); itC != snps.end(); ++itC)
  {
    abscor = fabs( data.computeCorrelation( *itC, aSNP) ); // computes the correlation between the two snps
    if (abscor > max_correlation)
    {
      bestCorrelatedSNP = *itC;
      max_correlation = abscor;
    } 
  }  
  return max_correlation;
}
*/

void GA::makeClusers(vector<set<snp_index_t> > &tab)
{
  //data.calculateIndividualTests();          // działamy tylko na snpach, których p-val < 0.1
  vector<snp_index_t> realSNPs;             // snps from the real model
  modelReader(parameter.causalModelFilename, realSNPs);
  goodSNPs.reset();
  unsigned int n = 0;
  long double abscorr;
  snp_index_t aSNP;
  
  while ( n < data.getSnpNo())
  {
    aSNP = n;
    ++n;
    //cout << n;
    for (vector<snp_index_t>::iterator it = realSNPs.begin(); it != realSNPs.end(); ++it) // casual SNPs
    {
      abscorr = fabs( data.computeCorrelation( *it, aSNP) ); // computes the correlation between the two snps   
      
      if (abscorr > 0.5)
      {
        tab[ cs[*it] ].insert(aSNP);
        //cout << "inserted snp:" << aSNP << endl;
        //char c; cout << "press.."; cin >> c;
      }  
    }  
  }
  int realSNPs_id[] = {11, 802, 1601, 2399, 3201, 4000, 4799, 5599, 6402, 7201, 8001, 8801, 9601, 10399, 11199, 12002,
                    12801, 13603, 14401, 15205, 16004, 16803, 17602, 18402, 19196, 20001, 20801, 21599, 22402, 23201};
  cout << tab.size() << endl;                    

 for (unsigned int i = 0; i < tab.size(); ++i)// vector<set<snp_index_t> >::iterator it = tab.begin(); it != tab.end(); ++it)
  {
    cout << realSNPs_id[i] << ">> {";
    for (set<snp_index_t>::iterator it = tab[i].begin(); it != tab[i].end(); ++it)
    {
      if (it == tab[i].begin())
        cout << *it;
      else
        cout << ", " << *it;
    }    
    cout << "}" << endl;
  }
}


void GA::piMassCluserPOWER(const string &fileName, const string &outFileName, map<snp_index_t, int> &mapSNPCausal_ind, vector< multiset<long double> > &tabCausalPost)
{
  makeVectCluster(mapSNPid_label);
  int intVal;
  double dVal;
  string line, s;  
  ofstream outFile;
  outFile.open((outFileName + "piMass_PowerFRD.txt").c_str(), ios::app); 
  ofstream modelsFile;
  modelsFile.open((outFileName + "piMassModels.txt").c_str(), ios::app);
  
  stringstream sp_sort;
  stringstream ssModels;
  map<int, unsigned int> snp;  // liczebność poszczególnych snp'ów
  
  ifstream  aFile;
  stringstream ss;
  map<snp_index_t, long double> Pmi_Y;
  multimap<long double, snp_index_t> Pmi_Ysort;
  aFile.open( (fileName + "/pref_#" + int2str(parameter.in_values_int) + ".mcmc.txt").c_str(), ios::in ); // pref_#1.mcmc.txt
  aFile.clear();    
  cout << "a file: " << (fileName + "/pref_#" + int2str(parameter.in_values_int) + ".mcmc.txt").c_str() << " is opened" << endl;
  getline(aFile, line);
  stringstream sInfo;
  while (!aFile.eof())
  {
    aFile >> s;                     // rs
    aFile >> intVal;                // chr
    aFile >> intVal;                // position
    aFile >> dVal;                  // postc
    --intVal;
    Pmi_Y.insert(pair<snp_index_t, long double> (intVal, dVal));
    aFile >> dVal;                  // postb
    //Pmi_Y.insert(pair<snp_index_t, long double> (intVal, dVal));
    aFile >> dVal;                  // beta
    aFile >> dVal;                  // betarb
  }  // odczyt id i poster z jednego pliku
  aFile.close();
  cout << "Pmi_Y.size = " << Pmi_Y.size() << endl;
  for (map<snp_index_t, long double>::iterator itMap = Pmi_Y.begin(); itMap != Pmi_Y.end(); ++itMap)
  {
    Pmi_Ysort.insert(pair<long double, snp_index_t>( (*itMap).second, (*itMap).first));
  }
  
  set<snp_index_t> mySnps;    
  bool isBad = false;
  //snp_index_t badSNP;  // najlepszy snp z posteriori <= 0.5
  multimap<long double, snp_index_t>::reverse_iterator itM;
  for (itM = Pmi_Ysort.rbegin(); itM != Pmi_Ysort.rend() && (*itM).first >  0.5; itM++)
  {
    mySnps.insert((*itM).second);
    cout << setw(15) << setprecision(11) << (*itM).first << ": " << setw(7) << (*itM).second << " " << setw(10) << data.getSNP( (*itM).second ).getSnpId() << endl;
  }
  if (itM != Pmi_Ysort.rend())
  {  
//     badSNP = (*itM).second;
    isBad = true;
  }  
  // w mySnp mamy model wynikowy piMasa
  
  TPOWER_FDR powerFDR;
  powerFDR.POWER = 0.0;
  powerFDR.FDR = 0.0;
  powerFDR.FDcount = 0;
  powerFDR.badSNP = 0;       // the best SNP of posteriori <= 0.5
  powerFDR.posteriorBad = 0.0;
  powerFDR.posteriorBadCluster = 0.0;  
  
  if (parameter.causalModelFilename.length() > 0)
  {
    vector<snp_index_t> realSNPs;             // snps from real model
    modelReader(parameter.causalModelFilename, realSNPs);  
    
    if (realSNPs.size() > 0 && isBad == true)
    {
      set<snp_index_t> trueSNPs;
      //calculatePOWER_FDR_clust(mySnps, realSNPs, POWER, FDR, FDcount, trueSNPs, tabClust, Pmi_Y);
      cout << "causal model: OK" << endl; 
      Model m(data);
      vector<snp_index_t> vM(mySnps.begin(), mySnps.end());
      if (vM.size() > 0)
      {
        m.addManySNP(vM);
      }
      if (m.computeRegression())
        m.computeMSC();
      else
      {
        m.computeMSCfalseRegression();
        cerr << "FALSE_regression" << endl;
      }
      modelsFile << m << endl;
      modelsFile.close();

      // calculatePOWER_FDR_clust_max(realSNPs, powerFDR, Pmi_Y, recognizedSNPs_clusterMax, mapSNPCausal_ind, tabCausalPost);      
      //calculatePOWER_FDR_clust_sum(realSNPs, powerFDR, Pmi_Y, recognizedSNPs, mapSNPCausal_ind, tabCausalPost);
      calculatePOWER_FDR_clust_max(realSNPs, powerFDR, Pmi_Y, recognizedSNPs_piMass, mapSNPCausal_ind, tabCausalPost);
      cout << "POWER: " <<  powerFDR.POWER <<  ", FDR: " << powerFDR.FDR << ", FDcount: " << powerFDR.FDcount << endl;
      cout << "piMassModel: " << m << endl;
      sp_sort << powerFDR.POWER << "\t" << powerFDR.FDR << "\t" << powerFDR.FDcount << "\t";
      sp_sort << powerFDR.badSNP << "\t";
      sp_sort << powerFDR.posteriorBad << "\t";
      sp_sort << powerFDR.posteriorBadCluster << "\t";
      sp_sort <<  m.getMSC() << endl;
    }
    else
    {
//          ssModels << "causalModel:msc: 0.0, size: 0, []" << endl;
      cout << "causal model: 0 SNPs" << endl;      
      sp_sort << powerFDR.POWER << "\t" << powerFDR.FDR << "\t" << powerFDR.FDcount << "\t";
      sp_sort << powerFDR.badSNP << "\t";
      sp_sort << powerFDR.posteriorBad << "\t";
      sp_sort << powerFDR.posteriorBadCluster << "\t";
      sp_sort <<  0.0 << endl;
    }
  }
  else
  {
    cout << "casual snps: N/A" << endl;
    //causalModel:msc: 265.75508305, size: 3, [123932, 717206, 774597]
    //      ssModels << "causalModel: msc: 0.0, size: 0, []" << endl;
    /*    sp_sort << "\tPOWER: " << "NA";
     *    sp_sort << "\tFPD: " << "NA" << "\t";    
     *    sp_sort << "badSNP: " << "NA";
     *    sp_sort << "\tposteriori: " << "NA";
     *    sp_sort << "\tpost_cluster: " << "NA" << endl;
     */    
    sp_sort << powerFDR.POWER << "\t" << powerFDR.FDR << "\t" << powerFDR.FDcount << "\t";
    sp_sort << powerFDR.badSNP << "\t";
    sp_sort << powerFDR.posteriorBad << "\t";
    sp_sort << powerFDR.posteriorBadCluster << "\t";
    
    vector<snp_index_t> modelSNPs_max;
    bool clusterSum = true;
    calculate_clusters(Pmi_Y, modelSNPs_max, !clusterSum, "piMass_cluster_Max");
    
    Model *mx = new Model(data);
    if (modelSNPs_max.size() > 0)
    {
      mx->addManySNP(modelSNPs_max);
    }
    if (mx->computeRegression())
      mx->computeMSC();
    else
      mx->computeMSCfalseRegression();
    cout << "cluster_Pimass_model: " << *mx << endl;
    modelsFile << *mx << endl;
    modelsFile.close();
    sp_sort <<  mx->getMSC() << endl;
    //      sp_sort << "cluster_max: " << *mx << endl;
    delete mx;    
    
  }
  outFile << sp_sort.str();// << endl;    
  outFile.flush();
  outFile.close();
  ofstream file;      
  file.open((outFileName + "piMassModels_result_snp.txt").c_str(), ios::out);  
  for (map<int, unsigned int>::iterator it = snp.begin(); it != snp.end(); ++it)
    file << (*it).first << "\t" << (*it).second << endl;
  file.close();  
}


/**
 * @brief Calculates the statistics (POWER, FDR, False Discovery count) of the given SNPs. 
 * @param mySnp       - SNPs to test
 * @param realSNPs    - the causal SNPs
 * @param powerFDR    - returns the POWER, FDR and FD count (struct)
 * @param Pmi_Y       - the map of Pmi_Y value of the snps
 **/
/*
void GA::calculatePOWER_FDR_clust(set<snp_index_t> &mySnps, vector<snp_index_t> &causalSNPs, TPOWER_FDR &powerFDR, map<snp_index_t, long double> &mapSNPid_Pmi_Y, 
                                snp_index_t badSNP, map<snp_index_t, int> &recognizedSNPs, map<snp_index_t, int> &mapSNPCausal_ind,
                                vector< multiset<long double> > &tabCausalPost)
{
  bool TEST = false;
  
  assert(causalSNPs.size() > 0);
  assert(mySnps.size() > 0);
  
  int sumTP = 0;
  int sumFP = 0;
  set<snp_index_t> snp2check(mySnps.begin(), mySnps.end());
  // calcuates mapSNPid_Pmi_Y for the clusters
  int label;
  snp_index_t aSNP;
  
  for (map<snp_index_t, long double>::iterator it = mapSNPid_Pmi_Y.begin(); it != mapSNPid_Pmi_Y.end(); ++it)
  {
    aSNP = (*it).first;
    label = mapSNPid_label[ aSNP ];
    mapLabel_PmiY[label] += mapSNPid_Pmi_Y[aSNP];
  }
  
  stringstream ssLabel_PmiY;
  if (TEST) 
  {  
    ssLabel_PmiY << "label\tPmi_Y" << endl;
    for (map<int, long double>::iterator it =  mapLabel_PmiY.begin(); it != mapLabel_PmiY.end(); ++it)
    {
      ssLabel_PmiY << (*it).first << "\t" << (*it).second << endl;
    }
    ssLabel_PmiY << "snp_id\tlabel" << endl;
    for (set<snp_index_t>::iterator it = mySnps.begin(); it != mySnps.end(); ++it)
      ssLabel_PmiY << *it << "\t" << mapSNPid_label[*it] << endl;
  }
  
  ofstream file;
  
  if (TEST) 
  {
    file.open((parameter.out_file_name + "Label_PmiY.txt").c_str(), ios::trunc);
    file << ssLabel_PmiY.str() << endl;
    file.close();
  }  
  // the number of labels of mySnps
  stringstream ssmySnpsLabel;
  set<int> mySnpsLabel;
  for (set<snp_index_t>::iterator it = mySnps.begin(); it != mySnps.end(); ++it)
  {
    mySnpsLabel.insert( mapSNPid_label[*it] );  // calculates the number of cluster in the causal model
    if (TEST) ssmySnpsLabel << *it << " label " << mapSNPid_label[*it] << endl;
  }  
  
  if (TEST) 
  {
    for (set<int>::iterator it = mySnpsLabel.begin(); it != mySnpsLabel.end(); ++it)
    {
      ssmySnpsLabel << *it << " ";
    }
    file.open((parameter.out_file_name + "mySnpsLabel.txt").c_str(), ios::trunc);
    file << ssmySnpsLabel.str() << endl;
    file.close();
  }
  
  stringstream ssPower;
  
  bool isInCausal;
  bool isInResult;
  int ind;

  if (TEST) 
  {
    stringstream ssClusters;  
    for (map<int, long double>::iterator it =  mapLabel_PmiY.begin(); it != mapLabel_PmiY.end(); ++it)
    {
      ind = mapLabel_ind[ (*it).first ];
      ssClusters << "[" << (*it).first << "] (" << ind << ")=> {";
      for (set<snp_index_t>::iterator itSet = vectClust[ind].begin(); itSet != vectClust[ind].end(); ++itSet)
      {
        if (itSet == vectClust[ind].begin())
          ssClusters << *itSet;
        else
          ssClusters << " " << *itSet;
      }
      ssClusters << "}" << endl;
    }
    file.open((parameter.out_file_name + "myClusters.txt").c_str(), ios::trunc);
    file << ssClusters.str() << endl;
    file.close();
  }
  if (TEST) ssPower << snp2check << endl;
  snp_index_t recognizedSNP;
  for (map<int, long double>::iterator it =  mapLabel_PmiY.begin(); it != mapLabel_PmiY.end(); ++it)
  { // for each cluster
    isInCausal = false;
    isInResult = false;
    if ( (*it).second > 0.5 )  // if a posterior of the cluster is > 0.5
    {
      ind = mapLabel_ind[ (*it).first ];
      if (TEST) ssPower << "label: " << (*it).first << ", ind: " << ind << endl;
      for (set<snp_index_t>::iterator itSet = vectClust[ind].begin(); itSet != vectClust[ind].end() && isInCausal == false; ++itSet)
      { // for each snp in the cluster
        if (TEST) ssPower << "snp: " << *itSet;  
        vector<snp_index_t>::iterator itFind = find(causalSNPs.begin(), causalSNPs.end(), (*itSet));
        if (itFind != causalSNPs.end())
        { 
          isInCausal = true;
          if (TEST) ssPower << " " << "+c" << endl << "-------------" << endl;
          recognizedSNP = *itFind;
        }
        else
          if (TEST) ssPower << endl;
      } 
      if (isInCausal == true)
      {
        for (set<snp_index_t>::iterator itSet = vectClust[ind].begin(); itSet != vectClust[ind].end(); ++itSet)
        {
          if (TEST) ssPower << "GAsnp: " << *itSet;  
          set<snp_index_t>::iterator itFind = snp2check.find(*itSet);
          if (itFind != snp2check.end())
          { 
            isInResult = true;
            //snp2check.erase(*itSet);
            if (TEST) ssPower << " " << "+r" << endl;
          }
          else
            if (TEST) ssPower << " " << "-" << endl;
        }
      }
      if (isInCausal == true && isInResult == true)
      {
         ++sumTP;
         ++recognizedSNPs[recognizedSNP];
         for (set<snp_index_t>::iterator itSet = vectClust[ind].begin(); itSet != vectClust[ind].end(); ++itSet)
         {
           snp2check.erase(*itSet);
         }         
         int ind = mapSNPCausal_ind[recognizedSNP];
         long double Pmi = mapSNPid_Pmi_Y[recognizedSNP];
         tabCausalPost[ind].insert(Pmi);
         
      }
      if (TEST) ssPower << endl;
      if (TEST) ssPower << snp2check << endl;
    }
  }
  
  if (TEST) ssPower << endl << "reset: " << endl;
  // to co zostało w snp2check jest FP, liczymy klastry - the rest of snps is FP, calculates the number of the clusters
  set<int> rest;
  for (set<snp_index_t>::iterator it = snp2check.begin(); it != snp2check.end(); ++it)
  {
    rest.insert( mapSNPid_label[*it] );
    if (TEST) ssPower << "rest snp: " << *it << ", label: " << mapSNPid_label[*it] << endl;
  }  
  sumFP += rest.size();    
  
  
  set<int> causalLabels;
  for (vector<snp_index_t>::iterator it = causalSNPs.begin(); it != causalSNPs.end(); ++it)
  {
    causalLabels.insert( mapSNPid_label[*it] );  // calculates the number of cluster in the causal model
    //if (TEST) ssmySnpsLabel << *it << " label " << mapSNPid_label[*it] << endl;
  }
  
  powerFDR.POWER = (sumTP + 0.0) / causalLabels.size();
  if (sumFP + sumTP > 0)
  {
    powerFDR.FDR = (sumFP + 0.0) / (sumFP + sumTP);
  }
  else
  {
//    cerr << "Something wrong: sumFP + sumTP = 0!" << endl;
    powerFDR.FDR = 0;
  }  
  //powerFDR.FDR = (sumFP + 0.0) / mySnpsLabel.size();
  powerFDR.FDcount = sumFP;
  
  // looks for the best snp of the worst (posteriori)
  ssPower << "orig bad: " << badSNP << endl;
  if (snp2check.size() > 0)  // looks for bed in the set of rest 
  { // for bed
    badSNP = *(snp2check.begin());
    for (set<snp_index_t>::iterator it = snp2check.begin(); it != snp2check.end(); ++it)
    {
      if (it == snp2check.begin() ||  mapSNPid_Pmi_Y[*it] < mapSNPid_Pmi_Y[badSNP])
        badSNP = *it;
    }
  } 
  ssPower << "new bad: " << badSNP << endl;
  powerFDR.badSNP = badSNP;
  powerFDR.posteriorBad = mapSNPid_Pmi_Y[badSNP];
  label = mapSNPid_label[ badSNP ];
  powerFDR.posteriorBadCluster = mapLabel_PmiY[label];

  if (TEST) 
  {
    file.open((parameter.out_file_name + "_power.txt").c_str(), ios::trunc);
    file << ssPower.str() << endl;
    file.close();
  }  
}
*/

/** !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 * @brief Calculates the statistics (POWER, FDR, False Discovery count) of the given SNPs. 
 * @param mySnp       - SNPs to test
 * @param realSNPs    - the causal SNPs
 * @param powerFDR    - returns the POWER, FDR and FD count (struct)
 * @param Pmi_Y       - the map of Pmi_Y value of the snps
 **/
/*
void GA::calculatePOWER_FDR_clust(set<snp_index_t> &mySnps, vector<snp_index_t> &causalSNPs, TPOWER_FDR &powerFDR, map<snp_index_t, long double> &mapSNPid_Pmi_Y,
                                map<snp_index_t, int> &mapSNPCausal_ind, vector< multiset<long double> > &tabCausalPost_b)
{
  bool TEST = false;
  
  assert(causalSNPs.size() > 0);
  assert(mySnps.size() > 0);

  int sumTP = 0;
  int sumFP = 0;
  set<snp_index_t> snp2check(mySnps.begin(), mySnps.end());
  // calcuates mapSNPid_Pmi_Y for the clusters
  int label;
  snp_index_t aSNP;
  
  for (map<snp_index_t, long double>::iterator it = mapSNPid_Pmi_Y.begin(); it != mapSNPid_Pmi_Y.end(); ++it)
  {
    aSNP = (*it).first;
    label = mapSNPid_label[ aSNP ];
    mapLabel_PmiY[label] += mapSNPid_Pmi_Y[aSNP];
    
  }
  
  stringstream ssLabel_PmiY;
  if (TEST) 
  {  
    ssLabel_PmiY << "label\tPmi_Y" << endl;
    for (map<int, long double>::iterator it =  mapLabel_PmiY.begin(); it != mapLabel_PmiY.end(); ++it)
    {
      ssLabel_PmiY << (*it).first << "\t" << (*it).second << endl;
    }
    ssLabel_PmiY << "snp_id\tlabel" << endl;
    for (set<snp_index_t>::iterator it = mySnps.begin(); it != mySnps.end(); ++it)
      ssLabel_PmiY << *it << "\t" << mapSNPid_label[*it] << endl;
  }
  
  ofstream file;
  
  if (TEST) 
  {
    file.open((parameter.out_file_name + "Label_PmiY.txt").c_str(), ios::trunc);
    file << ssLabel_PmiY.str() << endl;
    file.close();
  }  
  // the number of labels of mySnps
  stringstream ssmySnpsLabel;
  set<int> mySnpsLabel;
  for (set<snp_index_t>::iterator it = mySnps.begin(); it != mySnps.end(); ++it)
  {
    mySnpsLabel.insert( mapSNPid_label[*it] );  // calculates the number of cluster in the causal model
    if (TEST) ssmySnpsLabel << *it << " label " << mapSNPid_label[*it] << endl;
  }  
  
  if (TEST) 
  {
    for (set<int>::iterator it = mySnpsLabel.begin(); it != mySnpsLabel.end(); ++it)
    {
      ssmySnpsLabel << *it << " ";
    }
    file.open((parameter.out_file_name + "mySnpsLabel.txt").c_str(), ios::trunc);
    file << ssmySnpsLabel.str() << endl;
    file.close();
  }
  
  stringstream ssPower;
  
  bool isInCausal,
       isInResult;
  int ind;

  if (TEST) 
  {
    stringstream ssClusters;  
    for (map<int, long double>::iterator it =  mapLabel_PmiY.begin(); it != mapLabel_PmiY.end(); ++it)
    {
      ind = mapLabel_ind[ (*it).first ];
      ssClusters << "[" << (*it).first << "] (" << ind << ")=> {";
      for (set<snp_index_t>::iterator itSet = vectClust[ind].begin(); itSet != vectClust[ind].end(); ++itSet)
      {
        ssClusters << *itSet << " ";
      }
      ssClusters << "}" << endl;
    }
    file.open((parameter.out_file_name + "myClusters.txt").c_str(), ios::trunc);
    file << ssClusters.str() << endl;
    file.close();
  }
  
  snp_index_t recognizedSNP;
  for (map<int, long double>::iterator it =  mapLabel_PmiY.begin(); it != mapLabel_PmiY.end(); ++it)
  { // for each cluster
    isInCausal = false;
    isInResult = false;
    if ( (*it).second > 0.5 )  // if a posterior of the cluster is > 0.5
    {
      ind = mapLabel_ind[ (*it).first ];
      if (TEST) ssPower << "label: " << (*it).first << ", ind: " << ind << endl;
      for (set<snp_index_t>::iterator itSet = vectClust[ind].begin(); itSet != vectClust[ind].end() && isInCausal == false; ++itSet)
      { // for each snp in the cluster
        if (TEST) ssPower << "snp: " << *itSet;  
        vector<snp_index_t>::iterator itFind = find(causalSNPs.begin(), causalSNPs.end(), (*itSet));
        if (itFind != causalSNPs.end())
        { 
          isInCausal = true;
          if (TEST) ssPower << " " << "+c" << endl;
          recognizedSNP = *itFind;
        }
        else
          if (TEST) ssPower << endl;
      } 
      if (isInCausal == true)
      {
        for (set<snp_index_t>::iterator itSet = vectClust[ind].begin(); itSet != vectClust[ind].end(); ++itSet)
        {
          if (TEST) ssPower << "snp: " << *itSet;  
          set<snp_index_t>::iterator itFind = snp2check.find(*itSet);
          if (itFind != snp2check.end())
          { 
            isInResult = true;
//            snp2check.erase(*itSet);
            if (TEST) ssPower << " " << "+r" << endl;
          }
          else
            if (TEST) ssPower << endl;
        }
      }
      if (isInCausal == true && isInResult == true)
      {
         for (set<snp_index_t>::iterator itSet = vectClust[ind].begin(); itSet != vectClust[ind].end(); ++itSet)
         {
           snp2check.erase(*itSet);
         }         
         ++sumTP;
         //++recognizedSNPs[recognizedSNP];
         int ind = mapSNPCausal_ind[recognizedSNP];
         long double Pmi = mapSNPid_Pmi_Y[recognizedSNP];
         tabCausalPost_b[ind].insert(Pmi);
      }
      ssPower << endl;
    }
  }
  
  if (TEST) ssPower << endl << "reset: " << endl;
  // to co zostało w snp2check jest FP, liczymy klastry - the rest of snps is FP, calculates the number of the clusters
  set<int> rest;
  for (set<snp_index_t>::iterator it = snp2check.begin(); it != snp2check.end(); ++it)
  {
    rest.insert( mapSNPid_label[*it] );
    if (TEST) ssPower << "snp: " << *it << ", label: " << mapSNPid_label[*it] << endl;
  }  
  sumFP += rest.size();    
  
  
  set<int> causalLabels;
  for (vector<snp_index_t>::iterator it = causalSNPs.begin(); it != causalSNPs.end(); ++it)
  {
    causalLabels.insert( mapSNPid_label[*it] );  // calculates the number of cluster in the causal model
    //if (TEST) ssmySnpsLabel << *it << " label " << mapSNPid_label[*it] << endl;
  }
  
  powerFDR.POWER = (sumTP + 0.0) / causalLabels.size();
  if (sumFP + sumTP > 0)
  {
    powerFDR.FDR = (sumFP + 0.0) / (sumFP + sumTP);
  }
  else
  {
//    cerr << "Something wrong: sumFP + sumTP = 0!" << endl;
    powerFDR.FDR = 0;
  }  
  //powerFDR.FDR = (sumFP + 0.0) / mySnpsLabel.size();
  powerFDR.FDcount = sumFP;
  
  if (TEST) 
  {
    file.open((parameter.out_file_name + "_power.txt").c_str(), ios::trunc);
    file << ssPower.str() << endl;
    file.close();
  }  
}
*/


/** MAX
 * @brief Calculates the statistics (POWER, FDR, False Discovery count) of the given SNPs. 
 * Oblicza POWER, FDR
 * 1. Oblicza prawd. poster. dla klastra jako max prawd. poster. jego SNPów
 * 2. Jeżeli w klastrze z prawd. post. > 0.5 znajduje się SNP przyczynowy, to mamy TP
 * 3. w przeciwnym razie FP
 * @param mySnp       - SNPs to test
 * @param realSNPs    - the causal SNPs
 * @param powerFDR    - returns the POWER, FDR and FD count (struct)
 * @param Pmi_Y       - the map of Pmi_Y value of the snps
 **/
void GA::calculatePOWER_FDR_clust_max(vector<snp_index_t> &causalSNPs, TPOWER_FDR &powerFDR, map<snp_index_t, long double> &mapSNPid_Pmi_Y, 
                                map<snp_index_t, int> &recognizedSNPs, map<snp_index_t, int> &mapSNPCausal_ind, vector< multiset<long double> > &tabCausalPost)
{
//  cout << "POWER, FDR - cluster MAX" << endl;
  bool TEST = false;//true;
  assert(causalSNPs.size() > 0);
  int sumTP = 0;
  int sumFP = 0;
  // calcuates mapSNPid_Pmi_Y for the clusters
  int label;
  snp_index_t aSNP;
  
  for (map<snp_index_t, long double>::iterator it = mapSNPid_Pmi_Y.begin(); it != mapSNPid_Pmi_Y.end(); ++it)
  {
    aSNP = (*it).first;
    label = mapSNPid_label[ aSNP ];
    if (mapSNPid_Pmi_Y[aSNP] > mapLabel_PmiY[label])
      mapLabel_PmiY[label] = mapSNPid_Pmi_Y[aSNP];
  }
  
  for (map<int, long double>::iterator it = mapLabel_PmiY.begin(); it != mapLabel_PmiY.end(); ++it)
  {
    if (it->second > 0.5)
    {
      if ( mapLabel_count.find( it->first) == mapLabel_count.end() )
        mapLabel_count[ it->first ] = 1;
      else
        ++mapLabel_count[ it->first ]; 
    }
  }

  stringstream ssLabel_PmiY;
  if (TEST) 
  {  
    ssLabel_PmiY << "label\tPmi_Y" << endl;
    for (map<int, long double>::iterator it =  mapLabel_PmiY.begin(); it != mapLabel_PmiY.end(); ++it)
    {
      ssLabel_PmiY << (*it).first << "\t" << (*it).second << endl;
    }
    ssLabel_PmiY << "snp_id\tlabel" << endl;
  }
  
  ofstream file;
  if (TEST) 
  {
    file.open((parameter.out_file_name + "Label_PmiY.txt").c_str(), ios::trunc);
    file << ssLabel_PmiY.str() << endl;
    file.close();
  }  
  stringstream ssPower;
  bool isInCausal;
  int ind;

  if (TEST) 
  {
    stringstream ssClusters;  
    for (map<int, long double>::iterator it =  mapLabel_PmiY.begin(); it != mapLabel_PmiY.end(); ++it)
    {
      ind = mapLabel_ind[ (*it).first ];
      ssClusters << "[" << (*it).first << "] (" << ind << ")=> {";
      for (set<snp_index_t>::iterator itSet = vectClust[ind].begin(); itSet != vectClust[ind].end(); ++itSet)
      {
        if (itSet == vectClust[ind].begin())
          ssClusters << *itSet;
        else
          ssClusters << " " << *itSet;
      }
      ssClusters << "}" << endl;
    }
    file.open((parameter.out_file_name + "myClusters.txt").c_str(), ios::trunc);
    file << ssClusters.str() << endl;
    file.close();
  }
  snp_index_t recognizedSNP;
  map<int, long double>::iterator itBad = mapLabel_PmiY.end();
  for (map<int, long double>::iterator it =  mapLabel_PmiY.begin(); it != mapLabel_PmiY.end(); ++it)
  { // for each cluster
    isInCausal = false;
    if ( (*it).second > 0.5 )  // if a posterior of the cluster is > 0.5
    {
      ind = mapLabel_ind[ (*it).first ];
      if (TEST) ssPower << "label: " << (*it).first << ", ind: " << ind << ", Pmi_Y: " << (*it).second << endl;
      for (set<snp_index_t>::iterator itSet = vectClust[ind].begin(); itSet != vectClust[ind].end() && isInCausal == false; ++itSet)
      { // for each snp in the cluster
        if (TEST) ssPower << "snp: " << *itSet;  
        vector<snp_index_t>::iterator itFind = find(causalSNPs.begin(), causalSNPs.end(), (*itSet));
        if (itFind != causalSNPs.end())
        { 
          isInCausal = true;
          ++sumTP;
          if (TEST) ssPower << " " << "+c " << sumTP << endl << "-------------" << endl;
          recognizedSNP = *itFind;
          
          ++recognizedSNPs[recognizedSNP];
          int indC = mapSNPCausal_ind[recognizedSNP];
          long double Pmi = mapSNPid_Pmi_Y[recognizedSNP];
          tabCausalPost[indC].insert(Pmi);
        }
        else
          if (TEST) ssPower << endl;
      } 
      
      if (isInCausal == false)
      {
        ++sumFP;
        if (TEST) ssPower << sumFP << " !!!" << endl;
        if (itBad == mapLabel_PmiY.end())
          itBad = it;
        else 
          if ((*itBad).second < (*it).second)
            itBad = it;
      }
      if (TEST) ssPower << endl;
    }
    else
    {
      if (itBad == mapLabel_PmiY.end())
        itBad = it;
      else
        if ( (*itBad).second < (*it).second )
          itBad = it;
    }
  }
  if (TEST) ssPower << endl << "reset: " << endl;
  set<int> causalLabels;
  for (vector<snp_index_t>::iterator it = causalSNPs.begin(); it != causalSNPs.end(); ++it)
  {
    causalLabels.insert( mapSNPid_label[*it] );  // calculates the number of cluster in the causal model
  }
  
  powerFDR.POWER = (sumTP + 0.0) / causalLabels.size();
  if (sumFP + sumTP > 0)
  {
    powerFDR.FDR = (sumFP + 0.0) / (sumFP + sumTP);
  }
  else
  {
//    cerr << "Something wrong: sumFP + sumTP = 0!" << endl;
    powerFDR.FDR = 0;
  }  
  powerFDR.FDcount = sumFP;
  
  // looks for the best snp of the worst (posteriori)
  int indBad = mapLabel_ind[ (*itBad).first ];
  snp_index_t badSNP;
  if (TEST) ssPower << "label_Bad: " << (*itBad).first << ", ind: " << ind << endl;
  if (vectClust[indBad].size() > 0)
  {
    badSNP = *(vectClust[indBad].begin());
  }
  else
  {
    cerr << "vectClust[ind].size() is equal 0!" << endl;
    exit(-1);
  }
  for (set<snp_index_t>::iterator it = vectClust[indBad].begin(); it != vectClust[indBad].end(); ++it)
  {
    if (mapSNPid_Pmi_Y[badSNP] < mapSNPid_Pmi_Y[*it])
      badSNP = *it;
  }
  ssPower << "badSNP: " << badSNP << endl;
  powerFDR.badSNP = badSNP;
  powerFDR.posteriorBad = mapSNPid_Pmi_Y[badSNP];
  label = mapSNPid_label[ badSNP ];
  powerFDR.posteriorBadCluster = mapLabel_PmiY[label];
  ssPower << "badSNPposterBad: " << powerFDR.posteriorBad << endl;
  ssPower << "badSNPposterClustBad: " << powerFDR.posteriorBadCluster << endl;
  ssPower << "POWER: " << powerFDR.POWER << endl;
  ssPower << "FDR: " << powerFDR.FDR << endl;
  ssPower << "TPcount: " << sumTP << endl;
  ssPower << "FDcount: " << powerFDR.FDcount << endl;
  
  if (TEST) 
  {
    file.open((parameter.out_file_name + "_" + int2str(parameter.in_values_int) + "_power.txt").c_str(), ios::trunc);
    file << ssPower.str() << endl;
    file.close();
  }  
}

/**
 * @brief Calculates the statistics (POWER, FDR, False Discovery count) of the given SNPs. 
 * Oblicza POWER, FDR
 * 1. Oblicza prawd. poster. dla klastra jako suma prawd. poster. jego SNPów
 * 2. Jeżeli w klastrze z prawd. post. > 0.5 znajduje się SNP przyczynowy, to mamy TP
 * 3. w przeciwnym razie FP
 * @param mySnp       - SNPs to test
 * @param realSNPs    - the causal SNPs
 * @param powerFDR    - returns the POWER, FDR and FD count (struct)
 * @param Pmi_Y       - the map of Pmi_Y value of the snps
 **/
void GA::calculatePOWER_FDR_clust_sum(vector<snp_index_t> &causalSNPs, TPOWER_FDR &powerFDR, map<snp_index_t, long double> &mapSNPid_Pmi_Y, 
                                map<snp_index_t, int> &recognizedSNPs, map<snp_index_t, int> &mapSNPCausal_ind, vector< multiset<long double> > &tabCausalPost)
{
  cout << "POWER, FDR - cluster SUM" << endl;
  bool TEST = false;
  assert(causalSNPs.size() > 0);
  int sumTP = 0;
  int sumFP = 0;
  // calcuates mapSNPid_Pmi_Y for the clusters
  int label;
  snp_index_t aSNP;
  
  for (map<snp_index_t, long double>::iterator it = mapSNPid_Pmi_Y.begin(); it != mapSNPid_Pmi_Y.end(); ++it)
  {
    aSNP = (*it).first;
    label = mapSNPid_label[ aSNP ];
    mapLabel_PmiY[label] += mapSNPid_Pmi_Y[aSNP];
  }
  
  for (map<int, long double>::iterator it = mapLabel_PmiY.begin(); it != mapLabel_PmiY.end(); ++it)
  {
    if (it->second > 0.5)
    {
      if ( mapLabel_count.find( it->first) == mapLabel_count.end() )
        mapLabel_count[ it->first ] = 1;
      else
        ++mapLabel_count[ it->first ]; 
    }
  }
  
  stringstream ssLabel_PmiY;
  if (TEST) 
  {  
    ssLabel_PmiY << "label\tPmi_Y" << endl;
    for (map<int, long double>::iterator it =  mapLabel_PmiY.begin(); it != mapLabel_PmiY.end(); ++it)
    {
      ssLabel_PmiY << (*it).first << "\t" << (*it).second << endl;
    }
    ssLabel_PmiY << "snp_id\tlabel" << endl;
  }
  
  ofstream file;
  if (TEST) 
  {
    file.open((parameter.out_file_name + "Label_PmiY.txt").c_str(), ios::trunc);
    file << ssLabel_PmiY.str() << endl;
    file.close();
  }  
  stringstream ssPower;
  bool isInCausal;
  int ind;

  if (TEST) 
  {
    stringstream ssClusters;  
    for (map<int, long double>::iterator it =  mapLabel_PmiY.begin(); it != mapLabel_PmiY.end(); ++it)
    {
      ind = mapLabel_ind[ (*it).first ];
      ssClusters << "[" << (*it).first << "] (" << ind << ")=> {";
      for (set<snp_index_t>::iterator itSet = vectClust[ind].begin(); itSet != vectClust[ind].end(); ++itSet)
      {
        if (itSet == vectClust[ind].begin())
          ssClusters << *itSet;
        else
          ssClusters << " " << *itSet;
      }
      ssClusters << "}" << endl;
    }
    file.open((parameter.out_file_name + "myClusters.txt").c_str(), ios::trunc);
    file << ssClusters.str() << endl;
    file.close();
  }
  snp_index_t recognizedSNP;
  map<int, long double>::iterator itBad = mapLabel_PmiY.end();
  for (map<int, long double>::iterator it =  mapLabel_PmiY.begin(); it != mapLabel_PmiY.end(); ++it)
  { // for each cluster
    isInCausal = false;
    if ( (*it).second > 0.5 )  // if a posterior of the cluster is > 0.5
    {
      ind = mapLabel_ind[ (*it).first ];
      if (TEST) ssPower << "label: " << (*it).first << ", ind: " << ind << ", Pmi_Y: " << (*it).second << endl;
      for (set<snp_index_t>::iterator itSet = vectClust[ind].begin(); itSet != vectClust[ind].end() && isInCausal == false; ++itSet)
      { // for each snp in the cluster
        if (TEST) ssPower << "snp: " << *itSet;  
        vector<snp_index_t>::iterator itFind = find(causalSNPs.begin(), causalSNPs.end(), (*itSet));
        if (itFind != causalSNPs.end())
        { 
          isInCausal = true;
          ++sumTP;
          if (TEST) ssPower << " " << "+c " << sumTP << endl << "-------------" << endl;
          recognizedSNP = *itFind;
          
          ++recognizedSNPs[recognizedSNP];
          int indC = mapSNPCausal_ind[recognizedSNP];
          long double Pmi = mapSNPid_Pmi_Y[recognizedSNP];
          tabCausalPost[indC].insert(Pmi);
        }
        else
          if (TEST) ssPower << endl;
      } 
      
      if (isInCausal == false)
      {
        ++sumFP;
        if (TEST) ssPower << sumFP << " !!!" << endl;
        if (itBad == mapLabel_PmiY.end())
          itBad = it;
        else 
          if ((*itBad).second < (*it).second)
            itBad = it;
      }
      if (TEST) ssPower << endl;
    }
    else
    {
      if (itBad == mapLabel_PmiY.end())
        itBad = it;
      else
        if ( (*itBad).second < (*it).second )
          itBad = it;
    }
  }
  
  if (TEST) ssPower << endl << "reset: " << endl;
  set<int> causalLabels;
  for (vector<snp_index_t>::iterator it = causalSNPs.begin(); it != causalSNPs.end(); ++it)
  {
    causalLabels.insert( mapSNPid_label[*it] );  // calculates the number of cluster in the causal model
  }
  
  powerFDR.POWER = (sumTP + 0.0) / causalLabels.size();
  if (sumFP + sumTP > 0)
  {
    powerFDR.FDR = (sumFP + 0.0) / (sumFP + sumTP);
  }
  else
  {
//    cerr << "Something wrong: sumFP + sumTP = 0!" << endl;
    powerFDR.FDR = 0;
  }  
  powerFDR.FDcount = sumFP;
  
  // looks for the best snp of the worst (posteriori)
  int indBad = mapLabel_ind[ (*itBad).first ];
  snp_index_t badSNP;
  if (TEST) ssPower << "label_Bad: " << (*itBad).first << ", ind: " << ind << endl;
  if (vectClust[indBad].size() > 0)
  {
    badSNP = *(vectClust[indBad].begin());
  }
  else
  {
    cerr << "vectClust[ind].size() is equal 0!" << endl;
    exit(-1);
  }
  for (set<snp_index_t>::iterator it = vectClust[indBad].begin(); it != vectClust[indBad].end(); ++it)
  {
    if (mapSNPid_Pmi_Y[badSNP] < mapSNPid_Pmi_Y[*it])
      badSNP = *it;
  }
  ssPower << "badSNP: " << badSNP << endl;
  powerFDR.badSNP = badSNP;
  powerFDR.posteriorBad = mapSNPid_Pmi_Y[badSNP];
  label = mapSNPid_label[ badSNP ];
  powerFDR.posteriorBadCluster = mapLabel_PmiY[label];
  ssPower << "badSNPposterBad: " << powerFDR.posteriorBad << endl;
  ssPower << "badSNPposterClustBad: " << powerFDR.posteriorBadCluster << endl;
  ssPower << "POWER: " << powerFDR.POWER << endl;
  ssPower << "FDR: " << powerFDR.FDR << endl;
  ssPower << "TPcount: " << sumTP << endl;
  ssPower << "FDcount: " << powerFDR.FDcount << endl;
  
  cout << "file_BEFORE: " << (parameter.out_file_name + "_" + int2str(parameter.in_values_int) + "_power.txt").c_str() << endl;
  if (TEST) 
  {
    file.open((parameter.out_file_name + "_" + int2str(parameter.in_values_int) + "_power.txt").c_str(), ios::trunc);
    file << ssPower.str() << endl;
    file.close();
  }  
}


///////// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

/**
 * @brief Calculates the statistics 
 * Oblicza 
 * 1. dla clusterSum = true - oblicza prawd. poster. dla klastra jako suma prawd. poster. jego SNPów
 *    dla clusterSum = false - oblicza prawd. poster. dla klastra jako wartość maksymalna prawd. poster. jego SNPów
 * @param mySnp       - SNPs to test
 * @param realSNPs    - the causal SNPs
 * @param powerFDR    - returns the POWER, FDR and FD count (struct)
 * @param Pmi_Y       - the map of Pmi_Y value of the snps
 * UWAGA clusterSum == otwiera plik *_recognized w trybie append
 **/
void GA::calculate_clusters(map<snp_index_t, long double> &mapSNPid_Pmi_Y, vector<snp_index_t> &modelSNPs, bool clusterSum, string method_Name)
{
  set<int>clusterSet;
  modelSNPs.clear();
  bool TEST = false;
  int label;
  snp_index_t aSNP;

  for (map<snp_index_t, long double>::iterator it = mapSNPid_Pmi_Y.begin(); it != mapSNPid_Pmi_Y.end(); ++it)
  {
    aSNP = (*it).first;
    label = mapSNPid_label[ aSNP ];
    mapLabel_PmiY[label] = 0.0;
  }
  // calcuates mapSNPid_Pmi_Y for the clusters    // Obliczamy prawd. a posteriori dla klastra metodą: SUM (suma prawd. a posteriori SNPów z tego klastra lub jako wartość MAX)
  if (clusterSum == true)
  {
    for (map<snp_index_t, long double>::iterator it = mapSNPid_Pmi_Y.begin(); it != mapSNPid_Pmi_Y.end(); ++it)
    {
      aSNP = (*it).first;
      label = mapSNPid_label[ aSNP ];
      mapLabel_PmiY[label] += mapSNPid_Pmi_Y[aSNP];
    }
  }
  else
  {
    for (map<snp_index_t, long double>::iterator it = mapSNPid_Pmi_Y.begin(); it != mapSNPid_Pmi_Y.end(); ++it)
    {
      aSNP = (*it).first;
      label = mapSNPid_label[ aSNP ];
      if (mapSNPid_Pmi_Y[aSNP] > mapLabel_PmiY[label])
        mapLabel_PmiY[label] = mapSNPid_Pmi_Y[aSNP];
    }
  }
  int ind;  
  stringstream ssClusters;  
  stringstream ssSNPs_ID;
//  ssClusters << method_Name << endl;
//  ssSNPs_ID << method_Name << endl;
  // data.getSNP( (*itM).second ).getSnpId()
  ofstream file;
  for (map<int, long double>::iterator it = mapLabel_PmiY.begin(); it != mapLabel_PmiY.end(); ++it)
  { // Counts the labels of clusters for many runs of GA.  // zlicza klastry dla wielu uruchomień GA
    if (it->second > 0.5)
    {
      if ( mapLabel_count.find( it->first) == mapLabel_count.end() )
        mapLabel_count[ it->first ] = 1;
      else
        ++mapLabel_count[ it->first ]; 
      
      ind = mapLabel_ind[ (*it).first ];
      ssClusters << (*it).first << "\t>>\t";
      ssSNPs_ID << (*it).first << "\t>>\t";

      for (set<snp_index_t>::iterator itS = vectClust[ind].begin(); itS != vectClust[ind].end(); ++itS)
      {
        if (itS == vectClust[ind].begin())
        {
          aSNP = *itS;
          ssClusters << *itS;
          ssSNPs_ID << data.getSNP( *itS ).getSnpId();
        }
        else
        {
          if (mapSNPid_Pmi_Y[ *itS ] > mapSNPid_Pmi_Y[ aSNP ])
            aSNP = *itS;
          ssClusters << " " << *itS;
          ssSNPs_ID << " " << data.getSNP( *itS ).getSnpId();
        }
      
      }
      if (vectClust[ind].size() > 0)
        modelSNPs.push_back(aSNP);
      ssClusters << endl;
      ssSNPs_ID << endl;
    }
  }  
  Model model(data);
  if (modelSNPs.size() > 0)
    model.addManySNP(modelSNPs);
  if (model.computeRegression())
    model.computeMSC();
  else
    model.computeMSCfalseRegression();
  if (isCluster == true)
    saveRecognizedSNPinfo(modelSNPs, method_Name, mapSNPid_Pmi_Y);
  
  stringstream ssLabel_PmiY;
  if (TEST) 
  {  
    ssLabel_PmiY << "cluster label\tPmi_Y" << endl;
    for (map<int, long double>::iterator it =  mapLabel_PmiY.begin(); it != mapLabel_PmiY.end(); ++it)
    {
      ssLabel_PmiY << (*it).first << "\t" << (*it).second << endl;
    }
  }
  if (TEST) 
  {
    file.open((parameter.out_file_name +  "_Label_PmiY.txt").c_str(), ios::trunc);
    file << ssLabel_PmiY.str() << endl;
    file.close();
  }  
  stringstream ssPower;
  if (TEST) 
  {
    stringstream ssClusters;  
    for (map<int, long double>::iterator it =  mapLabel_PmiY.begin(); it != mapLabel_PmiY.end(); ++it)
    {
      ind = mapLabel_ind[ (*it).first ];
      ssClusters << ind << ")\t" <<  (*it).first << "\t => {";
      for (set<snp_index_t>::iterator itSet = vectClust[ind].begin(); itSet != vectClust[ind].end(); ++itSet)
      {
        if (itSet == vectClust[ind].begin())
          ssClusters << *itSet;
        else
          ssClusters << " " << *itSet;
      }
      ssClusters << "}" << endl;
    }
    file.open((parameter.out_file_name + "_myClusters.txt").c_str(), ios::trunc);
    file << ssClusters.str() << endl;
    file.close();
  }
}


/**
 * @brief Local improvement of a given model
 * @param model - improvement model
 * @param threshold - minimal correlation value
 * @param correlationRange - correlation for snp is compute for snps from range [snp - correlationRange, snp + correlationRange]
 */
void GA::localImprovement(Model &model, double threshold, int correlationRange) 
{
  const bool TEST = false;
  
  if (model.computeRegression())
    model.computeMSC();
  else
    model.computeMSCfalseRegression();
  if (TEST) cout << "A model to improve: " << model << endl;
  vector <snp_index_t> snps = model.getModelSnps();  // snps of the model
  set<snp_index_t> newSnps;   // snps of a new model
  
  vector<snp_index_t> v; 
  int label,
      ind,
      randVal;
  snp_index_t aSNP;
  set<snp_index_t>::iterator itSNP;            
  int newModelNo = 3 * snps.size();  // the number of the new models
  Model theBestModel(data);
  theBestModel = model;  // zamienić na vector i zmienną msc !!!!!!!!!!!!!!
  Model *newModel;
  for (int i = 0; i < newModelNo; ++i)
  {
    newModel = new Model(data);
    newSnps.clear();
    for (unsigned int s = 0; s < snps.size(); ++s)  // mades a new model
    { // I exchange a snp aSNP for a randomly choosen snps from the cluster
      aSNP = snps[s];
      label = mapSNPid_labelLO[ aSNP ];
      ind = mapLabel_ind[ label ];
      randVal = random() % vectClust[ind].size();
      if (TEST) cout << "ind:" << ind << " randVal: " << randVal << endl << aSNP << " -> cluster " << vectClust[ind] << endl;
      for (itSNP = vectClust[ind].begin(); itSNP != vectClust[ind].end() && randVal > 0; ++itSNP, --randVal)
      {
        ;
        //if (TEST) cout << "looking... " << *itSNP << ", randVal: " << randVal << endl;
      }
      assert(itSNP != vectClust[ind].end());
      if (TEST) cout << "a new snp  " << *itSNP << endl;
      newSnps.insert(*itSNP);
    }
    v.assign(newSnps.begin(), newSnps.end());
    if (v.size() > 0)
      newModel->addManySNP(v);
    if (newModel->computeRegression())
      newModel->computeMSC();
    else
      newModel->computeMSCfalseRegression();
    if (newModel->getMSC() < theBestModel.getMSC() || i == 0)
      theBestModel = *newModel;
    // to write a model to the pool  
    toPool(newModel, 'L');
    writePoolofSize(p_30K, 30000);
    writePoolofSize(p_100K, 100000);
    if (TEST) cout << "a new snps  " << newSnps << endl;
    if (TEST) cout << "a new model " << *newModel << endl << endl;
    delete newModel;
  }
  if (TEST) cout << "przed: " << model <<  endl << endl;
  if (theBestModel.getMSC() < model.getMSC())
    model = theBestModel;
  if (TEST) cout << "po: " << model <<  endl << endl;
  if (TEST) cout << "A model to return " << model << endl;
}

void GA::makeVectCluster(map<snp_index_t, int>& mSNPid_label)
{
  clusterLabel.clear();
  mapLabel_ind.clear();
  vectClust.clear();
  
  
  if (mSNPid_label.size() == 0)
  {
    cerr << "No clusters (makeVectCluster(...))" << endl;
    exit(-1);
  }
  set<int>labels;
  for (map<snp_index_t, int>::iterator it =  mSNPid_label.begin(); it != mSNPid_label.end(); ++it)
  {
    labels.insert((*it).second);
  }
  unsigned int i = 0;
  
  clusterLabel.resize(labels.size());
  for (set<int>::iterator it = labels.begin(); it != labels.end(); ++it, ++i)
  {
    mapLabel_ind.insert(pair<int, int>(*it, i));
    clusterLabel[i] = *it;
  }
  vectClust.resize(labels.size()); 
  for (map<snp_index_t, int>::iterator it =  mSNPid_label.begin(); it != mSNPid_label.end(); ++it)
  {
    vectClust[ mapLabel_ind[(*it).second] ].insert( (*it).first );
  }
  bool writeToFile = false;
  if (writeToFile == true)
  {
    stringstream ssClusters;
    for (i = 0; i < vectClust.size(); ++i)
    {
      ssClusters << "[" << clusterLabel[i] << "]-> {";
                 //<< vectClust[i] << endl;
      for (set<snp_index_t>::iterator it = vectClust[i].begin(); it != vectClust[i].end(); ++it)
      {
         if (it == vectClust[i].begin())
           ssClusters << data.getSNP( *it ).getSnpId();
         else
           ssClusters << " " << data.getSNP( *it ).getSnpId();
      }
      ssClusters << "}" << endl;
    }
    ofstream file;
    file.open("myClusters.txt", ios::trunc);
    file << ssClusters.str() << endl;
    file.close();
  }
}

/**
 * @brief Reads a file of clusters (of format: SNP_id, label) and prepares 2 maps: SNP_id -> label, and label - > Pmi_Y (the sum of the posteriors)
 * @param mapSNPid_label - reference to the map SNP_id -> label
 * 
 * readClusterLabel(mapSNPid_label, clusterFile_POWER, SNP_Names_Id, &mapLabel_PmiY);
 */ 
void GA::readClusterLabel(map<snp_index_t, int>& mapSNPid_label, const string& fileName, map<string, int>& SNP_Names_Id, map<int, long double>* mapLabel_PmiY)
{
  stringstream ssError;
  ifstream aFile;
  string SNPname;
  int label;
  set<int> labels;
  aFile.open( fileName.c_str(), ios::in );
  aFile.clear();
  int badName = 0;
  if ( aFile.is_open())
  {

    while (!aFile.eof())
    {
      aFile >> SNPname;
      aFile >> label;
      if (SNP_Names_Id.find(SNPname) == SNP_Names_Id.end())
      {
        cerr << "The set of SNPs from file \"" << fileName << "\"  and from the data set is not the same.\r\nThe SNP \"" 
             << SNPname << "\" is not in the data set." << endl;
        ssError << SNPname << endl;
        ++badName;
      }
      else
      {
        if (mapLabel_PmiY != 0)
        {
          labels.insert(label);
          mapLabel_PmiY->insert(pair<int, long double>(label, 0.0L));
        }
        mapSNPid_label.insert(pair<snp_index_t, int>(SNP_Names_Id[SNPname], label));
      }
    }
    aFile.close();
  }
  else
  {
    printLOG( "No clusters file" );
    return;
  }
  if (mapLabel_PmiY != 0)
  {
    for(set<int>::iterator it = labels.begin(); it != labels.end(); ++it)
    {
      mapLabel_PmiY->insert(pair<int, long double>(*it, 0.0L));
    }
    
  }
  ofstream errFile;
  errFile.open((parameter.out_file_name +  fileName + "_missing.txt").c_str());
  errFile << "bad names: " << badName << endl << "missing SNPs:\r\n " << ssError.str();
  errFile.close();
  isCluster = true;
}

/** @brief creates map snp_i -> index
 * 
 */
void GA::initCausalPost( map<snp_index_t, int> &mapSNPCausal_ind )
{
  if (parameter.causalModelFilename.length() > 0)  // The simulation. We have the causalModel
  {
    vector<snp_index_t> realSNPs;                // snps from the casual model
    modelReader(parameter.causalModelFilename, realSNPs);
    if (realSNPs.size() > 0)
    {
      int i = 0;
      for (vector<snp_index_t>::iterator it = realSNPs.begin(); it != realSNPs.end(); ++it, ++i)
        mapSNPCausal_ind.insert( pair<snp_index_t, int>(*it, i) );
    }
    else
    {
      cerr << "initCausalPost(): the size of real model is 0!" << endl;
      exit(-1);
    }
  }
  else
  {
    cerr << "initCausalPost(): no real model! causalModelFilename: " << parameter.causalModelFilename << endl;
    //exit(-1);
  }
}

void writePosterior(string fileName, map<snp_index_t, int> &mapSNPCausal_ind, vector< multiset<long double> > &tabCausalPost, int size)
{
  bool TEST = false;
  stringstream ssClusters;
  int id = 1;
  for (map<snp_index_t, int>::iterator it = mapSNPCausal_ind.begin(); it != mapSNPCausal_ind.end(); ++it)
  {
    if (TEST) cout << "%" << (*it).first << endl;
    ssClusters << "%" << (*it).first << endl;
    ssClusters << "I(" << id << ", : ) = [";
    int i = 0;
    int ind = (*it).second;
    for (multiset<long double>::iterator itC = tabCausalPost[ind].begin(); itC != tabCausalPost[ind].end(); ++itC, ++i)
    {
      ssClusters << *itC << " ";
      if (TEST) cout << *itC << " ";
    }
    if (TEST) cout << endl << "i = " << i << endl;
    for (; i < size; ++i)
      ssClusters << 0 << " ";
    ssClusters << "];" << endl;
    ++id;
  }
  ofstream file;
  file.open(fileName.c_str(), ios::trunc);
  file << ssClusters.str() << endl;
  file.close();
}

/**
 * 
 */
void GA::makeCausalClasters(string fileName, vector<snp_index_t> causalSNPs)
{
  long double abscor;
  map<snp_index_t, int> mapSNPid_ind;
  vector< set<snp_index_t> > tab_Corr;
  assert(causalSNPs.size() > 0);
  tab_Corr.resize(causalSNPs.size());
  for (unsigned int i = 0; i < causalSNPs.size(); ++i)  // makes a map SNP_id -> index
  {
    mapSNPid_ind.insert(pair<snp_index_t, int> (causalSNPs[i], i));
  }
  int ind;
  bool isCorrelated;
  stringstream ssLabels;
  for (unsigned int j = 0; j < data.getSnpNo(); ++j)
  {
    isCorrelated = false;
    for (unsigned int i = 0; i < causalSNPs.size(); ++i)
    {
      abscor= fabs( data.computeCorrelation( causalSNPs[i], j ) ); // compute correlation between model SNP and SNP j 
      if (abscor > 0.5)
      {
        ind = mapSNPid_ind[ causalSNPs[i] ];
        tab_Corr[ind].insert(j);
        if (isCorrelated == true)
        {
          cerr << causalSNPs[i] << " <---> " <<  j << ":" << abscor << endl;
        }
        isCorrelated = true;
        ssLabels <<  data.getSNP(j).getSnpId() << "\t" << i << endl;
      }  
    }  
    if (isCorrelated == false)
    {
      ssLabels <<  data.getSNP(j).getSnpId() << "\t" << causalSNPs.size() << endl;
    }
  }
  stringstream ssClusters;
  for (unsigned int i = 0; i < causalSNPs.size(); ++i)  // makes a map SNP_id -> index
  {
     ssClusters << "[" << causalSNPs[i] << "] -> {";
     ind = mapSNPid_ind[ causalSNPs[i] ];
     for (set<snp_index_t>::iterator it = tab_Corr[ind].begin(); it != tab_Corr[ind].end(); ++it)
     {
       if (it == tab_Corr[ind].begin())
         ssClusters << *it;
       else
         ssClusters << ", " << *it;
     }
     ssClusters << "}" << endl;
  }
  ofstream file;
  file.open((fileName + "_causalLabelAG").c_str(), ios::trunc);
  file << ssLabels.str() << endl;
  file.close();
  
  file.open((fileName + "_causalClastersAG").c_str(), ios::trunc);
  file << ssClusters.str() << endl;
  file.close();
}

/**
 * 
 */
void GA::saveLabelCount(const string &fileName)
{

  stringstream ss;
  ss << "label\tcount\tlabel\tSNPs" << endl;
  for (map<int, int>::iterator it = mapLabel_count.begin(); it != mapLabel_count.end(); ++it)
  {
    if (it->second > 0)
    {
      ss << (*it).first << "\t" << (*it).second << endl;
      int ind = mapLabel_ind[it->first];
      ss << clusterLabel[ind] << "\t{";
      for (set<snp_index_t>::iterator itSet = vectClust[ind].begin(); itSet != vectClust[ind].end(); ++itSet)
      {
        if (itSet == vectClust[ind].begin())
          ss << *itSet;
        else
          ss << ", " << *itSet;
      } 
      ss << "}" << endl;
    }   
  }
  ofstream  outFile;
  outFile.exceptions ( ofstream::eofbit | ofstream::failbit | ofstream::badbit ); // checks if poolfile can be written
  try
  {
    outFile.open((parameter.out_file_name + fileName).c_str(),  fstream::out | fstream::trunc );
    outFile << ss.str() << endl;
    outFile.flush();
    outFile.close();
  }
  catch (ofstream::failure e)
  {
    cerr << "Could not write a piMass reco file" <<endl;
    exit(-1);
  }
}

/**
 * @brief Initilizes the map (SNPs -> count) with 0
 * @param map_SNP2count - the map to initialize
 */
void GA::setCausalCount(map<snp_index_t, int> &map_SNP2count)
{
  vector<snp_index_t> causalSNPs;             // snps from the real model
  modelReader(parameter.causalModelFilename, causalSNPs);  
  for (unsigned int i = 0; i < causalSNPs.size(); ++i)
  {
     map_SNP2count[ causalSNPs[i] ] = 0;    
  }
}


/**
 * @brief Computes correlation for snps of causal model
 * @param threshold - threshold value of correlation. 
 * @return vector of snps with a correlation above threshold
 */
void GA::stronglyCorrelatedSnpsCluster(const vector<snp_index_t> & tabSNPs, const double& threshold)
{
  stringstream ss;
  
//  int tabCl[] = // 20SNPs {46, 448, 884, 1275, 1642, 2013, 2416, 2749, 3043, 3352, 3712, 4096, 4448, 4848, 5185, 5570, 5960, 6378, 6764, 7254};
  //{25, 314, 592, 846, 1107, 1388, 1618, 1822, 2093, 2366, 2592, 2818, 3006, 3238, 3452, 3675, 3926, 4193, 4414, 4691, 4921, 5153, 5393, 5675, 5939, 6199, 6484, 6726, 7003, 7280};
//  {242029, 302459};

//  vector<snp_index_t> causalSNPs;             // snps from real model
//  modelReader(parameter.causalModelFilename, causalSNPs);  

  stringstream ssAll;                 // to save the output
  ssAll << "causalSNP\tclusterSNP\tcorrelation" << endl;
  stringstream ssClusters;
  multimap<double, int> StrongCor; // to sort the correlated SNPs

  double abscor;                   // abs|Correlation| of two SNPs
  unsigned int j;                           
  
  
  //  697067,  599750
  
  vector<snp_index_t> tabSpr;
  tabSpr.push_back(130633);
  tabSpr.push_back(267614);
  tabSpr.push_back(267621);
  tabSpr.push_back(267637);

  for (unsigned int ind = 0; ind < tabSNPs.size(); ++ind)
  {
    cout << "calculation for " << tabSNPs[ind] << " <<>> " << data.getSNP( tabSNPs[ind] ).getSnpId() << endl;
    for (j = 0; j < tabSpr.size(); j++)
    {
      abscor = fabs( data.computeCorrelation( tabSNPs[ind], tabSpr[j] ) ); // compute correlation between model SNP and SNP j
      if (abscor  >= threshold) // add if correlation is big enough
      {
//        ssAll << causalSNP[ind] << "\t" << j << "\t" << abscor << endl; data.getSNP( (*itM).second ).getSnpId()
        ssAll << data.getSNP( tabSNPs[ind] ).getSnpId() << "\t" << data.getSNP( tabSpr[j] ).getSnpId() << "\t" << abscor << "\t" << endl; // tabCl[ind] << endl;
        cout << data.getSNP( tabSNPs[ind] ).getSnpId() << "\t" << data.getSNP( tabSpr[j] ).getSnpId() << "\t" << abscor << "\t" << endl; // tabCl[ind] << endl;
/*        if (j == rangeFrom)
          ssClusters << endl << data.getSNP( tabSNPs[ind] ).getSnpId() << "\t" << tabCl[ind] << endl;
        ssClusters << data.getSNP( j ).getSnpId() << "\t" << tabCl[ind] << endl;
*/        
      }
    }
    ssAll << endl;
    cout << endl;
  }
  /*
  unsigned int rangeFrom = 0,
               rangeTo = data.getSnpNo();
  // search for strongly correlated SNPs
  for (unsigned int ind = 0; ind < tabSNPs.size(); ++ind)
  {
    cout << "calculation for " << tabSNPs[ind] << " <<>> " << data.getSNP( tabSNPs[ind] ).getSnpId() << endl;
    for (j = rangeFrom; j < rangeTo; j++)
    {
      abscor = fabs( data.computeCorrelation( tabSNPs[ind], j ) ); // compute correlation between model SNP and SNP j
      if (abscor  >= threshold) // add if correlation is big enough
      {
//        ssAll << causalSNP[ind] << "\t" << j << "\t" << abscor << endl; data.getSNP( (*itM).second ).getSnpId()
        ssAll << data.getSNP( tabSNPs[ind] ).getSnpId() << "\t" << data.getSNP( j ).getSnpId() << "\t" << abscor << "\t" << endl; // tabCl[ind] << endl;
  //      if (j == rangeFrom)
     //    ssClusters << endl << data.getSNP( tabSNPs[ind] ).getSnpId() << "\t" << tabCl[ind] << endl;
    //    ssClusters << data.getSNP( j ).getSnpId() << "\t" << tabCl[ind] << endl;
        
      }
    }
    ssAll << endl;
  }
  */
  /*
  abscor = fabs( data.computeCorrelation( 13, 15 ) ); 
  ssAll << endl << data.getSNP( 13 ).getSnpId() << "\t" << data.getSNP( 15 ).getSnpId() << "\t" << abscor << endl;
  abscor = fabs( data.computeCorrelation( 13, 16 ) );
  ssAll << data.getSNP( 13 ).getSnpId() << "\t" << data.getSNP( 16 ).getSnpId() << "\t" << abscor << endl;
  abscor = fabs( data.computeCorrelation( 7207, 7279 ) );
  ssAll << endl << data.getSNP( 7207 ).getSnpId() << "\t" << data.getSNP( 7279 ).getSnpId() << "\t" << abscor << endl;
  abscor = fabs( data.computeCorrelation( 15616, 15622 ) );
  ssAll << endl << data.getSNP( 15616 ).getSnpId() << "\t" << data.getSNP( 15622 ).getSnpId() << "\t" << abscor << endl;
*/

  
  ofstream  outFile;
  outFile.exceptions ( ofstream::eofbit | ofstream::failbit | ofstream::badbit ); // checks if poolfile can be written
  try
  {
    outFile.open((parameter.out_file_name + "_corrAG_Example.txt").c_str(),  fstream::out | fstream::trunc );
    outFile << ssAll.str() << endl;
    outFile.flush();
    outFile.close();
  }
  catch (ofstream::failure e)
  {
    cerr << "Could not write a _corrAG file" <<endl;
    exit(-1);
  }
  
  try
  {
    outFile.open((parameter.out_file_name + "_cluAG_30.clu").c_str(),  fstream::out | fstream::trunc );
    outFile << ssClusters.str() << endl;
    outFile.flush();
    outFile.close();
  }
  catch (ofstream::failure e)
  {
    cerr << "Could not write a _cluAG file" <<endl;
    exit(-1);
  }
}  

//SNP_A-1828353 SNP_A-2157434 SNP_A-1794641 SNP_A-1815281 SNP_A-2202441 SNP_A-2309459 SNP_A-2160092 SNP_A-1850477 SNP_A-2289125 SNP_A-1829559 SNP_A-2208065 SNP_A-1786242


void GA::saveRecognizedSNPinfo(const vector<snp_index_t> &v, const string &method_Name, map<snp_index_t, long double> &mapSNPid_Pmi_Y)
{
  stringstream ssClusters;  
  stringstream ssSNPs_ID;  
  
  TSNP_Info aSNP_info;
  set<TSNP_Info> aRegion;
  stringstream ssInfo;
  
  ssInfo << setw(7) << "Nr" << "\t" << setw(15) << "SNP Id" << "\t" << setw(4) << "Chr" << "\t" << setw(12) << "Pos" << "\t" << setw(12) << "posteriori" << "\t" 
         << endl << endl;  
  double sum, max;
/*   //orig
  stringstream ssInfo_orig;
  ssInfo_orig << "Original" << endl;
  for (unsigned int i = 0; i < v.size(); ++i)
  {
    max = 0.0;
    sum = 0.0;
    int cluster_label = mapSNPid_label[ v[i] ];
    int ind = mapLabel_ind[ cluster_label ];
    aRegion.clear();
    for (set<snp_index_t>::iterator itS = vectClust[ind].begin(); itS != vectClust[ind].end(); ++itS)
    {
      aSNP_info.SNP = *itS;
      aSNP_info.Chr = (data.getSNP(*itS))->getChromosome();
      aSNP_info.pos = (data.getSNP(*itS)).getBasePairPosition();
      aRegion.insert(aSNP_info);
      if (mapSNPid_Pmi_Y[ *itS ] > max)
        max = mapSNPid_Pmi_Y[ *itS ];
      sum += mapSNPid_Pmi_Y[ *itS ];
    }
    for (set<TSNP_Info>::iterator it = aRegion.begin(); it != aRegion.end(); ++it)
    {
      ssInfo_orig << setw(7) << it->SNP << "\t" << setw(17) << data.getSNP( it->SNP ).getSnpId() << "\t" << setw(4) << it->Chr << "\t" << setw(12) << it->pos << "\t" 
             << setw(12) << mapSNPid_Pmi_Y[ it->SNP ] << endl;
    }
    ssInfo_orig << setw(7) << cluster_label << "\t" << setw(17) << "cluster(sum, max)" << "\t" << setw(4) << "->" << "\t"
           << setw(12) << sum << "\t" << setw(12) << max << endl << endl;
  }
*/  
  set<TSNP_Info> *aRegionPtr;
  set<TRegionSet_info> setNewRegion;
  for (unsigned int i = 0; i < v.size(); ++i)
  {
    int cluster_label = mapSNPid_label[ v[i] ];
    int ind = mapLabel_ind[ cluster_label ];
    aRegion.clear();
    aRegionPtr = new set<TSNP_Info>;
    for (set<snp_index_t>::iterator itS = vectClust[ind].begin(); itS != vectClust[ind].end(); ++itS)
    {
      aSNP_info.SNP = *itS;
      aSNP_info.Chr = (data.getSNP(*itS)).getChromosome();
      aSNP_info.pos = (data.getSNP(*itS)).getBasePairPosition();
      aRegion.insert(aSNP_info);
      aRegionPtr->insert(aSNP_info);
      if (mapSNPid_Pmi_Y[ *itS ] > max)
        max = mapSNPid_Pmi_Y[ *itS ];
      sum += mapSNPid_Pmi_Y[ *itS ];
    }
    TRegionSet_info anItem;
    anItem.s = aRegionPtr;
    setNewRegion.insert(anItem);
  }
  
  for (set<TRegionSet_info>::iterator itS = setNewRegion.begin(); itS != setNewRegion.end(); ++itS)
  {
    TRegionSet_info s = *itS;
    set<TSNP_Info> *aSNP = s.s;
    sum = 0.0;
    max = 0.0;
//    int cluster_label = mapSNPid_label[ aSNP->begin()->SNP ];
    for (set<TSNP_Info>::iterator it = aSNP->begin(); it != aSNP->end(); ++it)
    {
      ssInfo << setw(7) << it->SNP << "\t" << setw(17) << data.getSNP( it->SNP ).getSnpId() << "\t" << setw(4) << it->Chr << "\t" << setw(12) << it->pos << "\t" 
             << setw(12) << mapSNPid_Pmi_Y[ it->SNP ] << endl;  
      if (mapSNPid_Pmi_Y[ it->SNP ] > max)
        max = mapSNPid_Pmi_Y[ it->SNP ];
      sum += mapSNPid_Pmi_Y[ it->SNP ];             
    }        
    if (method_Name.compare("cluster_SUM") == 0 || method_Name.compare("cluster_MAX") == 0 || method_Name.compare("piMass_cluster_Max") == 0)
    {
      ssInfo << setw(7) << " " << "\t" << setw(17) << "cluster(sum, max)" << "\t" << setw(4) << "-->" << "\t"
             << setw(12) << sum << "\t" << setw(12) << max << endl << endl;
    }
    else
      ssInfo << endl;
    delete s.s;
  }
 
  Model model(data);
  if (v.size() > 0)
    model.addManySNP(v);
  if (model.computeRegression())
    model.computeMSC();
  else
  {
    model.computeMSCfalseRegression();
  }  
  ofstream file;
  file.open((parameter.out_file_name + "_Report.txt").c_str(), ios::app);
  file << method_Name << endl << endl;
  file << ssInfo.str() 
       << model << endl;
//       << "+++++++++++++++++++++++++++++++++" << endl;
//  file << ssInfo_orig.str() << endl;
  file << "________________________________________________________________________________" << endl << endl;
  file.close(); 
}


/**
 * @brief writes pool of size 30K or 100K to the file
 * 
 */
void GA::writePoolofSize(bool &flag, int maxSize)
{ 
  if (maxSize == 30000)
  {
    if (flag && pool.size() >= 30000)
    {
      stringstream ss;
      const time_t time_end = time(NULL);
      ss << "GA time: " << (time_end > time_start? sec2time(time_end - time_start) : sec2time(24 * 3600 + time_end - time_start));  
      writePoolToFile(ss, "_30K");
      flag = false;
      
    }
  }
  else
  {
    if (flag && pool.size() >= 100000)
    {
      stringstream ss;
      const time_t time_end = time(NULL);
      ss << "GA time: " << (time_end > time_start? sec2time(time_end - time_start) : sec2time(24 * 3600 + time_end - time_start));  
      writePoolToFile(ss, "_100K");
      flag = false;
      printLOG( ("100K" + ss.str() + "\n\r") );
    }
  }
}

/**
 * @brieb wyznacza SNPy, 
 * @param th - próg
*/
void GA::regionStrategy(const string method_Name, Model &model, multimap<long double, snp_index_t>& Pmi_Ysort, map<snp_index_t, long double> &mapSNPid_Pmi_Y, double th, set<TRegionSet_info> &setNewRegion)
{
//  cout << "th: " << th << endl;
//  cout << "Pmi_Ysort.size(): " << Pmi_Ysort.size() << endl;
//  cout << "mapSNPid_Pmi_Y.size(): " << mapSNPid_Pmi_Y.size() << endl;  
  
  const bool TEST = false;
  map<snp_index_t, unsigned int> mapSNP2ind;   // odwzorowanie SNP -> indeks w tablicy tabCorr oraz tablicy tabCorr[]
  vector<vector<double> > tabCorr;
  vector<snp_index_t> SNPs;                    // SNPy, które biorą udział w operacji, odwzorowanie index tablicy -> SNP
  vector< set<snp_index_t> > tabRegionSNPs;    // tablica zawierająca SNPy, które należą do rejonu
  vector<double> SNPs_posteriori;              // posteriori rejonu
  snp_index_t aSNP;
  double absCorr;
  unsigned int i = 0;
  multimap<long double, snp_index_t>::reverse_iterator itM;
  vector<double> vv;
  vector<snp_index_t> v_all;
  set<snp_index_t>s;
  ofstream  outFile;
  if (TEST) outFile.open((parameter.out_file_name + "_region_mapSNP2ind.txt").c_str(),  fstream::out | fstream::trunc );
  if (TEST) outFile << "aSNP" << "\t" << "mapSNP2ind" << "\t" << "mapSNPid_Pmi_Y" << endl;
  for (itM = Pmi_Ysort.rbegin(); itM != Pmi_Ysort.rend() && (*itM).first > th; ++itM)
  {
    aSNP = itM->second;
    mapSNP2ind.insert(pair<snp_index_t, unsigned int> (aSNP, i));
    tabCorr.push_back(vv);
    SNPs.push_back(aSNP);
    SNPs_posteriori.push_back(mapSNPid_Pmi_Y[aSNP]);
    tabRegionSNPs.push_back(s);
    if (TEST) outFile << aSNP << "\t" << mapSNP2ind[aSNP] << "\t" << mapSNPid_Pmi_Y[aSNP] << endl;
    if (i != 0)
      tabCorr[i].resize(i);
    for (unsigned int k = 0; k < i; k++)
    {
      absCorr = fabs( data.computeCorrelation( SNPs[i], SNPs[k] ) ); // compute correlation between model SNP and SNP j     
      tabCorr[ mapSNP2ind[aSNP] ][k] = absCorr;
    }
    ++i;
  }
  outFile.close();
  if (TEST)
  {
    outFile.exceptions ( ofstream::eofbit | ofstream::failbit | ofstream::badbit ); // checks if poolfile can be written
    try
    {
      outFile.open((parameter.out_file_name + "_regionStrategy.txt").c_str(),  fstream::out | fstream::trunc );
      
      for (i = 0; i < tabCorr.size(); ++i)
      {
        aSNP = SNPs[i];
        outFile << i << "\t" << aSNP << "\t" << SNPs_posteriori[i];
        for (unsigned int k = 0; k < i; ++k)
        {
          outFile << "\t" << tabCorr[i][k];
        }
        outFile << endl;
      }
      // the last line
      outFile << "no\tSNP_ind \tpost. \t";
      
      for (unsigned int k = 0; k < i; ++k)
      {
        outFile << "\t" << k;
      }
      outFile << endl;    
      
      outFile.flush();
      outFile.close();
    } 
    catch (ofstream::failure e)
    {
      cerr << "Could not write a _regionStrategy file" <<endl;
      exit(-1);
    }
    outFile.open((parameter.out_file_name + "_regionStrategy.txt").c_str(),  fstream::out | fstream::app );
    outFile << endl;
    outFile << "SNP_ind\tposter.\tposter(region)\tregion" << endl;
  }
  // for each SNP which is correlated each other we sum up a posteriori probability
  for (i = 0; i < tabCorr.size(); ++i)
  {
    aSNP = SNPs[i];
    for (unsigned int k = 0; k < i; ++k)
    {
      if (tabCorr[i][k] > 0.5  && (data.getSNP( SNPs[i] )).getChromosome() == (data.getSNP( SNPs[k] )).getChromosome() )  
      {                                                                                                     // 3 conditions: 1. correlations > 0.5, 2. the same getChromosome, 
        if ( fabs((data.getSNP(SNPs[i])).getBasePairPosition() - (data.getSNP(SNPs[k])).getBasePairPosition()) < 10000000) // 3. distance
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
  
  vector<snp_index_t> v;
  v_all.clear();
  TSNP_Info aSNP_info;
  set<TSNP_Info> aRegion;
  set<TSNP_Info> *aRegionPtr;// = new set<TSNP_Info>;
  stringstream ssInfo;

  ssInfo << setw(7) << "Nr" << "\t" << setw(15) << "SNP Id" << "\t" << setw(4) << "Chr" << "\t" << setw(12) << "Pos" << "\t" << setw(12) << "posteriori" << "\t" 
         << setw(12) << "posteriori (region)" << endl << endl;
//  set<TRegionSet_info> setNewRegion;
  for (i = 0; i < SNPs.size(); ++i)
  {
    if (SNPs_posteriori[i] > 0.5)  // If a posteriori probability of a region is greater than 0.5 we report this region.
    { 
      if (TEST) outFile << SNPs[i] << "\t" << mapSNPid_Pmi_Y[ SNPs[i] ] << "\t" << SNPs_posteriori[i] << "\t";
      if (TEST) outFile << SNPs[i];
      v.push_back(SNPs[i]);
      for (set<snp_index_t>::iterator it = tabRegionSNPs[i].begin(); it != tabRegionSNPs[i].end(); ++it)
      {
          if (TEST) outFile << " " << SNPs[*it];
          if ( mapSNPid_Pmi_Y[ SNPs[*it] ] > mapSNPid_Pmi_Y[ v[v.size() -1] ] )
            v[v.size() -1] = SNPs[*it];    // wybór najlepszego SNPa z danego rejonu, żeby reprezentował ten rejon w modelu
      }
      v_all.push_back(SNPs[i]);
      if (TEST)  
      {
        outFile << endl;  
        outFile << "\t" << data.getSNP( SNPs[i] ).getSnpId();
        for (set<snp_index_t>::iterator it = tabRegionSNPs[i].begin(); it != tabRegionSNPs[i].end(); ++it)
          outFile << " " << data.getSNP( SNPs[*it] ).getSnpId(); 
        outFile << endl;
      }
      aRegionPtr = new set<TSNP_Info>;
      aRegion.clear();
      aSNP_info.SNP = SNPs[i];
      aSNP_info.Chr = (data.getSNP(SNPs[i])).getChromosome();
      aSNP_info.pos = (data.getSNP(SNPs[i])).getBasePairPosition();
      aRegion.insert(aSNP_info);
      aRegionPtr->insert(aSNP_info);
      for (set<snp_index_t>::iterator it = tabRegionSNPs[i].begin(); it != tabRegionSNPs[i].end(); ++it)
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
      /*
      for (set<TSNP_Info>::iterator it = aRegion.begin(); it != aRegion.end(); ++it)
      {
        SNPs_posteriori[ mapSNP2ind[it->SNP] ] = -1;
      }
      */
    }
  }    

  for (set<TRegionSet_info>::iterator itS = setNewRegion.begin(); itS != setNewRegion.end(); ++itS)
  {
    TRegionSet_info s = *itS;
    set<TSNP_Info> *aSNP = s.s;
    for (set<TSNP_Info>::iterator it = aSNP->begin(); it != aSNP->end(); ++it)
    {
      ssInfo << setw(7) << it->SNP << "\t" << setw(15) << data.getSNP( it->SNP ).getSnpId() << "\t" << setw(4) << it->Chr << "\t" << setw(12) << it->pos << "\t" 
             << setw(12) << mapSNPid_Pmi_Y[ it->SNP ] << "\t" << setw(12) << SNPs_posteriori_orig[ mapSNP2ind[it->SNP] ] << endl;  
    }        
    ssInfo << endl;
  }
  
  
  ///// a compact raport
  SNPs_posteriori = SNPs_posteriori_orig;

  // zwolnić pamięć dla aRegionPtr
  for (set<TRegionSet_info>::iterator itS = setNewRegion.begin(); itS != setNewRegion.end(); ++itS)
  {
    TRegionSet_info s = *itS;
    delete s.s;
    s.s = 0;
  }  
  
  
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
      for (set<snp_index_t>::iterator it = tabRegionSNPs[i].begin(); it != tabRegionSNPs[i].end(); ++it)
      {
        aSNP_info.SNP = SNPs[*it];
        aSNP_info.Chr = (data.getSNP(SNPs[*it])).getChromosome();
        aSNP_info.pos = (data.getSNP(SNPs[*it])).getBasePairPosition();
//        cout << i << "\taSNP_info: SNP:" << aSNP_info.SNP << "\t" << aSNP_info.Chr << "\t" << aSNP_info.pos << endl;
        aRegion.insert(aSNP_info);
        aRegionPtr->insert(aSNP_info);

      }
      TRegionSet_info anItem;
      anItem.s = aRegionPtr;
      setNewRegion.insert(anItem);
      for (set<TSNP_Info>::iterator it = aRegion.begin(); it != aRegion.end(); ++it)
      {
        SNPs_posteriori[ mapSNP2ind[it->SNP] ] = -1;
//        cout << "SNP: " << it->SNP << ", SNPs_posteriori[ mapSNP2ind[it->SNP] ]: " << SNPs_posteriori[ mapSNP2ind[it->SNP] ] << endl;
      }
    }
  }
  
  stringstream ssInfo_Compact;
  ssInfo_Compact << setw(7) << "Nr" << "\t" << setw(15) << "SNP Id" << "\t" << setw(4) << "Chr" << "\t" << setw(12) << "Pos" << "\t" << setw(12) << "posteriori" << "\t" 
         << setw(12) << "posteriori (region)" << endl << endl;

  for (set<TRegionSet_info>::iterator itS = setNewRegion.begin(); itS != setNewRegion.end(); ++itS)
  {
    TRegionSet_info s = *itS;
    set<TSNP_Info> *aSNP = s.s;
    for (set<TSNP_Info>::iterator it = aSNP->begin(); it != aSNP->end(); ++it)
    {
      ssInfo_Compact << setw(7) << it->SNP << "\t" << setw(15) << data.getSNP( it->SNP ).getSnpId() << "\t" << setw(4) << it->Chr << "\t" << setw(12) << it->pos << "\t" 
                     << setw(12) << mapSNPid_Pmi_Y[ it->SNP ] << "\t" << setw(12) << SNPs_posteriori_orig[ mapSNP2ind[it->SNP] ] << endl;  
    }        
    ssInfo_Compact << endl;
  }

  Model m(data);
  if (v.size() > 0)
    m.addManySNP(v);
  if (m.computeRegression())
    m.computeMSC();
  else
    m.computeMSCfalseRegression();
  if (v_all.size() > 0)
    model.addManySNP(v_all);
  if (model.computeRegression())
    model.computeMSC();
  else
  {
    model.computeMSCfalseRegression();
  }  
  if (TEST) outFile << model << endl;
  if (TEST) outFile.close();
  ofstream file;
  file.open((parameter.out_file_name + "_Report.txt").c_str(), ios::app);
  file << method_Name << " (compact report)"<< endl << endl;
  file << ssInfo_Compact.str()
       << m << endl;
  file << "________________________________________________________________________________" << endl << endl;  
  
  file << method_Name << " (full report)"<< endl << endl;
  file << ssInfo.str()
       << m << endl;
  file << "________________________________________________________________________________" << endl << endl;  

  
  file.close();            
  model = m;
}


void GA::calculatePOWER_FDR_clustRegion(vector<snp_index_t> &causalSNPs, TPOWER_FDR &powerFDR, set<TRegionSet_info> &setNewRegion, map<snp_index_t, int>& recognizedSNPs)
{
  assert(causalSNPs.size() > 0);
  bool isTP;
  set<int> TP_Set;
  int FP = 0;
  set<int> causal_Set(causalSNPs.begin(), causalSNPs.end());
  for (set<TRegionSet_info>::iterator itS = setNewRegion.begin(); itS != setNewRegion.end(); ++itS)
  {
    TRegionSet_info s = *itS;
    set<TSNP_Info> *aInfo = s.s;
    isTP = false;
    for (set<TSNP_Info>::iterator it = aInfo->begin(); it != aInfo->end(); ++it)
    {
      if (causal_Set.find(it->SNP) != causal_Set.end())
      {
        recognizedSNPs[it->SNP]++;
        TP_Set.insert(it->SNP);
        isTP = true;
        break;
      }
    }
    if (isTP == false)
      ++FP;
  }
  int sumTP = TP_Set.size();
  int sumFP = FP;
  
  powerFDR.POWER = (sumTP + 0.0) / causalSNPs.size();
  if (sumFP + sumTP > 0)
  {
    powerFDR.FDR = (sumFP + 0.0) / (sumFP + sumTP);
  }
  else
  {
    powerFDR.FDR = 0;
  }  
  powerFDR.FDcount = sumFP;
}


void GA::test_region()
{
/* 
  157462       rs7666489    4      4497328  1.53703e-07     0.835508
 194054       rs6553419    4    170282434  1.26595e-07     0.835497
 203492      rs10866471    5     10272168  1.66906e-05     0.835584
 203495       rs4701866    5     10276371   1.8096e-07      0.83554
 203496       rs4702680    5     10276423  1.86269e-06     0.835584
 203497      rs12653646    5     10277631  3.87043e-06     0.835584
 203500       rs1355096    5     10279025  1.79843e-07     0.835524
 203504        rs906686    5     10281825  2.11606e-07     0.835543
 203506       rs7720661    5     10283270  2.35151e-07     0.835543
 203517       rs2445872    5     10316692  5.82796e-07     0.835589
 203520      rs10074613    5     10322935  2.21168e-07     0.835559
 203521      rs10053516    5     10323505  2.65332e-07     0.835561
 450882       rs2862465   11     42336068  1.59971e-06     0.835561
 511769       rs2776982   13     28979384   5.2529e-07     0.835523
 511771       rs2482093   13     28998007   6.2172e-07       0.8355
 652626       rs6044502   20     16832497  4.33657e-07     0.835495
 652628       rs2208146   20     16845102  1.14033e-06      1.51687
 652630       rs2208150   20     16851467  4.04198e-07     0.835547
 652655       rs6075165   20     16986263   4.9917e-07     0.835556
 652656       rs2011424   20     16987629  4.10683e-07      0.83556
 652657       rs6080544   20     16987959  4.59708e-07     0.835546
 652666      rs11700038   20     17000919  5.30952e-07     0.835472
 697044       rs6635091   23    134665706  2.40151e-07     0.835524
 697063       rs7882434   23    134842709  2.53974e-05     0.835494
 697067       rs4829804   23    134858912     0.835469     0.835577
 697074        rs903143   23    134976149  1.46358e-06     0.835556
 697075       rs5929713   23    134976949  3.58562e-06     0.835554
 697076       rs5975658   23    134977630  9.23227e-07      0.83552
 697077       rs6635218   23    134979824  2.59903e-06     0.835547
 697078       rs6633796   23    134981217  2.86623e-06     0.835521
 697079       rs5929715   23    134985417  1.83803e-05     0.835562
 697080       rs6528335   23    134988715  2.64183e-06     0.835524
 697082       rs5975665   23    134998159  3.25724e-06     0.835548
 697083       rs5930859   23    135000657  1.54988e-05     0.835559
*/ 

  int tab[] = {
               157462, 194054, 203492, 203495, 203496, 203497, 203500, 203504, 203506, 203517, 
               203520, 203521, 450882, 511769, 511771, 652626, 652628, 652630, 652655, 652656, 
               652657, 652666, 697044, 697063, 697067, 697074, 697075, 697076, 697077, 697078, 
               697079, 697080, 697082, 697083
              };
  const int tab_size = sizeof(tab)/sizeof(tab[0]);
  vector<snp_index_t> SNPs(tab, tab + tab_size);                    // SNPy, które biorą udział w operacji, odwzorowanie index tablicy -> SNP
  map<snp_index_t, unsigned int> mapSNP2ind;   // odwzorowanie SNP -> indeks w tablicy tabCorr oraz tablicy tabCorr[]
  vector<vector<double> > tabCorr(SNPs.size());
  vector< set<snp_index_t> > tabRegionSNPs;    // tablica zawierająca SNPy, które należą do rejonu
  snp_index_t aSNP;
  double absCorr;

  ofstream  outFile;
  outFile.open((parameter.out_file_name + "_region_TEST_mapSNP2ind.txt").c_str(),  fstream::out | fstream::trunc );
  outFile << "mapSNP2ind" << "\t" << "aSNP" << endl;
  for (  unsigned int i = 0; i < SNPs.size(); ++i)
  {
    aSNP = SNPs[i];
    mapSNP2ind.insert(pair<snp_index_t, unsigned int> (aSNP, i));
    outFile << mapSNP2ind[aSNP] << "\t" << aSNP << "\t" << endl;
    if (i != 0)
      tabCorr[i].resize(i);
    for (unsigned int k = 0; k < i; k++)
    {
      absCorr = fabs( data.computeCorrelation( SNPs[i], SNPs[k] ) ); // compute correlation between model SNP and SNP j     
      tabCorr[ mapSNP2ind[aSNP] ][k] = absCorr;
    }
  }
  outFile.close();
  outFile.exceptions ( ofstream::eofbit | ofstream::failbit | ofstream::badbit ); // checks if poolfile can be written
  try
  {
    outFile.open((parameter.out_file_name + "_regionStrategy_TEST.txt").c_str(),  fstream::out | fstream::trunc );
    
    for (unsigned i = 0; i < tabCorr.size(); ++i)
    {
      aSNP = SNPs[i];
      outFile << i << "\t" << aSNP; // << "\t" << SNPs_posteriori[i];
      for (unsigned int k = 0; k < i; ++k)
      {
        outFile << "\t" << tabCorr[i][k];
      }
      outFile << endl;
    }
    // the last line
    outFile << "no\tSNP_ind \tpost. \t";
    
    for (unsigned int k = 0; k < tabCorr.size(); ++k)
    {
      outFile << "\t" << k;
    }
    outFile << endl;    
    
    outFile.flush();
    outFile.close();
  } 
  catch (ofstream::failure e)
  {
    cerr << "Could not write a _regionStrategy file" <<endl;
    exit(-1);
  }
  outFile.open((parameter.out_file_name + "_regionStrategy_TEST.txt").c_str(),  fstream::out | fstream::app );
  outFile << endl;
  outFile << "SNP_ind\tposter.\tposter(region)\tregion" << endl;
  outFile.close();
  
  outFile.open((parameter.out_file_name + "_correlation_TEST.txt").c_str(),  fstream::out | fstream::trunc );
  outFile << "no\tSNP_Id\tSNP_ind \t SNP_ind\t correlation" << endl;
  for (  unsigned int i = 0; i < SNPs.size(); ++i)
  {
    absCorr = fabs( data.computeCorrelation( SNPs[i], 697067 ) ); // compute correlation between model SNP and SNP j     
    outFile << i << "\t" << data.getSNP( SNPs[i] ).getSnpId() << "\t" << SNPs[i] << "\t" << 697067 << "\t" << absCorr << endl;
    cout << i << "\t" << data.getSNP( SNPs[i] ).getSnpId() << "\t" << SNPs[i] << "\t" << 697067 << "\t" << absCorr << endl;
    
  }
  outFile.close();
}


  