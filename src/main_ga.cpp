#include "Parameter.hpp"
#include "Helpfull.hpp"
#include "Log.hpp"
#include "MData.hpp"
#include "Model.hpp"
#include "GA.hpp"
#include <iomanip>
#include <getopt.h>
#include "Exception.hpp"

/** global variable for logfile */
ofstream LOG;

/** Global parameter holder */
Parameter parameter;


string timeStemple(const time_t t)
{
  int hour, min, sec;
  sec = t;
  hour = sec / 3600;
  sec = sec % 3600;
  min = sec / 60;
  sec = sec % 60;
  return int2strPadWith( hour, 2, '0' ) + ":" + int2strPadWith( min, 2, '0' ) + ":" + int2strPadWith(   sec, 2, '0' );  
}

void readSortFile(string fileName, string &gaTime, set<snp_index_t>& gaSNPs, double &msc);      // odczytuje pliki sord GA
void readStepwisetFile(string fileName, string &gaTime, set<snp_index_t>& gaSNPs, double &msc); // odczytuje pliku Stepwise 
void saveMapToFile(const string &fileName, map<snp_index_t, int>& map_SNP2count);
void savePOWER_FDR_ToFile(const string &fileName, vector <TPOWER_FDR>& vPOWER_FDR, vector <TPOWER_FDR>& vPOWER_FDR_ClGA, string *times, double *msc);

void runGA(int outNo, const string & modelsFileName);
void runPiMassConvert();
void runPiMassResults(const string & inputFileDir);
void runPoolReader(const string & pool_filename);
void runClusterFromPosterior();
void runGA2nd(int artNo);
void runMOSGWA(int outNo);
void generateY();

void print_usage (char const* program_name)
{
  cerr <<  "Usage:" << program_name << " [options] conf_fileName" << endl;
  cerr << " -h --help\t Displays this usage information." << endl;
  cerr << " -p --pool pool_filename conf_fileName\t Reads pool data from pool_filename" << endl;
  cerr << " -x --extractsPiMass conf_fileName\t Extracts some information from piMass results files" << endl;
  cerr << " -c --convertsPiMass conf_fileName\t Converts input data files to piMass format" << endl;
  cerr << " -n --number no conf_fileName\t Runs GA_MOSGWA for the trait \'no\', e.g. -n 5 for the 5th trait'. DOES NOT take into account the trait from conf_fileName. It is for multirun" << endl;
  cerr << " -i --initPop initPopulationfilenName conf_filename\t Reads the initial population from file initPopulationfilenName and runs GA_MOSGWA" << endl;
  cerr << " -o --outNo no conf_fileName\t Sets the number wchich is used in output file names" << endl << endl;
  cerr << " -l --clusters conf_fileName\tDON'T USE! Calculates clusters models from *sort files" << endl;
  cerr << " -a --art no\tDON'T USE!. Calculates statistics of no-th simulation for the second publication, " << endl;
  cerr << " -s --stepwise\t Calculates Stepwise method (original MOSGWA)" << endl;
  cerr << " -g --generateY conf_fileName\t Generates Y file." << endl;
  //exit (exit_code);
}

/**; Application entry point */
int main ( const int argc, const char *argv[] )
{
  int next_option;
  const char* const short_options = "hp:x:cn:i:o:la:sg";
  const struct option long_options[] =
    {
      {"help",           0, NULL, 'h'},
      {"pool",           1, NULL, 'p'},  // reads pool and calculates results
      {"extractsPiMass", 1, NULL, 'x'},  // extracts same results from piMass outpu files
      {"convertsPiMass", 0, NULL, 'c'},  // converts PLink input files into piMass input files
      {"number",         1, NULL, 'n'},  // runs the GA for Y of number n
      {"initPop",        1, NULL, 'i'},  // reads the initial population from a file and runs the Ga
      {"outNo",          1, NULL, 'o'},  // saves results file using the number outNo
      {"clusters",       0, NULL, 'l'},  // calculates the cluster model or cluster models
      {"art",            1, NULL, '2'},  // calculates the statistics for 2nd article
      {"stepwise",       0, NULL, 's'},  // calculates the stepwise model
      {"generateY",      0, NULL, 'g'},  // generates Y
      {NULL,             0, NULL, 0}
    };

  const char* pool_filename = NULL;
//  string inputFileDir;
//  int verbose = 0;
//  bool piMassConvert = false;
//  bool poolReader = false;
//  bool initPopulation = false;
//  bool piMassResult = false;
  int no = -1; 
  int outNo = -1;
  int artNo;
  stringstream ss_no, ss_outNo; 
  string initPopulationfilenName = "";
  string inputFileDir;
  stringstream ssF;

  enum Operations {opRun = 0, opPoolReader, opPiMassResults, opPiMassConvert, opRunNinputFile, opInitPopulation, opCluster, op2nd, opMosgwa, opGenerateY};
  Operations operation = opRun;
  
  bool setArgs = false;
  do {
    next_option = getopt_long(argc, (char* const*) argv, short_options, long_options, NULL);
    switch (next_option) {
    case 'h':
      print_usage (argv[0]);
      return 0;
      break;
    case 'p':
      setArgs = true;
      pool_filename = optarg;
      operation = opPoolReader;
      cout << "1st: operation = opPoolReader: " << operation << endl;
      break;
    case 'x': // extracts PiMass results
      setArgs = true;
      ssF << optarg;
      cout << ssF.str() << endl;
      ssF >> inputFileDir;//;
      operation = opPiMassResults;
      break;
    case '?':
      print_usage ( argv[0]);
      return 0;
      break;
    case 'c':
      setArgs = true;
      operation = opPiMassConvert;
    break;
    case 'n':
      ss_no << optarg;
      ss_no >> no;
      cout << "n: " << no << endl;
      operation = opRunNinputFile;
      break;
    case 'o':
      ss_outNo << optarg;
      ss_outNo >> outNo;
      
//      cout << "parameter.out_file_name: " << parameter.out_file_name  << endl;
      break;
    case 'i':
      setArgs = true;
      initPopulationfilenName = optarg;
      operation = opInitPopulation;
      break;
    case 'l':
      setArgs = true;
      operation = opCluster;
      break;
    case 'a':
      ss_no << optarg;
      ss_no >> artNo;
      cout << "artNo: " << artNo << endl;
      operation = op2nd;
      setArgs = true;
      break;
      case 's':
      operation = opMosgwa;
      setArgs = true;
      break;
      case 'g':
        operation = opGenerateY;
        setArgs = true;
        break;
    default: // opRun - run GA
      ;
    } // switch
  } while (next_option != -1);  

  cout << "argc: " << argc << endl;
  for (int i = 0; i < argc; i++)
    cout << i << "), argv: " << argv[i] << endl;
  //cout << "optind: " << optind << endl;
  //cout << "argv[optind]: " << argv[optind] << endl;
//  cout << "poolReader: " << poolReader << endl;
  cout << "no: " << no << endl;  
  cout << "outNo: " << outNo << endl;
  //cout << "setArgs: " << setArgs << endl;
  printStartScreen();
  
  if (setArgs == true || no >= 0 || outNo >= 0)
  {
    int argcc = 2;
    const char* argvv[2] = {argv[0] /*"./GA_MOSGWA"*/, argv[optind]};
//    cout << "argcc: " << argcc << ", AG argv[1]: " << argvv[1] << ":: " << argv[optind]<< endl;
    parameter.setParameters( argcc, argvv );
    
    if (no >= 0)
    {
//      stringstream ss;   // After calculations I increase the value of trait_name_in_yvm 
//      ss << no;
      //parameter.in_values_name = ss.str(); // declare( "data", "trait_name_in_yvm", in_values_name );

      parameter.in_values_int = no;
      cout << "parameter.in_values_int: " << parameter.in_values_int << endl; 
      parameter.out_file_name = parameter.out_file_name + "_#" + int2str(parameter.in_values_int);
    }
    if (outNo >= 0)
    {
      parameter.out_file_name = parameter.out_file_name + "_" + int2str(outNo);
      cout << "parameter.out_file_name: " << parameter.out_file_name  << endl;
      
    }
  }  
  else
  {
    parameter.setParameters( argc, argv );
  }
  
  string logFileName(parameter.out_file_name + /*"_#" + int2str(parameter.in_values_int) + (no >= 0? "_" + int2str(no): "") + */ ".log" );
  try 
  {
    printLOG( "Start: open log to file \"" + logFileName + "\"" );
    LOG.open( logFileName.c_str(), fstream::out );
  } 
  catch ( ofstream::failure e ) 
  {
    cerr << "Could not open logfile \"" + logFileName + "\"" <<endl;
  }
  // checks if logfile can be written.
  LOG.exceptions ( ofstream::eofbit | ofstream::failbit | ofstream::badbit );

  string realModelFilename = parameter.causalModelFilename;
  string modelsFileName = initPopulationfilenName;
  cout << "operation: " << operation << endl;
  switch (operation)
  {
    case opPoolReader:
      runPoolReader(pool_filename);
      break;
    case opPiMassResults:
      runPiMassResults(inputFileDir);
      break;
    case opPiMassConvert:
      runPiMassConvert();
      break;
    case opCluster:
      runClusterFromPosterior();
      break;
 //   case opRunNinputFile:
 //     break;
 //   case opInitPopulation:
      //break;
    case op2nd:
      //runGA(artNo);   // 
      runGA2nd(artNo);
      break;
    case opMosgwa:
      runMOSGWA(outNo);
      break;
    case opGenerateY:
        generateY();
      break;
    default:
      runGA(outNo, modelsFileName);
  }
  try 
  {
    printLOG("End");
    LOG.close();
  } 
  catch ( ofstream::failure e ) 
  {
    cerr << "Could not close logfile \"" + logFileName + "\"" << endl;
  }
}

/*
 * Converts Plink data files into piMass data files
 */
void runPiMassConvert()
{
//  const int gensNo = 100; // liczba genów, dla których preparujemy pliki Y dla piMassa
//  string gensNamesTab[gensNo] = {"hmm25278-S", "hmm26651-S", "hmm32074-S", "hmm34610-S"};
//  int gensIDNo[gensNo] = {41133, 41810, 44935, 45851};
  //string gensNamesTab[gensNo] = {"0", "1"};
  //int gensIDNo[gensNo] = {1, 2};
  
  const string in_file_name = parameter.in_files_plink + "_Plink.txt";
  cout << "in_file_name: " << in_file_name << endl;
  
  MData data;
  data.calculateIndividualTests();    
  
  
  
  // podmienia genotyp na ten, który jest w pliku Matlaba   
//  data.writeBEDfromMatlab(in_file_name);
//  cout << "BED and mgt files form tex genotype" << endl;
//  exit(0);

  
  
  for (unsigned int i = 1; i <= 1; ++i)
  {    
    
    //parameter.in_values_name = gensNamesTab[i];
    //parameter.in_values_int = gensIDNo[i];
    // ANOTHER WAY
//    parameter.in_values_int = i; 
    //parameter.in_values_int = i + 1;  // "trait_position_in_yvm",

//    cout << "pring Geno data: " << endl;
//    cout << endl << "print individuals: " << endl;
//    cout << i << "parameter.in_values_int: " << parameter.in_values_int << endl; 
// 07.07    data.writePiMassFormatFiles(parameter.in_values_int); // ! tylko dla tego argumentu (tu 1), generowany jest plik genotyp .mgt!!!!!!
//    parameter.in_values_int = i + 1;  // "trait_position_in_yvm",
  }
}

/*
 * Extracts same information from piMass fies
 * None of parameters are used
 */
void runPiMassResults(const string & inputFileDir)
{
  unsigned int modelsNo_ = 1;
  unsigned int maxNoProgressIter_ = parameter.maxNoProgressIter; // 1000;
  double pCross_ = parameter.pCross;                             // 0.65;
  double pMutation_ = parameter.pMutation;                       // 0.05;
  unsigned int tournamentSize_ = parameter.tournamentSize;       // 2;
  double correlationThreshold_ = parameter.correlationThreshold; // 0.7;
  int correlationRange_ = parameter.correlationRange;            // 500  
  double regionMinCorrelation_ = parameter.regionMinCorrelation;
  int B_ = 1;
  string modelsFileName = "";

  stringstream dirName(inputFileDir);    // katalog z wynikiami
                                               //            sim_0M_piMass.2014.06.24
  stringstream outName("piMass_PowerFRD.txt");
  string outFileName = dirName.str();// + "max_"; //"Real_2013.04.25/hmm26651-S/piMass";
  
  cout << "inputFileDir: " << inputFileDir.c_str() << endl;
  cout << "dirName: " << dirName.str() << endl;
  
  map<snp_index_t, int> recognizedSNPs;   // liczebność snp'ow w i-iteracjach
  stringstream sp_sort;
  sp_sort << "Power\tFDR\tFD_count\tbadSNP\tposterior of bad\tposterior of Cluster of Bad\tmBIC2";
  ofstream outFile;
  
//  outFile.open("piMass_cluster_prior.txt",  fstream::out | fstream::trunc );  // wyczyszczenie pliku z prawd a poster dla klastrów
//  outFile.close();                                                                      // ścieżka na stałe w calculatePOWER_FDR_clust 
  
  outFile.open((outFileName + "piMass_PowerFRD.txt").c_str(), ios::trunc);
  outFile << sp_sort.str() << endl;
  outFile.close();
  outFile.open((outFileName + "piMassModels.txt").c_str(), ios::trunc);
  outFile.close();
  cout << "start" << endl;
  // UWAGA dla tych samych Y GA ga przed petla
  //GA ga(modelsNo_, maxNoProgressIter_, pCross_, pMutation_, tournamentSize_, B_, modelsFileName, correlationThreshold_, correlationRange_, true);
  stringstream fileName;
  fileName << "output";
  
  vector< multiset<long double> > tabCausalPost;  // zawiera prawd. poster. odkrytych snpów przyczynowych
  map<snp_index_t, int> mapSNPCausal_ind; // zawiera odwzorowanie snp_id -> ind w tabCausalPost
  
  ofstream modelsFile;
  modelsFile.open((outFileName + "piMassModels.txt").c_str(), ios::trunc);
  modelsFile.close();
  int start = 1, end = 6;
  string orig_out_file_name = parameter.out_file_name;
  for (int i = start; i <= end; ++i)
  {
    parameter.out_file_name = orig_out_file_name + "_#" + int2str(i);
    // UWAGA!!! Sprawdź, czy trzeba ustawić numer cechy!
    // parameter.in_values_int = i;
    cout << "parameter.in_values_int: " << parameter.in_values_int << endl;   
    stringstream command;
    
//    command << "tar -xvzf " << dirName.str() << "output_POPRES_CH6_20SNPs_#" << i << ".tar.gz " << endl; // output_CH6_0.tar.gz output_hmm25278-S_2.tar.gz output_CH6_clusters_29.tar.gz output_POPRES_CH6_20SNPs_#5.tar.gz
    command << "tar -xvzf " << dirName.str() << "output_#" << i << ".tar.gz " << endl;
    //                                           output_CH6_0M_#45.tar.gz                  output_CH6_50SNPs_#11.tar.gz
    cout << command.str() << endl;               //output_CH6_0M_#16.tar.gz
    system(command.str().c_str());
    stringstream com_del;
    com_del << "rm output/pref_#" << i << ".gamma.txt";
    system(com_del.str().c_str());
    cout << endl << "rozpakowany" << endl;
    // UWAGA dla roznych Y GA ga.. tu!
    GA ga(modelsNo_, maxNoProgressIter_, pCross_, pMutation_, tournamentSize_, B_, modelsFileName, correlationThreshold_, correlationRange_, regionMinCorrelation_, true);
    
    
    if (i == start)
    {
      ga.setRecognisedSNPs();
      ga.initCausalPost( mapSNPCausal_ind );
      tabCausalPost.resize( mapSNPCausal_ind.size() );
    } 
    cout << "before cluster Power" << endl;
    ga.piMassCluserPOWER(fileName.str(), outFileName, mapSNPCausal_ind, tabCausalPost);
    //ga.calculatePOWER_FDR_clust_max(realSNPs, powerFDR, Pmi_Y, recognizedSNPs, mapSNPCausal_ind, tabCausalPost);
    //calculatePOWER_FDR_clust(mySnps, realSNPs, powerFDR, Pmi_Y, badSNP, recognizedSNPs);
//!!! GA uses this:    calculatePOWER_FDR_clust_max(realSNPs, powerFDR, Pmi_Y, recognizedSNPs, mapSNPCausal_ind, tabCausalPost);
    // calculatePOWER_FDR_clust_max(realSNPs, powerFDR, Pmi_Y, recognizedSNPs, mapSNPCausal_ind, tabCausalPost);
    cout << "piMass Extract" << endl;
    stringstream fileName;
    fileName << "output";
    //!!!! ga.piMassExtract(fileName.str(), outFileName);//, recognizedSNPs);
  }  
  writePosterior((parameter.out_file_name + "_post_piMass.txt").c_str(), mapSNPCausal_ind, tabCausalPost, end);
  /*
  stringstream ssClust;
  long double sum;
  for (int i = 0; i < tabClust.size(); ++i)
  {
    ssClust << "[" << realSNPs__[i] << "]" << endl;
    sum = 0.0L;
    for (set<snp_index_t>::iterator it = tabClust[i].begin(); it != tabClust[i].end(); ++it)
    {
      ssClust << *it << "\t" << setprecision(15) <<  clustPmi_Y[*it]/end << endl;
      sum += clustPmi_Y[*it];
    }
    ssClust << "average a post. cluster: " << setprecision(15) << sum/(tabClust[i].size()*end) << endl;
    ssClust << "sum a post. cluster: " << setprecision(15) << sum/(end) << endl;
    ssClust << endl;
  }
  */
  ofstream  recoFile;
  /*
  recoFile.exceptions ( ofstream::eofbit | ofstream::failbit | ofstream::badbit ); // checks if poolfile can be written
  try
  {
    recoFile.open((outFileName + "piMass_CLUSTERS.txt").c_str(),  fstream::out | fstream::trunc );
    cout << "Zapisane w " << (outFileName + "piMass_CLUSTERS.txt").c_str() << endl;
    //recoFile << ssClust.str() << endl;
    recoFile.flush();
    recoFile.close();
  }
  catch (ofstream::failure e)
  {
    cerr << "Could not write a piMass reco file" <<endl;
    exit(-1);
  }
  */
  recognizedSNPs = GA::getRecognisedSNPs_piMass();
  
  stringstream ss;
  for (map<snp_index_t, int>::iterator it = recognizedSNPs.begin(); it != recognizedSNPs.end(); ++it)
  {
    ss << (*it).first << "\t" << (*it).second << endl;
  }
  
  recoFile.exceptions ( ofstream::eofbit | ofstream::failbit | ofstream::badbit ); // checks if poolfile can be written
  try
  {
    recoFile.open((outFileName + "piMassModels_reco.txt").c_str(),  fstream::out | fstream::trunc );
    recoFile << ss.str() << endl;
    recoFile.flush();
    recoFile.close();
  }
  catch (ofstream::failure e)
  {
    cerr << "Could not write a piMass reco file" <<endl;
    exit(-1);
  }
}  

/**
 * @brief Runs Genetics Algorithm as a model selection method.
 * @param outNo - it is number of output. It is used to run GA in a loop and save all the results in different files.
 * @param modelsFileName - it is the file name which contains information about an initial population. 
 *                         If a modelsFileName is empty GA creates a new initial population
 * 
 */
void runGA(int outNo, const string &modelsFileName)
{
  cerr << "runGA" << endl;
  unsigned int modelsNo_ = parameter.modelsNo;
  unsigned int maxNoProgressIter_ = parameter.maxNoProgressIter; 
  double pCross_ = parameter.pCross;                             
  double pMutation_ = parameter.pMutation;                       
  unsigned int tournamentSize_ = parameter.tournamentSize;       
  double correlationThreshold_ = parameter.correlationThreshold; 
  int correlationRange_ = parameter.correlationRange;           
  double regionMinCorrelation_ = parameter.regionMinCorrelation;
  int B_ = parameter.B;
  string old_out_file_name;
  
  cout <<   "modelsNo_ " <<  parameter.modelsNo << endl
       << "maxNoProgressIter_ = " << parameter.maxNoProgressIter << endl
  << "pCross_ = " << parameter.pCross << endl
  << "pMutation_ = " << parameter.pMutation << endl                      
  << "tournamentSize_ = " << parameter.tournamentSize << endl
  << "correlationThreshold_ = " << parameter.correlationThreshold << endl
  << "correlationRange_ = " << parameter.correlationRange << endl
  << "regionMinCorrelation_ = " <<  parameter.regionMinCorrelation << endl
  << "B_ = " << parameter.B << endl
  << "causalModelFilename = \"" << parameter.causalModelFilename << "\"" << endl;
  
/*
[genetic_algorithm]
B = 10
modelsNo = 10
maxNoProgressIter = 1000
pCross = 0.9
pMutation = 0.05
tournamentSize = 2
correlationThreshold = 0.5
correlationRange = 50
regionMinCorrelation = 0.001
# WARNING!!! An empty path "" runs GA on real data. Set a path of causal model to run a simulation.
causalModelFilename = "../../data/CH6/CH6_20_SNPs.mod.orig"

 */  
//   return 0;
  
  
  if (outNo >= 0)  
  {
    old_out_file_name = parameter.out_file_name;
  }

  const time_t time_start = time(NULL);
  stringstream ss;
  timespec ts;
  clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &ts);
  time_t now;
  time(&now);
  srand(now);

  map<snp_index_t, int> mapSNPCausal_ind; 
  vector< multiset<long double> > tabCausalPost; 
  vector< multiset<long double> > tabCausalPost_b; 

  GA ga(modelsNo_, maxNoProgressIter_, pCross_, pMutation_, tournamentSize_, B_, modelsFileName, correlationThreshold_, correlationRange_, regionMinCorrelation_);
  ga.run();

  ss.str("");
  ss.clear();
  ss << "GA time: ";
  ga.writePoolToFile(ss);
  
  ga.initCausalPost( mapSNPCausal_ind );
  tabCausalPost.resize( mapSNPCausal_ind.size() );
  tabCausalPost_b.resize( mapSNPCausal_ind.size() );
  
  stringstream ssModels;
  ssModels << "";
  ga.computePosteriorProbability(ssModels, mapSNPCausal_ind, tabCausalPost, tabCausalPost_b);//, minPosterior);
  writePosterior((parameter.out_file_name + "_post_12.txt").c_str(), mapSNPCausal_ind, tabCausalPost, 1);

  const time_t time_end = time(NULL);          //get current calendar time
/*
  cout << "start at " << timeStemple(time_start) << endl;
  cout << "start at " << time_start << endl;
  cout << "end at " << timeStemple(time_end) << endl;
  cout << "end at " << time_end << endl;
*/  
  ss << (time_end > time_start? sec2time(time_end - time_start) : sec2time(24 * 3600 + time_end - time_start));
  cout << ss.str() << endl;
  ga.writePoolToFile(ss);

  ofstream  Pmi_YsortFile;
  Pmi_YsortFile.exceptions ( ofstream::eofbit | ofstream::failbit | ofstream::badbit ); 
  try
  {
    Pmi_YsortFile.open( ( parameter.out_file_name + "_PjMi_YsortFile" /*+ int2str(parameter.in_values_int)*/ + ".txt" ).c_str(),  fstream::out | fstream::trunc );
    Pmi_YsortFile << ss.str() << endl << ssModels.str() << endl;
    Pmi_YsortFile.flush();
    Pmi_YsortFile.close();
  }
  catch (ofstream::failure e)
  {
    cerr << "Could not write PMi_Y-sorted File" <<endl;
  }
  cout << "Posterior probalibities of SNPs are in the file: " 
       <<  ( parameter.out_file_name + "_PjMi_YsortFile" /*+ int2str(parameter.in_values_int)*/ + ".txt" ).c_str() << endl;
  cout << "A report is in the file: " <<  ( parameter.out_file_name + "_recognized_SNPs" /*+ int2str(parameter.in_values_int)*/ + ".txt" ).c_str() << endl;  
  printLOG(ss.str());   
  if (outNo >= 0)  
  {
    parameter.out_file_name = old_out_file_name;
  }
}


/*
 * Reads the pool and calculates results
 * UWAGI 
 */
void runPoolReader(const string & pool_filename)
{
  cout << "POOL READER" << endl;
  unsigned int modelsNo_ = 1;
  unsigned int maxNoProgressIter_ = parameter.maxNoProgressIter; //1000;
  double pCross_ = parameter.pCross;                             // 0.65;
  double pMutation_ = parameter.pMutation;                       // 0.05;
  unsigned int tournamentSize_ = parameter.tournamentSize;       //2;
  double correlationThreshold_ = parameter.correlationThreshold; //0.5;
  int correlationRange_ = parameter.correlationRange;            // 500  
  int B = 1;
  string modelsFileName = "";
  int real_modelsNo = parameter.modelsNo;
  double regionMinCorrelation_ = parameter.regionMinCorrelation;
  
//  multiset<long double> *tab_Prob = new multiset<long double>[60];
//  multiset<long double> *tab_Prob_b = new multiset<long double>[60];
  
  //map<snp_index_t, int> recognizedSNPs;   // liczebność snp'ow w i-iteracjach
  
  map<snp_index_t, int> recognizedSNPs_Region;
  map<snp_index_t, int> recognizedSNPs_bestGA;
  map<snp_index_t, int> recognizedSNPs_mosgwa;
  map<snp_index_t, int> recognizedSNPs_posterioriModel;
  map<snp_index_t, int> recognizedSNPs_clusterMax;  
  map<snp_index_t, int> recognizedSNPs_clusterSum;  

  
  /*
  map<int, int> cs;
  int realSNPs__[] = {11, 802, 1601, 2399, 3201, 4000, 4799, 5599, 6402, 7201, 8001, 8801, 9601, 10399, 11199, 12002,
                    12801, 13603, 14401, 15205, 16004, 16803, 17602, 18402, 19196, 20001, 20801, 21599, 22402, 23201};
  int realSNPsize = sizeof(realSNPs__)/sizeof(realSNPs__[0]);  
  for (int i = 0; i < realSNPsize; ++i)
    cs.insert(pair<snp_index_t, int>(realSNPs__[i], i));  // snp_id -> indeks w tablicy, gdzie jest klaster dla tego snp
  vector<multiset<snp_index_t> > tabClust;
  tabClust.resize(30);
  map<snp_index_t, long double> clustPmi_Y;  // suma post. dla snpów z klastra
  */
//  long double minPosterior = 0.01;
  
  vector< multiset<long double> > tabCausalPost;  // zawiera prawd. poster. odkrytych snpów przyczynowych
  vector< multiset<long double> > tabCausalPost_b;  // zawiera prawd. poster. odkrytych snpów przyczynowych
  map<snp_index_t, int> mapSNPCausal_ind; // zawiera odwzorowanie snp_id -> ind w tabCausalPost
  int start = 1,
       end = 100;
  // UWAGA dla tych samych Y GA poza petla
  // GA ga(modelsNo_, maxNoProgressIter_, pCross_, pMutation_, tournamentSize_, B, modelsFileName, correlationThreshold_, correlationRange_);
  for (int i = start; i <= end; ++i)  // nazwa pliku bez rozszerzenia
  {
    // UWAGA dla każdej iteracji ustawione są inne Y
    //GA ga(modelsNo_, maxNoProgressIter_l, pCross_, pMutation_, tournamentSize_, B, modelsFileName, correlationThreshold_, correlationRange_);
    parameter.in_values_int = i;
    cout << "parameter.in_values_int: " << parameter.in_values_int << endl; 

    const bool statisticsOnly = true;
//   char c; cout << "statisticsOnly: " << statisticsOnly << ", press a key.."; cin >> c;
    string old_out_file_name = parameter.out_file_name;
    stringstream sOutFileName;
    sOutFileName << parameter.out_file_name << "_#" << i;
    parameter.out_file_name = sOutFileName.str();
    
    
    GA ga(modelsNo_, maxNoProgressIter_, pCross_, pMutation_, tournamentSize_, B, modelsFileName, correlationThreshold_, correlationRange_, regionMinCorrelation_, statisticsOnly);
    // dla różnych Y GA w petli
    if (i == start)
    {
      ga.setRecognisedSNPs();
      ga.initCausalPost( mapSNPCausal_ind );
      tabCausalPost.resize( mapSNPCausal_ind.size() );
      tabCausalPost_b.resize( mapSNPCausal_ind.size() );
      //ga.makeClusers(tabClust);
    }
    
    stringstream ssTar;
    //ssTar << "tar -xzf Imp_hmm25278-S_sort_max_" << i << ".tar.gz" << endl;  //Imp_hmm25278-S_sort_max_95.tar.gz
    //system(ssTar.str().c_str());
    stringstream ssPool;
   
//    ssPool << pool_filename << i << "_pool.txt";  // aktualne
    ssPool << pool_filename << i << "_pool_30K.txt";  // aktualne
    //ssPool << pool_filename << i << "_pool_#hmm34610-S.txt";  //hmm32074-S  // out_1_pool_#hmm26651-S.txt
    cout << ssPool.str() << endl;
    stringstream sGA_Time;
//    char chr;
    //cout << "ga.poolReader() in. Press a key.."; cin >> chr;
    ga.poolReader(ssPool.str(), sGA_Time, real_modelsNo); 
    stringstream ssModels;
    stringstream ssRm;
    
    ga.computePosteriorProbability(ssModels, mapSNPCausal_ind, tabCausalPost, tabCausalPost_b);//, minPosterior);
    
    parameter.out_file_name = old_out_file_name;  
    stringstream sortFileName;
    ofstream  Pmi_YsortFile;
    try
    {
      sortFileName << parameter.out_file_name << "_PjMi_YsortFile_PR_#" << i << ".txt";
      cout << "sortFileName: " << sortFileName.str() << endl; 
      Pmi_YsortFile.exceptions ( ofstream::eofbit | ofstream::failbit | ofstream::badbit ); // checks if poolfile can be written
      Pmi_YsortFile.open( sortFileName.str().c_str(),  fstream::out | fstream::trunc );
      if ( Pmi_YsortFile.is_open() == false)
      {
        cerr << " Can not open file: " << sortFileName.str() << endl;
        exit(-1);
      }
      Pmi_YsortFile << sGA_Time.str() << ssModels.str() << endl;
      Pmi_YsortFile.flush();
      Pmi_YsortFile.close();
      cerr << sortFileName.str() << " closed!" << endl;
    }
    catch (ofstream::failure e)
    {
      cerr << "Could not write PMi_Y-sorted File in the main() - poolReader:" << sortFileName.str() << endl;
      exit(-1);
    }
    cout << "Posterior probalibities of models are in the file: " <<  sortFileName.str().c_str() << endl;
    //fName << "Real_2013.04.25/hmm34610-S/out_" << i << "_PjMi_YsortFile_#hmm34610-S.txt"; 
//    ga.calculateClusterPosterior(sortFileName.str(), 0.01);
    old_out_file_name = parameter.out_file_name;
    //stringstream sOutFileName;
    sOutFileName << parameter.out_file_name << "_#PR#_" << i;
    parameter.out_file_name = sOutFileName.str();
    cout << "out_file_name: " << parameter.out_file_name << endl;
    //ga.writePoolToFile(sGA_Time);
    parameter.out_file_name = old_out_file_name;
    
    // zapis rozpoznanych klastrów
//    if (i == end-1)
//     ga.saveLabelCount(string("_clusters.txt"));  
  } 
  
//  writePosterior((parameter.out_file_name + "_post_12.txt").c_str(), mapSNPCausal_ind, tabCausalPost, end);
//  writePosterior((parameter.out_file_name + "_post_12_b.txt").c_str(), mapSNPCausal_ind, tabCausalPost_b, end);
  /*
  stringstream ssClust;
  long double sum;
  for (int i = 0; i < tabClust.size(); ++i)
  {
    ssClust << "[" << realSNPs__[i] << "]" << endl;
    sum = 0.0L;
    for (set<snp_index_t>::iterator it = tabClust[i].begin(); it != tabClust[i].end(); ++it)
    {
      ssClust << *it << "\t" << setprecision(15) <<  clustPmi_Y[*it]/end << endl;
      sum += clustPmi_Y[*it];
    }
    ssClust << "average a post. cluster: " << setprecision(15) << sum/(tabClust[i].size()*end) << endl;
    ssClust << "sum a post. cluster: " << setprecision(15) << sum/(end) << endl;
    ssClust << endl;
  }
  */
  // UWAGA zapis dla modelu testowego !!!!!!!!!!!!!!!!!!!!!!
  /*
  ofstream fileProp;
  fileProp.open((parameter.out_file_name + "_posteriori.txt").c_str(),  fstream::out | fstream::trunc );
  for (int i = 0; i < 30; ++i)
  {
    fileProp << "";
  //  cout << realSNPs__[i];
    if (tab_Prob[i].size() > 0)
      for (multiset<long double>::iterator itP = tab_Prob[i].begin(); itP != tab_Prob[i].end(); ++itP)
      {
        fileProp << "\t" << setprecision(10) << *itP;
        cout << "\t" << setprecision(10) << *itP;
      }  
    fileProp << endl;
    //cout << endl;  
  }
  fileProp.close();
  */
  /*
  ofstream fileProp_b;
  fileProp_b.open((parameter.out_file_name + "_posteriori_b.txt").c_str(),  fstream::out | fstream::trunc );
  for (int i = 0; i < 30; ++i)
  {
    fileProp_b << "";
    cout << realSNPs__[i];
    if (tab_Prob[i].size() > 0)
      for (multiset<long double>::iterator itP = tab_Prob_b[i].begin(); itP != tab_Prob_b[i].end(); ++itP)
      {
        fileProp_b << "\t" << setprecision(10) << *itP;
        cout << "\t" << setprecision(10) << *itP;
      }  
    fileProp_b << endl;
    cout << endl;  
  }
  fileProp_b.close();
  */
  
  // zapis rozpoznanych snpów
  
  recognizedSNPs_Region = GA::getRecognizedSNPs_Region();
  recognizedSNPs_bestGA = GA::getRecognizedSNPs_bestGA();
  recognizedSNPs_mosgwa = GA::getRecognisedSNPs_mosgwa();
  recognizedSNPs_posterioriModel = GA::getRecognisedSNPs_posterioriModel();
  recognizedSNPs_clusterMax = GA::getRecognisedSNPs_clusterMax();
  recognizedSNPs_clusterSum = GA::getRecognisedSNPs_clusterSum();
  
  map<snp_index_t, int>::iterator it_bestGA = recognizedSNPs_bestGA.begin();
  map<snp_index_t, int>::iterator it_mosgwa = recognizedSNPs_mosgwa.begin();
  map<snp_index_t, int>::iterator it_posterioriModel = recognizedSNPs_posterioriModel.begin();
  map<snp_index_t, int>::iterator it_clusterMax = recognizedSNPs_clusterMax.begin();
  map<snp_index_t, int>::iterator it_clusterSum = recognizedSNPs_clusterSum.begin();

  stringstream ss;
  ss << "SNP_nr\tRegion\tbestGA\tmosgwa\tposteriori\tcluster_Max\tcluster_Sum" << endl;
  for (map<snp_index_t, int>::iterator it = recognizedSNPs_Region.begin(); it != recognizedSNPs_Region.end(); ++it)
  {
    ss << (*it).first << "\t  \t" << it->second << "\t" << it_bestGA->second << "\t" << it_mosgwa->second << "\t" 
       << it_posterioriModel->second << "\t" << it_clusterMax->second << "\t" << it_clusterSum->second << endl;
    ++it_bestGA;   
    ++it_mosgwa;
    ++it_posterioriModel;
    ++it_clusterMax;
    ++it_clusterSum;
  }
  ofstream  recoFile;
  recoFile.exceptions ( ofstream::eofbit | ofstream::failbit | ofstream::badbit ); // checks if poolfile can be written
  try
  {
    recoFile.open((parameter.out_file_name + "GA_reco.txt").c_str(),  fstream::out | fstream::trunc );
    recoFile << ss.str() << endl;
    recoFile.flush();
    recoFile.close();
  }
  catch (ofstream::failure e)
  {
    cerr << "Could not write a piMass reco file" <<endl;
    exit(-1);
  }
  
  string label = "reco";
  string caption = "The number of detected SNPs during 100 simulations";
  
  ofstream file;
  file.open((parameter.out_file_name + "MOSGWA_reco.tex").c_str(), ios::out);
  file << "\\begin {table}" << endl;
  file << "\\caption {" << caption << "}" << endl;
  file << "{\\label {tab:" << label << "}}" << endl;
  file << "\\begin {tabular}{|l|l|l|l|l|l|l|}" << endl;
  file << "\\hline" << endl;

  it_bestGA = recognizedSNPs_bestGA.begin();
  it_mosgwa = recognizedSNPs_mosgwa.begin();
  it_posterioriModel = recognizedSNPs_posterioriModel.begin();
  it_clusterMax = recognizedSNPs_clusterMax.begin();
  it_clusterSum = recognizedSNPs_clusterSum.begin();
  
  file << "SNP nr & region & best & mosgwa & posteriori & cluster Max & cluster Sum \\\\ \\hline" << endl;  
  for (map<snp_index_t, int>::iterator it = recognizedSNPs_Region.begin(); it != recognizedSNPs_Region.end(); ++it)
  {
    file << (*it).first << " & " << it->second << " & " << it_bestGA->second << " & " << it_mosgwa->second << " & " 
       << it_posterioriModel->second << " & " << it_clusterMax->second << " & " << it_clusterSum->second 
       << "\\\\ \\hline" << endl;
    ++it_bestGA;   
    ++it_mosgwa;
    ++it_posterioriModel;
    ++it_clusterMax;
    ++it_clusterSum;   
  }
  file << "\\end {tabular}" << endl;
  file << "\\end {table}" << endl;

  file.close();
}

void runClusterFromPosterior()
{
  unsigned int modelsNo_ = parameter.modelsNo;                   // 5
  unsigned int maxNoProgressIter_ = parameter.maxNoProgressIter; //1000;
  double pCross_ = parameter.pCross;                             // 0.65;
  double pMutation_ = parameter.pMutation;                       // 0.05;
  unsigned int tournamentSize_ = parameter.tournamentSize;       //2;
  double correlationThreshold_ = parameter.correlationThreshold; //0.7;
  int correlationRange_ = parameter.correlationRange;            // 500
  double regionMinCorrelation_ = parameter.regionMinCorrelation;  // 0.1
//  unsigned int testsNo = parameter.testsNo;
  modelsNo_ = 1;    
  int B_ = 1;
  string modelsFileName = "";
  GA ga(modelsNo_, maxNoProgressIter_, pCross_, pMutation_, tournamentSize_, B_, modelsFileName, correlationThreshold_, correlationRange_, regionMinCorrelation_);
  
  /**
   * Oblicza korelacje pomiędzy snpami. 
   * @param clusterSNP - podajemy snp, ktorych korelacje liczymy
   * @param clusterSNPposterior - nie używwamy tutaj
   */
  // /*
//  int tab[] = {111904, 123932, 154723, 158315, 185298, 202224, 209224, 235362, 256793, 429043, 448797, 457295, 496170, 500290, 508753, 509116, 
//               577468, 583058, 595973, 616958, 655532, 665485, 671221, 698054, 707056, 707388, 707720, 717207}; // hmm25278
  /*int tab[] = {99190, 111904, 123932, 156203, 158315, 168299, 169303, 177247, 177892, 200554, 200770, 209224, 255224, 305928, 445024, 448797, 496170,
               498822, 502670, 506589, 509116, 557255, 562392, 578385, 587213, 668447, 696344, 696474, 697427, 707388, 710333, 716935, 717207, 723921,
               793522, 817812, 931130}; hmm26651
  */             
  /*int tab[] = {123932, 126323, 134932, 154723, 168299, 185298, 185299, 235966, 280336, 329593, 448797, 500290, 508753, 509116, 517132, 545690,
              573946, 577468, 583889, 616958, 617494, 630404, 651597, 665485, 671221, 697427, 717207, 841109};
  */ //hmm32074              
  int tab[] = {99190, 123932, 134950, 168299, 169303, 170530, 178841, 247198, 448797, 500290, 502670, 504497, 509116, 510473, 510501, 545690, 552011,
              695843, 697427, 710333, 715791, 717207, 883890};

  int size = sizeof(tab)/sizeof(tab[0]);
  vector<snp_index_t> clusterSNP(tab, tab + size);
  vector<long double> clusterSNPposterior(size, 1);
  //ga.calculateClusterPosterior(clusterSNP, clusterSNPposterior);
  //return; 
  // */
  
///* --------------- Clusters------------------------------------
 /* Obliczenie prawd. posteriori dla klastrów. Zaremować run i 
  * Wersje od 2013.05.01 automatycznie obliczają clastry
*/  
  for (int i = 1; i <= 1; ++i)
  {
    stringstream fName;
//    fName << "2012.12.28/hmm25278-S/out_" << i << "_PjMi_YsortFile_#hmm25278-S.txt";
//    fName << "2012.12.28/hmm32074-S/out_" << i << "_PjMi_YsortFile_#hmm32074-S.txt";
//    fName << "2012.12.28/hmm34610-S/out_" << i << "_PjMi_YsortFile_#hmm34610-S.txt";
//    fName << "2012.12.28/hmm26651-S/out_" << i << "_PjMi_YsortFile_#hmm26651-S.txt";
//      fName << "Test_2013.04.24/#0/TESTout_model_30SNPsb05_" << i << "_PjMi_YsortFile_#0.txt"; 
    //fName << "../../data/testrun/TESTout_model_30SNPsb05_PjMi_YsortFile_#" << i << ".txt"; // TESTout_model_30SNPsb05_PjMi_YsortFile_#1.txt
    fName << parameter.out_file_name + "_PjMi_YsortFile_#hmm25278-S.txt";
    ga.calculateClusterPosterior(fName.str(), 0.01);
  }  

}

void runGA2nd(int artNo)
{

/*
  cout << "runGA2nd RUNs" << endl;
  unsigned int modelsNo_ = 1;
  unsigned int maxNoProgressIter_ = parameter.maxNoProgressIter; //1000;
  double pCross_ = parameter.pCross;                             // 0.65;
  double pMutation_ = parameter.pMutation;                       // 0.05;
  unsigned int tournamentSize_ = parameter.tournamentSize;       //2;
  double correlationThreshold_ = parameter.correlationThreshold; //0.5;
  int correlationRange_ = parameter.correlationRange;            // 500  
  int B = 1;
  string modelsFileName = "";
  int start = 0,
       end = 100;
  //set<snp_index_t> gaSNPs;
  set<snp_index_t> mySnps;     
  string *gaTime = new string[end];
  double *msc = new double[end];
  map<snp_index_t, int> map_SNP2count_Old, 
                        map_SNP2count_ClustEA;
  vector<snp_index_t> causalSNPs;
  
  vector<TPOWER_FDR > vPOWER_FDR;
  vector<TPOWER_FDR > vPOWER_FDR_ClGA;
//  vector<double> vCausalMSC;
  vPOWER_FDR.resize(end);
  vPOWER_FDR_ClGA.resize(end);
//  vCausalMSC.resize(end);
  
  
  set<snp_index_t> trueSNPs;
  for (int i = start; i < end; ++i)  // nazwa pliku bez rozszerzenia
  {
    trueSNPs.clear();
    cout << "parameter.in_values_int: " << parameter.in_values_int << endl;     
    GA ga(modelsNo_, maxNoProgressIter_, pCross_, pMutation_, tournamentSize_, B, modelsFileName, correlationThreshold_, correlationRange_);
    stringstream ss;   // After reading the data I increase the value of trait_name_in_yvm 
    ss << i;
    //parameter.in_values_name = ss.str(); // declare( "data", "trait_name_in_yvm", in_values_name );
    parameter.in_values_int = i + 1;
    if (i == start)
    {
      ga.setCausalCount(map_SNP2count_Old);
      ga.setCausalCount(map_SNP2count_ClustEA);
      ga.modelReader(parameter.causalModelFilename, causalSNPs);
    }


    gaTime[i] = "";
    mySnps.clear();
    stringstream ssFileName;
    ssFileName << "results_" << int2strPadWith( artNo, 2, '0' ) << "/out_CH6_20_max_" << i << "_PjMi_sortFile_#" << i << ".txt"; // out_CH6_20_max_28_PjMi_YsortFile_#28.txt
    //ssFileName << "results_stw/out_CH6_20_max_" << i << "_Stepwise_#" << i << ".txt"; //out_CH6_20_max_0_Stepwise_#0.txt out_CH6_20_max_10_Stepwise_#10.txt
    string fileName = ssFileName.str();
    //readSortFile(fileName, gaTime[i], mySnps, msc[i]);
    readStepwisetFile(fileName, gaTime[i], mySnps, msc[i]);
    cout << i << ") GA 2nd, size:" << mySnps.size() << ", " << mySnps << endl;
    cout << gaTime[i] << endl;    

    ga.calculatePOWER_FDR(mySnps, causalSNPs, vPOWER_FDR[i].POWER, vPOWER_FDR[i].FDR, vPOWER_FDR[i].FDcount, trueSNPs);
    cout << "trueSNPs: " << trueSNPs << endl;
    for (set<snp_index_t>::iterator it = trueSNPs.begin(); it != trueSNPs.end(); ++it)
      ++map_SNP2count_Old[*it];  
    
    trueSNPs.clear();
//    ga.calculatePOWER_FDR_clustGA(mySnps, causalSNPs, vPOWER_FDR_ClGA[i].POWER, vPOWER_FDR_ClGA[i].FDR, vPOWER_FDR_ClGA[i].FDcount, trueSNPs);
    ga.calculatePOWER_FDR_clustGA_2ndArticle(mySnps, causalSNPs, vPOWER_FDR_ClGA[i].POWER, vPOWER_FDR_ClGA[i].FDR, vPOWER_FDR_ClGA[i].FDcount, trueSNPs);
    cout << "trueSNPs_cl: " << trueSNPs << endl;
    cout << "POWER " << vPOWER_FDR[i].POWER << ", FDR " << vPOWER_FDR[i].FDR << ", FD count " << vPOWER_FDR[i].FDcount
         << ", POWER_cl " << vPOWER_FDR_ClGA[i].POWER << ", FDR_cl " << vPOWER_FDR_ClGA[i].FDR << ", FD_cl count " << vPOWER_FDR_ClGA[i].FDcount << endl;
    for (set<snp_index_t>::iterator it = trueSNPs.begin(); it != trueSNPs.end(); ++it)
      ++map_SNP2count_ClustEA[*it];  
    
  }
  saveMapToFile((parameter.out_file_name + "GA_old_reco" + int2strPadWith(artNo, 2, '0') + ".txt"), map_SNP2count_Old);
  saveMapToFile((parameter.out_file_name + "GA_Clust_EA_reco" + int2strPadWith(artNo, 2, '0') + ".txt"), map_SNP2count_ClustEA);
  savePOWER_FDR_ToFile((parameter.out_file_name + "GA_POWEWR_FDR "+ int2strPadWith(artNo, 2, '0') + ".txt"), vPOWER_FDR, vPOWER_FDR_ClGA, gaTime, msc);
  if (gaTime != 0)
    delete [] gaTime;
  if (msc != 0)
    delete [] msc;
*/  
}

void readSortFile(string fileName, string &gaTime, set<snp_index_t>& gaSNPs, double &msc)
{
  int snp;
  string s;
  int size;
  char chr;
  string line;
  ifstream  poolFile;
//  unsigned int id;
  //double msc;
  stringstream ss;
  poolFile.open( fileName.c_str(), ios::in );
  cout << "Try open the file: " << fileName.c_str() << endl;
  bool isTime = true;
  stringstream sGATime;
  if ( poolFile.is_open() )
  {
    cout << "the file is opened" << endl;
    while ( ! poolFile.eof() ) 
    {
      if (isTime == true)
      {
        chr = poolFile.peek();
        if (chr == 'G')
        {
          poolFile >> s;
          //sGATime << s << " ";
          poolFile >> s;
          //sGATime << s << " ";
          poolFile >> s;
          sGATime << s;
          cout << "time: " << sGATime.str() << endl;
        }   
        isTime = false;
      }  
      do 
      {
        poolFile >> line;
      } while (line.compare("bestGA:") != 0);
      poolFile >> line;     // ga "msc:" 
      poolFile >> msc;      // ga msc
      if (poolFile.fail() == true)
      {
        cerr << "Something gone wrong. Incorect file format. Wrong msc value: " << msc << endl;
        exit(-1);
      }
      poolFile >> chr;   // ,
      poolFile >> line;  // "size:"
      poolFile >> size;  // size
      //cout << line << " " << size << endl;
      if (poolFile.fail() == true)
      {
        cerr << "Something gone wrong. Incorect file format. Wrong the model size: " << size << endl;
        exit(-1);
      }      
      poolFile >> chr;  // ,
      poolFile >> chr;  // " " 
      poolFile >> chr;  // [
      chr = poolFile.peek();   // check for zero model "[]"
      while (chr != ']')
      {
         poolFile >> snp;
         if (poolFile.fail() == true)
         {
           cerr << "Something gone wrong. Incorect file format, snp:" << snp << endl;
           exit(-1);
         }
         gaSNPs.insert(snp);
         chr = poolFile.peek();
         if (chr != ']')
           poolFile >> s;  //  ", "
      }
      break;
    }
    poolFile.close();
    cout << "GA model: " << gaSNPs << endl;
    gaTime = sGATime.str();
  }  
  else
  {
    cerr << "file: " << fileName << " is not open" << endl;
    exit (-1);
  }
}


void readStepwisetFile(string fileName, string &gaTime, set<snp_index_t>& gaSNPs, double &msc)
{
  int snp;
  string s;
  int size;
  char chr;
  string line;
  ifstream  poolFile;
//  unsigned int id;
  //double msc;
  stringstream ss;
  poolFile.open( fileName.c_str(), ios::in );
  cout << "Try open the file: " << fileName.c_str() << endl;
  //bool isTime = true;
  stringstream sGATime;
  int ind = 0;
  if ( poolFile.is_open() )
  {
    cout << "the file is opened" << endl;
    while ( ! poolFile.eof() ) 
    {
      //poolFile >> s;
      //sGATime << s << " ";
      //poolFile >> s;
      //sGATime << s << " ";
      poolFile >> s;
      sGATime << s;
      cout << "timeSTW: " << sGATime.str() << endl;
      poolFile >> line;     // ga "msc:" 
      poolFile >> msc;      // ga msc
      if (poolFile.fail() == true)
      {
        cerr << "Something gone wrong. Incorect file format. Wrong msc value: " << msc << endl;
        exit(-1);
      }
      poolFile >> chr;   // ,
      poolFile >> line;  // "size:"
      poolFile >> size;  // size
      //cout << line << " " << size << endl;
      if (poolFile.fail() == true)
      {
        cerr << "Something gone wrong. Incorect file format. Wrong the model size: " << size << endl;
        exit(-1);
      }      
      cout << "msc: " << msc << ", size: " << size << endl;
      poolFile >> chr;  // ,
      cout << ", " << chr << endl;
      poolFile >> chr;  // " " 
      cout << " " << chr << endl;
//      poolFile >> chr;  // [
      cout << " " << chr << endl;
      
      chr = poolFile.peek();   // check for zero model "[]"
      cout << chr << endl;
      while (chr != ']')
      {
         poolFile >> snp;
         if (poolFile.fail() == true)
         {
           cerr << "Something gone wrong. Incorect file format, snp:" << snp << ", ind: " << ind << endl;
           exit(-1);
         }
         ++ind;
         gaSNPs.insert(snp);
         chr = poolFile.peek();
         if (chr != ']')
           poolFile >> s;  //  ", "
      }
      break;
    }
    poolFile.close();
    cout << "STW model: " << gaSNPs << endl;
    gaTime = sGATime.str();
  }  
  else
  {
    cerr << "file: " << fileName << " is not open" << endl;
    exit (-1);
  }
}

/**
 * 
 */
void saveMapToFile(const string &fileName, map<snp_index_t, int>& map_SNP2count)
{
  ofstream  recoFile;
  recoFile.exceptions ( ofstream::eofbit | ofstream::failbit | ofstream::badbit ); // checks if poolfile can be written
  try
  {
    recoFile.open(fileName.c_str(),  fstream::out | fstream::trunc );
    for (map<snp_index_t, int>::iterator it = map_SNP2count.begin(); it != map_SNP2count.end(); ++it)
      recoFile << it->first << "\t" << it->second << endl;
    recoFile.close();
  }
  catch (ofstream::failure e)
  {
    cerr << "Could not write a file: " << fileName << endl;
    exit(-1);
  }  
}


void savePOWER_FDR_ToFile(const string &fileName, vector <TPOWER_FDR>& vPOWER_FDR, vector <TPOWER_FDR>& vPOWER_FDR_ClGA, string *times, double *msc)
{
  ofstream  recoFile;
  recoFile.exceptions ( ofstream::eofbit | ofstream::failbit | ofstream::badbit ); // checks if poolfile can be written
  try
  {
    recoFile.open(fileName.c_str(),  fstream::out | fstream::trunc );
    recoFile << "time\tPOWER\tFDR\tFD_count\tPOWER_cl\tFDR_cl\tFD_count_cl\tmsc" << endl;
    for (unsigned int i = 0; i < vPOWER_FDR.size(); ++i)
    {
      recoFile << times[i] << "\t"  << vPOWER_FDR[i].POWER << "\t" << vPOWER_FDR[i].FDR << "\t" << vPOWER_FDR[i].FDcount << "\t" 
               << vPOWER_FDR_ClGA[i].POWER << "\t" << vPOWER_FDR_ClGA[i].FDR << "\t" << vPOWER_FDR_ClGA[i].FDcount << "\t" << msc[i] 
               << endl;
    }  
    recoFile.close();
  }
  catch (ofstream::failure e)
  {
    cerr << "Could not write a file: " << fileName << endl;
    exit(-1);
  }    
}

void runMOSGWA(int outNo)
{
 // string logFileName(parameter.out_file_name + /*"_#" + int2str(parameter.in_values_int) + (no >= 0? "_" + int2str(no): "") + */ ".log" );
  /*
  try 
  {
    printLOG( "Start: open log to file \"" + logFileName + "\"" );
    LOG.open( logFileName.c_str(), fstream::out );
  } 
  catch ( ofstream::failure e ) 
  {
    cerr << "Could not open logfile \"" + logFileName + "\"" <<endl;
  }
  LOG.exceptions ( ofstream::eofbit | ofstream::failbit | ofstream::badbit );
*/  
  try {
    
    // Read in data by generating MData object1
    MData data;
    // calculate single marker test (need for modelselection)
    data.calculateIndividualTests();
    
    // complete model selection process
    //OLD data.selectModel();
    
    Model model0( data );
    Model firstmodel(data);
    model0.computeRegression();
    data.setLL0M( model0.getMJC() );
    Model *modelin=&model0;
    data.selectModel(
      modelin,
      parameter.PValueBorder,
      parameter.maximalModelSize,
      Parameter::selectionCriterium_mBIC_firstRound
    );
    modelin->printModel( "First round result" );
    data.selectModel(
      modelin,
      5000,
      parameter.maximalModelSize,
      parameter.selectionCriterium
    );
    
    
  } catch ( const Exception e ) {
    printLOG( e.what() );
  }
/*  
  try {
    printLOG("end");
    LOG.close();
  } catch ( ofstream::failure e ) {
    cerr << "Could not close logfile \"" + logFileName + "\"" << endl;
  }  
*/  
/*  
  
  if (outNo >= 0)  // przy uruchomienie dla tych samych danych wejściowych, "no" oznacza numer symulacji
  {
    stringstream ss;   // After calculations I increase the value of trait_name_in_yvm 
    ss << outNo;
//    parameter.out_file_name = parameter.out_file_name + "_" + ss.str(); // declare( "data", "trait_name_in_yvm", in_values_name );
    cout << "outNo: " << outNo << endl;
  }
  
  stringstream ss;
  timespec ts;
  clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &ts);
  time_t now;
  time(&now);
  cout << "time: " << now << endl;  
  srand(now);
 
  MData data;
  data.calculateIndividualTests();                                                              
  Model m(data);
//  data.selectModel(&m);
  cerr << "NOT WORK!!!" << endl; exit(1);
  
  timespec te;
  int hour, min,nano, sec;
  clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &te);
  time_diff(hour, min, sec, nano, ts, te);  
  cout << "diff: ";
  if (hour > 0)
    cout << hour << ":";
  if (hour > 0 || min > 0)
    cout << min << ":";
  cout << sec << "." << nano << endl;
  ss.str("");
  ss.clear();
  //ss << "MOSGWA time: ";
  if (hour > 0)
    ss << hour << ":";
  if (hour > 0 || min > 0)
    ss << min << ":";
  ss << sec << "," << nano << "\t";
  ss << m << endl;
  
  ofstream  recoFile;
  recoFile.exceptions ( ofstream::eofbit | ofstream::failbit | ofstream::badbit ); // checks if poolfile can be written
  try
  {
    recoFile.open(( parameter.out_file_name + "_Stepwise_#" + int2str(parameter.in_values_int) + ".txt" ).c_str(),  fstream::out | fstream::trunc );
    recoFile << ss.str() << endl;
    recoFile.close();
  }
  catch (ofstream::failure e)
  {
    cerr << "Could not write a file: " << ( parameter.out_file_name + "_Stepwise_#" + int2str(parameter.in_values_int) + ".txt" ).c_str() << endl;
    exit(-1);
  }   
*/  
}

void generateY()
{
  MData data ;
  
  Model firstmodel(data);//an  more or less empty model
  data.calculateIndividualTests();  
//  exit(0);
  snp_index_t tabSNPs [] =   
  
  /*
  {  // 10 SNPs CH6 23171 SNPs
    39, 2392, 5227, 7595, 9949, 12804, 15683, 18022, 20878, 23170
  };
  */
  /*
  { // 20 SNPs CH6 23171 SNPs
       39,   979,  2392,  3313,  4743,  5722,  7107,  8075,  9473, 10868, 
    12325, 13738, 15167, 16122, 17542, 18487, 19900, 20878, 22331, 23170
  };
  */
  
  { // 30 SNPs CH6 23171 SNPs
       39,   506,  1424,  2392,  2851,  3806,  4743,  5227,  6180,  7107, 
     8075,  8515,  9473, 10438, 10868, 11854, 12804, 13269, 14230, 15167, 
    16122, 17096, 18022, 18487, 19451, 20390, 20878, 21810, 22760, 23170
  };
  
  
  /*
  // 50 SNPs CH6 23171 SNPs
  {    39,   506,   979,  1424,  1915,  2392,  2851,  3313,  3806,  4268, 
     4743,  5227,  5722,  6180,  6670,  7107,  7595,  8075,  8515,  9003, 
     9473,  9949, 10438, 10868, 11382, 11854, 12325, 12804, 13269, 13738, 
    14230, 14701, 15167, 15683, 16122, 16579, 17096, 17542, 18022, 18487, 
    18973, 19451, 19900, 20390, 20878, 21325, 21810, 22331, 22760, 23170};
*/
    
    
    
    
    
  // POPRES_CH6 23159 SNPs!!!
  /* 
  {  // 10 SNPs
    37, 2389, 5224, 7591, 9945, 12799, 15674, 18013, 20866, 23158
  };
  */
  /*
  {
       37,   976,  2389,  3310,  4740,  5719,  7104,  8071,  9469, 10864,    // 20 SNPs
    12320, 13730, 15159, 16113, 17533, 18477, 19889, 20866, 22319, 23158
  };
  */
     // 30 SNPs
     /*
  {
       37,   504,  1421,  2389,  2848,  3803,  4740,  5224,  6177,  7104,  
     8071,  8511,  9469, 10434, 10864, 11849, 12799, 13264, 14222, 15159,
    16113, 17087, 18013, 18477, 19440, 20378, 20866, 21798, 22748, 23158

  };
  */
  /*
     {37,   504,   976,  1421,  1912,  2389,  2848,  3310,  3803,  4265,  4740,  5224,  5719,  6177,  6667,  7104,  7591,  8071,  8511,  8999, 
    9469,  9945, 10434, 10864, 11377, 11849, 12320, 12799, 13264, 13730, 14222, 14693, 15159, 15674, 16113, 16570, 17087, 17533, 18013, 18477, 
    18963, 19440, 19889, 20378, 20866, 21313, 21798, 22319, 22748, 23158};   // for model of 50 SNPs
*/    
  const int tabSNPsSize = sizeof(tabSNPs)/sizeof(tabSNPs[0]);
  vector<snp_index_t> sel(tabSNPsSize, 1); 
  for (int i = 0; i < tabSNPsSize; ++i)
    sel[i] = tabSNPs[i];
  
  gsl_vector * v = gsl_vector_alloc (sel.size() + 1);
  
  const int dataSize = tabSNPsSize;
  
  double betas[dataSize + 1] = 
  
    {0.0}; // for a zero model 
  
 /*
 {
    0.0,   // 10 SNPs,  h2 = 0.12907
    0.13, 0.139, 0.148, 0.157, 0.166, 0.175, 0.184, 0.193, 0.202, 0.211
  };
  */
  /*
  {
    0.0,       // 20 SNPs, h2 = 0.18795
    0.05, 0.06, 0.07, 0.08, 0.09, 0.1, 0.11, 0.12, 0.13, 0.14,
    0.15, 0.16, 0.17, 0.18, 0.19, 0.2, 0.21, 0.22, 0.23, 0.24
  };
  */
  /*
  {    // 30 SNPs, h2 = 0.12213

    0.00, 
    0.05, 0.053, 0.056, 0.059, 0.062, 0.065, 0.068, 0.071, 0.074, 0.077, 
    0.08, 0.083, 0.086, 0.089, 0.092, 0.095, 0.098, 0.101, 0.104, 0.107,
    0.11, 0.113, 0.116, 0.119, 0.122, 0.125, 0.128, 0.131, 0.134, 0.137
  };
*/
  
/*
  {    // 30 SNPs, h2 = 0.40280
    0.00, 
    0.05, 0.06, 0.07, 0.08, 0.09, 0.1, 0.11, 0.12, 0.13, 0.14,
    0.15, 0.16, 0.17, 0.18, 0.19, 0.2, 0.21, 0.22, 0.23, 0.24,
    0.25, 0.26, 0.27, 0.28, 0.29, 0.3, 0.31, 0.32, 0.33, 0.34
  };
 */ 


 /* 
  { // 30 SNPs, h2 = 0.75186
    0.00, 
    0.3, 0.31, 0.32, 0.33, 0.34, 0.35, 0.36, 0.37, 0.38, 0.39, 
    0.4, 0.41, 0.42, 0.43, 0.44, 0.45, 0.46, 0.47, 0.48, 0.49,
    0.5, 0.51, 0.52, 0.53, 0.54, 0.55, 0.56, 0.57, 0.58, 0.59
  };
  */
  
  

  
  // for a model of 50 SNPs, h2 = 0.71624
  /*
   {0.0,    
    0.1000, 0.1025, 0.105,  0.1075, 0.110,  0.1125, 0.115,  0.1175, 0.120,  0.1225, 
    0.1250, 0.1275, 0.130,  0.1325, 0.135,  0.1375, 0.140,  0.1425, 0.145,  0.1475, 
    0.1500, 0.1525, 0.155,  0.1575, 0.160,  0.1625, 0.165,  0.1675, 0.170,  0.1725, 
    0.1750, 0.1775, 0.180,  0.1825, 0.185,  0.1875, 0.190,  0.1925, 0.195,  0.1975, 
    0.2000, 0.2500, 0.300,  0.4000, 0.500,  0.6000, 0.700,  0.8000, 0.900,  1.0000};
 */
  int betaSize = sizeof(betas)/sizeof(betas[0]);
  gsl_vector_set (v, 0, betas[0]); //da ist der intercept
  for (int i = 1; i < betaSize; i++) 
    gsl_vector_set (v, i, betas[i]); //sollte das nicht i=0 sein
  firstmodel.addManySNP(sel);
  firstmodel.setBeta(v);
    
  firstmodel.printModel(" ");
//  cout << "Ycontinous()" << endl;
  firstmodel.Ycontinous();//this is for continous variables
//  cout << "Ycontinous() - OK" << endl;
  printLOG("<Hallo das richtige wird"); 
}
