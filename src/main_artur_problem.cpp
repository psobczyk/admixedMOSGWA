#include "Parameter.hpp"
#include "Helpfull.hpp"
#include "Log.hpp"
#include "MData.hpp"
#include "Model.hpp"
#include <gsl/gsl_vector.h>
//nur für die Zufallspermutation für SNP Auswahl 
     #include <gsl/gsl_rng.h>
     #include <gsl/gsl_randist.h>
     #include <gsl/gsl_permutation.h>
//hier endet der Einschub 
//include <GA.hpp>
/** global variable for Logfile */
ofstream LOG;

/** Global parameter holder */
Parameter parameter;

/** Application entry point */
int main ( const int argc, const char *argv[] ) {

	// print logo on screen
	printStartScreen();

	// set parameters ... BB working on this
 //#include "HeadRD.cc" the *cc have to be set here!
        cout<<"!!!";
	parameter.setParameters( argc, argv );

	// init logging
	string logFileName( parameter.out_file_name + ".logg" );

 try  { printLOG( "Start: open log to file \"" + logFileName + "\"" );
      LOG.open( logFileName.c_str(), fstream::out ); }
 catch ( ofstream::failure e ) { cerr << "Could not open Logfile \"" + logFileName + "\"" <<endl; }
  
 // checks if logfile can be written.
 LOG.exceptions ( ofstream::eofbit | ofstream::failbit | ofstream::badbit );
 
 // Read in data by generating MData object
  MData data ;
 
 // imputate data if nessessary
 if ( ! parameter.imp_is_imputated )
 {     data.imputateMissingValues(); //data.writeBEDfile(); das ist ein schlechte file für plink
       //data.writeBEDfilePlink();//ist dann nicht offiziell impputiert
       data.writeBEDfilePlink();
 }  

 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 /* -----------------------------------------------------------------------------------------------------------
  * 
  *                                        I've got problem with this model
  * 
  */
 
 data.calculateIndividualTests();
 
 /*snp_index_t snpsTable[] = {42, 71, 89, 151, 172, 192, 225, 228, 277, 305, 317, 351, 363, 411, 483, 486, 501, 541, 555, 597, 625, 630, 644, 675, 713, 735, 765, 807, 813, 814, 883, 891, 957, 1004, 1009, 1059, 1149, 1152, 1160, 1197, 1273, 1291};
 
 int snpsTableSize = sizeof(snpsTable)/sizeof(snpsTable[0]);
 
 vector<snp_index_t> snps;
 snps.assign(snpsTable, snpsTable + snpsTableSize);
 */
 vector<snp_index_t> snps= {42, 71, 89, 151, 172, 192, 225, 228, 277, 305, 317, 351, 363, 411, 483, 486, 501, 541, 555, 597, 625, 630, 644, 675, 713, 735, 765, 807, 813, 814, 883, 891, 957, 1004, 1009, 1059, 1149, 1152, 1160, 1197, 1273, 1291};
 Model model(data);
 model.addManySNP(snps);
cout<<"addManySNPs"<<endl; 
 if (model.computeRegression())
 {  
    model.computeMSC();
    cout << endl << endl << "model = " << model << endl;
 }   
 else
 {
   cout << "Can't compute regression!" << endl;
 }   
 
 cout << endl << "I remove first snp from the model" << endl;
 model.removeSNPfromModel(2);
 if (model.computeRegression())
 {  
    model.computeMSC();
    cout << "model = " << model << endl;
 }   
 else
 {
   cout << "Model without first snp. Can't compute regression!" << endl;
 }
 
 
 cout << endl << "I add removed snp to the model" << endl;
 model.addSNPtoModel(snps[2]);
 if (model.computeRegression())
 {  
    model.computeMSC();
    cout << "model = " << model << endl;
 }   
 else
 {
   cout << "Model after add removed snp. Can't compute regression for this model:" << endl;
   cout << model << endl;
   cout << endl << "There is snp 42 is in this model at the end" << endl;
   cout << "There is old value of msc because I can't compute regression for this model. Function model.computeRegression() returns false" 
        << endl << endl;
 }
 return 0;
} 

