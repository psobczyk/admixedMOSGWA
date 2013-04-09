#include "Parameter.hpp"
#include "Helpfull.hpp"
#include "Log.hpp"
#include "MData.hpp"
#include "Model.hpp"
#include <gsl/gsl_vector.h>
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
	parameter.setParameters( argc, argv );

	// init logging
	string logFileName( parameter.out_file_name + ".logg" );

 try  { printLOG( "Start: open log to file \"" + logFileName + "\"" );
      LOG.open( logFileName.c_str(), fstream::out ); }
 catch ( ofstream::failure e ) { cerr << "Could not open Logfile \"" + logFileName + "\"" <<endl; }
  
 // checks if logfile can be written.
 LOG.exceptions ( ofstream::eofbit | ofstream::failbit | ofstream::badbit );
 
 // Read in data by generating MData object
 MData data;
 
 // imputate data if nessessary
 if ( ! parameter.imp_is_imputated )
 { data.imputateMissingValues();data.writeBEDfile();} 
//data.printHyper();
// data.printGenoData();
 	//~ 
 //data.printIndivduals();
 	//~ 
 //data.printY();
 Model newmodel(data);
 vector<snp_index_t> sel(6, 1);
 sel[1]=2;
 sel[2]=31;
 sel[3]=142;
 sel[5]=199;
  gsl_vector * v = gsl_vector_alloc (6);
 // for (int i = 1; i < 6; i++) gsl_vector_set (v, i, 1/(double)i);
 for (int i = 1; i < 6; i++) gsl_vector_set (v, i,50);
  newmodel.addManySNP(sel);
//v[0]=intercept
newmodel.setBeta(v);
newmodel.printModel(" ");
 
newmodel.expXbeta();//this produces the Y
//newmodel.Ycontinous();//this is for continous variables
newmodel.Ybinary();// call expXbeta()

newmodel.printYvec();
vector<bool> t;
newmodel.getY(t);
data.printHyper();
data.printGWASselect(newmodel);
//from original main
// calculate single marker test (need for modelselection)                                     
data.calculateIndividualTests();                                                              
                                                                                                       
// complete model selection process                                                           
data.selectModel(); 
// these lines are from the original main
 try { printLOG("End");
      LOG.close();
 } catch ( ofstream::failure e ) {
      cerr << "Could not close Logfile \"" + logFileName + "\"" <<endl;
 }
}
