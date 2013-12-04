/********************************************************************************
 *	This file is part of the MOSGWA program code.				*
 *	Copyright ©2011–2013, Erich Dolejsi, Bernhard Bodenstorfer.		*
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

#define CPP11

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
	//for random number generator
const gsl_rng_type * T;
       gsl_rng * r;
       /* create a generator chosen by the
          environment variable GSL_RNG_TYPE */

       gsl_rng_env_setup();

       T = gsl_rng_default;
       r = gsl_rng_alloc (T);

	// print logo on screen
	printStartScreen();
printLOG( "HUBBART " );

	// set parameters ... BB working on this
	parameter.setParameters( argc, argv );
if(1==parameter.test) 
	//this creates a new version of yvm files therefor no yvm file should be loaded
	parameter.y_value_extra_file=false;

	// init logging
	string logFileName( parameter.out_file_name + ".logg" );

 try  {
         LOG.open( logFileName.c_str(),  fstream::out ); 
	 printLOG( "Start: open log to file \"" + logFileName + "\"" );
 }
 catch ( ofstream::failure e ) { cerr << "Could not open Logfile \"" + logFileName + "\"" <<endl;exit(1); }
  
 // checks if logfile can be written.
try { LOG.exceptions ( ofstream::eofbit | ofstream::failbit | ofstream::badbit );}
 catch ( ofstream::failure e ) { cerr << "Could not written Logfile \"" + logFileName + "\"" <<endl;exit(1); }

 // Read in data by generating MData object
  MData data ;
  //data.printIndividuals();
 
 // imputate data if nessessary
 if ( ! parameter.imp_is_imputated )
 {     data.imputateMissingValues(); //data.writeBEDfile(); das ist ein schlechte file für plink
       data.writeBEDfile();
 }  
//data.printHyper();
// data.printGenoData();
 	//~ 
 //data.printIndivduals();
 	//~ 
 //data.printY();
 //sollte für Simulation nichts tun
 Model firstmodel(data);//an  more or less empty model
 cerr<<"HUBBART data generator"<<parameter.test<<endl;
 if(3==parameter.test)
 {   	data.printHyper(); 
 }
 if(1==parameter.test) //creates a new Y depending on the list in SNPNames and parameter.binary which says
	                   //in the first case when you WANT a binary trait:
					   //binary=true
					   //and by a 
					   //quantitaive trait
					   //          binary=true
 {

	 
 vector<string> SNPNames ={"SNP_A-2257716","SNP_A-4242616","SNP_A-4254650","SNP_A-2240856","SNP_A-4194567","SNP_A-4278454","SNP_A-2001840","SNP_A-2039683","SNP_A-4269736","SNP_A-2311261","SNP_A-4202772","SNP_A-2054595","SNP_A-2294949","SNP_A-2304163","SNP_A-1961168","SNP_A-1895202","SNP_A-1847062","SNP_A-2187886","SNP_A-2165981","SNP_A-1868008","SNP_A-2126190","SNP_A-2066489","SNP_A-2276551","SNP_A-1920013","SNP_A-1802784","SNP_A-2190800","SNP_A-2162708","SNP_A-1853541","SNP_A-2143131","SNP_A-1903532","SNP_A-2078012","SNP_A-2026715","SNP_A-2234206","SNP_A-1964012","SNP_A-1900553","SNP_A-1966868","SNP_A-1968294","SNP_A-2213530","SNP_A-4247841","SNP_A-1927256","SNP_A-2302761","SNP_A-1815474","SNP_A-1960766","SNP_A-4263654","SNP_A-1977563","SNP_A-2177101","SNP_A-2269286","SNP_A-4280828","SNP_A-4201310","SNP_A-1788695","SNP_A-1832806","SNP_A-1978963","SNP_A-1909904","SNP_A-2271864","SNP_A-2192793","SNP_A-2133703","SNP_A-1906069","SNP_A-2228652","SNP_A-2299116","SNP_A-1812008","SNP_A-1841448","SNP_A-2276504","SNP_A-4240556","SNP_A-2237144","SNP_A-2256716","SNP_A-2121010","SNP_A-2225802","SNP_A-4206207","SNP_A-2210551","SNP_A-2212437","SNP_A-1990968","SNP_A-2278113","SNP_A-4205113","SNP_A-2276279","SNP_A-4242878","SNP_A-4298897","SNP_A-1896334","SNP_A-4280493","SNP_A-2058509","SNP_A-1861258","SNP_A-2102342","SNP_A-4291482","SNP_A-4290594","SNP_A-1962295","SNP_A-1852439","SNP_A-2153198","SNP_A-2205950","SNP_A-1962762","SNP_A-2306175","SNP_A-2275814","SNP_A-1962957","SNP_A-1941808","SNP_A-1917665","SNP_A-1963138","SNP_A-1803104","SNP_A-1804392","SNP_A-2201952","SNP_A-4199461","SNP_A-4272756","SNP_A-2235469","SNP_A-4244004","SNP_A-1901999","SNP_A-4262020","SNP_A-1841438","SNP_A-4253701","SNP_A-1785943","SNP_A-1818110","SNP_A-2174727","SNP_A-1964016","SNP_A-1917711","SNP_A-2135420","SNP_A-2273048","SNP_A-2209499","SNP_A-1964525","SNP_A-2149934","SNP_A-4234178","SNP_A-2107496","SNP_A-2048911","SNP_A-2133535","SNP_A-4204992","SNP_A-2033711","SNP_A-4223190","SNP_A-2089932","SNP_A-1868288","SNP_A-2143974","SNP_A-2273557","SNP_A-1965818","SNP_A-4216892","SNP_A-2073233","SNP_A-1877984","SNP_A-4226486","SNP_A-1871062","SNP_A-1871612","SNP_A-2157764","SNP_A-1789556","SNP_A-1846426","SNP_A-1872383","SNP_A-1941423","SNP_A-4300704","SNP_A-1811033","SNP_A-4202326","SNP_A-2118781","SNP_A-2285005","SNP_A-1966990","SNP_A-2209508","SNP_A-1947508","SNP_A-1877839","SNP_A-4275747","SNP_A-4216428","SNP_A-4276514","SNP_A-1818582","SNP_A-1919535","SNP_A-2062373","SNP_A-1925223","SNP_A-1883155","SNP_A-1818240","SNP_A-2211321","SNP_A-2211213","SNP_A-1812628","SNP_A-1944230","SNP_A-2227457","SNP_A-1840753","SNP_A-4198256","SNP_A-1831064","SNP_A-2214139","SNP_A-1817363","SNP_A-1783095","SNP_A-4294924","SNP_A-2296070","SNP_A-1932145","SNP_A-4209427","SNP_A-2252448","SNP_A-4235502","SNP_A-2035822","SNP_A-4293482","SNP_A-4195755","SNP_A-2175876","SNP_A-2269158","SNP_A-4246383","SNP_A-4247365","SNP_A-2140880","SNP_A-1822945","SNP_A-4240320","SNP_A-2008939","SNP_A-2069242","SNP_A-4209636","SNP_A-1899591","SNP_A-2043416","SNP_A-1974205","SNP_A-1974220","SNP_A-1974221","SNP_A-1942565","SNP_A-1974242","SNP_A-4302781","SNP_A-1922571","SNP_A-4206309","SNP_A-2060855","SNP_A-4238504","SNP_A-2167124","SNP_A-1974502"};
vector< unsigned int> index;
data.findSNPIndex(SNPNames, index);
 firstmodel.addManySNP(index);
 vector<unsigned int> zwischen;
 firstmodel.printStronglyCorrelatedSnps2(0,0.7,zwischen,40);
  firstmodel.printStronglyCorrelatedSnps(0.1); 
 // return(0);
//automatic setting for beta

  gsl_vector * beta = gsl_vector_alloc (index.size()+1); //the intercept is 0!
 gsl_vector_set (beta, 0,-1.2); //da ist der intercept -1.8 für 0.6 bis
//16.04.2012 12:41:39: 0.500157 
// const double up=0.2,lo=0.2; //vorher 1 und 0.3
/*      double step=(up  -lo)/(index.size()-1); //1.5:0.5 in 14 Stufen war 1 große Studie 
  for(int i=1;i<=index.size();++i)
  { gsl_vector_set(beta,i,up-(i-1)*step);
  //DEBUG 
  //cout<< 1.8,1.5,1.0-(i-1)*step<<endl;
  }
 */
		 double be=0.2;
for(int i=1;i<=10;++i)
    gsl_vector_set(beta,i,be); //uniform beta for all 10 variables
for(int i=11;i<=200;i++)
    gsl_vector_set(beta,i, gsl_ran_gaussian (r, be/10));
 firstmodel.setBeta(beta);
 //firstmodel.computeRegression();
 firstmodel.printModel("");
gsl_rng_free (r);
 //this should be tested 
 //newmodel.expXbeta();//this produces the Y for the binary case
 if(false==parameter.binary)
 {firstmodel.Ycontinous();}//this is for continous variables
 else if(true==parameter.binary)
 {firstmodel.Ybinary();// call expXbeta()
firstmodel.printModel("generierter Prototyp");	 return(0);
 vector<bool> YN(data.getSnpNo(),false); //all set to false 0
 firstmodel.getYvec(YN); //I want the Yvec in any case
 // zur Sicherheit wird auch das Individual gesetzt
 data.setYfromMOSGWA(YN);
 for(int i=0;i<data.getIdvNo();i++)
     data.setY(i,YN[i]); 
 }
  	 firstmodel.printModel("generierter Prototyp");
 }

//generates the inputfiles for  HLasso and GWASelect
if(2==parameter.test){
 Model newmodel(data);
 newmodel.printYvec(false); //false or  true/
 //achtung kennt gar kein contiuierliche Variablen
data.printmyY(); 
cerr<<"yvm_local for checking against Yvecout only"<<endl;
vector<bool> t;
printLOG("vor getYvec(t)");
newmodel.getYvec(t); //hier wird kein YVec irgen
     //printYforGWASselect
     //should print a Yout file


//newmodel.printYvec(true);
data.printYforGWASselect();
//printLOG("vor printHyper() ");
data.printHyper();
data.printGWASselect(newmodel);
//printLOG("vor printSNPs()");
}
if(4==parameter.test)
{//für prostata die genotypen
	// nehme die ersten 100 Zeilen aus /dokumente/prostata/merge/prost/ergebnisse1_IT.txt 
	// und generiere ein MATLAB file
vector<string> SNPNames= {"SNP_A-2257716","SNP_A-4242616","SNP_A-4254650","SNP_A-2240856","SNP_A-4194567","SNP_A-4278454","SNP_A-2001840","SNP_A-2039683","SNP_A-4269736","SNP_A-2311261"};
vector< unsigned int> index;
data.findSNPIndex(SNPNames, index);
 firstmodel.addManySNP(index);
 vector<unsigned int> zwischen;
 //automatic setting for beta

  gsl_vector * beta = gsl_vector_alloc (index.size()+1); //the intercept is 0!
 gsl_vector_set (beta, 0,-1.05); //da ist der intercept -1.8 für 0.6 bis
//16.04.2012 12:41:39: 0.500157 
// const double up=0.2,lo=0.2; //vorher 1 und 0.3
/*      double step=(up  -lo)/(index.size()-1); //1.5:0.5 in 14 Stufen war 1 große Studie 
  for(int i=1;i<=index.size();++i)
  { gsl_vector_set(beta,i,up-(i-1)*step);
  //DEBUG 
  //cout<< 1.8,1.5,1.0-(i-1)*step<<endl;
  }
 */
		 double be=0.2;
for(int i=1;i<=10;++i)
    gsl_vector_set(beta,i,be); //uniform beta for all 10 variables

 firstmodel.setBeta(beta);
 //firstmodel.computeRegression();
 firstmodel.printModel("");
gsl_rng_free (r);
 //this should be tested 
 //newmodel.expXbeta();//this produces the Y for the binary case
 if(false==parameter.binary)
 {firstmodel.Ycontinous();}//this is for continous variables
 else if(true==parameter.binary)
 {firstmodel.Ybinary();// call expXbeta()
firstmodel.printModel("generierter Prototyp");	 return(0);
 vector<bool> YN(data.getSnpNo(),false); //all set to false 0
 firstmodel.getYvec(YN); //I want the Yvec in any case
 // zur Sicherheit wird auch das Individual gesetzt
 data.setYfromMOSGWA(YN);
 for(int i=0;i<data.getIdvNo();i++)
     data.setY(i,YN[i]); 
 }
  	 firstmodel.printModel("generiertes");
}
if (0==parameter.test)//eigentliche Studien
{// calculate single marker test (need for modelselection)
	Model model(data);
	// calculate single marker test (need for modelselection)
	data.calculateIndividualTests();
// compute the log-likelihood of the 0 Model
	Model model0( data );
	model0.computeRegression();
	data.setLL0M( model0.getMJC() );
 }
// these lines are from the original main
 try { printLOG("End");
      LOG.close();
 } catch ( ofstream::failure e ) {
      cerr << "Could not close Logfile \"" + logFileName + "\"";}
}
