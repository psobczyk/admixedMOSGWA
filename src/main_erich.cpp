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
#include <boost/algorithm/string.hpp> //für mehrfache Modelle
#include "TIP"
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
        cerr<<"srct-Version"<<TIP<<endl;//the TIP file will created automatically
	// set parameters ... BB working on this
	parameter.setParameters( argc, argv );
	for(int i=0;i<parameter.SNPs.size();i++)
	cerr<<parameter.SNPs[i]<<endl;
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
 
 // imputate data if nessessary
 if ( ! parameter.imp_is_imputated )
 {     data.imputateMissingValues(); //data.writeBEDfile(); das ist ein schlechtes file für plink
       data.writeBEDfile();
       exit(0); //after imputation we stop
 }  
//data.printHyper();
// data.printGenoData();
 	//~ 
 //data.printIndivduals();
 	//~ 
 //data.printY();
 //sollte für Simulation nichts tun
 Model firstmodel(data);//an  more or less empty model
 cerr<<"PARAMETER TEST 0 is for starting modelselection="<<parameter.test<<endl;
 if(3==parameter.test)
 {   	 
	 
  vector< unsigned int> index;
  index.push_back(10);
  index.push_back(50);
  index.push_back(100);
  index.push_back(150);
  index.push_back(230);
  firstmodel.addManySNP(index);

//  vector<unsigned int> zwischen;
// firstmodel.printStronglyCorrelatedSnps2(2,0.2,zwischen,80);
// for(int i=1;i<zwischen.size();++i)
//	 cout<<"zwischen="<<zwischen[i]<<endl;
 
   data.calculateIndividualTests();//needed
  Model model0(data);//needed for check against
  data.setLL0M( model0.getMJC() ); //needed to set
int typeNr=0;
  
  firstmodel.computeMSC(typeNr);
  firstmodel.printModel("with 18 :23 000");
  Model back(data); 
 //  cerr<< firstmodel.makeBackwardStepED( back )<<endl;
  //  cerr<< back.makeBackwardStepED( back )<<endl;
 // back.printModel("Backward");
 // firstmodel.printModelInMatlab();

firstmodel.printModelInMatlab ("vorher");  
firstmodel.replaceModelSNPbyNearCorrelated(typeNr);
firstmodel.printModelInMatlab ("nachher");  

  
  /*DEBUG
   * for(int i=0;i<index.size();++i)
  	cout<<index[i]<<endl;
  	*/
   
  
   }
 if(1==parameter.test) //creates a new Y depending on the list in SNPNames and parameter.binary which says
	                   //in the first case when you WANT a binary trait:
					   //binary=true
					   //and by a 
					   //quantitaive trait
					   //          binary=true
 {


//#if defined (CPP11)
//this is old for checking 2012-08-26
//vector<string> SNPNames = {"SNP_A-2131632",
//"SNP_A-2271680", "SNP_A-4303951", "SNP_A-2110108", "SNP_A-2109718", "SNP_A-2066310", "SNP_A-1833123", "SNP_A-1790769", "SNP_A-4265550", "SNP_A-1843509", "SNP_A-1944795", "SNP_A-1914343", "SNP_A-4201081", "SNP_A-2079135", "SNP_A-1837563", "SNP_A-2114091", "SNP_A-4203993", "SNP_A-2216749"};

	//für Cluster 26-08 
/*	vector<string> SNPNames = {

"SNP_A-2164045",
"SNP_A-1870851",
"SNP_A-2059748",
"SNP_A-1963031",
"SNP_A-1792750",
"SNP_A-2276054",
"SNP_A-1830349",
"SNP_A-2151686",
"SNP_A-2064321",
"SNP_A-1934354",
"SNP_A-4301792",
"SNP_A-1980261",
"SNP_A-4301792",
"SNP_A-2296765",
"SNP_A-1815922",
"SNP_A-1813070",
"SNP_A-4238604",
"SNP_A-2251276",
"SNP_A-1872692"
};
*/
vector <string> SNPNames(parameter.SNPs);
for(int i=0;i<SNPNames.size();i++)
	cerr<<SNPNames[i]<<endl;
vector< unsigned int> index;
data.findSNPIndex(SNPNames, index);
 firstmodel.addManySNP(index);
 vector<unsigned int> zwischen;
// firstmodel.printStronglyCorrelatedSnps2(0,0.7,zwischen,40);
 firstmodel.printStronglyCorrelatedSnps(0.1);  
printLOG("es wird nur printStronglyCorrelatedSnps(0.1) für ein 400 Fenster gerechnet");
 
//automatic setting for beta
for(int i=0;i<index.size();i++)
	cerr<<index[i]<<endl;
cerr<<"index.size()"<<index.size()<<endl;
  gsl_vector * beta = gsl_vector_alloc (index.size()+1); //the intercept is 0!
 gsl_vector_set (beta, 0,parameter.inter);//da ist der intercept -1.8 für 0.6 bis

//16.04.2012 12:41:39: 0.500157 
 const double up=parameter.betaup,lo=parameter.betadown; //vorher 1 und 0.3
 printLOG("up="+ double2str(parameter.betaup) + " lo=" +  double2str(parameter.betadown));
 double step=(up  -lo)/(index.size()-1); //1.5:0.5 in 14 Stufen war 1 große Studie 
  for(int i=1;i<=index.size();++i)
  { gsl_vector_set(beta,i,up-(i-1)*step);
  //DEBUG 
  //cout<< 1.8,1.5,1.0-(i-1)*step<<endl;
  }
//nur für das 2te Szenario mit den verdeckten SNPs, siehe He und Lin
//printLOG("15");
 // gsl_vector_set(beta,15,-0.3002833755);
// printLOG("18"); 
//  gsl_vector_set(beta,17,-0.4110974071); //siehe makebeta
// printLOG("alles ok"); 
 firstmodel.setBeta(beta);
 //firstmodel.computeRegression();
 firstmodel.printModel("");

//generierter Prototyp

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
 //achtung kennt gar kein continuierliche Variablen
data.printmyY(); 
cerr<<"yvm_local for checking against Yvecout only"<<endl;
vector<bool> t;
printLOG("vor getYvec(t)");
//newmodel.getYvec(t); //hier wird kein YVec irgen
     //printYforGWASselect
     //should print a Yout file


//newmodel.printYvec(true);
//data.printYforGWASselect(); //braucht GWASelect nicht !! zumindest komnetiere ich das mal weg
//printLOG("vor printHyper() ");
data.printHyper();
data.printGWASselect(newmodel);
//printLOG("vor printSNPs()");
}
if(4==parameter.test)
{//für prostata die genotypen
	// nehme die ersten 100 Zeilen aus /dokumente/prostata/merge/prost/ergebnisse1_IT.txt 
	// und generiere ein MATLAB file
vector<string> index = {"rs6718629","rs9987306","rs2005384","rs17109239","rs11160443","rs6982012","rs16866998","rs12965947","rs11199434","rs7195562","rs17167034","rs10800436","rs16956904","rs7538501","rs3933923","rs2826778","rs11923244","rs346238","rs935228","rs1656966","rs3763980","rs10833508","rs9829039","rs9786140","rs8000933","rs3759347","rs7431999","rs1681595","rs9786916","rs16839682","rs2581602","rs7731850","rs7924025","rs11018411","rs8071149","rs1520589","rs7867776","rs765557","rs9786111","rs9786824","rs9785897","rs11096435","rs9786773","rs1276034","rs1997342","rs12641380","rs7815779","rs11777846","rs2051288","rs17048832","rs981844","rs537160","rs2260290","rs11736598","rs9327566","rs17333948","rs4574650","rs10423551","rs4737999","rs36103066","rs17036882","rs11636327","rs7637934","rs12136793","rs11164351","rs10903481","rs10760798","rs1175645","rs4782841","rs10840478","rs10782716","rs17149878","rs8076441","rs1036951","rs4765528","rs8133006","rs10867227","rs1332470","rs2392025","rs2406986","rs1490855","rs16947962","rs2009111","rs10959903","rs6849738","rs1940012","rs407894","rs17509794","rs4753969","rs7229918","rs10745408","rs225712","rs4143570","rs1114970","rs723859","rs1381710","rs9348156","rs838267","rs11876341"};


data.printSelectedSNPsInMatlab(index,string("dieErsten100"));
}
if(5==parameter.test)
{// calculate single marker test (need for modelselection)
        Model model0( data );
	model0.computeRegression();
	data.setLL0M( model0.getMJC() );
	Model model(data);
	Model Gmodel(data);
	Model smaller(data);
        Model Gsmaller(data);
	Model osmaller(data);
	Model omodel(data);
	Model Lmodel(data);
	Model Lsmaller(data);	
/*	
		// calculate single marker test (need for modelselection)
	data.calculateIndividualTests();
	vector <string> SNPNames(parameter.SNPs);
       for(int i=0;i<SNPNames.size();i++)
	cerr<<"'"<<SNPNames[i]<<"'"<<endl;
        vector< unsigned int> index;
       data.findSNPIndex(SNPNames, index);
       firstmodel.addManySNP(index);
firstmodel.computeRegression();
firstmodel.computeMSC(0);
firstmodel.printModel("HLmodel","_HLorig");
firstmodel.printStronglyCorrelatedSnps( 0.999, int2str(parameter.in_values_int) + "fullmodel" );
firstmodel.saveguardbackwardstep(smaller);
firstmodel.printModel("HLmodelback","_HLback");
firstmodel.printStronglyCorrelatedSnps( 0.999, int2str(parameter.in_values_int) + "backmodel" );
////////////GWASelect
vector <string> SNPNamesG(parameter.SNPs1);
       for(int i=0;i<SNPNamesG.size();i++)
	cerr<<"'"<<SNPNamesG[i]<<"'"<<endl;
        vector< unsigned int> indexG;
       data.findSNPIndex(SNPNamesG, indexG);
       Gmodel.addManySNP(indexG);
Gmodel.computeRegression();
Gmodel.computeMSC(0);
Gmodel.printModel("Gmodel","_Gorig");
Gmodel.printStronglyCorrelatedSnps( 0.999, int2str(parameter.in_values_int) + "Gfullmodel" );
Gmodel.saveguardbackwardstep(Gsmaller);
Gmodel.printModel("Gmodelback","_Gback");
Gmodel.printStronglyCorrelatedSnps( 0.999, int2str(parameter.in_values_int) + "Gbackmodel" );

//////////////ORIGINALES  MODEL
vector<string> oSNPs={"SNP_A-2020880", "SNP_A-2231834", "SNP_A-2070276", "SNP_A-1815642", "SNP_A-1962975", "SNP_A-2077185", "SNP_A-2175818", "SNP_A-4262563", "SNP_A-1919646", "SNP_A-4260995", "SNP_A-2218280", "SNP_A-1923094", "SNP_A-2232157", "SNP_A-1978011", "SNP_A-2275264", "SNP_A-4218005", "SNP_A-2094426", "SNP_A-4286805", "SNP_A-2220513", "SNP_A-4297940", "SNP_A-1985005", "SNP_A-1869537", "SNP_A-4211904", "SNP_A-1941585"};

	vector <string> oSNPNames(oSNPs);

     for(int i=0;i<oSNPNames.size();i++)
	cerr<<"'"<<oSNPNames[i]<<"'"<<endl;
       vector< unsigned int> oindex;
       data.findSNPIndex(oSNPNames, oindex);
       omodel.addManySNP(oindex);
     
       omodel.computeRegression();
       omodel.computeMSC(0);
       omodel.printStronglyCorrelatedSnps( 0.999, int2str(parameter.in_values_int) + "orig_fullmodel" );
 omodel.printModel("orig_model","_of");


       omodel.saveguardbackwardstep(osmaller);
       omodel.printStronglyCorrelatedSnps( 0.999, int2str(parameter.in_values_int) + "orig_backmodel" );
        omodel.printModel("orig_backmodel","_ob");

//////////////////ORIGINALES MODEL aber ohne die entfernten SNPS aber mit den in SNPs in starken LD Fenster 40 und r2=0.8
//orig snps
//"SNP_A-202088","SNP_A-181564","SNP_A-217581","SNP_A-426256","SNP_A-426099","SNP_A-192309","SNP_A-227526","SNP_A-421800","SNP_A-222051","SNP_A-429794","SNP_A-198500","SNP_A-194158",


*/
vector<string> LDSNPs={"SNP_A-2020880","SNP_A-1815642","SNP_A-2175818","SNP_A-4262563","SNP_A-4260995","SNP_A-1923094","SNP_A-2275264","SNP_A-4218005","SNP_A-2220513","SNP_A-4297940","SNP_A-1985005","SNP_A-1941585", "SNP_A-2104965","SNP_A-2225019","SNP_A-2207150","SNP_A-2089059","SNP_A-2107431","SNP_A-1865752","SNP_A-2086212","SNP_A-2141521","SNP_A-1800559","SNP_A-2301126","SNP_A-1801278","SNP_A-2160387","SNP_A-1930022","SNP_A-2302790","SNP_A-4289013","SNP_A-2134130","SNP_A-2163978","SNP_A-1780694","SNP_A-4277866","SNP_A-2062455","SNP_A-2084317","SNP_A-2033911","SNP_A-1956074","SNP_A-1913873","SNP_A-2217489","SNP_A-2073899","SNP_A-1870534","SNP_A-4291857","SNP_A-2038959","SNP_A-2155630","SNP_A-1953574","SNP_A-1859780","SNP_A-4221344","SNP_A-2197439","SNP_A-2010616","SNP_A-2010617","SNP_A-4215015","SNP_A-2010618","SNP_A-2121761","SNP_A-2079004","SNP_A-2203549","SNP_A-1986190","SNP_A-1825815","SNP_A-1786036","SNP_A-1986192","SNP_A-2214585","SNP_A-2125952","SNP_A-4211683","SNP_A-4211684","SNP_A-4289986","SNP_A-2135688","SNP_A-1987242","SNP_A-4265870"};
	vector <string> LSNPNames(LDSNPs);

     for(int i=0;i<LSNPNames.size();i++)
	cerr<<"'"<<LSNPNames[i]<<"'"<<endl;
       vector<unsigned int> Lindex;
       data.findSNPIndex(LSNPNames, Lindex);
       Lmodel.addManySNP(Lindex);
     
       Lmodel.computeRegression();
       Lmodel.computeMSC(0);
       Lmodel.printStronglyCorrelatedSnps( 0.999, int2str(parameter.in_values_int) + "LD_fullmodel" );
 Lmodel.printModel("LDorig_model",0,"_Lf");


       Lmodel.saveguardbackwardstep(Lsmaller);
       Lmodel.printStronglyCorrelatedSnps( 0.999, int2str(parameter.in_values_int) + "LD_backmodel" );
       Lmodel.printModel("LDorig_backmodel",0,"_Lb");

}
if (6==parameter.test)
{      printLOG("calculate single marker tests only");
	// calculate single marker test (needed for modelselection)
	Model model(data);
	/// calculate single marker test (needed in modelselection for selction of the order new SNP will be added to the model)
	data.calculateIndividualTests();
        //exit(0); //check
}
if (7==parameter.test)
{ //model file names Models.txt wird gelesen und
         string STRING;
	ifstream infile;
	infile.open (parameter.models_file,ios::out);
       if(	infile.is_open())
       {
     
	int i=0;
	string ST;
	getline(infile,ST);
        std::vector<std::string> strs; 
	boost::split(strs, ST, boost::is_any_of("\t "));
        vector<int> len(strs.size());	
        vector<vector<string>> models(strs.size());

	for(int j=0;j<strs.size();j++)
	{
	 cerr<<strs[j]<<endl;
	 len[j]=stoi(strs[j]);
	}
//	return(0);
	vector<vector<string>> mod(strs.size());
	vector<string>line;
        for(int b=0;b<len.size();b++)
	for(int c=0; c<len[b];c++)
	{ //cerr<<"b="<<b<<endl;
          //cerr<<"c="<<c<<endl;
          if(!infile.eof())
	  {
           getline(infile,STRING); // Saves the line in STRING.
	   models[b].push_back(STRING);
	   //cout<<STRING<<endl; // Prints our STRING.
	  }
	  else	; //nothing
	}     
	infile.close();
	for(int i=0;i<models.size();i++)
	{for(int j=0; j<models[i].size();j++)
	cerr<<models[i][j]<<";";
	cerr<<endl;
	}
cerr<<i<<endl;
Model model(data);
	/// calculate single marker test (need for modelselection)
	data.calculateIndividualTests();
	 Model model0( data );
	model0.computeRegression();
	data.setLL0M( model0.getMJC() );

for(int s=0;s<models.size();s++)
 {vector <string> SNPNames(models[s]);
  vector< unsigned int> index;
  data.findSNPIndex(SNPNames, index);
  firstmodel.addManySNP(index);
   
 // Model model0( data );
 //	model0.computeRegression();
 //	data.setLL0M( model0.getMJC() );
  firstmodel.computeRegression();
  firstmodel.computeMSC(0); 
  firstmodel.printModel(int2str(s),0,'_'+int2str(s));
  firstmodel=model;
  index.clear();
}
	return(0);

       } //the file open
}
if (9==parameter.test) //9 generiert nur Xmatrix
{         string STRING;
	ifstream infile;
	infile.open (parameter.models_file,ios::out);
       if(	infile.is_open())
       {
     
	int i=0;
	string ST;
	getline(infile,ST);//erste Zeile enthält die Längen der Modelle
        std::vector<std::string> strs; 
	boost::split(strs, ST, boost::is_any_of("\t "));
        vector<int> len(strs.size());	
        vector<vector<string>> models(strs.size());
//the first line are the lebgth of the models
	for(int j=0;j<strs.size();j++)
	{
	 len[j]=stoi(strs[j]);
	}
	vector<vector<string>> mod(strs.size());
	vector<string>line;
        for(int b=0;b<len.size();b++)
	for(int c=0; c<len[b];c++)
	{ 
          if(!infile.eof())
	  {
           getline(infile,STRING); // Saves the line in STRING.
	   models[b].push_back(STRING);
	  }
	  else	; //nothing:
	}     
	infile.close();
for(int s=0;s<models.size();s++)
 {
	 vector <string> SNPNames(models[s]);
         vector<  unsigned int> index;
         data.findSNPIndex(SNPNames, index);
 
         firstmodel.addManySNP(index);
         firstmodel.printXmat(int2str(s));
	 firstmodel.printModelInMatlab ("selected");
  index.clear();
}

}
} //if 9

if (10==parameter.test)
{
data.printALLmat();}//end 10

if (8==parameter.test)//speziell wenn man glm hat zB.
	//suche zumindest nach Modellen der Größe 35 in parameter.models_file
{ //model file names Models.txt wird gelesen und
         string STRING;
	ifstream infile;
	infile.open (parameter.models_file,ios::out);
       if(	infile.is_open())
       {
     
	int i=0;
	string ST;
	getline(infile,ST);
        std::vector<std::string> strs; 
	boost::split(strs, ST, boost::is_any_of("\t "));
        vector<int> len(strs.size());
	vector<int> dif(strs.size());
               vector<vector<string>> models(strs.size());
			   //für was ist das gut?
int mindiff=100000; //much to big
int minindex=0; //the first index will indeed a possible minimum
	for(int j=0;j<strs.size();j++)
	{
	 //cerr<<strs[j]<<endl;
	 len[j]=stoi(strs[j]);
	 //suche zuminest nach Modellen der Größe 35 
	 int mm=abs(len[j]-35); //this should be search for a min difference to 35
	 if (mindiff>mm)
	 {mindiff=mm;
	         minindex=j;}
	          
	}
	//cerr<<"11"<<endl;
	vector<vector<string>> mod(strs.size());
	vector<string>line;
        //cerr<<len.size()<<endl;
        for(int b=0;b<len.size();b++)
	for(int c=0; c<len[b];c++)
	{           if(!infile.eof())
	  {
           getline(infile,STRING); // Saves the line in STRING.
	   models[b].push_back(STRING);
	   //cout<<STRING<<endl; // Prints our STRING.
	  }
	  else	; //nothing
	}     
	infile.close();
Model model(data);
	/// calculate single marker test (need for modelselection)
	data.calculateIndividualTests();
//for(int s=0;s<models.size();s++)
int s=minindex; 
{
	 vector <string> SNPNames(models[s]);
  vector< unsigned int> index;
  data.findSNPIndex(SNPNames, index);
  firstmodel.addManySNP(index);
  Model model0( data );
  Model smaller(data );
  
  model0.computeRegression();
  data.setLL0M( model0.getMJC() );
  firstmodel.computeRegression();
  firstmodel.computeMSC(0); 
  firstmodel.saveguardbackwardstep(smaller);
  firstmodel.printModel(int2str(s),0,'_'+int2str(s));
  firstmodel=model;
  index.clear();
}
	return(0);

       } //the file open
}

if (0==parameter.test)//eigentliche Studien
{// calculate single marker test (needed for modelselection)
	Model model(data);
	/// calculate single marker test (need for modelselection)
	data.calculateIndividualTests();
//	exit(0); //check 
//	data.printSNPs();vector<string> SNPNames ={ "SNP_A-2110024", "SNP_A-2237238", "SNP_A-4241915", "SNP_A-2047003"};


	vector< unsigned int> index;
		
        Model model0( data );
	model0.computeRegression();
	data.setLL0M( model0.getMJC() );
Model *modelin=&model0;

modelin->printModel("we are before select model",3);//3 ist mBIC
data.selectModel(modelin,parameter.PValueBorder,parameter.expected_causal_snps1,35,3);//5 parameter takes 3 

modelin->printModel("erstes Endergebnis");
data.selectModel(modelin,5000,parameter.ms_ExpectedCausalSNPs);//3 parameter

//modelin->printModel("zweites Endergebnis");
//now search again with the weak criteria
//data.selectModel(modelin,parameter.PValueBorder,parameter.expected_causal_snps1,modelin->getModelSize()+10,3);
//modelin->printModel("drittes Endergebnis",3);
//data.selectModel(modelin,5000,parameter.ms_ExpectedCausalSNPs);
 }
// these lines are from the original main
 try { printLOG("End");
      LOG.close();
 } catch ( ofstream::failure e ) {
      cerr << "Could not close Logfile \"" + logFileName + "\"";}
}
