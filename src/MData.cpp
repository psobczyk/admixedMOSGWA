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

#include "MData.hpp"
#include <cmath>	//in c++11 for nan(...)
#include "Model.hpp"
#include "GenotypeFreq.hpp"
#include "PermSort.hpp"
#include "io/PlinkInput.hpp"
#include "io/Hdf5Input.hpp"
#include <sstream>
#include <map>
#include <memory>
#include <cmath>	// for nan(...)
#include <cfloat>	// for maximal double
#include <omp.h>
#include <hdf5.h>

using namespace linalg;
using namespace io;

////++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//  class MData
////++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

const Vector MData::getXcolumn ( const size_t snp ) { return xMat.columnVector( snp ); }

const string& MData::getCovMatElementName( const size_t cov ) const {
	return covNames.at( cov );
}

const Vector MData::getCovariateColumn ( const size_t cov ) { return covMat.columnVector( cov ); }

const Vector MData::getY () { return yVec; }

void MData::setLL0M ( const double ll ) {
	loglikelihood0Model_ = ll;
}

void MData::setY ( const size_t index, const double value ) {
	yVec.set(index,value);
	// REMARK<BB>: Does not alter value stored in Individual.
}

/** Default Constructor: reads the input-files, sets parameters, deallambda.hpps with missing phenotypes */
MData::MData ( io::Input *input ) : xMat( 0, 0 ), yVec( 0 ), covMat( 0, 0 ) {

	const bool allocateInput = NULL == input;
	if ( allocateInput ) {
		if ( parameter.in_file_hdf5.empty() ) {
			input = new PlinkInput( parameter.in_files_plink.c_str() );
		if (0==parameter.nSNPKriterium)
                parameter.nSNPKriterium=getSnpNo();
		} else {
			input = new Hdf5Input( parameter.in_file_hdf5.c_str(), parameter.cov_extra_file );
		}
	}

	const size_t
		snps = input->countSnps(),
		idvs = input->countIndividuals(),
		covs = input->countCovariates();
	        //ED setting default Value for nSNPKriterium when not set 
	        if(0==parameter.nSNPKriterium)
                  parameter.nSNPKriterium=snps;
	xMat.exactSize( idvs, snps );
	yVec.exactSize( idvs );
	covMat.exactSize( idvs, covs );
	const SNP * snpArray = input->getSnps();
	for ( size_t snp = 0; snp < snps; ++snp ) {
		snpList.push_back( snpArray[ snp ] );
		Vector xVec = xMat.columnVector( snp );
		input->retrieveGenotypeVector( snp, xVec );
	}
	const Individual * idvArray = input->getIndividuals();
	for ( size_t idv = 0; idv < idvs; ++idv ) {
		idvList.push_back( idvArray[ idv ] );
	}
	input->retrievePhenotypeVector( parameter.in_values_int, yVec );
	const std::string * covariates = input->getCovariates();
	for ( size_t cov = 0; cov < covs; ++cov ) {
		covNames.push_back( covariates[ cov ] );
		Vector covVec = covMat.columnVector( cov );
		input->retrieveCovariateVector( cov, covVec );
	}
	Y_name_ = input->getTraits()[parameter.in_values_int];
	if ( allocateInput ) {
		delete input;
		input = NULL;
	}

	checkYValues();
}

MData::~MData () {
}

size_t MData::getSnpNo () const {
	return snpList.size();
}

size_t MData::getIdvNo () const {
	return idvList.size();
}

size_t MData::getCovNo () const {
	return covNames.size();
}

string MData::getFID ( const size_t index ) const {
	const Individual individual = idvList.at( index );
	return individual.getFamilyID();
}

string MData::getID ( const size_t index ) const {
	const Individual individual = idvList.at( index );
	return individual.getIndividualID();
}

const SNP * MData::getSNP ( const size_t snp ) const {
	return &( snpList.at( snp ) );
}

size_t MData::getOrderedSNP ( const size_t snp ) const {
	return snp_order_.getId( snp );
}

size_t MData::getCaseNo () const {
	return caseNo_;
}

size_t MData::getContNo () const {
	return contNo_;
}

double MData::getLL0M () const {
	return loglikelihood0Model_;
}

double MData::getSnp_order_Value ( const size_t index ) const {
	return snp_order_.getValue( index );
}

void MData::setSingleMarkerTestAt ( const size_t index, const double value ) {
	snpList.at( index ).setSingleMarkerTest( value );
}

void MData::fillSnp_order_Vec ( const size_t snpNo, size_t* SNPList, double* TestStat ) {
	snp_order_.fillVec( snpNo, SNPList, TestStat );
}

/** Remark: case and control counts and other values derived from the input data are not updated here! */
void MData::removeIndividual ( const size_t idv ) {
	idvList.erase( idvList.begin() + idv );
	xMat.removeRow( idv );
	yVec.removeDimension( idv );
	covMat.removeRow( idv );
	covNames.erase( covNames.begin() + idv );
}



/** Check Y-Values and store them in a vector.
 MISSING: to check for missing values and remove then the individual, and the genotype-info */
void MData::checkYValues()
{
	
	int noOfRemovedIndividuals=0;
	parameter.affection_status_phenotype = true;	// intitial setting, Y-values are tested if the are just (0,1) (or missing)
	
	caseNo_=0;
	contNo_=0;
	
	for ( size_t idv = 0; idv < getIdvNo(); ++idv ) {
		const double indPheno = yVec.get( idv );
		// do not consider individuals with missing phenotype
		if ( indPheno == parameter.missing_phenotype_code ) {
			removeIndividual( idv );
			// the current position was removed,
			// so a diffrent individual is now at the i-th position
			--idv;
			++noOfRemovedIndividuals;
		}
		else // no missing phenotype, determine if phenotype is affection (case-control) or quantitative
		{
			if ( indPheno == parameter.case_value) // was set 1 check affection status
			{
                         //yVec.set(idv,1); //internally a case ist 1

				caseNo_++;
			}
			else if ( indPheno == parameter.control_value)  //was set 0
			{
			 //yVec.set(idv,0); //internally a control ist 0	
				contNo_++;

			}
			else // some other number than 0/1 detected, phenotype must be qua
			{ cerr<<"i="<<idv<<" indPheno="<<indPheno<<endl;
				parameter.affection_status_phenotype = false;
			}
		}
	
	}
	
	vector<Individual*>::iterator it_ind;
	// fill Y_vector_ with phenotype, test of cases 
	
	
	if ( parameter.affection_status_phenotype ) {
		if ( caseNo_ + contNo_ == getIdvNo() ) {
			printLOG("Affection status phenotype detected: "+int2str(caseNo_)+" Cases, " + int2str(contNo_) + " Controls, " + int2str(noOfRemovedIndividuals) + " missing");
		}
		else
		{
			printLOG( "Error in affection status: " + int2str(caseNo_) + " Cases, " + int2str(contNo_) + " Controls, " + int2str( getIdvNo() ) + " Total");	
			exit(7);
		}	
	}
	else	
	{
		printLOG( "Quantitive phenotype detected: " + int2str( getIdvNo() ) +" Individuals, " +int2str(noOfRemovedIndividuals) +" Individuals with missing phenotype" );	
	}	

}


void MData::setYfromMOSGWA( const vector<bool>& Y ) {
	for( size_t idv = 0; idv < getIdvNo(); ++idv ) {
		yVec.set( idv, Y[idv] ? 1.0 : 0.0 );
	}
}



/** TESTING */
void MData::printSNPs () {
	for (
		vector<SNP>::const_iterator itsnps = snpList.begin();
		itsnps < snpList.end();
		++itsnps
	) {
		cout << *itsnps << endl;
	}
}



/** TESTING */
void MData::printIndividuals () {
	for (
		vector<Individual>::const_iterator itidvs = idvList.begin();
		itidvs < idvList.end();
		++itidvs
	) {
		cout << *itidvs << endl;
	}
}

void MData::printIndividuals ( ostream &s ) {
	const size_t idvs = getIdvNo();
	for ( size_t idv = 0; idv < idvs; ++idv ) {
		const Individual individual = idvList.at( idv );
		s
			<< individual.getFamilyID() << " "
			<< individual.getIndividualID() << " "
			<< yVec.get( idv ) << endl;
	}
}

//WARNING this is not used for generation of a new yvm file with generated data
void MData::printmyY()
{ cout<<parameter.out_file_name + ".yvm_local";
	ofstream Y;
	Y.exceptions ( ofstream::eofbit | ofstream::failbit | ofstream::badbit );
	try {
		Y.open( ( parameter.out_file_name + ".yvm" ).c_str(), fstream::out );
	        Y<<"FID IID G"<<endl;  
		                     
                 printIndividuals (Y); //ED version
                Y.close();
	} catch ( ofstream::failure e/*xception*/ ) {
		cerr << "Could not write plink  --pheno  file" << ( parameter.out_file_name + ".yvm_local" ).c_str()<< endl;
	}
}


/** printY printY to a file */
void MData::printYforGWASselect () {
	ofstream Y;
	size_t counter =yVec.countDimensions();
	//	cout << counter<<endl;

	Y.exceptions ( ofstream::eofbit | ofstream::failbit | ofstream::badbit );
	try {

		Y.open( ( parameter.out_file_name + "Yextra" ).c_str(), fstream::out );
	} catch ( ofstream::failure e/*xception*/ ) {
		cerr << "Could not open Y for GWASselect file:" << ( parameter.out_file_name + "Yextra" ).c_str()<< endl;
	}

		for ( size_t dim = 0; dim < counter; ++dim ) 
			Y << yVec.get( dim ) << endl;
	try {	Y.close();
		}
	catch ( ofstream::failure e/*xception*/ ) {
		cerr << "Could not closing Y for GWASselect file:" << ( parameter.out_file_name + "Yextra" ).c_str()<< endl;
	}
}


/** TESTING */
/** caco is case control 0 control 1 case*/
void MData::printGenoData ( ostream& hyper, const vector<bool>& sel, const bool caco ) const {
	if ( sel.size() != getIdvNo() ) {
		cerr << "printGenoData should have sel.size() == getIdvNo() !" << endl;
	} else {
		for ( size_t snp = 0; snp < getSnpNo(); ++snp ) {
			hyper << ( 0 < snp ? " " : "" );
			snpList.at( snp ).getSnpId();
			const Vector xVec = const_cast<MData*>( this )->getXcolumn( snp );
			for ( size_t idv = 0; idv < getIdvNo(); ++idv ) {
				if ( caco == sel[idv] ) {
					hyper << " ";
					hyper << xVec.get( idv );
				}
			}
			hyper << endl;
		}
	}
}

void MData::printGenoData () const {
	for ( size_t snp = 0; snp < getSnpNo(); ++snp ) {
		const Vector xVec = const_cast<MData*>( this )->getXcolumn( snp );
		for ( size_t idv = 0; idv < getIdvNo(); ++idv )  {
			cout << ( 0 < idv ? " " : "" );
			cout << xVec.get( idv );
		}
		cout << endl;
	}
}

/** printGenoData prints the whole GenoData Matrix in 
a ofstream hyper 
Example: see printHyper
*/
void MData::printGenoData ( ostream& hyper ) const {
	for ( size_t idv = 0; idv < getIdvNo(); ++idv ) {
		for ( size_t snp = 0; snp < getSnpNo(); ++snp ) {
			const Vector genome = const_cast<MData*>( this )->getXcolumn( snp );
			cout << ( 0 < snp ? " " : "" );
			hyper << genome.get( idv );
		}
		hyper << endl;
	}
}


/** wer braucht die SNP id noch*/
void MData::printSNPId ( ostream& hyper ) const {
	for ( size_t snp = 0; snp < getSnpNo(); ++snp ) {
		hyper << ( 0 < snp ? " " : "" );
		snpList.at( snp ).getSnpId();
	}
	hyper << endl;
}

void  MData::findSNPIndex(vector<string>& SNPNames, vector<unsigned int>& index) const
{ //first sort the vector
 sort(SNPNames.begin(),SNPNames.end());	 //ERICH at first sort the input
 //create a vector of strings from

vector< int> tindex(SNPNames.size(),-999); //bring index to final size -999 is the not available SNP
 vector<string> allSNPs(getSnpNo()); //this is the extension of the vector the SNP names 
                                    // have probably different size.

	//this is maybe a large loop 500 000 is possible
	for( size_t snp = 0; snp < getSnpNo(); ++snp ) {
		allSNPs[snp] = snpList.at( snp ).getSnpId();
	}
 std::vector<int> permutation;
 sortingPermutation(allSNPs, permutation);
	for ( size_t i=0; i < SNPNames.size(); ++i ) {
		for ( size_t j=0; j < allSNPs.size(); ++j ) {
			if ( 0 == SNPNames[i].compare( allSNPs[ permutation[j] ] ) ) {
				tindex[i] = permutation[j];	//permutation beginnt mit 1
	  			break;
			}
		}
	}
	for ( size_t i = 0; i < tindex.size(); ++i ) {
		if ( -999 == tindex[i] ) {
			cout	<< "SNP name "
				<<SNPNames[i]
				<< " At position "
				<< i
				<< " does not exist nicht in this dataset"
				<< endl;
		}
	}
vector<int>::iterator iter,nend;

nend=remove(tindex.begin(),tindex.end(),-999);
for(iter=tindex.begin();iter<nend;++iter)
{if(false) cout<<"Position:"<<setw(6)<<*iter<<" is "<<allSNPs[*iter]<<endl; //DEBUG
 index.push_back(*iter);}

 //remove duplicates
 //that model works
 sort(index.begin(),index.end() );

 index.erase( unique( index.begin(), index.end() ), index.end() );
/*for_each(index.begin(),index.end(),
		boost::lambda::if_then(boost::lambda::_1==-999,
		cout<< boost::lambda::_1 <<"is not available\n"));*/
}

void MData::printHyper() const {
bool	create=false;

ifstream ifile((parameter.singlefile /*out_file_name*/ + "HyperLasso.in").c_str()); //check if I could open the file 
if (!ifile) create=true;
ifile.close(); //close the file
ofstream Hyper;
if(create)
{
Hyper.exceptions ( ofstream::eofbit | ofstream::failbit | ofstream::badbit ); // checks if Logfile can be writt
try
 {Hyper.open( ( parameter.singlefile /*out_file_name*/ + "HyperLasso.in" ).c_str(), fstream::out );

  printSNPId(Hyper);
  printGenoData(Hyper);

  Hyper.close();
 }
catch (ofstream::failure e)//why e?
 {
 cerr << "Could not write Hyper-Lasso file" <<endl;
 }
} //no create needed
}


void MData::printGWASselect(Model & newmodel) const {
//#pragma omp parallel shared(GWASselectcases,GWASselectcontrol,sel) private(i) //has to be alterd to fit
bool createcas=false, createcont=false;
vector<bool> sel;
	newmodel.getYvec(sel);
ifstream ifile(( parameter.out_file_name + "GWASelect.cont" ).c_str());
   if (!ifile)
	createcont=true;
ifile.close(); //close the file
if(createcont)	
{ofstream GWASselectcontrol;
// can be written
GWASselectcontrol.exceptions ( ofstream::eofbit | ofstream::failbit | ofstream::badbit ); // checks if Logfile can be written
/*omp parallel section*/
//#pragma omp section
try
{GWASselectcontrol.open( ( parameter.out_file_name + "GWASelect.cont" ).c_str(), fstream::out ); //changed file names 

  printGenoData(GWASselectcontrol,sel,false);

  GWASselectcontrol.close();
 }
catch (ofstream::failure e)//why e?
 {
 cerr << "Could not write GWASselect file" <<endl;
 }	
}

ifstream ifile1(( parameter.out_file_name + "GWASelect.cas" ).c_str());
   if (!ifile1)
	createcas=true;
ifile1.close(); //close the file
if (createcas)
{//#pragma omp section
ofstream GWASselectcases;
GWASselectcases.exceptions ( ofstream::eofbit | ofstream::failbit | ofstream::badbit ); // checks if Logfile
try
 {GWASselectcases.open( ( parameter.out_file_name + "GWASelect.cas" ).c_str(), fstream::out ); //changed file names 

  printGenoData(GWASselectcases,sel,true);

  GWASselectcases.close();
 }

catch (ofstream::failure e)//why e?
 {
 cerr << "Could not write GWASselect file" <<endl;
 }	
}
}


void MData::printSelectedSNPsInR ( vector<string> SNPList ) const {
	ofstream	SNPL;

	SNPL.exceptions ( ofstream::eofbit | ofstream::failbit | ofstream::badbit ); // checks if Logfile can be writt
	try
	{
		SNPL.open( ( parameter.out_file_name + "_SNPList.txt" ).c_str(), fstream::out );
		vector<string>::iterator it_SNPList;
		
		for (it_SNPList = SNPList.begin(); it_SNPList < SNPList.end(); it_SNPList++)
		{
			size_t snp;
			for (
				snp = 0;
				snp < getSnpNo() && snpList.at( snp ).getSnpId() != *it_SNPList;
				++snp
			);

			if ( snp < getSnpNo() ) {
				SNPL << *it_SNPList <<" <- c(";
				const Vector xVec = const_cast<MData*>( this )->getXcolumn( snp );
				for ( size_t idv = 0; idv < getIdvNo(); ++idv ) {
					SNPL << ( 0 < idv ? "," : "" );
					SNPL << xVec.get( idv );
				}
				SNPL<< ")"<< endl;
			}
			else
			{
				cout << "SNP " << *it_SNPList << "not found";	
			}
		}
		
		for ( size_t cov = 0; cov < getCovNo(); ++cov ) {
			SNPL << getCovMatElementName( cov ) << " <- c(";
			const Vector covVec = const_cast<MData*>( this )->getCovariateColumn( cov );
			for ( size_t idv = 0; idv < getIdvNo(); ++idv ) {
				SNPL << ( 0 < idv ? "," : "" );
				SNPL << covVec.get( idv );
			}
			SNPL<< ")"<< endl;	
		}

		SNPL << "Y" <<" <- c(";
		for ( size_t idv = 0; idv < getIdvNo(); ++idv ) {
			SNPL << ( 0 < idv ? "," : "" );
			SNPL << yVec.get( idv );
		}
		SNPL<< ")"<< endl;

		SNPL << "Intercept" <<" <- c(";
		for ( size_t idv = 0; idv < getIdvNo(); ++idv ) {
			SNPL << ( 0 < idv ? "," : "" );
			SNPL << 1;
		}
		SNPL<< ")"<< endl;

		//~ SNPL << "summary(lm(Y~Dummy1 + Dummy2 + Dummy3 +";
		SNPL <<"X <- cbind( Intercept";
		for (it_SNPList = SNPList.begin(); it_SNPList < SNPList.end(); it_SNPList++)
		{
			//~ if ( it_SNPList != SNPList.begin() )
			// TODO: and how about covariates?
			{SNPL << ",";}	// TODO<BB> Why {}?
			SNPL << *it_SNPList;
		}
		SNPL << ")"<<endl;
		
		SNPL.close();
		printLOG( "Written R-File \"" + parameter.out_file_name + "_SNPList.txt\"." );
	}
	catch (ofstream::failure e)
	{
		cerr << "Could not write R-File" <<endl;
	}	
}

//everything in h5

bool MData::printALLmat(const string& extra){ } 
//	char mat [getIdvNo()][getModelSize()];
//const int64_t  IndN=getIdvNo(),SnpN=getSnpNo();
//const int64_t sA= IndN*SnpN;
//
//	char * mat= (char *)malloc(sA);
//	// genotype-data
//	for(int64_t i=0;i<SnpN;++i){
//	   const Vector genotypes = getX().columnVector( i );//getX unbekannt
//	   for (int64_t   j = 0; j <IndN ; ++j ) {
//		mat[i*IndN+j]=
//		 genotypes.get(j );}
//    //init of genotype matrix
//}	
////from h5 compress
//    hid_t    file_id, dataset_id, dataspace_id; /* identifiers */
//    hid_t    plist_id; 
//
//    size_t   nelmts;
//    unsigned flags, filter_info;
//    H5Z_filter_t filter_type;
//
//    herr_t   status;
//    hsize_t  dims[2];
//    hsize_t cdims[2];
// 
//    int      idx;
//    int      i,j, numfilt;
//   // int      buf[DIM0][DIM1];
//  //  int      rbuf [DIM0][DIM1];
//  /* Uncomment these variables to use SZIP compression 
//    unsigned szip_options_mask;
//    unsigned szip_pixels_per_block;
//    */
//
//    /* Create a file.  */
//string filename=	parameter.models_file+extra+".h5";
//cerr<<filename<<endl;
//file_id = H5Fcreate (filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
//
//
//    /* Create dataset "Compressed Data" in the group using absolute name. 
//	   getIdvNo()][getModelSize()]*/
//    dims[1] = getIdvNo();//verdreht!
//    dims[0] = getSnpNo();
//    dataspace_id = H5Screate_simple (2, dims, NULL); //2 ist der Rang
//
//    plist_id  = H5Pcreate (H5P_DATASET_CREATE);
//
//    /* Dataset must be chunked for compression */
//    cdims[1] = 20;//verdreht!
//    cdims[0] = getSnpNo();
//    status = H5Pset_chunk (plist_id, 2, cdims);
//
//    /* Set ZLIB / DEFLATE Compression using compression level 6.
//     * To use SZIP Compression comment out these lines. 
//    */ 
//    status = H5Pset_deflate (plist_id, 6); 
//
//    /* Uncomment these lines to set SZIP Compression 
//    szip_options_mask = H5_SZIP_NN_OPTION_MASK;
//    szip_pixels_per_block = 16;
//    status = H5Pset_szip (plist_id, szip_options_mask, szip_pixels_per_block);
//    */
//    /*nun nicht H5T_STD_I32BE sondern I8BE)*/
//    dataset_id = H5Dcreate2 (file_id, "X", H5T_STD_I8BE, 
//                            dataspace_id, H5P_DEFAULT, plist_id, H5P_DEFAULT); 
//
//    status = H5Dwrite (dataset_id, H5T_NATIVE_CHAR, H5S_ALL, H5S_ALL, H5P_DEFAULT, mat);
//
//    status = H5Sclose (dataspace_id);
//    status = H5Dclose (dataset_id);
//    status = H5Pclose (plist_id); //für die cunchs 
////%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
////%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//
//
//    status = H5Fclose (file_id);
//
//
//}


























/** printSelectedSNPsInMatlab create an *Octave.h5 file with the selected SNPs (as a vector of strings)  and their X and of course the phenotype Y
 * additionaly it writes a */ 
void MData::printSelectedSNPsInMatlab ( vector<string> SNPList , string extra) const  {
	ofstream	SNPL;
         {        	ofstream	SNPL;
		 //create the same in hdf
		 hid_t file,fid,dataset,space,/*dset,memtype,*/ filetype,props;
                 herr_t status;
                 hsize_t dim[]={SNPList.size(),getIdvNo()};//transponiert
                 hsize_t di[]={SNPList.size()};
		 double  zwischen[getIdvNo()*SNPList.size()];
		 double  Y[getIdvNo()];
		 vector<string>::iterator it_SNPList;
		 int i;
printLOG( parameter.out_file_name + extra + "Octave.h5");
file=H5Fcreate((parameter.out_file_name + extra + "Octave.h5").c_str(),H5F_ACC_TRUNC,H5P_DEFAULT,H5P_DEFAULT); //standart
const char *S[SNPList.size()];
	filetype = H5Tcopy (H5T_C_S1);
        status = H5Tset_size (filetype,H5T_VARIABLE);
	/* create the datatype properties                                 */
  	props = H5Pcreate (H5P_DATASET_CREATE);
	space = H5Screate_simple (1, di, NULL);
        /*
         * Create the dataset and write the variable-length string data to
         * it.
         */
	dataset = H5Dcreate2( file, "SNPnames", filetype, space, H5P_DEFAULT, props, H5P_DEFAULT );
       //ED funktioniert nicht; dataset = H5Dcreate( file, "SNPnames", filetype, space, props );
	for (it_SNPList = SNPList.begin(),i=0; it_SNPList < SNPList.end(); it_SNPList++,i++)
		{
                        S[i]=it_SNPList->c_str();
				
		}//that should have done the job
unsigned int jj=0;
vector<unsigned int> index;
findSNPIndex(SNPList, index);
//sortd according the names1
unsigned int ii=0;
	for (jj=0; jj < index.size(); ++jj)
	{
		const Vector tmp = const_cast<MData*>( this )->getXcolumn( index[jj] );
		for(ii=0;ii<getIdvNo();++ii)
		       	zwischen[jj*getIdvNo()+ii] = tmp.get(ii);
	}	
        for ( size_t idv = 0; idv < getIdvNo(); ++idv ) {
		Y[idv] = yVec.get( idv );
	}

      	status = H5Dwrite (dataset,filetype , H5S_ALL, H5S_ALL, H5P_DEFAULT, S);
	status = H5Tclose(filetype);	
	status = H5Pclose(props);
	status = H5Sclose(space);
	status = H5Dclose(dataset);
	fid=H5Screate_simple(2,dim,NULL);
//dataset = H5Dcreate (file, "X",filetype, space, H5P_DEFAULT,props,H5P_DEFAULT);	
	dataset = H5Dcreate2( file, "X", H5T_NATIVE_DOUBLE, fid, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );
//ED funktioniert nicht ;dataset=H5Dcreate( file, "X", H5T_NATIVE_DOUBLE, fid, H5P_DEFAULT );
	status=H5Dwrite(dataset,H5T_NATIVE_DOUBLE,H5S_ALL,H5S_ALL,H5P_DEFAULT,zwischen);
	status=H5Dclose (dataset);
        status=H5Sclose (fid);
dim[0]=getIdvNo();
dim[1]=1;
fid=H5Screate_simple(2,dim,NULL);
	dataset=H5Dcreate2( file, "Y", H5T_NATIVE_DOUBLE, fid, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );
//ED funktioniert nicht; dataset=H5Dcreate( file, "Y", H5T_NATIVE_DOUBLE, fid, H5P_DEFAULT );
	status=H5Dwrite(dataset,H5T_NATIVE_DOUBLE,H5S_ALL,H5S_ALL,H5P_DEFAULT,Y);
status=H5Dclose (dataset);
        status=H5Sclose (fid);
status=H5Fclose(file);

//write SNP names to dataset

	 }	 
	SNPL.exceptions ( ofstream::eofbit | ofstream::failbit | ofstream::badbit ); // checks if file can be written
	try
	{
		SNPL.open( ( parameter.out_file_name + "_SNPList.m" ).c_str(),  fstream::out );
		vector<string>::iterator it_SNPList;
                  
		SNPL << "SNPnames = { ";
for (it_SNPList = SNPList.begin(); it_SNPList < SNPList.end(); it_SNPList++)
		{
                        SNPL << "'" <<*it_SNPList<<"' ";
				
		}
                SNPL << "}"<<endl;

			SNPL << "X= [";
		for (it_SNPList = SNPList.begin(); it_SNPList < SNPList.end(); it_SNPList++)
		{
			size_t snp;
			for (
				snp = 0;
				snp < getSnpNo() && snpList.at( snp ).getSnpId() != *it_SNPList;
				++snp
			);
			
			if ( snp < getSnpNo() ) {
				const Vector xVec = const_cast<MData*>( this )->getXcolumn( snp );
				for ( size_t idv = 0; idv < getIdvNo(); ++idv ) {
					SNPL << ( 0 < idv ? " " : "" );
					SNPL << xVec.get( idv );
				}
				SNPL<< ";"<< endl;
			}
			else
			{
				cout << "SNP " << *it_SNPList << "not found";	
			}
		}
			SNPL << "]'"<<endl;
		
		for ( size_t cov = 0; cov < getCovNo(); ++cov ) {
			SNPL << getCovMatElementName( cov ) << " = [";
			const Vector covVec = const_cast<MData*>( this )->getCovariateColumn( cov );
			for ( size_t idv = 0; idv < getIdvNo(); ++idv ) {
				SNPL << ( 0 < idv ? ";" : "" );
				SNPL << covVec.get( idv );
			}
			SNPL<< "]'"<< endl;	
		}

		SNPL << "Y" <<" = [";
		for ( size_t idv = 0; idv < getIdvNo(); ++idv ) {
			SNPL << ( 0 < idv ? ";" : "" );
			SNPL << yVec.get( idv );
		}
		SNPL << "]" << endl;
		SNPL.close();
		printLOG( "Written MATLAB-File \"" + parameter.out_file_name + "_SNPList.m \"" );
	}
	catch (ofstream::failure e)
	{
		cerr << "Could not write MATLAB-File" <<endl;
	}	
}

void MData::writeBEDfile() {
	printLOG( "Writing genotype bitfile to \"" + parameter.in_files_plink + ".bed\"" );
	ofstream bed( ( parameter.in_files_plink + ".bed" ).c_str(), ios::out | ios::binary );
  
	// Magic numbers for .bed file
	const char magic[3] = {
		0x6c,
		char( parameter.imp_is_imputated ? 0x9a : 0x1b ),
		0x01	// SNP-major true
	};
	bed.write( magic, sizeof( magic ) );
  
	for ( size_t snp = 0; snp < getSnpNo(); ++snp ) {
		const Vector xVec = const_cast<MData*>( this )->getXcolumn( snp );
		unsigned int count = 0;
		char accumulator = 0;
		for ( size_t idv = 0; idv < getIdvNo(); ++idv ) {
			const double x = xVec.get( idv );
			unsigned int pattern;
			if ( -1.0 == x ) {
				pattern = 0x0;
			} else if ( 0.0 == x ) {
				pattern = 0x1;
			} else if ( 1.0 == x ) {
				pattern = 0x3;
			} else if ( ::isnan( x ) ) {
				pattern = 0x2;
			} else {
				throw Exception(
					"Cannot convert genotype matrix containing a value of %f at (%u,%u)"
					" to BED file."
					" Only values -1, 0, +1 and NaN are allowed.",
					x,
					idv,
					snp
				);
			}
			accumulator |= pattern << count;
			if ( 8 >= ( count+= 2 ) ) {
				bed.write( &accumulator, sizeof( accumulator ) );
			}
		}
	}

	bed.close();
}


void MData::writeBEDfilePlink()
{	
	bool orig = parameter.imp_is_imputated;
	parameter.imp_is_imputated = false;
	writeBEDfile();
	parameter.imp_is_imputated = orig;
}



/** computes correlation between two snps, utilise the structur of the data */
double MData::computeCorrelation ( const size_t locus1, const size_t locus2 ) const {
	const Vector
		v1 = const_cast<MData*>( this )->getXcolumn( locus1 ),
		v2 = const_cast<MData*>( this )->getXcolumn( locus2 );
	double
		sum1 = 0.0,
		sum2 = 0.0;
	const size_t idvs = getIdvNo();
	size_t n = 0;

	for ( size_t idv = 0; idv < idvs; ++idv ) {
		const double
			x1 = v1.get( idv ),
			x2 = v2.get( idv );
		if ( !::isnan( x1 ) && !::isnan( x2 ) ) {
			n++; // count for average;
			sum1 += x1;
			sum2 += x2;
		}
	}
	
	if  (n==0) // snps could not be compared because not matching values found (just missing)
	{				
		//~ printLOG
		//~ (
			//~ "Could not compute Correlation between SNPs \""+
			//~ (*snps_.at(getSnpVPFromGenoMatIt(locus1))).getSnpId()+
			//~ "\" and \""+
			//~ (*snps_.at(getSnpVPFromGenoMatIt(locus2))).getSnpId()+
			//~ "\""
		//~ );
		return 0.0; 
	}
	const double
		avg1 = sum1 / n,
		avg2 = sum2 / n;

	double sumCov = sum1 = sum2 = 0.0;
	
	for ( size_t idv = 0; idv < idvs; ++idv ) {
		const double
			x1 = v1.get( idv ),
			x2 = v2.get( idv );
		if ( !::isnan( x1 ) && !::isnan( x2 ) ) {
			const double
				diff1 = x1 - avg1,
				diff2 = x2 - avg2;
			sumCov += diff1 * diff2;
			sum1 += diff1 * diff1;
			sum2 += diff2 * diff2;
		}
	}
	
	if ( 0.0 == sum1 ) {
		//~ printLOG(
			//~ "Variance for SNP \""+
			//~ (*snps_.at(getSnpVPFromGenoMatIt(locus1))).getSnpId()+
			//~ "\" is 0!"
		//~ );
		return 0.0;
	}
		
	if ( 0.0 == sum2 ) {
		//~ printLOG(
			//~ "Variance for SNP \""+
			//~ (*snps_.at(getSnpVPFromGenoMatIt(locus2))).getSnpId()+
			//~ "\" is 0!"
		//~ );
		return 0.0;
	}

	return sumCov / sqrt( sum1 * sum2 );

}





/** according to Piotr's code

	Suppose x_ij, the genotype of SNP j for the i-th individual, is missing. We search for the
	"parameter.imp_Best_SNPs_No" SNPs with strongest correlation to SNP j fulfilling two conditions: They are in
	a neighborhood of "parameter.imp_Neighbours_No" SNPs upstream or downstream of SNP j and their values
	for the i-th individual are not missing. If we find individuals who have exactly
	the same values as the i-th individual on these "parameter.imp_Best_SNPs_No" SNPs, then we predict the value
	of x_ij as the most frequent value of SNP j among these individuals. If we cannot
	find individuals fulfilling the above mentioned condition, then the most frequent
	value of the jth SNP among all individuals is imputed. */
void MData::imputateMissingValues () {
	printLOG( "Begin Imputation of Missing Values." );
	int 	j, k, l; // loop variables
	int 	Predicted, PredictedAPriori;
	AutoVector pattern( 2 * parameter.imp_Neighbours_No );
	// auto_ptr used to deallocate upon exception
	auto_ptr<double> correlationCoeffs( new double[ 2 * parameter.imp_Neighbours_No + 1 ] );
	auto_ptr<size_t>
		neighboursOrder( new size_t[ 2 * parameter.imp_Neighbours_No + 1 ] ),
		bestSNPs( new size_t[ 2 * parameter.imp_Neighbours_No ] );
	size_t Frequency[3];
	bool 	Conformity, UsePriorPrediction;
	
	int 	MissingEntries=0; 	// Counter for missing entries
	int 	MissingHelp; 		// to check for a SNP there are missing entries
	int 	MissingSNPs =0; 	// Counter for SNPs with missing entries

	// Create Matrix
	AutoMatrix correlationMatrixT( getSnpNo(), parameter.imp_Neighbours_No );
	int progressCounter = 0;
	for ( size_t snp = 0; snp < getSnpNo(); ++snp ) {
		Vector xVec = getXcolumn( snp );

		MissingHelp=MissingEntries;
	
		// print progress (ED) this print only 1 time
		if ( 200 == progressCounter ) {
			printf( "\rDone %2.2f%%...", snp * 100.0 / getSnpNo() );
			fflush(stdout);
			progressCounter = 0;
		} else {
			++progressCounter;
		}
 
		// for the snp-th SNP determine predecessors and successors:

		//  predecessors do not exist
		for ( j = 0; j < parameter.imp_Neighbours_No - snp; ++j ) {
			correlationCoeffs.get()[j] = -10.0;
		}

		//  predecessors allready known
		for ( ; j < parameter.imp_Neighbours_No; ++j ) {
			correlationCoeffs.get()[j] = correlationMatrixT.get( snp, j );
		}
		// the snp-th
		correlationCoeffs.get()[ parameter.imp_Neighbours_No ] = -10.0;
               
		// compute successors
		// compute some safeguard
		if ( getSnpNo() <= parameter.imp_Neighbours_No+1 ) {
			printLOG( "Imputation should be used with this >>consider_n_neighbours<< setting, maybe there are only few snps" );
			exit(22);
		}
		for (
			j = parameter.imp_Neighbours_No + 1;
			j <= 2 * parameter.imp_Neighbours_No
			&&
			snp + j - parameter.imp_Neighbours_No < getSnpNo();   
			++j
		) {
			// WARNING could be less then 0
			correlationCoeffs.get()[j] = fabs( computeCorrelation( snp, j ) );
			correlationMatrixT.set(
				snp + j - parameter.imp_Neighbours_No,
				2 * parameter.imp_Neighbours_No - j,
				correlationCoeffs.get()[j]
			);
		}

		//  successors do not exist
		for ( j = 2 * parameter.imp_Neighbours_No; snp + j - parameter.imp_Neighbours_No >= getSnpNo(); --j ) {
			correlationCoeffs.get()[j] = -10.0;
		}

		// for sorting
		for ( j = 0; j <= 2 * parameter.imp_Neighbours_No; ++j ) {
			neighboursOrder.get()[j] = j;
		}

		// sort the neighbours w.r.t ascending correlation coefficients 
		SortVec pSort(
			2 * parameter.imp_Neighbours_No + 1,
			neighboursOrder.get(),
			correlationCoeffs.get()
		);

		// determine the abs. frequency of genotypes for the snp-th SNP 
		Frequency[0] = Frequency[1] = Frequency[2] = 0;
		for ( size_t idv = 0; idv < getIdvNo(); ++idv ) {
			const double x = xVec.get( idv );
			if ( -1.0 == x ) {
				++Frequency[0];
			} else if ( 0.0 == x ) {
				++Frequency[1];
			} else if ( 1.0 == x ) {
				++Frequency[2];
			} else if ( ::isnan( x ) ) {
				// treat later
			} else {
				throw Exception(
					"Cannot imputate genotype matrix containing a value of %f at (%u,%u)."
					" Only values -1, 0, +1 and NaN are allowed.",
					x,
					idv,
					snp
				);
			}
		}

		// use highest frequency as predictor
		if (
			Frequency[0] > Frequency[1]
			&&
			Frequency[0] > Frequency[2]
		) {
			PredictedAPriori = 0;
		} else if (
			Frequency[1] >= Frequency[0]
			&&
			Frequency[1] >= Frequency[2]
		) {
			PredictedAPriori = 1;
		} else {
			PredictedAPriori = 2;
		}
		
		// check for all individuals, whether the snp-th SNP is missing
		for ( size_t idv = 0; idv < getIdvNo(); ++idv ) {
			const double x = xVec.get( idv );
			if ( ::isnan( x ) ) {
				MissingEntries++;
				UsePriorPrediction = false;
				
				// search non-missing entries until parameter.imp_Best_SNPs_No SNPs are reached, sorted by correlation coefficent

				for (
					k = 0, l = 2 * parameter.imp_Neighbours_No;
					k < parameter.imp_Best_SNPs_No;
					++k, --l
				) {
					while( true ) {
						const double x1 = getXcolumn( snp + pSort.getId( l ) - parameter.imp_Neighbours_No ).get( j );
						if ( -1.0 != x1 && 0.0 != x1 && 1.0 != x1 ) {
							--l;
						} else {
							break;
						}
					}
					
					// no strongly correlated SNPs with missing SNP are found, so the most-frequent is used for prediction 	
					if ( pSort.getValue(l) < 0.0)
					{
						UsePriorPrediction = true;
						break;
					}
					
					bestSNPs.get()[k] = pSort.getId(l); // position of stronly correlated SNP
					// genotype of stronly correlated SNP
					pattern.set( k, getXcolumn( snp + bestSNPs.get()[k] - parameter.imp_Neighbours_No ).get( j ) );
				}
				
				Predicted = PredictedAPriori;
				// strongly correlated SNPs with missing SNP are found
				if (!UsePriorPrediction)
				{
					Frequency[0] = Frequency[1] = Frequency[2] = 0;
					// check for each individual if pattern is matched
					for (l = 0; l < getIdvNo(); l ++)
					{
						Conformity = true;
						for ( k = 0; Conformity && k < parameter.imp_Best_SNPs_No; ++k  )
							if (
								getXcolumn( snp + bestSNPs.get()[k] - parameter.imp_Neighbours_No ).get( l )
								!=
								pattern.get( k )
							) {
								Conformity = false;
							}
						if ( Conformity ) {
							const double x2 = xVec.get( l );
							if ( -1.0 == x2 ) {
								++Frequency[0];
							} else if ( 0.0 == x2 ) {
								++Frequency[1];
							} else if ( 1.0 == x2 ) {
								++Frequency[2];
							}
						}
					}
					if (
						Frequency[0] > Frequency[PredictedAPriori]
						&&
						Frequency[0] > Frequency[1]
						&&
						Frequency[0] > Frequency[2]
					) {
						Predicted = 0;
					} else if (
						Frequency[1] > Frequency[PredictedAPriori]
						&&
						Frequency[1] >= Frequency[0]
						&&
						Frequency[1] >= Frequency[2]
					) {
						Predicted = 1;
					} else if (
						Frequency[2] > Frequency[PredictedAPriori]
						&&
						Frequency[2] >= Frequency[0]
						&&
						Frequency[2] > Frequency[1]
					) {
						Predicted = 2;
					}
				}
				
				xVec.set( j, Predicted - 1.0 );
			}
		}
		if (MissingHelp != MissingEntries)
		{
			++MissingSNPs;
		}
	}
	
	printLOG("Imputation finished: "+int2str(MissingEntries)+ " missing values in "+int2str(MissingSNPs)+" SNPs replaced");
	parameter.imp_is_imputated = true;
}



void MData::calculateIndividualTests()
{
	
	printLOG("Start Individual Tests");
	
	Model SingleSNP( *this );	// create Model with current MData
	auto_ptr<size_t> snpArray( new size_t[ getSnpNo() ] );	// to store information for SortVec snp_order_
	auto_ptr<double> TestStat( new double[ getSnpNo() ] ); // to store information for SortVec snp_order_
// return allways  1 when asking for omp_get_num_threads()
	cerr<<"!parameter.affection_status_phenotype="<<!parameter.affection_status_phenotype<<endl;
	int ont=omp_get_num_threads();
	//cerr<<"omp_get_num_threads()"<<ont<<endl;
	//if(0==parameter.affection_status_phenotype)
	//{cerr<<"reel"<<endl;	omp_set_num_threads(1);}
	//else
	//{cerr<<"1"<<endl;   omp_set_num_threads(1);}
cerr<<"omp_get_num_threads()"<<omp_get_num_threads()<<endl;
//
	#pragma omp parallel for
	for ( size_t snp = 0; snp < getSnpNo(); ++snp ) {
		// for sorting, store positon 
		snpArray.get()[snp] = snp;
		// compute p-value of single marker test, and store for sorting
		TestStat.get()[snp] = SingleSNP.computeSingleRegressorTest( snp );
	}
	omp_set_num_threads( 4 );
	for ( size_t snp = 0; snp < getSnpNo(); ++snp ) {
		snpList.at( snp ).setSingleMarkerTest( TestStat.get()[snp] );
	}
	
	// sort the SNPs w.r.t ascending p-values in snp_order_
	snp_order_.fillVec( getSnpNo(), snpArray.get(), TestStat.get() );

	// output the SNP order an the p-values in a file
	ofstream IT;
	IT.open( ( parameter.out_file_name + "_IT.txt" ).c_str(), ios::out );
	
	IT << "SNP_no. \t SNP_name \t Chr \t Pos \t p-value" << endl;
	for ( size_t snp = 0; snp < getSnpNo(); ++snp ) {
		const size_t snpi = snp_order_.getId( snp );
		const SNP& snpo = snpList.at( snpi );
		IT 	<< snpi << "\t"
			<< snpo.getSnpId() <<"\t"
			<< snpo.getChromosome() << "\t"
			<< snpo.getBasePairPosition() << "\t"
			<< snp_order_.getValue( snp ) << endl;
	}
	
	IT.close();
	
	printLOG( "Individual Tests finished, written to \"" + parameter.out_file_name + "_IT.txt\"" );

}
/**calculate PValuePorder needs parameter.ms_MaximalPValueForwardStep but set in the conf-file
 *and of course the MData variable.
 determines the SNPs with p-Value < parameter.ms_MaximalPValueForwardStep
 these are tested in the Forward Step
 */
size_t MData::calculatePValueBorder () const {
	size_t PValueBorder;
	for (
		PValueBorder = getSnpNo() - 1;
		0 < PValueBorder
		&&
		snp_order_.getValue( PValueBorder ) > parameter.ms_MaximalPValueForwardStep;
		--PValueBorder
	);
//checking this will bring something when PValueBorder will be very small
	PValueBorder = min( 100u, (unsigned int) PValueBorder );
	if ( PValueBorder < 100 ) {
		PValueBorder = min( 100u, (unsigned int) getSnpNo()-1 );	//this is only for very small SNP files. REMARK<BB>: how about 0 SNPs?
	}
	return PValueBorder;
}
bool MData::selectModel( Model *currentModel, size_t PValueBorder, int ExpectedCausalSNPs /*no effect */, int maxModel, int criterium ) {
	int JJ=0;
	printLOG("Model Selection started: ");
	bool 	stop = false;
//	int     removedSNP=-1;
	int *startIndex; //start at the begin 
	int dummy=0;
	     startIndex=&dummy;
       // int PValueBorder =calculatePValueBorder();
	//or take the setting from the conf file
	PValueBorder = min( getSnpNo() - 1, PValueBorder );	// REMARK<BB>: how about 0 SNPs?
//	 int PValueBorder=parameter.PValueBorder;
PValueBorder=min(getSnpNo()-1,PValueBorder);
	cerr<<"PValueBorder"<< PValueBorder<<endl;
	// compute the log-likelihood of the 0 Model
	Model model0( *this );
	
	model0.computeRegression();
        model0.printYvec(true);
	setLL0M( model0.getMJC() );	// REMARK<BB>: functionality seems redundant with that in main.cpp
	Model model1( *this ), model2( *this );
        //Mod
	// Using pointers to avoid expensive copying of Model objects
	Model
        //	*currentModel,// = &model0, //,original 
		*forwardModel = &model1,
		*backwardModel = &model2;
       
//	schreibt es in die Variable für mBIC2 auch wenn man mBIC haben will
if (currentModel->getModelSize()>0)
{if (currentModel->computeMSC(criterium)>0.000001)
	//only here a backwardstep makes sense
	currentModel->saveguardbackwardstep( *backwardModel,criterium);
}

	currentModel->makeMultiForwardStep(PValueBorder,criterium,startIndex);
	currentModel->computeRegression();

	printLOG( "Start stepwise selection" );
double bestMSC= currentModel->getMSC();
//currentModel->replaceModelSNPbyNearCorrelated(0);
//exit(0);
bool improvment=false;
	while ( !stop ) {
         improvment=currentModel->replaceModelSNPbyNearFromCAT( *startIndex, PValueBorder, criterium );
	 improvment=currentModel->saveguardbackwardstep( *backwardModel, criterium );
         
	 //ED break because of computational limitations
	 if (currentModel->getModelSize()>min(parameter.maximalModelSize,maxModel))
		 break;
        /* linear case normal forward step
         */
          if (!parameter.affection_status_phenotype)//quantitative
		improvment = currentModel->makeForwardStepLinear( forwardModel, JJ, &bestMSC, PValueBorder, startIndex );
          else if (parameter.affection_status_phenotype)/*PRÄSELECTION nur bis Revision 274*/
        	improvment = currentModel->makeForwardStepLogistic( JJ, &bestMSC, PValueBorder, startIndex, criterium );	  
          else 
		  cerr<<" not linear nor logistic, this should not happen"<<endl;
	stop = currentModel->finalizeModelSelection( *backwardModel, JJ, improvment, PValueBorder, startIndex, criterium );
}//while
size_t reference=350;	// REMARK<BB>: Where does 350 come from? Also mind 0 SNPs case below.
	if( parameter.affection_status_phenotype)
      {if (350>getSnpNo()-1) reference=getSnpNo()-1;
	    
	while(currentModel->selectModel(*backwardModel,max(reference,PValueBorder),maxModel,criterium)) //minimum 100 or PValueBorder//with 1000 instead 100 it takes 2' on a model with 17 SNPS
		{
		 //currentModel->replaceModelSNPbyNearCorrelated1();//should
		 //improvment=currentModel->replaceModelSNPbyNearFromCAT(*startIndex , PValueBorder, criterium);
	         //if(improvment)	 currentModel->saveguardbackwardstep( *backwardModel,criterium);
		}
	} 
                  improvment=currentModel->replaceModelSNPbyNearFromCAT(*startIndex , PValueBorder, criterium);
	         if(improvment)	 currentModel->saveguardbackwardstep( *backwardModel,criterium);

//when all fails	

//currentModel->replaceModelSNPbyNearCorrelated(); bringt nichts
//currentModel=backwardModel;
currentModel->printStronglyCorrelatedSnps( 0.999, int2str(parameter.in_values_int) + "the_result" );

currentModel->printModel("finalModel",criterium);

return true ; 
}

///////////////////////////////////////////////////7
//selectModel ohne Argument
bool MData::selectModel()
{
       int JJ=0;
	printLOG("OLD Model Selection started: ");
	Model SelectedModel( *this ); 
	Model Forward( *this );
	Model Backward( *this );
	
	bool 	stop = false;
	int *startIndex; //start at the begin 
	int dummy=0;
	     startIndex=&dummy;
       // int PValueBorder =calculatePValueBorder();
	//or take the setting from the conf file
	size_t PValueBorder = parameter.PValueBorder;
	size_t eins=1;
	PValueBorder = max(min( getSnpNo()-1, PValueBorder ),eins);	// REMARK<BB>: how about 0 SNPs?
//if(getSnpNo()-1>350)
//PValueBorder=max(350,PValueBorder);
	cerr<<"PValueBorder"<< PValueBorder<<endl;
	// compute the log-likelihood of the 0 Model
	Model model0( *this );
	
	model0.computeRegression();
        model0.printYvec(true);
	setLL0M( model0.getMJC() );
	Model model1( *this ), model2( *this );
        //Mod
	// Using pointers to avoid expensive copying of Model objects
	Model
		*currentModel = &model0,
		*forwardModel = &model1,
		*backwardModel = &model2;
int backup_for_causal_snps=parameter.ms_ExpectedCausalSNPs;
if(parameter.expected_causal_snps1>parameter.ms_ExpectedCausalSNPs)
parameter.ms_ExpectedCausalSNPs=parameter.expected_causal_snps1;
	currentModel->makeMultiForwardStep(PValueBorder,1,startIndex);
if(parameter.expected_causal_snps1>parameter.ms_ExpectedCausalSNPs)
	parameter.ms_ExpectedCausalSNPs=backup_for_causal_snps;
//Achtung hier soll wieder alles beim alten sein
	currentModel->computeRegression();
       // best=currentModel;
	printLOG( "Start stepwise selection " );
        currentModel->computeMSC(0);
        //currentModel->printModel("check1");
double bestMSC= currentModel->getMSC();
bool improvment=false;
/*  This is the backward step
 *  one take the full model and remove every SNP.
 *  Then one look for that model with  lowest MSC 
 *  then we have a new start model.
 *
 *  We repeat this until the 1 model is selected.
 */
	while ( !stop ) {
         improvment=currentModel->replaceModelSNPbyNearFromCAT(*startIndex , PValueBorder);
	 improvment=currentModel->saveguardbackwardstep( *backwardModel);	 
	 currentModel=backwardModel; 
        /* linear case normal forward step
         */
         if (!parameter.affection_status_phenotype)//quantitative
	      improvment=currentModel->makeForwardStepLinear( forwardModel,  JJ,  &bestMSC,
		                                	 PValueBorder, startIndex);
          else if (parameter.affection_status_phenotype)/*PRÄSELECTION nur bis Revision 274*/
              improvment=currentModel->makeForwardStepLogistic(JJ, &bestMSC,  PValueBorder, startIndex);	  
          else cerr<<" not linear nor logistic, this should not happen"<<endl;
	  stop=currentModel->finalizeModelSelection( *backwardModel, JJ,  improvment,  PValueBorder,  startIndex);
	}//while
	size_t reference = 350;	// REMARK<BB>: Where does 350 come from?
	if( parameter.affection_status_phenotype)
	{
		if ( 350 > getSnpNo()-1 ) reference=getSnpNo()-1;	// REMARK<BB>: Beware 0 SNPs case.
	    
	while(currentModel->getModelSize()&&currentModel->selectModel(*backwardModel,max(reference,PValueBorder)))
		 //minimum 100 or PValueBorder//with 1000 instead 100 it takes 2' on a model with 17 SNPS
		{
		 improvment=currentModel->replaceModelSNPbyNearFromCAT(*startIndex , PValueBorder);
		 currentModel->saveguardbackwardstep( *backwardModel);
		}
	} 

currentModel->printStronglyCorrelatedSnps( 0.999, int2str(parameter.in_values_int) + "the_result" );

currentModel->printModel("finalModel");
cerr<<parameter.ms_ExpectedCausalSNPs<<endl;

return true;
}



/** just for TESTING !memory leak! */
void MData::readInSNPOrder ( const string& filename ) {

	ifstream	ITF;
	string 		helper;
	int 		Pos,SNPChr,SNPPos;
	string 		SNPId;
	double 		pValue;
	
	int 		i=0;
	auto_ptr<size_t> snpArray( new size_t[ getSnpNo() ] );
	auto_ptr<double> testStat( new double[ getSnpNo() ] );
	
	vector <bool> CheckSNP ( getSnpNo(), false);
	vector <bool>::iterator it_CheckSNP;
	ITF.open((filename).c_str(), ios::in);	
	if (ITF.is_open())
	{
		getline(ITF,helper ); //remove fist line
		while(! ITF.eof())
		{
			ITF >> Pos;    // get information per line
			ITF >> SNPId;
			ITF >> SNPChr; // not needed
			ITF >> SNPPos; // not needed
			ITF >> pValue;	
	
			if (getSNP(Pos)->getSnpId() != SNPId) // check if position match
			{
				printLOG("SNP Position of : \""+ filename +"\" and the Data do not match.");
				exit(1);
			}
			else
			{
				snpList.at(Pos).setSingleMarkerTest( pValue );
				CheckSNP.at(Pos)=true;
				snpArray.get()[i]=Pos;
				testStat.get()[i]=pValue;
			}
			i++;
		}
		
		for (it_CheckSNP = CheckSNP.begin(); it_CheckSNP < CheckSNP.end(); it_CheckSNP++) // check if Single Marker Test for every SNP was in file
		{
			if ( !*it_CheckSNP  )
			{
					printLOG("Not all SNPs in file \""+ filename +"\".");
					exit(2);	
			}
		}
		// set the SNP order in Model (a bit redundant, since sortet)
		snp_order_.fillVec( getSnpNo(), snpArray.get(), testStat.get() );
	}
	else
	{
		printLOG("Could-Not open file: \""+ filename +"\" as SNP-Orderfile.");
		exit(3);
	}
}

/** Artur new code:
   * @brief writes ordered snps to file. File name is set up in parameters class.
*/
void MData::printSnpOrder()
{
  // output the SNP order an the p-values in a file
  ofstream IT;
  IT.open( ( parameter.out_file_name + "_IT.txt" ).c_str(), ios::out );
  
  IT << "SNP_no. \t SNP_name \t Chr \t Pos \t p-value" << endl;
	for ( size_t snp = 0; snp < getSnpNo(); ++snp ) {
		const size_t snpi = snp_order_.getId( snp );
		const SNP& snpo = snpList.at( snpi );
		IT	<< snpi << "\t" 
			<< snpo.getSnpId() <<"\t"
			<< snpo.getChromosome() << "\t"
			<< snpo.getBasePairPosition() << "\t"
			<< snp_order_.getValue( snp ) << endl;
	}
  IT.close();
  printLOG( "Individual Tests finished, written to \"" + parameter.out_file_name + "_IT.txt\"" );
}
