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
#include <cmath>	// for nan(...)
#include <cfloat>	// for maximal double
#include "Model.hpp"
#include "GenotypeFreq.hpp"
#include "io/Hdf5Input.hpp"
#include <omp.h>
#include "PermSort.hpp"
#include <hdf5.h> //for hdf5
bool selec=false ;//false; //for debug  and linearisation
using namespace linalg;
using namespace io;

////++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//  class MData
////++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

const Vector MData::getXcolumn ( const size_t dim ) { return xMat.columnVector( dim ); }
Matrix MData::getCovariableMatrix () { return covMat; }
Matrix MData::getDummyCovariableMatrix () { return dummyCovMat; }
const Vector MData::getY () { return yVec; }

int MData::data2Geno ( const bool one, const bool two ) const {
	if ( one ) { 
		if (two) {//cerr<<"(o= 2";	// 11
		 	return 2;	// Homozygote "2"/"2"/
		} else {//cerr<<"(o= 9";	// 10
			return 9;
		}
	} else {
		if (two) {//cerr<<"(o= 1";	// 01
			return 1;	// Heterozygote
		} else {//cerr<<"(o= 0";	// 00
			return 0;	// Homozygote "1"/"1"
		}
	}
}

// getters for elements of the Genotype-Matrix with various parameters

/** get binary information for snp of the idv-individum */
int MData::getGenoMatElement ( const int snp, const int idv ) const {
	return data2Geno(genoMat_[snp]->one_[idv],genoMat_[snp]->two_[idv]);
}

/** @returns the binary information to {0,1,2} "risk" alleles or 9 for missing */
int MData::getGenoMatElement ( OneSnpAllInd* const & locus, const int idv ) const {
	// get binary information from the SNP pointed by locus of the idv-individum
	bool one = locus->one_.at(idv);
	bool two = locus->two_.at(idv);

	return data2Geno(one, two); // convert the binary information to {0,1,2}
}

int MData::getGenoMatElement(vector<SNP*>::iterator const& it_snp, vector<Individual*>::iterator const& it_ind) // slow!!! 
{	
	int snpdist = distance(snps_.begin(), it_snp); //  distance computes the number of vector-elements between two pointers 
	int inddist = distance(individuals_.begin(), it_ind); 
	
	return getGenoMatElement(snpdist, inddist);

}

/** @param entry = 0,1,2 */
void MData::setGenoMatElement ( const int snp, const int idv, const int entry ) {
	bool one=true; // default setting for to missing genotype
	bool two=false;
	
	// set boolens accoring to entry
	if (entry == 0)
	{
		one=false;
		two=false;
	}
	else if (entry == 1)
	{
		one=false;
		two=true;	
	}
	else if (entry == 2)
	{
		one=true;
		two=true;	
	}
	
	// set element in genotype-matrix
	genoMat_.at(snp)->one_.at(idv) = one;
	genoMat_.at(snp)->two_.at(idv) = two;
	
}

void MData::setLL0M ( const double ll ) {
	loglikelihood0Model_ = ll;
}

void MData::setY ( const size_t index, const double value) {
	yVec.set(index,value);
}

void MData::setY ( const size_t index, const int value) {
	yVec.set(index,value);
}	

/** Default Constructor: reads the input-files, sets parameters, deallambda.hpps with missing phenotypes */
MData::MData ( io::Input *input ) : xMat( 0, 0 ), covMat( 0, 0 ), dummyCovMat( 0, 0 ), yVec( 0 ) {

	if ( NULL == input && parameter.in_file_hdf5.empty() ) {
	////////////////////////////////////////////////////////////////////
	// Read SNP-Data (bim-File)
	ifstream 	BIM;
	
	BIM.open( parameter.in_files_plink_bim.c_str(), ios::in );
	BIM.clear();
	
	if ( BIM.is_open() ) {
		while( ! BIM.eof() ) {
			string chromosomeId, snpId;
			double geneticDistance;
			int basePairPosition;
			char allele1, allele2;
			BIM >> chromosomeId;
			BIM >> snpId;
			BIM >> geneticDistance;
			BIM >> basePairPosition;
			BIM >> allele1;
			BIM >> allele2;

			if ( BIM.eof() ) {
				break;
			}

			//SNP locus;
			SNP *locus = new SNP(
				chromosomeId,
				snpId,
				geneticDistance,
				basePairPosition,
				allele1,
				allele2
			);

			snps_.push_back(locus);
		}

		printLOG( "Read "+ int2str( getSnpNo() ) + " SNPs from file \"" + parameter.in_files_plink_bim+"\"" );
		BIM.close();
	} else {
		printLOG( "ERROR: Could not open .bim-file \"" + parameter.in_files_plink_bim + "\"" );
		exit(1);
	}
	
	////////////////////////////////////////////////////////////////////
	// Read Individual-Data (fam-File)	
	ifstream FAM;
	FAM.open( parameter.in_files_plink_fam.c_str(), ios::in );
	
	if ( FAM.is_open() ) {
		while( ! FAM.eof() ) {
			string familyId, individualId, paternalId, maternalId;
			unsigned int sexCode;
			double phenotype;
			FAM
				>> familyId
				>> individualId
				>> paternalId
				>> maternalId
				>> sexCode
				>> phenotype;

			if ( FAM.eof() ) {
				break;
			}

			Individual * person = new Individual(
				familyId, individualId, paternalId, maternalId,
				1 == sexCode ? Individual::MALE
					: 2 == sexCode ? Individual::FEMALE
						: Individual::MISSING,
				phenotype
			);

			individuals_.push_back(person);
			// REM: no checking of phenotype, no missings, no sexcode... see plinksource line 52806
		}
		printLOG( "Read " + int2str( getIdvNo() ) + " Individuals from file \"" + parameter.in_files_plink_fam + "\"" );
		FAM.close();
	} else {
		printLOG( "ERROR: Could not open .fam-file \"" + parameter.in_files_plink_fam + "\"" );
		exit(2);
	}
	
	////////////////////////////////////////////////////////////////////
	// Read binary genotype information file: (.bed -File)
	
	printLOG( "Reading genotype bitfile \"" + parameter.in_files_plink_bed + "\"" );
 
	ifstream BED;
	// Check File Format
	openBinaryFile( parameter.in_files_plink_bed, BED );
	
	// Allocate Memory
	genoMat_.resize( getSnpNo() );
	for ( int i = 0; i < getSnpNo(); ++i ) {
		OneSnpAllInd * newlocus = new OneSnpAllInd;
		newlocus->one_.resize( getIdvNo() ); // vector resize to allocate memmory
		newlocus->two_.resize( getIdvNo() );
		genoMat_.at(i)=newlocus;
    	} 
	// TODO<BB>: Improve error tolerance, e.g. if bed-File contains too many or too few entries.
    
	// Read data  
	OneSnpAllInd * snp;
	  
	// Outer loop for SNPs
	for ( int s = 0; s < getSnpNo(); ++s ) {
		snp = genoMat_.at(s);
	      
		// Inner loop for individuals
	      
		for ( int indx = 0; indx < getIdvNo(); ) {

			char ch;
			BED.read( &ch, 1 );	// TODO<BB>: Use easier system function?
			if ( ! BED ) {
				printLOG("ERROR: Problem with the BED file...has the FAM/BIM file been changed?");
				exit(4);		
			}

			bitset<8> b = ch;
		  	for (
				int c = 0;
				c < 7 && indx < getIdvNo();
				++indx
			) {
				snp->one_.at( indx ) = b[ c++ ];
				snp->two_.at( indx ) = b[ c++ ];
			}
		}	  
	}
	
    	char ch;
    	BED.read( &ch, 1 );
    	if ( BED )  {
		printLOG( "ERROR: Problem with the BED file... has the FAM/BIM file been changed?" );
	}

	BED.clear();
	BED.close();

	printLOG("Successfully read-in binary file");	
	
	////////////////////////////////////////////////////////////////////
	// Read Y-values file: (.yvm - File)

	if ( parameter.y_value_extra_file ) {
		readInYValues();
	}
	checkYValues(); //.yvm - File

	// Migrate data to modern Matrix/Vector API data
	// TODO<BB>: Known BUG:
	// Covariable and dummy covariable matrix should not be read
	// before checkYValues has determined the final number of individuals.
	// So far, their sizes do not match the requirements, which should yield a GSL error!
	if (parameter.cov_extra_file)
               readCovariablesFile();
	const size_t
		idvs = getIdvNo(),
		snps = getSnpNo();
	xMat.exactSize( idvs, snps );
	yVec.exactSize( idvs );
	for ( size_t idv = 0; idv < idvs; ++idv ) {
		for ( size_t snp = 0; snp < snps; ++snp ) {
			int geno = getGenoMatElement( snp, idv );
			double x;
			switch( geno ) {
				case 0: x = -1.0; break;
				case 1: x = 0.0; break;
				case 2: x = +1.0; break;
				default: x = nan( "missing" );
			}
			xMat.set( idv, snp, x );
		}
		yVec.set( idv, getYvalue( idv ) );
	}
	} else {
		const bool hdf5Input = NULL == input;
		if ( hdf5Input ) {
			input = new Hdf5Input( parameter.in_file_hdf5.c_str() );
		}
		const size_t
			snps = input->countSnps(),
			idvs = input->countIndividuals();
		xMat.exactSize( idvs, snps );
		yVec.exactSize( idvs );
		for ( size_t i = 0; i < snps; ++i ) {
			SNP *snp = new SNP( input->getSnp( i ) );
			snps_.push_back( snp );
		}
		for ( size_t j = 0; j < idvs; ++j ) {
			Individual *idv = new Individual( input->getIndividual( j ) );
			individuals_.push_back( idv );
		}
		for ( size_t i = 0; i < snps; ++i ) {
			const Vector genotypes = input->getGenotypeVector( i );
			Vector xVec = xMat.columnVector( i );
			xVec.copy( genotypes );
		}
		const Vector phenotypes = input->getPhenotypeVector();
		yVec.copy( phenotypes );
		if ( hdf5Input ) {
			delete input;
			input = NULL;
		}
	}
}



/** Destructor */
MData::~MData () {

	// clear SNPs
	for(
		vector<SNP*>::iterator it_snp = snps_.begin();
		it_snp < snps_.end();
		++it_snp
	) {
		delete ( *it_snp );	// essential memory release
	}
	snps_.clear();

	// clear Individuals
	for(
		vector<Individual*>::iterator it_ind = individuals_.begin();
		it_ind < individuals_.end();
		++it_ind
	) {
		delete ( *it_ind );	// essential memory release
	}
	individuals_.clear();

	// clear Binary Data	
	for (
		vector<OneSnpAllInd*>::iterator it_osai= genoMat_.begin();
		it_osai < genoMat_.end();
		++it_osai
	) {
		(*it_osai)->one_.clear();	// TODO<BB>: not essential?
		(*it_osai)->two_.clear(); 
		delete ( *it_osai );	// essential memory release
	}
	genoMat_.clear();
}

size_t MData::getSnpNo () const {
	return snps_.size();
}

size_t MData::getIdvNo () const {
	return individuals_.size();
}

string MData::getFID ( const size_t index ) const {
	return individuals_.at( index )->getFamilyID();
}

string MData::getID ( const size_t index ) const {
	return individuals_.at( index )->getIndividualID();
}

double MData::getYvalue ( const size_t idv ) const {
	return yVec.get( idv );
}

SNP* MData::getSNP ( const size_t snp ) const {
	return snps_.at( snp );
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
	(snps_.at(index))->setSingleMarkerTest(value);
}

void MData::fillSnp_order_Vec ( const size_t snpNo, int*  SNPList, double* TestStat ) {
	snp_order_.fillVec( snpNo, SNPList, TestStat );
}

/** from Plink */
void MData::openBinaryFile ( const string filename, ifstream & BIT ) {
	// no checking if file is open!!!!
	BIT.open( filename.c_str(), ios::in | ios::binary );

  	// 1) Check for magic number
  	// 2) else check for 0.99 SNP/Ind coding
  	// 3) else print warning that file is too old
  
  	char ch;
  	BIT.read( &ch, 1 );
 
  	bool v1_bfile = false;

	// If v1.00 file format
	// Magic numbers for .bed file: 00110110 11011000 = v1.00 bed file
	// Achtu/ng immer von hinten lesen!!! 
	// in Datei steht dann immer 01101100 00011011 = 0x6c1b
  	if ( (char) 0x6c == ch ) {

 		// Next number 00011011
 		BIT.read( &ch, 1 );
		if (
			( (char) 0x1b == ch )
	 		or	// Next number 10011010 = 0x9a for imputated
			( (char) 0x9a == ch )
	   	) {
			if ( (char) 0x9a == ch ) {
				parameter.imp_is_imputated = true;
				printLOG( "Detected that binary PED file is already imputated" );
			}

			// Read SNP/Ind major coding
 			BIT.read( &ch, 1 );
			// if bfile_SNP_major
			if ( 0x1 & ch ) {
				printLOG( "Detected that binary PED file is v1.00 SNP-major mode" );
			} else {
				printLOG( "ERROR: Detected that binary PED file is v1.00 individual-major mode\n Only SNP-major mode is supported" );
				exit(3);
			}
			v1_bfile = true;
		}
	}

	// Reset file if < v1
	if ( ! v1_bfile ) {
		printLOG( "ERROR, old BED file <v1.00 : ...\n" );
		printLOG( "(try in plink --make-bed from PED )\n" );
		exit(4);
	}
}



/** remove the genotype of the idv-th individual from all SNPs */
void MData::removeIndividual ( const int idv ) {
	vector<OneSnpAllInd*>::iterator it_osai;	
	for(it_osai = genoMat_.begin(); it_osai < genoMat_.end(); it_osai++ )
	{
		(*it_osai)->one_.erase((*it_osai)->one_.begin() + idv);
		(*it_osai)->two_.erase((*it_osai)->two_.begin() + idv);
	}

	// remove the idv-th individual from all individuals_
	delete ( individuals_.at(idv) );
	individuals_.erase(individuals_.begin() + idv );
	// TODO: reflect removal in the linalg::AutoMatrix. Won't be fixed, because reading goes to the io package.
}



/** Check Y-Values and store them in a vector.
 MISSING: to check for missing values and remove then the individual, and the genotype-info */
void MData::checkYValues()
{
	
	int noOfRemovedIndividuals=0;
	double indPheno;
	parameter.affection_status_phenotype = true;	// intitial setting, Y-values are tested if the are just (0,1) (or missing)
	
	caseNo_=0;
	contNo_=0;
	
	for ( int i = 0; i < getIdvNo(); ++i ) {
		indPheno = individuals_.at(i)->getPhenotype();
	//	if(i<200)cerr<<"i="<<i<<" indPheno="<<indPheno<<endl;
		// do not consider individuals with missing phenotype
		if ( indPheno == parameter.missing_phenotype_code ) {
			removeIndividual(i); // remove all information concerning the i-th individual
			i--; // the current position was removed, so a diffrent individual is now at the i-th position
			noOfRemovedIndividuals++;
		}
		else // no missing phenotype, determine if phenotype is affection (case-control) or quantitative
		{
			const size_t dim = yVec.countDimensions();
			yVec.upSize( 1 + dim );
			yVec.set( dim, indPheno );
	                
			if ( indPheno == parameter.case_value) // was set 1 check affection status
			{
				caseNo_++;
			}
			else if ( indPheno == parameter.control_value)  //was set 0
			{
				contNo_++;
			}
			else // some other number than 0/1 detected, phenotype must be qua
			{ //cerr<<"i="<<i<<" indPheno="<<indPheno<<endl;
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



/** Read in Y-Values from a yvm-file 
 *  1. Row  the Headerline  first 2 column  
 *  FamilyId Individual ID
 *  Warning interactive Part when no column selected in conf file!*/
void MData::readInYValues()
{
	
	ifstream		YVM;
	string			headline;
	string			buffer;
	vector <string>	headers;
cerr<<"YVM file:"<<	parameter.in_files_values_yvm<<endl;
	YVM.open( parameter.in_files_values_yvm.c_str(), ios::in );
	
	if (YVM.is_open())
	{
		getline(YVM,headline,'\n'); // get first line of file
		stringstream ss(headline);
		 
		while( ss >> buffer) // read in headers and store them in vector
		{
			headers.push_back(buffer);
		}
		
		int headersize = headers.size(); // to compare int (otherwise int and unsigned int)
		
		if ( headersize <= 2) // no information about y-values
		{
			printLOG( "File \"" + parameter.in_files_values_yvm + "\" does not contain enough information to obtain Y-Value." );	
			exit(4);
		}
		
		if( parameter.in_values_name.empty() ) {
			while( parameter.in_values_int > headersize-2 or parameter.in_values_int < 1 ) {
				printLOG("\"" + int2str( parameter.in_values_int )+ "\" gives no valid Position for Y-Value.");
				cout << "Please chose a number between 1 and " << headersize-2 << endl;
				cin >> parameter.in_values_int;
			}
		}
		else
		{
			bool found = false; // bool if name is found
			vector <string>::iterator it_headers;
			for( it_headers= headers.begin()+2 ; it_headers < headers.end(); it_headers++)
			{
				if ( *it_headers == parameter.in_values_name ) // checks if name is found
				{
					// set found position, 1 indexed
					parameter.in_values_int = ( it_headers - headers.begin() - 1 );
					found = true;
				}
			}
			if ( !found )  // name is not found, user can choose out of possible traits
			{
				printLOG( "\"" + parameter.in_values_name + "\" gives no valid name for Y-Value." );
				cout << "Please chose a number according to the following traits: \n";
				for(int i = 1; i <=  headersize - 2 ; i++)
				{
					cout << i << " " << headers.at(i+1)<< endl; // state traits
				}
				cin >> parameter.in_values_int;
				// check whether input is valid now
				while ( parameter.in_values_int > headersize-2 or parameter.in_values_int < 1 ) {
					printLOG( "\"" + int2str( parameter.in_values_int )+ "\" gives no valid Position for Y-Value." );
					cout << "Please chose a number between 1 and " << headersize - 2 << endl;
					cin >> parameter.in_values_int;
				}	
			}
			
				
		
		}
		//searchPos =  parameter.in_values_int;
		
		parameter.y_value_name = headers.at( parameter.in_values_int + 1 );	// 1  - indexed and starting with 2 -> searchPos + 1
		printLOG( "Using trait \"" + parameter.y_value_name + "\" as Target" );
	        //parameter.out_file_name=parameter.out_file_name+parameter.y_value_name; this is to late
		//because the *logg file needs this name to!
		// Start Read-in:
		map<string,Individual*> uid; // map to search for the individuals
		map<string,Individual*>::iterator it_uid; // map to search for the individuals
		for ( int i = 0; i < getIdvNo(); ++i ) {
			uid.insert(make_pair (individuals_.at(i)->getFamilyID()+"_"+individuals_.at(i)->getIndividualID(),individuals_.at(i))  );  
			// set all invididuals to missing phenotype in order to check if the phenotype is really set
			individuals_.at(i)->setPhenotype( parameter.missing_phenotype_code );
		}
		
		
		while ( !YVM.eof() )
		{
			string fid;
			string iid;
			double d_yvalue;
			string st_yvalue;
			int i = 1;
			getline(YVM,headline,'\n');
			stringstream ss(headline);
			
			ss >> fid;
			ss >> iid;
			
			
			// search for the i-th entry, treat entries as strings
			while ( i <= parameter.in_values_int ) {
				ss >> st_yvalue; 
				i++;
			}
			
			// converte string entry to a double (if not possible, y=0)
			stringstream ys(st_yvalue);
			ys >> d_yvalue;
		//	cerr<<d_yvalue<<";";//
//person will become a Individual 			
		    it_uid = uid.find(fid+"_"+iid);
			if (it_uid != uid.end() ) // person was found
			{
				Individual * person = it_uid -> second;
				(person)->setPhenotype(d_yvalue);
			}
	
		}
		
	}
	else
	{
		printLOG( "Error: Could not open file \"" + parameter.in_files_values_yvm + "\"" );
		exit(-5);
	}
	
}
void MData::setYfromMOSGWA(const vector<bool> Y )
{
	for(int i=0;i < getIdvNo(); ++i )
	{
     individuals_.at(i)->setPhenotype( Y[i] );
    }
}






/** TESTING */
void MData::printSNPs () {
	for (
		vector<SNP*>::iterator itsnps = snps_.begin();
		itsnps < snps_.end();
		++itsnps
	) {
		cout << (*itsnps)->getChromosome() << " "
			 << (*itsnps)->getSnpId()	<< " " 
			 << (*itsnps)->getGeneticDistance()	<< " "  			 
			 << (*itsnps)->getBasePairPosition()	<< " " 
			 << (*itsnps)->getAllele1()	<< " " 
			 << (*itsnps)->getAllele2()	<< " "  
			 << (*itsnps)->getSingleMarkerTest()	<< " " 
			 << endl;
	}
}



/** TESTING */
void MData::printIndividuals () {

	vector<Individual*>::iterator itidvs;
	for (itidvs = individuals_.begin(); itidvs < individuals_.end();  itidvs++)
	{
		cout << (*itidvs) << endl;
	}
	
}

void MData::printIndividuals (ofstream &Y) {

	vector<Individual*>::iterator itidvs;
	for (itidvs = individuals_.begin(); itidvs < individuals_.end();  itidvs++)
	{
		Y << (*itidvs)->getFamilyID() << " "
	
			 << (*itidvs)->getIndividualID() << " " 
  
		//	 << (*itidvs)->paternalID_	<< " " 
		 
		//	 << (*itidvs)->maternalID_	<< " " 

		//	 << (*itidvs)->sexCode_	<< " " 

		// here pheno 1-num pheno;
		//here should be the new phenotype 
			 << (*itidvs)->getPhenotype() << " " 
			 << endl;
            // cout<<(*itidvs)->phenotype_;
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
			Y  << yVec.get( dim ) << endl;
	try {	Y.close();
		}
	catch ( ofstream::failure e/*xception*/ ) {
		cerr << "Could not closing Y for GWASselect file:" << ( parameter.out_file_name + "Yextra" ).c_str()<< endl;
	}
}


/** TESTING */
/** caco is case control 0 control 1 case*/
void MData::printGenoData (ofstream &Hyper, vector <bool> sel,bool caco) const {
	if(sel.size()!=getIdvNo())
	{cerr<<"printGenoData should have sel.size() == getIdvNo() !"<<endl;
	}
	else	
	{
	for ( int i= 0; i < getSnpNo(); ++i ) {
		Hyper <<(snps_[i])->getSnpId()<<' ';
		const Vector xVec = const_cast<MData*>( this )->getXcolumn( i );
		for ( int j = 0; j < getIdvNo(); ++j ) {
			if(sel[j]==caco)
			Hyper << xVec.get( j ) << " ";
		}
		Hyper << endl;
	}    
	}//from else
}
void MData::printGenoData () const {
	for ( int i= 0; i < getSnpNo(); ++i ) {
		const Vector xVec = const_cast<MData*>( this )->getXcolumn( i );
		for ( int j = 0; j < getIdvNo(); ++j )  {
			cout << xVec.get( j ) << " ";
		}
		cout << endl;
	}
}

/** printGenoData prints the whole GenoData Matrix in 
a ofstream  Hyper 
Example: see  printHyper
*/
void MData::printGenoData (ofstream &Hyper) const {
	for ( int i= 0; i < getIdvNo() ; ++i ) //transponiert
	{
		for ( int j = 0; j <getSnpNo() ; ++j ) {
			const Vector genome = const_cast<MData*>( this )->getXcolumn( j );
			Hyper << genome.get( i ) << " ";
		}
		Hyper << endl;
	}
}




/** wer braucht die SNP id noch*/
void MData::printSNPId(ofstream &Hyper) const {
for(int i=0; i< getSnpNo();++i)
    Hyper<<(snps_[i])->getSnpId()<<' ';
Hyper<<endl;	
}
void  MData::findSNPIndex(vector<string>& SNPNames, vector<unsigned int>& index) const
{ //first sort the vector
 sort(SNPNames.begin(),SNPNames.end());	 //ERICH at first sort the input
 //create a vector of strings from

vector< int> tindex(SNPNames.size(),-999); //bring index to final size -999 is the not available SNP
 vector<string> allSNPs(getSnpNo()); //this is the extension of the vector the SNP names 
                                    // have probably different size.

 for(int i=0; i< getSnpNo();++i)//this is maybe a large loop 500 000 is possible
 {
         allSNPs[i]=(snps_[i])->getSnpId();
 }
 std::vector<int> permutation;
 sortingPermutation(allSNPs, permutation);
 unsigned int j=0;
for(int i=0;i<SNPNames.size();++i)
    for(j=0;j<allSNPs.size();++j)
        if(0==SNPNames[i].compare(allSNPs[permutation[j]]))
	{ tindex[i]=(permutation[j]); //permutation beginnt mit 1
	if(false)	cout<<tindex[i]<<endl; //DEBUG
	   //the J in for loop head should not be called!	
	  break;
	}

for(int i=0;i<tindex.size();++i)
	if(tindex[i]==-999)
	{cout<<"SNP name "<<SNPNames[i]<<" At position " << i<<" does not exist nicht in this dataset \n";}
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
}


void MData::printGWASselect(Model & newmodel) const {
//#pragma omp parallel shared(GWASselectcases,GWASselectcontrol,sel) private(i) //has to be alterd to fit
vector<bool> sel;
	newmodel.getYvec(sel);
if (false) //Debug code 
{cerr<<"printGWASselect sel size "<<sel.size()<<endl;
 for(int  i=0;i<sel.size();++i)
	cerr<<sel[i];
}        

ofstream GWASselectcontrol;
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

//#pragma omp section
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



void MData::printSelectedSNPsInR ( vector<string> SNPList ) const {
	ofstream	SNPL;
	int			i;
	
	SNPL.exceptions ( ofstream::eofbit | ofstream::failbit | ofstream::badbit ); // checks if Logfile can be writt
	try
	{
		SNPL.open( ( parameter.out_file_name + "_SNPList.txt" ).c_str(), fstream::out );
		vector<string>::iterator it_SNPList;
		
		for (it_SNPList = SNPList.begin(); it_SNPList < SNPList.end(); it_SNPList++)
		{
			//cout << "Search " << *it_SNPList<< endl;
			i=0;
			while( i<getSnpNo() && snps_.at(i)->getSnpId() !=  *it_SNPList)
			{
				i++;
			}
			
			if (i < getSnpNo() )
			{
				SNPL << *it_SNPList <<" <- c(";
				const Vector xVec = const_cast<MData*>( this )->getXcolumn( i );
				for ( int j = 0; j < getIdvNo(); ++j ) {
					SNPL << ( 0 == j ? "" : "," );
					SNPL << xVec.get( j );
				}
				SNPL<< ")"<< endl;
			}
			else
			{
				cout << "SNP " << *it_SNPList << "not found";	
			}
		}
		
		for ( int i = 0; i < parameter.covariables; ++i ) {
			SNPL << getCovMatElementName(i)<<" <- c(";
			SNPL<<  getCovMatElement(i, 0);
			for ( int j = 1; j < getIdvNo(); ++j ) {	// TODO<BB>: Why not start with j=0?
				SNPL <<"," <<   getCovMatElement(i, j);
			}
			SNPL<< ")"<< endl;	
		}
		
		for ( int i = 0; i < parameter.dummy_covariables; ++i ) {
			SNPL << getDummyCovMatElementName(i)<<" <- c(";
			SNPL<<  getDummyCovMatElement(i, 0);
			for ( int j = 1; j < getIdvNo(); ++j ) {	// TODO<BB>: Why not start with j=0?
				SNPL <<"," <<   getDummyCovMatElement(i, j);
			}
			SNPL<< ")"<< endl;	
		}
	

			SNPL << "Y" <<" <- c(";
			SNPL<< individuals_.at(0)->getPhenotype();
			for ( int j = 1; j < getIdvNo(); ++j ) {	// TODO<BB>: Why not start with j=0?
				SNPL <<","<< individuals_.at(j)->getPhenotype();
			}
			SNPL<< ")"<< endl;
		SNPL <<endl;
		
			SNPL << "Intercept" <<" <- c(1";
			for ( int j = 1; j < getIdvNo(); ++j ) {	// TODO<BB>: Why not start with j=0?
				SNPL <<","<<1;
			}
			SNPL<< ")"<< endl;
		
		//~ SNPL << "summary(lm(Y~Dummy1 + Dummy2 + Dummy3 +";
		SNPL <<"X <- cbind( Intercept";
		for (it_SNPList = SNPList.begin(); it_SNPList < SNPList.end(); it_SNPList++)
		{
			//~ if ( it_SNPList != SNPList.begin() )
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

/** printSelectedSNPsInMatlab create an *Octave.h5 file with the selected SNPs (as a vector of strings)  and their X and ofcourse the phenotype Y
 * additionaly it writes a */ 
void MData::printSelectedSNPsInMatlab ( vector<string> SNPList , string extra) const  {
	ofstream	SNPL;
         { //create the same in hdf
		 hid_t file,fid,dataset,space,/*dset,memtype,*/ filetype,props;
                 herr_t status;
                 hsize_t dim[]={SNPList.size(),getIdvNo()};//transponiert
                 hsize_t di[]={SNPList.size()};
		 double  zwischen[getIdvNo()*SNPList.size()];
		 double  Y[getIdvNo()];
		 vector<string>::iterator it_SNPList;
		 int i;
		  //zwischen=malloc(getSnpNo()*SNPList.size()*sizof(double));
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
unsigned int ii=0;
	for (jj=0; jj < index.size(); ++jj)
	{
		const Vector tmp = const_cast<MData*>( this )->getXcolumn( index[jj] );
	//	cerr<<index[jj]<<endl;
/*DEBUG*/	if (false)	cerr<<"neues SNP"<<jj<<"index"<<index[jj]<<endl;
		for(ii=0;ii<getIdvNo();++ii)
        	{    //cerr<<tmp.get(ii)<<endl;
		     //cerr<<	jj*getIdvNo()+ii<<" ";
		       	zwischen[jj*getIdvNo()+ii] = tmp.get(ii);
	        }
	}	
        for ( int j = 0; j < getIdvNo(); ++j )
       		
				Y[j]=individuals_.at(j)->getPhenotype();

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

//wrtite SNP names to dataset

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
			int i = 0;

			//cout << "Search " << *it_SNPList<< endl;
			while( i<getSnpNo() && snps_.at(i)->getSnpId() !=  *it_SNPList)
			{
				i++;
			}
			
			if (i < getSnpNo() )
			{
				//SNPL << *it_SNPList <<" = [";
				const Vector xVec = const_cast<MData*>( this )->getXcolumn( i );
				for ( int j = 0; j < getIdvNo(); ++j ) {
					SNPL << ( 0 == j ? "" : " " );
					SNPL << xVec.get( j );
				}
				SNPL<< ";"<< endl;
			}
			else
			{
				cout << "SNP " << *it_SNPList << "not found";	
			}
		}
			SNPL << "]'"<<endl;
		
		for ( int i = 0; i < parameter.covariables; ++i ) {
			SNPL << getCovMatElementName(i)<<" = [";
			SNPL<<  getCovMatElement(i, 0);
			for ( int j = 1; j < getIdvNo(); ++j ) {	// TODO<BB>: Why not start with j=0?
				SNPL <<";" <<   getCovMatElement(i, j);
			}
			SNPL<< "]'"<< endl;	
		
		}
		
		for ( int i = 0; i < parameter.dummy_covariables; ++i ) {
			SNPL << getDummyCovMatElementName(i)<<" = [";
			SNPL<<  getDummyCovMatElement(i, 0);
			for ( int j = 1; j < getIdvNo(); ++j ) {	// TODO<BB>: Why not start with j=0?
				SNPL <<";" <<   getDummyCovMatElement(i, j);
			}
			SNPL<< "]"<< endl;	
		}
	

			SNPL << "Y" <<" = [";
			SNPL<< individuals_.at(0)->getPhenotype();
			for ( int j = 1; j < getIdvNo(); ++j ) {	// TODO<BB>: Why not start with j=0?
				SNPL <<";"<< individuals_.at(j)->getPhenotype();
			}
			SNPL<< "]"<< endl; //here the ; remain 
		SNPL <<endl;
		SNPL.close();
		printLOG( "Written MATLAB-File \"" + parameter.out_file_name + "_SNPList.m \"" );
	}
	catch (ofstream::failure e)
	{
		cerr << "Could not write MATLAB-File" <<endl;
	}	
}

void MData::writeBEDfile()
{

	
  //from plink output.cpp in_files_plink_bed
  ofstream BED;
	printLOG( "Over_Writing genotype bitfile to \"" + parameter.in_files_plink_bed+"\"" );

	BED.open( ( parameter.in_files_plink_bed).c_str(), ios::out | ios::binary );
  

  bitset<8> b;
  char ch[1];


  // Magic numbers for .bed file: 00110110 11011000 = v1.00 bed file

  b.reset();
  b.set(2);  b.set(3);  b.set(5);  b.set(6);
  ch[0] = (char)b.to_ulong();
  BED.write(ch,1);
	
  b.reset();
	// set to 01011001
	if ( parameter.imp_is_imputated ) {
	  b.set(1);  b.set(3);  b.set(4); b.set(7);
	}
  else
  {
	b.set(0);  b.set(1);  b.set(3);  b.set(4);
  }
  ch[0] = (char)b.to_ulong();
  BED.write(ch,1);


  // BIT represents status of SNP-major (true) 

  b.reset();
  b.set(0);
  ch[0] = (char)b.to_ulong();
  BED.write(ch,1);
  
   
  
  vector<OneSnpAllInd*>::iterator it_osai = genoMat_.begin();

    // Outer loop over SNPs
  while ( it_osai != genoMat_.end() )
  {
	  vector<bool>::iterator it_1 = (*it_osai)->one_.begin();
	  vector<bool>::iterator it_2 = (*it_osai)->two_.begin();


	  // Inner loop over individuals
	  while ( it_1 != (*it_osai)->one_.end() )
	    {
	      bitset<8> b;
	      b.reset();
	      int c=0;      
	      
	      while (c<8 && it_1 != (*it_osai)->one_.end() )
		  {
			if ( *(it_1) )
			{
				b.set(c);
			}
			it_1++;
			c++;

			if ( *(it_2) )
			{
				b.set(c);
			}		  
			it_2++;
			c++;
		}
	      
	      char ch[1];
	      ch[0] = (char)b.to_ulong();
	      BED.write(ch,1);
	      
	    }

	  // next SNP
	  it_osai++;
  }
   
  BED.close();
  
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
void MData::imputateMissingValues()
{

	

	printLOG("Begin Imputation of Missing Values.");
	int 	i, j, k, l; // loop variables
	int 	x;
	int 	Predicted, PredictedAPriori;
	double*	CorrelationCoeffs = new double [ 2 * parameter.imp_Neighbours_No + 1 ];
	int*	NeighboursOrder = new int[ 2 * parameter.imp_Neighbours_No + 1 ];
	int*	BestSNPs = new int[ 2 * parameter.imp_Neighbours_No ];
	int*	Pattern = new int[ 2 * parameter.imp_Neighbours_No ];
	int		Frequency[3];
	bool 	Conformity, UsePriorPrediction;
	
	int 	MissingEntries=0; 	// Counter for missing entries
	int 	MissingHelp; 		// to check for a SNP there are missing entries
	int 	MissingSNPs =0; 	// Counter for SNPs with missing entries

	// Create Matrix
	vector<vector <float> > CorrelationMatrix ( getSnpNo(), vector <float> ( parameter.imp_Neighbours_No ) );
int num= 0;
	for (i = 0; i < getSnpNo(); i ++)
	{
		
		MissingHelp=MissingEntries;
		
	
		// print progress (ED) this print only 1 time
		if (num==200)
		{	
		printf ("\rDone %2.2f%%...", i*100.0 / (getSnpNo()));
		fflush (stdout);}
		if(200==num)
			num=0;
		else
			num++;
		
         	
		// for the i-th SNP determine predecessors and successors:

		//  predecessors do not exist
		for ( j = 0; j < parameter.imp_Neighbours_No - i; ++j ) {
			CorrelationCoeffs[j] = -10.0;
		}
		
		//  predecessors allready known
		for ( ; j < parameter.imp_Neighbours_No; ++j ) {
			CorrelationCoeffs[j] = double (CorrelationMatrix[i][j]);
		}
		// the i-th
		CorrelationCoeffs[ parameter.imp_Neighbours_No ] = -10.0;
               
		// compute successors
		//compute some saveguard 
		if (getSnpNo()<=parameter.imp_Neighbours_No+1)
		{printLOG("Imputation should be used with this >>consider_n_neighbours<< setting, maybe there are only few snps  ");exit(22);}
		for (
			j = parameter.imp_Neighbours_No + 1;
			j <= 2 * parameter.imp_Neighbours_No && i + j - parameter.imp_Neighbours_No < getSnpNo();   
			++j
		) {   //cerr<<"i="<<i<<"j="<<j<<endl; 
			CorrelationCoeffs[j] = fabs ( computeCorrelation(i, j) ); //WARNING could be j could be less then 0
			CorrelationMatrix[ i + j - parameter.imp_Neighbours_No ][ 2 * parameter.imp_Neighbours_No - j ] = (float) CorrelationCoeffs[j];
		}
		
		//  successors do not exist
		for ( j = 2 * parameter.imp_Neighbours_No; i + j - parameter.imp_Neighbours_No >= getSnpNo(); --j ) {
			CorrelationCoeffs[j] = -10.0;
		}
		
		// for sorting
		for ( j = 0; j <= 2 * parameter.imp_Neighbours_No; ++j ) {
			NeighboursOrder[j] = j;
		}
		
		// sort the neighbours w.r.t ascending correlation coefficients 
		SortVec pSort( 2 * parameter.imp_Neighbours_No + 1, NeighboursOrder, CorrelationCoeffs );
		
		// determine the abs. frequency of genotypes for the i-th SNP 
		Frequency[0] = Frequency[1] = Frequency[2] = 0;
int GME=0; //should be a local variable instead of 3 function calls
		for (j = 0; j < getIdvNo(); ++j )
                   { GME= getGenoMatElement(i,j);
			if ( GME == 0)
				++ Frequency[0] ;
			else if ( GME == 1)
				++ Frequency[1] ;
			else if ( GME == 2)
				++  Frequency[2] ;
                               }		
		// use highest frequency as predictor
		if (Frequency[0] > Frequency[1] && Frequency[0] > Frequency[2])
			PredictedAPriori = 0;
		else if (Frequency[1] >= Frequency[0] && Frequency[1] >= Frequency[2])
			PredictedAPriori = 1;
		else
			PredictedAPriori = 2;
		
		// check for all individuals (j) if the i-th SNP is missing
		for (j = 0; j < getIdvNo(); j ++)
		{
			if (x = getGenoMatElement(i,j), x != 0 && x != 1 && x != 2)
			{
				MissingEntries++;
				UsePriorPrediction = false;
				
				// search non-missing entries until parameter.imp_Best_SNPs_No SNPs are reached, sorted by correlation coefficent

				for (
					k = 0, l = 2 * parameter.imp_Neighbours_No;
					k < parameter.imp_Best_SNPs_No;
					++k, --l
				) {
					while (
						x = getGenoMatElement( i + pSort.getId( l ) - parameter.imp_Neighbours_No, j ),
						0 != x && 1 != x && 2 != x
					) {
						--l;
					}
					
					// no strongly correlated SNPs with missing SNP are found, so the most-frequent is used for prediction 	
					if ( pSort.getValue(l) < 0.0)
					{
						UsePriorPrediction = true;
						break;
					}
					
					BestSNPs[k] = pSort.getId(l); // position of stronly correlated SNP
					// genotype of stronly correlated SNP
					Pattern[k] = getGenoMatElement( i + BestSNPs[k] - parameter.imp_Neighbours_No, j );
				}
				
				Predicted = PredictedAPriori;
				// strongly correlated SNPs with missing SNP are found
				if (!UsePriorPrediction)
				{
					Frequency[0] = Frequency[1] = Frequency[2] = 0;
					// check for each individual if Pattern is matched
					for (l = 0; l < getIdvNo(); l ++)
					{
						Conformity = true;
						for ( k = 0; Conformity && k < parameter.imp_Best_SNPs_No; ++k  )
							if (
								getGenoMatElement( i + BestSNPs[k] - parameter.imp_Neighbours_No, l )
								!=
								Pattern[k]
							) {
								Conformity = false;
							}
						if (Conformity)
						{
							if (getGenoMatElement(i,l) == 0)
								Frequency[0] ++;
							else if (getGenoMatElement(i,l) == 1)
								Frequency[1] ++;
							else if (getGenoMatElement(i,l) == 2)
								Frequency[2] ++;
						}
					}
					if (Frequency[0] > Frequency[int(PredictedAPriori)] && Frequency[0] > Frequency[1] && Frequency[0] > Frequency[2])
						Predicted = 0;
					else if (Frequency[1] > Frequency[int(PredictedAPriori)] && Frequency[1] >= Frequency[0] && Frequency[1] >= Frequency[2])
						Predicted = 1;
					else if (Frequency[2] > Frequency[int(PredictedAPriori)] && Frequency[2] >= Frequency[0] && Frequency[2] > Frequency[1])
						Predicted = 2;
				}
				
				setGenoMatElement(i,j, Predicted);	
			}
		}
		if (MissingHelp != MissingEntries)
		{
			MissingSNPs++;	
		}
	}	
	
	// deallocate arrays
	delete [] CorrelationCoeffs;
	delete [] NeighboursOrder;
	delete [] BestSNPs;
	delete [] Pattern; 
		
	printLOG("Imputation finished: "+int2str(MissingEntries)+ " missing values in "+int2str(MissingSNPs)+" SNPs replaced");
	parameter.imp_is_imputated = true;
}



void MData::calculateIndividualTests()
{
	
	printLOG("Start Individual Tests");
	
	Model SingleSNP( *this );	// create Model with current MData
	int*	SNPList = new int[getSnpNo()]; // to store information for SortVec snp_order_
	double*	TestStat = new double[getSnpNo()]; // to store information for SortVec snp_order_
// return allways  1 when asking for omp_get_num_threads()
	cerr<<"!parameter.affection_status_phenotype="<<!parameter.affection_status_phenotype<<endl;
	cerr<<"omp_get_num_threads()"<<omp_get_num_threads()<<endl;
	if(0==parameter.affection_status_phenotype)
	{cerr<<"reel"<<endl;	omp_set_num_threads(1);}
	else
	{cerr<<"1"<<endl;   omp_set_num_threads(1);}
cerr<<"omp_get_num_threads()"<<omp_get_num_threads()<<endl;
//
       #pragma omp parallel for	
	for ( int i = 0; i < getSnpNo(); ++i ) {
		SNPList[i]=i; // for sorting, store positon 
		TestStat[i]= SingleSNP.computeSingleRegressorTest(i); // compute p-value of single marker test, and store for sorting
		//(snps_.at(i))->setSingleMarkerTest(TestStat[i]); // store single  marker test in SNP
		
	}
 omp_set_num_threads(4);
for ( int i = 0; i < getSnpNo(); ++i ) 
	(snps_.at(i))->setSingleMarkerTest(TestStat[i]);// store single  marker test in SNP
	
	snp_order_.fillVec(getSnpNo(),SNPList, TestStat); // sort the SNPs w.r.t ascending p-values in snp_order_
	
	
	// output the SNP order an the p-values in a file
	ofstream IT;
	IT.open( ( parameter.out_file_name + "_IT.txt" ).c_str(), ios::out );
	
	IT << "SNP_no. \t SNP_name \t Chr \t Pos \t p-value" << endl;
	for ( int i = 0; i < getSnpNo(); ++i ) {
		IT 	<< snp_order_.getId(i)<< "\t" 
			<< (snps_.at( snp_order_.getId(i))->getSnpId() )  <<"\t"
			<< (snps_.at( snp_order_.getId(i))->getChromosome())<< "\t"
			<< (snps_.at( snp_order_.getId(i))->getBasePairPosition())<< "\t"
			<< snp_order_.getValue(i)	<< endl;
	}
	
	IT.close();
	
	delete[] SNPList;
	delete[] TestStat;
	printLOG( "Individual Tests finished, written to \"" + parameter.out_file_name + "_IT.txt\"" );

}
/**calculate PValuePorder needs parameter.ms_MaximalPValueForwardStep but set in the conf-file
 *and of course the MData variable.
 determines the SNPs with p-Value < parameter.ms_MaximalPValueForwardStep
 these are tested in the Forward Step
 */
int MData::calculatePValueBorder() const
{	int PValueBorder = getSnpNo()-1;
	while (
		0 < PValueBorder
		&&
		snp_order_.getValue(PValueBorder) > parameter.ms_MaximalPValueForwardStep
	) {
		--PValueBorder;
	}
//checking this will bring something when PValueBorder will be very small
PValueBorder=min(100,PValueBorder);
if (PValueBorder<100) //hg tip 303 
PValueBorder = min( (size_t) 100, getSnpNo()-1 );	//this is only for very small SNP files. REMARK<BB>: how about 0 SNPs?
return PValueBorder;
}
int MData::calculatePValueBorder() 
{	int PValueBorder = getSnpNo()-1;
	while (
		0 < PValueBorder
		&&
		snp_order_.getValue(PValueBorder) > parameter.ms_MaximalPValueForwardStep
	) {
		--PValueBorder;
	}
//when returning 0 the forward step cycle though all SNPs that is like selecting PValueBorder
//parameter.ms_MaximalPValueForwardStep=1
//

//checking this will bring something when PValueBorder will be very small
PValueBorder=min(350,PValueBorder);
if (PValueBorder<350) //hg tip 302 
PValueBorder = min( (size_t) 350, getSnpNo()-1 );	//this is only for very small SNP files. REMARK<BB>: how about 0 SNPs?
return PValueBorder;

}
bool MData::selectModel ( Model *currentModel, size_t PValueBorder, int ExpectedCausalSNPs, int maxModel ) {
       int JJ=0;
	printLOG("Model Selection started: ");
        bool forgetreplaceonce=false;	
	bool 	stop = false;
	int 	addedSNP=-1;
//	int     removedSNP=-1;
	int *startIndex; //start at the begin 
	int dummy=0;
	     startIndex=&dummy;
       // int PValueBorder =calculatePValueBorder();
	//or take the setting from the conf file
	PValueBorder = min( getSnpNo() - 1, PValueBorder );	// REMARK<BB>: how about 0 SNPs?
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
        //	*currentModel,// = &model0, //,original 
		*forwardModel = &model1,
		*backwardModel = &model2;
	int backup_for_causal_snps=parameter.ms_ExpectedCausalSNPs;
       
	forgetreplaceonce=true;
	parameter.ms_ExpectedCausalSNPs= ExpectedCausalSNPs;
        cerr<<parameter.ms_ExpectedCausalSNPs<<endl;

	currentModel->makeMultiForwardStep(PValueBorder,1,startIndex);
//if(forgetreplaceonce)
//	parameter.ms_ExpectedCausalSNPs=backup_for_causal_snps;
//cerr<<parameter.ms_ExpectedCausalSNPs<<endl;
//Achtung hier soll wieder alles beim alten sein
	currentModel->computeRegression();
       // best=currentModel;
	// TODO<BB>: Apply a memory scheme to avoid duplicate calculation of Models
	// which have been calculated in a previous search step.
	// I.e. need a central registry of calculated Models.
	// But these with empty caches to sayum remove yum-packagekitve Memory.
	// Effectively a map from subsets of SNP-indices to MSC values.
        // currentModel->printModel("");
	printLOG( "Start stepwise selection" );
        currentModel->computeMSC(0);
        //currentModel->printModel("check1");
double bestMSC= currentModel->getMSC();
//currentModel->replaceModelSNPbyNearCorrelated(0);
//exit(0);
bool improvment=false;
double locMSC= DBL_MAX;
	while ( !stop ) {
         improvment=currentModel->replaceModelSNPbyNearFromCAT(*startIndex , PValueBorder);
	 improvment=currentModel->saveguardbackwardstep( *backwardModel);	 
         
	 //ED break because of computational limitations
	 if (currentModel->getModelSize()>min(parameter.maximalModelSize,maxModel))
		 break;
        /* linear case normal forward step
         */
          if (!parameter.affection_status_phenotype)//quantitative
	      improvment=currentModel->makeForwardStepLinear( forwardModel,  JJ,  &bestMSC,
		                                	 PValueBorder, startIndex);
          else if (parameter.affection_status_phenotype)/*PRÄSELECTION nur bis Revision 274*/
              improvment=currentModel->makeForwardStepLogistic(JJ, &bestMSC,  PValueBorder, startIndex);	  
          else 
		  cerr<<" not linear nor logistic, this should not happen"<<endl;
	  stop=currentModel->finalizeModelSelection( *backwardModel, JJ,  improvment,  PValueBorder,  startIndex);
}//while
cerr<<parameter.ms_ExpectedCausalSNPs<<endl;
size_t reference=350;	// REMARK<BB>: Where does 350 come from? Also mind 0 SNPs case below.
	if( parameter.affection_status_phenotype)
      {if (350>getSnpNo()-1) reference=getSnpNo()-1;
	    
	while(currentModel->selectModel(*backwardModel,max(reference,PValueBorder),maxModel)) //minimum 100 or PValueBorder//with 1000 instead 100 it takes 2' on a model with 17 SNPS
		{
		 //currentModel->replaceModelSNPbyNearCorrelated1();//should
		 improvment=currentModel->replaceModelSNPbyNearFromCAT(*startIndex , PValueBorder);
		 currentModel->saveguardbackwardstep( *backwardModel);
		}
	} 
//when all fails	

//currentModel->replaceModelSNPbyNearCorrelated(); bringt nichts
//currentModel=backwardModel;
currentModel->printStronglyCorrelatedSnps( 0.999, int2str(parameter.in_values_int) + "the_result" );

currentModel->printModel("finalModel");
//reset the .ms_ExpectedCausalSNPs
cerr<<parameter.ms_ExpectedCausalSNPs<<endl;

if(forgetreplaceonce)
	parameter.ms_ExpectedCausalSNPs=backup_for_causal_snps;
cerr<<parameter.ms_ExpectedCausalSNPs<<endl;

return true ; 
}

///////////////////////////////////////////////////7
//selctModel ohne Argument
bool MData::selectModel()
{
       int JJ=0;
	printLOG("Model Selection started: ");
	Model SelectedModel( *this ); 
	Model Forward( *this );
	Model Backward( *this );
	
	bool 	stop = false;
	int 	addedSNP=-1;
//	int     removedSNP=-1;
	int *startIndex; //start at the begin 
	int dummy=0;
	     startIndex=&dummy;
       // int PValueBorder =calculatePValueBorder();
	//or take the setting from the conf file
	size_t PValueBorder = parameter.PValueBorder;
	PValueBorder = min( getSnpNo()-1, PValueBorder );	// REMARK<BB>: how about 0 SNPs?
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
        /*this is the original multi forward step only used here*/
//vector<string> LDSNPs={"SNP_A-2020880","SNP_A-1815642","SNP_A-2175818","SNP_A-4262563","SNP_A-4260995","SNP_A-1923094","SNP_A-2275264","SNP_A-4218005","SNP_A-2220513","SNP_A-4297940","SNP_A-1985005","SNP_A-1941585", "SNP_A-2104965","SNP_A-2225019","SNP_A-2207150","SNP_A-2089059","SNP_A-2107431","SNP_A-1865752","SNP_A-2086212","SNP_A-2141521","SNP_A-1800559","SNP_A-2301126","SNP_A-1801278","SNP_A-2160387","SNP_A-1930022","SNP_A-2302790","SNP_A-4289013","SNP_A-2134130","SNP_A-2163978","SNP_A-1780694","SNP_A-4277866","SNP_A-2062455","SNP_A-2084317","SNP_A-2033911","SNP_A-1956074","SNP_A-1913873","SNP_A-2217489","SNP_A-2073899","SNP_A-1870534","SNP_A-4291857","SNP_A-2038959","SNP_A-2155630","SNP_A-1953574","SNP_A-1859780","SNP_A-4221344","SNP_A-2197439","SNP_A-2010616","SNP_A-2010617","SNP_A-4215015","SNP_A-2010618","SNP_A-2121761","SNP_A-2079004","SNP_A-2203549","SNP_A-1986190","SNP_A-1825815","SNP_A-1786036","SNP_A-1986192","SNP_A-2214585","SNP_A-2125952","SNP_A-4211683","SNP_A-4211684","SNP_A-4289986","SNP_A-2135688","SNP_A-1987242","SNP_A-4265870"};
//vector <string> LSNPNames(LDSNPs);

  // for(int i=0;i<LSNPNames.size();i++)
        //cerr<<"'"<<LSNPNames[i]<<"'"<<endl;
        //vector< unsigned int> Lindex;
        //findSNPIndex(LSNPNames, Lindex);
        //currentModel->addManySNP(Lindex);
     
        //currentModel->computeRegression();
        //currentModel->computeMSC(0);
        //currentModel->printStronglyCorrelatedSnps( 0.999, int2str(parameter.in_values_int) + "LD_fullmodel" );
        //currentModel->saveguardbackwardstep( *backwardModel);
int backup_for_causal_snps=parameter.ms_ExpectedCausalSNPs;
if(parameter.expected_causal_snps1>parameter.ms_ExpectedCausalSNPs)
parameter.ms_ExpectedCausalSNPs=parameter.expected_causal_snps1;
	currentModel->makeMultiForwardStep(PValueBorder,1,startIndex);
if(parameter.expected_causal_snps1>parameter.ms_ExpectedCausalSNPs)
	parameter.ms_ExpectedCausalSNPs=backup_for_causal_snps;
//Achtung hier soll wieder alles beim alten sein
	currentModel->computeRegression();
       // best=currentModel;
	// TODO<BB>: Apply a memory scheme to avoid duplicate calculation of Models
	// which have been calculated in a previous search step.
	// I.e. need a central registry of calculated Models.
	// But these with empty caches to sayum remove yum-packagekitve Memory.
	// Effectively a map from subsets of SNP-indices to MSC values.
  //      currentModel->printModel("");
	printLOG( "Start stepwise selectionseine Leistungen in osteuropäischen Ländern über die MPA in Ungarn " );
        currentModel->computeMSC(0);
        //currentModel->printModel("check1");
double bestMSC= currentModel->getMSC();
//currentModel->replaceModelSNPbyNearCorrelated(0);
//exit(0);
bool improvment=false;
double locMSC= DBL_MAX;
/** This is the backward step
 *  one take the full model and remove every SNP.
 *  Then one look for that model with  lowest MSC 
 *  then we have a new start model.
 *
 *  We repeat this until the 1 model is selected.
 */
//int breakfor=0;
	while ( !stop ) {
		// compute steps
	// currentModel->replaceModelSNPbyNearCorrelated1(); //at first search for an improvment locally
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
	    
	while(currentModel->getModelSize()&&currentModel->selectModel(*backwardModel,max(reference,PValueBorder))) //minimum 100 or PValueBorder//with 1000 instead 100 it takes 2' on a model with 17 SNPS
		{
		 //currentModel->replaceModelSNPbyNearCorrelated1();//should
		 improvment=currentModel->replaceModelSNPbyNearFromCAT(*startIndex , PValueBorder);
		 currentModel->saveguardbackwardstep( *backwardModel);
		}
	} 


//currentModel->replaceModelSNPbyNearCorrelated(); bringt nichts

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
	int*		SNPList = new int[getSnpNo()];
	double*		TestStat = new double[getSnpNo()];
	
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
				getSNP(Pos)->setSingleMarkerTest(pValue); // set Data
				CheckSNP.at(Pos)=true;
				SNPList[i]=Pos;
				TestStat[i]=pValue;
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
		snp_order_.fillVec(getSnpNo(),SNPList, TestStat); // set the SNP order in Model (a bit redundant, since sortet)
				
	}
	else
	{
		printLOG("Could-Not open file: \""+ filename +"\" as SNP-Orderfile.");
		exit(3);
	}


	
	//~ delete[] SNPList; //fehlt
	//~ delete[] TestStat;

}

/** Covariates */
void MData::readCovariablesFile()
{
	ifstream		COV;
	string			headline;
	string 			buffer;
	vector<string>	headers;	
	vector< vector<string> > TM;
	vector< string > temp;
	int 			i;
	double			acovarible;
	
	//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	// read in data form .cov file
	printLOG("Covariables will be read!");
        cerr<<parameter.cov_file_name.c_str()<<endl;
	COV.open( parameter.cov_file_name.c_str(), ios::in );
	
	parameter.covariables = 0; 				// the number of quantitive covariables
   	parameter.dummy_covariables = 0; 			// the number of qualitativ covariables
	int headersize;
	int lines=0;								// the number of lines with covarialbles read in 
		//printLOG("Covariables 1");
	if (COV.is_open())
	{
		
		// read in headers 
		getline(COV,headline,'\n'); // get first line of file
		stringstream SS(headline);
		
		while( SS >> buffer) // read in headers and store them in vector
		{
			headers.push_back(buffer);
		}
		
		headersize = headers.size();
		//cerr<<headersize;

		if ( headersize <= 2) // no information about y-valuses
		{
			printLOG( "File \"" + parameter.cov_file_name + "\" does not contain enough information to obtain covariables." );	
			exit(9);
		}
		
	//printLOG("Covariables 2");

		
		
		while (! COV.eof())
		{
			temp.clear();
			getline(COV,headline,'\n');
			stringstream SS(headline);
			
			while( SS >> buffer) // read in headers and store them in vector
			{
				//~ cout <<"buffer: "<< buffer << endl;
				temp.push_back(buffer);
			}
			//~ cout << "size: "<<  temp.size()<< endl;
			if ( temp.size() != 0 ) // do not check empty lines
			{
					TM.push_back(temp);
					lines++;
					if ( int( temp.size() ) != headersize ) 
					{
						printLOG( "Error in covariables file \""+ parameter.cov_file_name + "\": Number of Variables do not match." );
						exit(13);
					}
			}	
		}
	}
	else
	{
		printLOG( "Error: could not open covariables file \"" + parameter.cov_file_name + "\"" );
		exit(14);
	}
	
	//printLOG("Covariables 3");
	
	//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	// use data form .cov file
	
	//~ vector<bool> selectedVar(headersize, true); 
	vector<bool> isCategorical(headersize, false); 
	
	
	
	// test phase, missing selection of covariables, setting
		
		vector<bool> selectedVar(headersize, true); 
		//selectedVar.at(6)=true; //age
	
	
		
		//~ cout << "size :" << TM.size() << endl;
		//~ 
		//~ 
		//~ cout << "print"<< endl;
		//~ 
		//~ for( j= 0; j<6; j++)
		//~ {
				//~ for (i =0; i< 3; i++)
				//~ {
			//~ cout << TM[j][i]<< " ";
				//~ cout << TM[0][0] << endl;
				//~ cout << (*TM.at(j)).at(0);
				//~ //cout << temp.at(0);
				//~ }
			//~ cout << endl;
		//~ }
	
		// match covaribles to the right individual
		
		map<string, int> uid; 				// map to search for the individuals
		map<string, int>::iterator it_uid; 	// iterator of map to search for the individuals
		for ( i=0; i < getIdvNo(); i++ )
		{
			uid.insert(make_pair (individuals_.at(i)->getFamilyID()+"_"+individuals_.at(i)->getIndividualID(), i ));      // set pairs with identifier for individual + the position
		}
		
		//~ for(it_uid = uid.begin(); it_uid != uid.end();  it_uid++ )
		//~ {
			//~ cout << it_uid -> first << " " << it_uid -> second << endl;	
		//~ }
		
		
		// find the position of the read-in individuals
		int* PosList = new int[lines];
		
		for ( i=0; i < lines; i++ )
		{
			it_uid = uid.find( TM[i][0]+"_"+TM[i][1] ); // TM[i][0] gives the 1 entry of the line, the fid
			
			if (it_uid != uid.end() ) // indvidual was found
			{
				PosList[i] = it_uid -> second;
			}
			else // indvidual not found
			{
				PosList[i] = -1; // individual not (no longer) in MData
			}
			
		}
		
		// on the i-th position of int* PosList is now the positon of the individual (in MData) which was in the i-thl line of 
		
		//~ cout << "hier3"<<endl;
		//~ 
		//~ for ( i=0; i < lines; i++ )
		//~ {
			//~ cout << PosList[i]<< endl;
		//~ }
	covMat.upSize(lines,0);	
		for ( i=2; i< headersize; i++ ) // start with 2, since 0, 1 is fid, iid
		{
			if ( selectedVar.at(i))
			{
				if (isCategorical.at(i))
				{
				
				} else {	// variable is selected and quantitative
					Vector aQuantCovV = covMat.newColumn();

					// initialise all entries as missing
					aQuantCovV.fill( parameter.missing_phenotype_code );

					for ( int j = 0; j < lines; ++j ) {
						if (PosList[j] >= 0)  // if a position is returned
						{
							stringstream SS(TM[j][i]); // write to stringstream
							SS >> acovarible; // convert read-in string to a double
						 // cerr<<"j="<< j<<PosList[j]<<" "<<acovarible<<endl;
							aQuantCovV.set( PosList[j], acovarible );
						}	
					}
					++parameter.covariables;
					Cov_Names_.push_back( headers[i] ); 	// add the name of the covariable
				}
			}
		}
	
		

	
		//~ vector<double>::iterator it_osai;	
		//~ cout<<endl;
		//~ for(it_osai =headers.begin(); it_osai < headers.end(); it_osai++ )
		//~ {
			//~ cout <<	(*it_osai) <<" ";
		//~ }
		//~ cout<<endl;
		
	
	delete[] PosList;
	
	//~ for (j=0; j < getIdvNo(); j++)
	//~ {
		//~ cout << getCovMatElement(0, j)<<" ";
	//~ }
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
  for ( int i = 0; i < getSnpNo(); ++i ) {
    IT  << snp_order_.getId(i)<< "\t" 
        << (snps_.at( snp_order_.getId(i))->getSnpId() )  <<"\t"
        << (snps_.at( snp_order_.getId(i))->getChromosome())<< "\t"
        << (snps_.at( snp_order_.getId(i))->getBasePairPosition())<< "\t"
        << snp_order_.getValue(i) << endl;
  }
  IT.close();
  printLOG( "Individual Tests finished, written to \"" + parameter.out_file_name + "_IT.txt\"" );
}


