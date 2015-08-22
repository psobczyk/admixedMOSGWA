/********************************************************************************
 *	This file is part of the MOSGWA program code.				*
 *	Copyright ©2011–2015, Erich Dolejsi, Bernhard Bodenstorfer.		*
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
#include <cassert>
#include "Model.hpp"
#include "GenotypeFreq.hpp"
#include "PermSort.hpp"
#include "logging/Logger.hpp"
#include "io/PlinkInput.hpp"
#include "io/Hdf5Input.hpp"
#include <sstream>
#include <map>
#include <memory>
#include <limits>	// for nan(...)
#include <cmath>	// for nan(...)
#include <cfloat>	// for maximal double
#include <omp.h>
#include <hdf5.h>

using namespace linalg;
using namespace logging;
using namespace io;

////++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//  class MData
////++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void MData::getXcolumn ( const size_t snp, Vector& vector ) const {
	inputCache->retrieveGenotypeVector( snp, vector );
}

const string& MData::getCovMatElementName( const size_t cov ) const {
	return covNames.at( cov );
}

void MData::getCovariateColumn ( const size_t cov, Vector& vector ) const {
	input->retrieveCovariateVector( cov, vector );
}

void MData::getY ( Vector& vector ) const {
	input->retrievePhenotypeVector( parameter->in_values_int, vector );
}

void MData::setLL0M ( const double ll ) {
	loglikelihood0Model_ = ll;
}

/** Default Constructor: reads the input-files, sets parameters, deallambda.hpps with missing phenotypes */
MData::MData ( io::Input * const externalInput ) : input( externalInput ), allocateInput( NULL == externalInput ), covMat( 0, 0 ) {
	if ( allocateInput ) {
		if ( parameter->in_file_hdf5.empty() ) {
			input = new PlinkInput( parameter->in_files_plink.c_str() );
			if ( 0 == parameter->nSNPKriterium ) {
				parameter->nSNPKriterium=getSnpNo();
			}
		} else {
			input = new Hdf5Input( parameter->in_file_hdf5.c_str(), parameter->cov_extra_file );
		}
	}
	inputCache.reset( new CachedInput( *input, parameter->cache_limit ) );

	const size_t
		snps = input->countSnps(),
		idvs = input->countIndividuals(),
		covs = input->countCovariates();
	        //ED setting default Value for nSNPKriterium when not set 
	        if(0==parameter->nSNPKriterium)
                  parameter->nSNPKriterium=snps;
	covMat.exactSize( idvs, covs );
	singleMarkerTestResult.resize( snps, numeric_limits<double>::signaling_NaN() );
	const std::string * covariates = input->getCovariates();
	for ( size_t cov = 0; cov < covs; ++cov ) {
		covNames.push_back( covariates[ cov ] );
		Vector covVec = covMat.columnVector( cov );
		input->retrieveCovariateVector( cov, covVec );
	}
	Y_name_ = input->getTraits()[parameter->in_values_int];

	checkData();
}

MData::~MData () {
	if ( allocateInput ) {
		delete input;
		input = NULL;
	}
}

size_t MData::getSnpNo () const {
	return input->countSnps();
}

size_t MData::getIdvNo () const {
	return input->countIndividuals();
}

size_t MData::getCovNo () const {
	return covNames.size();
}

string MData::getFID ( const size_t index ) const {
	assert( index < getIdvNo() );
	const Individual * individual = input->getIndividuals() + index;
	return individual->getFamilyID();
}

string MData::getID ( const size_t index ) const {
	assert( index < getIdvNo() );
	const Individual * individual = input->getIndividuals() + index;
	return individual->getIndividualID();
}

const SNP & MData::getSNP ( const size_t snp ) const {
	assert ( snp < input->countSnps() );
	return input->getSnps()[snp];
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

double MData::getSingleMarkerTestAt ( const size_t index ) const {
	return singleMarkerTestResult.at( index );
}

void MData::setSingleMarkerTestAt ( const size_t index, const double value ) {
	singleMarkerTestResult.at( index ) = value;
}

void MData::fillSnp_order_Vec ( const size_t snpNo, size_t* SNPList, double* TestStat ) {
	snp_order_.fillVec( snpNo, SNPList, TestStat );
}

void MData::checkData () {
	const size_t idvs = getIdvNo();
	parameter->affection_status_phenotype = true;	// intitial setting, Y-values are tested if the are just (0,1) (or missing)
	
	AutoVector yVec( idvs );
	getY( yVec );
	for ( size_t idv = 0; idv < idvs; ++idv ) {
		const double indPheno = yVec.get( idv );
		// do not consider individuals with missing phenotype
		if ( indPheno == parameter->missing_phenotype_code ) {
			logger->warning(
				"missing individuals' phenotype for individual \"%s %s\"",
				getID( idv ).c_str(),
				getFID( idv ).c_str()
			);
		}
		else // no missing phenotype, determine if phenotype is affection (case-control) or quantitative
		{
			if ( indPheno == parameter->case_value)	{	// was set 1 check affection status
			} else if ( indPheno == parameter->control_value ) {	//was set 0
			} else {
				logger->info(
					"individual[%u]( \"%s\", \"%s\" ) has phenotype code %f"
					" which is neither case (%f) nor control (%f)."
					" Will use linear regression.",
					idv,
					getID( idv ).c_str(),
					getFID( idv ).c_str(),
					indPheno,
					parameter->case_value,
					parameter->control_value
				);
				parameter->affection_status_phenotype = false;
			}
		}
	}
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
		allSNPs[snp] = getSNP( snp ).getSnpId();
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
			// TODO<BB>: int should not be used for SNP indices, but size_t
			logger->warning(
				"SNP %s at position %d does not exist nicht in this dataset",
				SNPNames.at( i ).c_str(),
				i
			);
		}
	}
vector<int>::iterator iter,nend;

nend=remove(tindex.begin(),tindex.end(),-999);
	for ( iter = tindex.begin(); iter != nend; ++iter ) {
		logger->debug( "Position: %d is %s", *iter, allSNPs.at( *iter ).c_str() );
		index.push_back( *iter );
	}

 //remove duplicates
 //that model works
 sort(index.begin(),index.end() );

 index.erase( unique( index.begin(), index.end() ), index.end() );
/*for_each(index.begin(),index.end(),
		boost::lambda::if_then(boost::lambda::_1==-999,
		cout<< boost::lambda::_1 <<"is not available\n"));*/
}


void MData::printSelectedSNPsInR ( vector<string> SNPList ) const {
	const size_t idvs = getIdvNo();
	AutoVector vec( idvs );	// allocate outside loop
	ofstream	SNPL;

	SNPL.exceptions ( ofstream::eofbit | ofstream::failbit | ofstream::badbit ); // checks if Logfile can be writt
	try {
		SNPL.open( ( parameter->out_file_name + "_SNPList.txt" ).c_str(), fstream::out );
		vector<string>::iterator it_SNPList;
		
		for (it_SNPList = SNPList.begin(); it_SNPList < SNPList.end(); it_SNPList++)
		{
			size_t snp;
			for (
				snp = 0;
				snp < getSnpNo() && getSNP( snp ).getSnpId() != *it_SNPList;
				++snp
			);

			if ( snp < getSnpNo() ) {
				getXcolumn( snp, vec );
				SNPL << *it_SNPList <<" <- c(";
				for ( size_t idv = 0; idv < idvs; ++idv ) {
					SNPL << ( 0 < idv ? "," : "" );
					SNPL << vec.get( idv );
				}
				SNPL<< ")"<< endl;
			} else {
				logger->error( "SNP %s not found", it_SNPList->c_str() );
			}
		}
		
		for ( size_t cov = 0; cov < getCovNo(); ++cov ) {
			SNPL << getCovMatElementName( cov ) << " <- c(";
			getCovariateColumn( cov, vec );
			for ( size_t idv = 0; idv < idvs; ++idv ) {
				SNPL << ( 0 < idv ? "," : "" );
				SNPL << vec.get( idv );
			}
			SNPL<< ")"<< endl;	
		}

		SNPL << "Y" <<" <- c(";
		getY( vec );
		for ( size_t idv = 0; idv < idvs; ++idv ) {
			SNPL << ( 0 < idv ? "," : "" );
			SNPL << vec.get( idv );
		}
		SNPL<< ")"<< endl;

		SNPL << "Intercept" <<" <- c(";
		for ( size_t idv = 0; idv < idvs; ++idv ) {
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
		logger->info(
			"Written R-File \"%s_SNPList.txt\".",
			parameter->out_file_name.c_str()
		);
	} catch ( ofstream::failure e ) {
		logger->error( "Could not write R-File: %s", e.what() );
	}
}

/** printSelectedSNPsInMatlab create an *Octave.h5 file with the selected SNPs (as a vector of strings)  and their X and of course the phenotype Y
 * additionaly it writes a */ 
void MData::printSelectedSNPsInMatlab ( vector<string> SNPList , string extra) const  {
	const size_t idvs = getIdvNo();
	AutoVector vec( idvs );
	ofstream	SNPL;
         {        	ofstream	SNPL;
		 //create the same in hdf
		 hid_t file,fid,dataset,space,/*dset,memtype,*/ filetype,props;
                 herr_t status;
                 hsize_t dim[] = { SNPList.size(), idvs };	//transponiert
                 hsize_t di[]={SNPList.size()};
		 double  zwischen[getIdvNo()*SNPList.size()];
		 double  Y[idvs];
		 vector<string>::iterator it_SNPList;
		 int i;
		logger->info(
			"out_file=%s%sOctave.h5",
			parameter->out_file_name.c_str(),
			extra.c_str()
		);
		file = H5Fcreate(
			( parameter->out_file_name + extra + "Octave.h5" ).c_str(),
			H5F_ACC_TRUNC,
			H5P_DEFAULT,
			H5P_DEFAULT
		);
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
vector<unsigned int> index;
findSNPIndex(SNPList, index);
//sortd according the names1

	for ( size_t jj = 0; jj < index.size(); ++jj ) {
		getXcolumn( index[jj], vec );
		for( size_t ii = 0; ii < idvs; ++ii ) {
		       	zwischen[jj*getIdvNo()+ii] = vec.get(ii);
		}
	}
	// TODO: refactor with Vector Y to avoid additional copy
	getY( vec );
        for ( size_t idv = 0; idv < idvs; ++idv ) {
		Y[idv] = vec.get( idv );
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
dim[0]=idvs;
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
		SNPL.open( ( parameter->out_file_name + "_SNPList.m" ).c_str(),  fstream::out );
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
				snp < getSnpNo() && getSNP( snp ).getSnpId() != *it_SNPList;
				++snp
			);
			
			if ( snp < getSnpNo() ) {
				getXcolumn( snp, vec );
				for ( size_t idv = 0; idv < idvs; ++idv ) {
					SNPL << ( 0 < idv ? " " : "" );
					SNPL << vec.get( idv );
				}
				SNPL<< ";"<< endl;
			} else {
				logger->error( "SNP %s not found", it_SNPList->c_str() );
			}
		}
			SNPL << "]'"<<endl;
		
		for ( size_t cov = 0; cov < getCovNo(); ++cov ) {
			SNPL << getCovMatElementName( cov ) << " = [";
			getCovariateColumn( cov, vec );
			for ( size_t idv = 0; idv < idvs; ++idv ) {
				SNPL << ( 0 < idv ? ";" : "" );
				SNPL << vec.get( idv );
			}
			SNPL<< "]'"<< endl;	
		}

		SNPL << "Y" <<" = [";
		getY( vec );
		for ( size_t idv = 0; idv < idvs; ++idv ) {
			SNPL << ( 0 < idv ? ";" : "" );
			SNPL << vec.get( idv );
		}
		SNPL << "]" << endl;
		SNPL.close();
		logger->info(
			"Written MATLAB-File \"%sSNPList.m\"",
			parameter->out_file_name.c_str()
		);
	} catch ( ofstream::failure e ) {
		logger->error( "Could not write MATLAB-File: %s", e.what() );
	}
}


/** computes correlation between two snps, utilise the structur of the data */
double MData::computeCorrelation ( const size_t locus1, const size_t locus2 ) const {
	const size_t idvs = getIdvNo();
	AutoVector
		v1( idvs ),
		v2( idvs );
	getXcolumn( locus1, v1 ),
	getXcolumn( locus2, v2 );
	double
		sum1 = 0.0,
		sum2 = 0.0;
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


void MData::calculateIndividualTests()
{
	const size_t snps = getSnpNo();
	
	logger->info( "Start Individual Tests" );
	
	Model singleSNP( *this );	// create Model with current MData
	auto_ptr<size_t> snpArray( new size_t[ snps ] );
	auto_ptr<double> TestStat( new double[ snps ] );

	#pragma omp parallel for private(singleSNP)
	for ( size_t snp = 0; snp < snps; ++snp ) {
		if ( 0 == singleSNP.getModelSize() ) {
			singleSNP.addSNPtoModel( snp );
		} else {
			singleSNP.replaceSNPinModel( snp, 0 );
		}
		// for sorting, store position
		snpArray.get()[snp] = snp;
		// compute p-value of single marker test and store for sorting
		TestStat.get()[snp] = singleSNP.computeSingleRegressorTest();
	}
	omp_set_num_threads( 4 );
	for ( size_t snp = 0; snp < snps; ++snp ) {
		setSingleMarkerTestAt( snp, TestStat.get()[snp] );
	}
	
	// sort the SNPs w.r.t ascending p-values in snp_order_
	snp_order_.fillVec( getSnpNo(), snpArray.get(), TestStat.get() );

	// output the SNP order an the p-values in a file
	ofstream IT;
	IT.open( ( parameter->out_file_name + "_IT.txt" ).c_str(), ios::out );
	
	IT << "SNP_no. \t SNP_name \t Chr \t Pos \t p-value" << endl;
	for ( size_t snp = 0; snp < getSnpNo(); ++snp ) {
		const size_t snpi = snp_order_.getId( snp );
		const SNP& snpo = getSNP( snpi );
		IT 	<< snpi << "\t"
			<< snpo.getSnpId() <<"\t"
			<< snpo.getChromosome() << "\t"
			<< snpo.getBasePairPosition() << "\t"
			<< snp_order_.getValue( snp ) << endl;
	}
	
	IT.close();

	logger->info(
		"Individual Tests finished, written to \"%s_IT.txt\"",
		parameter->out_file_name.c_str()
	);
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
		snp_order_.getValue( PValueBorder ) > parameter->ms_MaximalPValueForwardStep;
		--PValueBorder
	);
//checking this will bring something when PValueBorder will be very small
	PValueBorder = min( 100u, (unsigned int) PValueBorder );
	if ( PValueBorder < 100 ) {
		PValueBorder = min( 100u, (unsigned int) getSnpNo()-1 );	//this is only for very small SNP files. REMARK<BB>: how about 0 SNPs?
	}
	return PValueBorder;
}
bool MData::selectModel (
	Model *currentModel,
	size_t PValueBorder,
	int maxModel,
	const int selectionCriterium
) {
	logger->info( "Model Selection started:" );
	bool 	stop = false;
//	int     removedSNP=-1;
	int *startIndex; //start at the begin 
	int dummy=0;
	     startIndex=&dummy;
       // int PValueBorder =calculatePValueBorder();
	//or take the setting from the conf file
	PValueBorder = min( getSnpNo() - 1, PValueBorder );	// REMARK<BB>: how about 0 SNPs?
//	 int PValueBorder=parameter.PValueBorder;
	PValueBorder = min( getSnpNo() - 1, PValueBorder );
	logger->debug( "PValueBorder", PValueBorder );
	// compute the log-likelihood of the 0 Model
	Model model0( *this );
	
	model0.computeRegression();
	model0.printYvec();
	setLL0M( model0.getMJC() );	// REMARK<BB>: functionality seems redundant with that in main.cpp
	Model model1( *this ), model2( *this );
        //Mod
	// Using pointers to avoid expensive copying of Model objects
	Model
        //	*currentModel,// = &model0, //,original 
		*forwardModel = &model1,
		*backwardModel = &model2;
       
//	schreibt es in die Variable für mBIC2 auch wenn man mBIC haben will
	if ( 0 < currentModel->getModelSize() ) {
		if ( 0.000001 < currentModel->computeMSC( selectionCriterium ) ) {
			//only here a backwardstep makes sense
			currentModel->saveguardbackwardstep( *backwardModel, selectionCriterium );
		}
	}

	currentModel->makeMultiForwardStep( PValueBorder, selectionCriterium, startIndex );
	currentModel->computeRegression();

	logger->info( "Start stepwise selection" );
double bestMSC= currentModel->getMSC();
//currentModel->replaceModelSNPbyNearCorrelated(0);
//exit(0);
	bool improvement = false;
	while ( !stop ) {
		improvement = currentModel->replaceModelSNPbyNearFromCAT( *startIndex, PValueBorder, selectionCriterium );
		improvement = currentModel->saveguardbackwardstep( *backwardModel, selectionCriterium );
         
	 //ED break because of computational limitations
		if ( currentModel->getModelSize() > min( parameter->maximalModelSize, maxModel ) ) {
			break;
		}
        /* linear case normal forward step
         */
		if ( !parameter->affection_status_phenotype ) {		//quantitative
			improvement = currentModel->makeForwardStepLinear( forwardModel, &bestMSC, PValueBorder, startIndex );
		} else {	/*PRÄSELECTION nur bis Revision 274*/
			improvement = currentModel->makeForwardStepLogistic( &bestMSC, PValueBorder, startIndex, selectionCriterium );
		}
		stop = currentModel->finalizeModelSelection(
			*backwardModel,
			improvement,
			PValueBorder,
			startIndex,
			selectionCriterium
		);
	}
	size_t reference = 350;		// REMARK<BB>: Where does 350 come from? Also mind 0 SNPs case below.
	if ( parameter->affection_status_phenotype ) {
		if ( 350 > getSnpNo() - 1 ) {
			reference = getSnpNo() - 1;
		}
	    
		while(
			currentModel->selectModel(
				*backwardModel,
				max( reference, PValueBorder ),
				maxModel,
				selectionCriterium
			)
		) {
			//minimum 100 or PValueBorder//with 1000 instead 100 it takes 2' on a model with 17 SNPS
		 //currentModel->replaceModelSNPbyNearCorrelated1();//should
		 //improvment=currentModel->replaceModelSNPbyNearFromCAT(*startIndex , PValueBorder, criterium);
	         //if(improvment)	 currentModel->saveguardbackwardstep( *backwardModel,criterium);
		}
	}
		improvement = currentModel->replaceModelSNPbyNearFromCAT(
			*startIndex , PValueBorder, selectionCriterium
		);
		if ( improvement ) {
			currentModel->saveguardbackwardstep( *backwardModel, selectionCriterium );
		}

//when all fails	

//currentModel->replaceModelSNPbyNearCorrelated(); bringt nichts
//currentModel=backwardModel;
	currentModel->printStronglyCorrelatedSnps(
		parameter->correlation_threshold,
		int2str(parameter->in_values_int) + "the_result"
	);
	currentModel->printModel( "finalModel", selectionCriterium );
	return true ; 
}

///////////////////////////////////////////////////7
//selectModel ohne Argument
bool MData::selectModel()
{
	logger->info( "OLD Model Selection started:" );

	Model SelectedModel( *this ); 
	Model Forward( *this );
	Model Backward( *this );
	
	bool 	stop = false;
	int *startIndex; //start at the begin 
	int dummy=0;
	     startIndex=&dummy;
       // int PValueBorder =calculatePValueBorder();
	//or take the setting from the conf file
	size_t PValueBorder = parameter->PValueBorder;
	size_t eins=1;
	PValueBorder = max(min( getSnpNo()-1, PValueBorder ),eins);	// REMARK<BB>: how about 0 SNPs?
	logger->debug( "PValueBorder", PValueBorder );
	// compute the log-likelihood of the 0 Model
	Model model0( *this );
	
	model0.computeRegression();
        model0.printYvec();
	setLL0M( model0.getMJC() );
	Model model1( *this ), model2( *this );
        //Mod
	// Using pointers to avoid expensive copying of Model objects
	Model
		*currentModel = &model0,
		*forwardModel = &model1,
		*backwardModel = &model2;
	currentModel->makeMultiForwardStep(PValueBorder,1,startIndex);
	currentModel->computeRegression();
       // best=currentModel;
	logger->info( "Start stepwise selection:" );
        currentModel->computeMSC(0);
        //currentModel->printModel("check1");
double bestMSC= currentModel->getMSC();
	bool improvement=false;
/*  This is the backward step
 *  one take the full model and remove every SNP.
 *  Then one look for that model with  lowest MSC 
 *  then we have a new start model.
 *
 *  We repeat this until the 1 model is selected.
 */
	while ( !stop ) {
		improvement = currentModel->replaceModelSNPbyNearFromCAT( *startIndex , PValueBorder );
		improvement = currentModel->saveguardbackwardstep( *backwardModel );
	 currentModel=backwardModel; 
        /* linear case normal forward step
         */
		if ( !parameter->affection_status_phenotype ) {
			improvement = currentModel->makeForwardStepLinear(
				forwardModel,
				&bestMSC,
				PValueBorder,
				startIndex
			);
		} else {
			improvement = currentModel->makeForwardStepLogistic( &bestMSC, PValueBorder, startIndex );
		}
		stop = currentModel->finalizeModelSelection(
			*backwardModel,
			improvement,
			PValueBorder,
			startIndex
		);
	}
	size_t reference = 350;	// REMARK<BB>: Where does 350 come from?
	if( parameter->affection_status_phenotype ) {
		if ( 350 > getSnpNo()-1 ) reference=getSnpNo()-1;	// REMARK<BB>: Beware 0 SNPs case.
		while (
			currentModel->getModelSize()
			&&
			currentModel->selectModel( *backwardModel, max( reference, PValueBorder ) )
		) {
			//minimum 100 or PValueBorder//with 1000 instead 100 it takes 2' on a model with 17 SNPS
			improvement = currentModel->replaceModelSNPbyNearFromCAT( *startIndex , PValueBorder );
			currentModel->saveguardbackwardstep( *backwardModel);
		}
	}
	currentModel->printStronglyCorrelatedSnps(
		parameter->correlation_threshold,
		int2str(parameter->in_values_int) + "the_result"
	);
	currentModel->printModel("finalModel");
	return true;
}


/** Artur new code:
   * @brief writes ordered snps to file. File name is set up in parameters class.
*/
void MData::printSnpOrder()
{
  // output the SNP order an the p-values in a file
  ofstream IT;
	IT.open( ( parameter->out_file_name + "_IT.txt" ).c_str(), ios::out );
  
  IT << "SNP_no. \t SNP_name \t Chr \t Pos \t p-value" << endl;
	for ( size_t snp = 0; snp < getSnpNo(); ++snp ) {
		const size_t snpi = snp_order_.getId( snp );
		const SNP& snpo = getSNP( snpi );
		IT	<< snpi << "\t" 
			<< snpo.getSnpId() <<"\t"
			<< snpo.getChromosome() << "\t"
			<< snpo.getBasePairPosition() << "\t"
			<< snp_order_.getValue( snp ) << endl;
	}
  IT.close();
	logger->info(
		"Individual Tests finished, written to \"%s_IT.txt\"",
		parameter->out_file_name.c_str()
	);
}
