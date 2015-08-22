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

#include "Model.hpp"
#include "GenotypeFreq.hpp"
#include "QRuncher.hpp"
#include "Exception.hpp"
#include "lookup/package.hpp"
#include "logging/Logger.hpp"
#include <fstream>
#include <cfloat>	// for maximal double
#include <cmath>	// for isinf
#include <cstdio>	//getpid rand number init
#include <cassert>

//for random number generator
#include <gsl/gsl_rng.h>

// for the various classes of random  numbers
#include <gsl/gsl_randist.h>

// for Probability Distributions
#include <gsl/gsl_cdf.h>

// for Linear Algebra
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
#include <algorithm> //min und max
#include <sys/time.h> //time

#include <hdf5.h> //for hdf5
#include "ScoreTestShortcut.hpp" //SortVec scoreTests ( const Model& model, const MData& mData )
using namespace std;
using namespace linalg;
using namespace lookup;
using namespace logging;
bool DEBUG=false,DEBUG2=false,DEBUG3=false;
////++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//  class Model
////++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



Model::Model ( const MData& mData )
:
	data_( &mData ),
	xMat( mData.getIdvNo(), 1 + mData.getCovNo() ),
	yVec( mData.getIdvNo() ),
	beta( 1 + mData.getCovNo() )
{
	const size_t covs = data_->getCovNo();
	size_t col = 0;
	xMat.columnVector( col++ ).fill( 1.0 );		// intercept is a column of 1
	for ( size_t cov = 0; cov < covs; ++cov ) {
		Vector xVec = xMat.columnVector( col++ );
		data_->getCovariateColumn( cov, xVec );
	}
	data_->getY( yVec );
	initializeModel();
}


void Model::initializeModel () {
	const size_t
		idvs = data_->getIdvNo(),
		msnps = getModelSize(),
		cols = getNoOfVariables();
	size_t col = 1 + data_->getCovNo();
	xMat.upSize( idvs, cols );
	for ( size_t modelSnp = 0; modelSnp < msnps; ++modelSnp ) {
		const size_t snp = modelSnps_.at( modelSnp );
		Vector xVec = xMat.columnVector( col++ );
		data_->getXcolumn( snp, xVec );
	}
	assert( cols == col );

	beta.upSize( cols );
	beta.fill( 0.0 );
	upToDateBetas_ = false;
}

Model& Model::operator= ( const Model & orig ) {
	if ( &orig != this ) {
		data_ = orig.data_;
		modelSnps_ = orig.modelSnps_; 
		xMat = orig.xMat;
		yVec = orig.yVec;
		beta = orig.beta;
		upToDateBetas_ = orig.upToDateBetas_;
		msc = orig.msc;
		modelJudgingCriterion_ = orig.modelJudgingCriterion_;
	}
	return *this;
}

/** TODO<BB>: Replace private variable {@link modelSnps_} by a {@link lookup::MutableMultiIndex}
* to take full advantage of the modern multi-index approach.
*/
ModelIndex Model::getIndex () const {
	return ModelIndex( modelSnps_ );
}

int Model::getModelSize () const {
	return modelSnps_.size();
}

int Model::getNoOfVariables () const {
	return 1 + data_->getCovNo() + getModelSize();
}

double Model::getMJC () const {
	return modelJudgingCriterion_;
}

double Model::getBeta ( const int i ) const {
	if ( i >= getNoOfVariables() ) {
		throw;
	} else {
		return beta.get( i );
	}	
}

string 	Model::getSNPId ( const size_t i ) const {
	const size_t size = getModelSize();
	if ( size <= i ) {
		logger->error(
			"Get requested SNP(%u) beyond modelsize %u you will see an SNP_false instead",
			i,
			size
		);
	//	throw; //simple but bad because this string is only for reporting!
	return "ERROR SNP";
	} else {	
		return data_->getSNP( modelSnps_.at(i) ).getSnpId();
	}
}

void Model::sortSNPsAccordingBetas () {
  vector<size_t> SNP( getNoOfVariables(), 0 );
  SortVec sbetas(getNoOfVariables());
  vector<double> betas(getNoOfVariables(),0);
 for (int i=0;i<getNoOfVariables();i++) //this assumes that getNoOfVariables	
 { betas[i]=getBeta(i);
   SNP[i] = modelSnps_[i];
 }
sbetas.fillVec(getNoOfVariables(), &SNP[0],&betas[0],false);
}

void Model::addSNPtoModel ( const size_t snp ) {
	const size_t
		idvs = data_->getIdvNo(),
		vars = getNoOfVariables();

	// TODO<BB>: I do not see: How is adding a SNP duplicate avoided?
	// only when you use addmanySNP 
	modelSnps_.push_back( snp );
	xMat.upSize( idvs, vars + 1 );
	beta.upSize( vars + 1 );
	
	Vector xVec = xMat.columnVector( vars );
	data_->getXcolumn( snp, xVec );
	beta.set( vars, 0.0 );
	upToDateBetas_= false; 
}

/** replaces replaceSNPinModel 
 * set snp at position  
 * position is the place where from 0 to  getNoOfVariables()
 * that says that even covariates could be replaced!!!!
 *
 * */

bool Model::replaceSNPinModel ( const size_t snp, const size_t position ) {
	if ( 0 == getModelSize() ) {
		logger->error(
			"Cannot replace SNP %u at position %u, since the model is empty",
			snp,
			position
		);
		return false;
	}
	if ( getModelSize() <= position ) {
		logger->error( "replaceSNPinModel position %u is not in Model", position );
		return false;
	}
	const size_t reset = 1 + data_->getCovNo() + position;		//position 0 is the first 
 //cout<<"reset="<<reset<<endl;
	modelSnps_.at( position ) = snp;
// cerr<<"reset="<<reset;
	// add the new column at the end
		Vector xVec = xMat.columnVector( reset );
		data_->getXcolumn( snp, xVec );
		beta.set( reset, 0.0 );		// 0.0 is a relatively good starting point for a new variable
		upToDateBetas_= false;		// but it is not the precise coefficient and others will change, too
return true;	
}

/** removes SNP from Model, SNP is relativ position at vector modelSnps_ this is covariate aware */
bool Model::removeSNPfromModel ( const size_t msnp ) {
	const size_t
		idvs = data_->getIdvNo();

	// check if SNP is a valid postition. 
	if ( 0 > msnp || getModelSize() <= msnp ) {
			return false;
	}

	modelSnps_.erase( modelSnps_.begin() + msnp );

	const size_t colsMinus1 = getNoOfVariables();
	for ( size_t col = 1 + data_->getCovNo() + msnp; col < colsMinus1; ++col ) {
		const Vector source = xMat.columnVector( col + 1 );
		Vector target = xMat.columnVector( col );
		target.copy( source );
		beta.set( col, beta.get( col + 1 ) );
	}
	xMat.upSize( idvs, colsMinus1 );	// reduce dimensions (without deallocation)
	beta.upSize( colsMinus1 );
	upToDateBetas_ = false;
	return true;
}

//**adds many SNP to the Model

void Model::addManySNP ( vector<size_t> selected ) {
	//tests
	if (
		0 == selected.size()
		||
		data_->getSnpNo () < selected.size()
	) {
		logger->error(
			"addManySNP selected SNP vector is of size 0 or is bigger than the number of SNP"
		);
		// TODO<BB>: do not use exit in library code
		exit( 1 );
	}
	//add further testing here
	// dubletts in selected
	// using < as comparison 
	
	//instead using sort use SortVec
	sort (selected.begin(), selected.end()); // sorts all
	
	if ( data_->getSnpNo () < selected.back() ) {
		logger->error( "addmanySNP index in selected is bigger as Number of SNP" );
		// TODO<BB>: throw Exception
		return;
	}
	vector<size_t> finalselect;
	//reserve enough room for adding all elements of selected
	finalselect.reserve( selected.size() );
	finalselect.push_back( selected[0] );
        /*remove doubletts */
	for( int i = 1; i < selected.size(); ++i ) {
		if( selected[i-1] < selected[i] ) {
			finalselect.push_back(selected[i]);
		}
	}
 
	// add snp in orderd sequence
  	for( int i = 0; i < finalselect.size(); ++i ) {
		addSNPtoModel(finalselect[i]);
	}
}

	/** set \beta could be used for generating Y
            the \betas should be  gsl_vector* Beta      
 */

void Model::setBeta ( const Vector& newBeta ) {
	beta.copy( newBeta );
	upToDateBetas_ = true;
}

/***
 *Ycontinous() 
 * the beta is set in advance
 * set an ran_gaussian this should be the continous!
 *
 */
void Model::Ycontinous () {
	const size_t idvs = data_->getIdvNo ();
	AutoVector tVec( idvs );
	tVec.gemv( 1.0, xMat, false, beta, 0.0 );
	for ( size_t idv = 0; idv < idvs; ++idv ) {
		logger->debug( "Y[%u] = %f", tVec.get( idv ) );
	}

//these initialisations are because of the random number generator
  const gsl_rng_type * T; 
 gsl_rng * r;//unknown
 gsl_rng_env_setup();
 T = gsl_rng_default;
 r = gsl_rng_alloc (T);
 gsl_rng_set (r, time (NULL));//every start another random number

 //here the file
string FID,ID;
ofstream Y;
	Y.exceptions ( ofstream::eofbit | ofstream::failbit | ofstream::badbit );
	try {
		Y.open( ( parameter->out_file_name + ".yvm" ).c_str(), fstream::out );
	} catch ( ofstream::failure e ) {
		logger->error(
			"Could not open additional pheno file %s for octave: %s",
			( parameter->out_file_name + ".yvm" ).c_str(),
			e.what()
		);
	}
// the header line is not that good when to use as  an input to octave
// but heres she is:
Y<<"FID IID ";
	for ( int j = 1; j <= parameter->replications; ++j ) {	//0 is not good 
		Y << j << " ";
	}
	Y << endl;
//this was the first line of yvm

    // the header line is not that good when to use as  an input to octave
    // but heres she is:
for (int i=0;i<idvs;++i)
   { FID=data_->getFID(i);
     ID=data_->getID(i);
     Y <<FID << " "<<ID<<" "; //hier werden die Variablen aufgerufen, das geschieht jede Zeile
		for( int j = 0; j < parameter->replications; ++j ) {
			Y << ( tVec.get( i ) + gsl_ran_gaussian( r, 1 ) ) << " ";
		}
		Y << endl;	//after all a newline
   }
		try {
			Y.close();
		} catch ( ofstream::failure e ) {
		logger->error(
			"Could not  additional pheno file %s for octave:",
			( parameter->out_file_name + ".yvm" ).c_str(),
			e.what()
		);
	}
 gsl_rng_free (r); 
}

void Model::printYvec () {
	ofstream Y;
	Y.exceptions ( ofstream::eofbit | ofstream::failbit | ofstream::badbit );
	try {
		Y.open( ( parameter->out_file_name + "Yvecout" ).c_str(), fstream::out );
		const size_t idvs = yVec.countDimensions();
		for ( size_t idv = 0; idv < idvs; ++idv ) {
			Y << yVec.get( idv ) << "\n";
		}
		Y.close();
	}
	catch ( ofstream::failure e ) {
		logger->error(
			"Could not write Yvecout file %s: %s",
			( parameter->out_file_name + "Yvecout" ).c_str(),
			e.what()
		);
	}
}

//REM 
void Model::printModel (
	const string& out,
	const int selectionCriterium,
	const string& filemodifier
) {
	stringstream ss; // to save output
	ofstream OUT; // output model to file
	const size_t
		covs = data_->getCovNo(),
		offset = 1 + covs;

	// generate output
	ss << out << endl;
	ss << "ModelSize: " << getModelSize()
		<< "\tMSC: " << computeMSC( selectionCriterium )
		<<"\tMJC: "<< getMJC() << endl;
	ss << "Nr \tSNPId       \tChr \tPos   \tbeta    \tp-Value Single Marker Test" <<endl;

	for ( int i = 0; i < getModelSize(); ++i ) {
		const SNP & snp = data_->getSNP(modelSnps_.at(i));
		ss << modelSnps_.at(i)<< "\t" 
				<< getSNPId(i) << "\t"
				<< snp.getChromosome()<<"\t"
				<< setw(10)<<snp.getBasePairPosition()<<"\t"   //the setw for formating 
				<< getBeta( offset + i ) << "\t"
				<< data_->getSingleMarkerTestAt( modelSnps_.at(i) )
				<< endl;
	}
        if(0==getModelSize())
		ss<<"\tIntercept \t \t \t is not available empty model"<<endl;
	else
	ss << "\tIntercept \t \t \t" << getBeta(0)<< endl;

	for ( int i = 0; i < covs; ++i ) {
		ss << "\t" << data_->getCovMatElementName(i) << "\t\t\t" << getBeta( 1 + i ) << endl;
	}

	// output to screen
	cout << ss.str();
	
	// output to file
	try {
		OUT.open( ( parameter->out_file_name +filemodifier+ ".mod" ).c_str(), ios::out );
		OUT << ss.str();
		OUT.close();
	} catch ( ofstream::failure e ) {
		logger->error(
			"Could not write Modelfile \"%s%s.mod\": %s",
			parameter->out_file_name.c_str(),
			filemodifier.c_str(),
			e.what()
		);
	}
}

/**
 *
 *printStronglyCorrelatedSnps  gives a number of correlated SNP with there
 *correlations from 1 to a threshold
 *this threshold is the parameter of the function
 */
void Model::printStronglyCorrelatedSnps ( const double threshold,  string extra ) const {
	
	stringstream 			ss; 			// to save output
	ofstream 				OUT; 			// output model to file
	multimap<double, int> 	StrongCor;		// to sort correlated SNPs
        const int fenster=400; //keine Fenster mehr left and right side of the fenster
        const int reserve=2*fenster+1; // if fenster=1 then there is the left and the right and the middle
	double                  zwischen[2*reserve];  //only here for hdf5
	const char                    *S[reserve]={NULL};
	int 					count=0;			// count strongly correlated SNPs for a given SNP
	double					abscor;			// |Correlation| for two SNPs
         //ED
	 //
	size_t snps = data_->getSnpNo();
//HDF5 
stringstream name;
hid_t file,fid,dataset,space,dset,/*memtype,*/ filetype;
herr_t status;
hsize_t dim[]={0,2};//dim für Vector der korrelierte
hsize_t di[]={0};
//if (0<getModelSize())
//{	
	file = H5Fcreate(
		( parameter->out_file_name + extra + "Corr.h5" ).c_str(),
		H5F_ACC_TRUNC,
		H5P_DEFAULT,
		H5P_DEFAULT
	);
//HDF5
	for ( size_t i = 0; i < getModelSize(); ++i ) {
		count =0;
	
		// search for strongly correlated SNPs
		// PARALLELISIERBAR kann GEWALTIG GRO
		//for (int j = 0; j < data_->getSnpNo(); ++j )
		//ternary operator :make
		//(logical expression)?if_true : if false
	
		size_t low = max( 0, (int)modelSnps_.at(i)-fenster );//unsigned is wrong  !
		const size_t high = min( snps, modelSnps_.at(i) + fenster );
		int      line=0;
		
		for ( size_t j = low; j < high; ++j ) {	//für 0U Literale füur unsigned int	
				if(false)
					abscor = data_->computeCorrelation( modelSnps_.at(i), j ) ;
		         	else//sometime one wants the signed version 
					abscor = fabs( data_->computeCorrelation( modelSnps_.at(i), j ) );
			// add if correlation is big enough
			 if ( abscor  >= threshold ) {
			 // if (fabs(abscor)>=threshold) {//for the signed version one is interesstet in the high correlated snps
				StrongCor.insert(pair<double, int>(abscor, j));
				S[line] = data_->getSNP(j).getSnpId().c_str();
				zwischen[line*2]=abscor;
				zwischen[line*2+1]=j;
				++line;
				++count;
			}
		 } 
//HDF5 
//here all the dataset generation and closing afterwards
        dim[0]=count;//dim[1]=2;

	fid=H5Screate_simple(2,dim,NULL);
	/**Create  dataset 
	 * take an name + number of selected SNP
	 * for the variable
	 */
	name.str("");//in the other case you get very long names!
	name<<"co"<<i; // 
	//some more H5Pcreate 

	 dataset = H5Dcreate2( file, name.str().c_str(), H5T_NATIVE_DOUBLE, fid, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );//name.str().c_str();
	// if (0>=dataset)
	// { /** write the data to the dataset
	// *  and then close nothing because it works 
	// */
	status=H5Dwrite(dataset,H5T_NATIVE_DOUBLE,H5S_ALL,H5S_ALL,H5P_DEFAULT,zwischen);
	status=H5Dclose (dataset);
        status=H5Sclose (fid);
        //nun die scheußlichen SNP name in die Datei schreiben
        //H5T_NATIVE_CHAR
        di[0]=count;

	filetype = H5Tcopy (H5T_C_S1);
        status = H5Tset_size (filetype,H5T_VARIABLE);
	/* create the datatype properties                                 */
	space = H5Screate_simple (1, di, NULL);
        /*
         * Create the dataset and write the variable-length string data to
         * it.
         */
        dataset = H5Dcreate2( file, (name.str()+"SNP").c_str(), filetype, space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );
      	status = H5Dwrite (dataset,filetype , H5S_ALL, H5S_ALL, H5P_DEFAULT, S);
	status = H5Tclose(filetype);	
	status = H5Sclose(space);
	status = H5Dclose(dataset);


	    // save strongly correlated SNPs for output 
		ss << getSNPId(i) << " has an absolute correlation of at least "<< threshold  <<" with the following "<< count << " SNPs:\n";
	
		for (multimap<double, int>::reverse_iterator it = StrongCor.rbegin(); it != StrongCor.rend(); it++)	
		{
			const int snpIndex = (*it).second;
			const SNP & snp = data_->getSNP( snpIndex );
			ss 	<< snpIndex
				<< "\t" << snp.getSnpId() 
				<< "\t" << snp.getChromosome() 
				<< "\t" << snp.getBasePairPosition() 
				<< "\t" << (*it).first << "\n";
		}
		StrongCor.clear();
	
	}
	
	// output to screen
	cout << ss.str();
	
	// output to file
	try {
        status = H5Fclose(file);//the file will closed!
		logger->info(
			"Close %s_%s_Corr.txt with status: %d",
			parameter->out_file_name.c_str(),
			extra.c_str(),
			status
		);
		OUT.open( ( parameter->out_file_name + "_" + extra + "_Corr.txt" ).c_str(), ios::out );
		OUT << ss.str();
		OUT.close();
	} catch( ofstream::failure e ) {
		logger->error(
			"Could not write the correlation file \"%s_%s_Corr.txt\": %s",
			parameter->out_file_name.c_str(),
			extra.c_str(),
			e.what()
		);
	}
}

/**
 * replaceModelSNPbyNearFromCAT 
 * if an SNP s enters the Model in the sequence of the accending p-values of the CAT test,
 * occasional the first variable wich enters improves the model  but there is an variable t situated near to s 
 *  which is even better than the selected variable.
 *  currentPosition could be set to 0 search all or the psoition from the forwardstep
 * */
// REMARK<BB> 0 should not be a reserved value, it is a valid SNP index.
bool Model::replaceModelSNPbyNearFromCAT (
	int currentPosition,
	const size_t PValueBorder,
	const int selectionCriterium
) {
	// current position is the position returned from fast forward,
	// you have to use this before reset
	// find the position of the SNP which has entered the Model;
	// the last SNP in the model will be checked here!
	if (0==getModelSize())
		return false;
	bool improve=false;
	int fenster=50; //be bold
	const int grace=5000;//war 500
	const size_t nSNP = data_->getSnpNo();
	logger->debug(
		"replaceModelSNPbyNearFromCAT currentPosition=%d window=%d grace=%d last model SNP=%d",
		 currentPosition,
		 fenster,
		 grace,
		 modelSnps_.at( getModelSize() - 1 )
	);
	double bestMSC=getMSC(); //the msc in the current model is the best one up to now
	Model model0( *data_ );
/*random permutation
 */
	const size_t N=getModelSize();
	//DEBUG cerr<<"N="<<getModelSize()<<endl; //DEBUG
	const gsl_rng_type *T;
	gsl_rng *r;
	gsl_permutation *a =gsl_permutation_alloc (N);
	gsl_permutation_init(a);
	//DEBUG gsl_permutation_fprintf(stdout,a," %u"); //DEBUG
    //gsl_permutation *q =gsl_permutation_alloc (N);
	gsl_rng_env_setup();
 T = gsl_rng_default;
 r = gsl_rng_alloc (T);
gsl_ran_shuffle (r, a->data, N, sizeof (size_t));//reference to data!!!

//gsl_permutation_fprintf (stdout, a, " %u");
//	for(int jo=getModelSize()-1; jo >=0; jo--)
for(int jo=0; jo<getModelSize(); jo++)
{ int j=gsl_permutation_get(a,jo);
	//DEBUG cerr<<"a["<<jo<<"]="<<j<<endl;
	for (
		size_t i = 0;
		i < min( max( PValueBorder + grace, 1000u ), nSNP );
		++i
	)	// safeguard against overrun of nSNP and only 500SNP in large Problems, with 1000 in most cases all relevant SNP will found
	        { //cerr<<(abs(ref-data_-> getOrderedSNP(i))<50)<<",";
		       	if (
				abs( (int) modelSnps_.at( j ) - (int) data_->getOrderedSNP( i ) ) < fenster
				&&
				data_->getOrderedSNP( i ) != modelSnps_.at( j )	//the original model should not be considered!
			) {
		          model0=*this;
				model0.replaceSNPinModel( data_->getOrderedSNP( i ), j );
       	                  double val=DBL_MAX;
                          model0.computeRegression(); //regression should be calculated!
				val = model0.computeMSC( selectionCriterium );
                    //  printLOG("ModelSNP="+ int2str(j) + "POS="+ int2str(data_-> getOrderedSNP(i))+"val=" +double2str(val));

                  if(val<bestMSC-0.0001)
		       	//saveguard against rouning errors in logistic regression
                      { double  alt =bestMSC;
			 bestMSC=val;
				logger->debug(
					"Better Model at position %u SNP=%u"
					" is replaced with %u oldMSC=%f newMSC=%f",
					j,
					modelSnps_.at(j),
					data_->getOrderedSNP( i ),
					alt,
					bestMSC
				);
				model0.printModel( "Replaced_inter", selectionCriterium );
				*this=model0;
				computeMSC( selectionCriterium );	// REMARK<BB>: compute again?
                         improve=true; 
		      }
                      }
		}
}//rand order

gsl_permutation_free (a);
gsl_rng_free (r);//
return improve;


}//replaceModelSNPbyNearFromCAT ends

bool Model::replaceModelSNPSCORE ( const int selectionCriterium ) {
 bool improve=false; //we don't know weather we could improve
 int fenster=50; //be bold
 const int grace=5000;
 double bestMSC=getMSC(); //the msc in the current model is the best one up to now
	Model model0( *data_ );
	SortVec score;
 	for(int j=getModelSize()-1; j >=0; j--)
	{//replace
		scoreTestWithOneSNPless( j, score );
		const size_t nscores = score.size();
		for ( size_t i = 0; i < nscores; ++i ) {
			model0 = *this;
                        model0.replaceSNPinModel (score.getId(i) ,  j );

	                double val=DBL_MAX;
                          model0.computeRegression(); //regression should be calculated!
			val= model0.computeMSC( selectionCriterium );
                
			if(bestMSC>0?val<0.9999*bestMSC:val<1.0002*bestMSC)
	                      { double  alt =bestMSC;
				 bestMSC=val;
				logger->debug(
					"Better Model at Modelposition %u"
					" SCOREPostition=%u"
					" SNP=%u"
					" is replaced with %u"
					" oldMSC=%u"
					" newMSC=%f",
					j,
					i,
					modelSnps_.at(j),
					score.getId(i),
					alt,
					bestMSC
				);
					model0.printModel( "Replaced_inter", selectionCriterium );
	                         *this=model0;   
	                         improve=true; 
			      }

		}
	}
return improve;
}
//--------------------------------------------
/** The code here uses modern tools and points into the direction of future code-consolidation.
* TODO<BB>: Store all calculated models' MSCs in a ResultStore.
* Or (cheaper) at least tick those off which are suboptimal.
*/


double Model::oraculateOptimalLinearForwardStep ( size_t *snp, size_t bound ) const {
	AutoVector xVec( data_->getIdvNo() );

	// Fast incremental linear regression calculator
	QRuncher qruncher( yVec );

	// Import the coefficient matrix into the QRuncher
	for ( size_t col = 0; col < xMat.countColumns(); ++col ) {
		// TODO: Take care of return value against adding linearly dependent columns
		qruncher.pushColumn( const_cast<AutoMatrix&>( xMat ).columnVector( col ) );
	}

	// To quickly search whether SNP is already "in"
	const ModelIndex modelIndex = getIndex();

	// Search best to add
	double bestRSS = DBL_MAX;
	size_t bestSNP;
	for ( size_t snpCol = 0; snpCol < bound; ++snpCol ) {
		const size_t orderedSnpIndex = data_->getOrderedSNP( snpCol );
		if ( modelIndex.contains( orderedSnpIndex ) ) continue;

		// Prepare new column
		data_->getXcolumn( orderedSnpIndex, xVec );
		qruncher.pushColumn( xVec );

		const double rss = qruncher.calculateRSS();
		// TODO<BB>: Here we should store the RSS in a ResultStore to avoid duplicate calculation
		if ( rss < bestRSS ) {
			bestRSS = rss;
			bestSNP = orderedSnpIndex;
		}
		qruncher.popColumn();
	}
	if ( bestRSS < DBL_MAX ) {
		*snp = bestSNP;
	}
	return bestRSS;
}

/** The code here uses modern tools and points into the direction of future code-consolidation.
* TODO<BB>: Store all calculated models' MSCs in a ResultStore.
* Or (cheaper) at least tick those off which are suboptimal.
* @see oraculateOptimalLinearForwardStep( size_t, size_t* )
*/
double Model::oraculateOptimalLinearBackwardStep( size_t *snp ) const {
	// Fast linear regression calculator with Bernhard's backward shortcut algorithm
	QRuncher qruncher( yVec );

	// Import the coefficient matrix into the QRuncher
	for ( size_t col = 0; col < xMat.countColumns(); ++col ) {
		// TODO: Take care of return value against adding linearly dependent columns
		qruncher.pushColumn( const_cast<AutoMatrix&>( xMat ).columnVector( col ) );
	}

	// To quickly search whether SNP is already "in"
	const ModelIndex modelIndex( modelSnps_ );

	// Ignore first (fixed) columns not corresponding to SNPs
	// TODO<BB>: Handle this differently once xMat will be pre-transformed by the fixed columns' Householder vectors
	const size_t colOffset = getNoOfVariables() - getModelSize();	// = 1 + covariables

	// Search best to remove
	double bestRSS = DBL_MAX;
	size_t bestCOL;
	
	for ( size_t col = 0; col < getModelSize(); ++col ) {
		const double rss = qruncher.calculateSkipColumnRSS( col + colOffset );
		logger->debug( "Oraculate column %u yields rss %f", col, rss );
		if ( rss < bestRSS ) {
			bestRSS = rss;
			bestCOL = col;
		}
	}
	if ( bestRSS < DBL_MAX ) {
		*snp = modelSnps_.at( bestCOL );
		logger->debug( "Oraculate SNP %u at column %u", *snp, bestCOL );
	}
	return bestRSS;
}

/** finalizeModelSelection()
 *   finalization work call some function and quit
 */
bool Model::finalizeModelSelection (
	Model &backwardModel,
	bool improvement,
	const size_t PValueBorder,
	int *startIndex,
	vector<size_t> score,
	const int selectionCriterium
) {
	if ( !improvement ) {
		printModel( "no improvement", selectionCriterium );
		*startIndex = 0;
		if ( !parameter->affection_status_phenotype ) {
			logger->info( "finalise with score is not implemented for continuous traits" );
			backwardModel=*this;
			improvement = saveguardbackwardstep( backwardModel, selectionCriterium);
			// not use makeMFFL in this case 
		} else {
			makeMFFL(
				max( parameter->reset, PValueBorder ),
				startIndex,
				score,
				selectionCriterium
			);
		}
		//don't set to an explicit value because of memory leaks when the number of variables is very small!
		backwardModel = *this;
		improvement = saveguardbackwardstep( backwardModel, selectionCriterium );
		if ( improvement ) {
			logger->info( "finalizeModelSelection" );
			*this=backwardModel;	// REMARK<BB>: Erich had this line after the return … ?
			return true;
		}

		printModel( "final model", selectionCriterium );
		return true;
	} else {
		logger->debug( "after forward step" );
		return false;
	}
}

bool Model::finalizeModelSelection (
	Model &backwardModel,
	bool improvement,
	const size_t PValueBorder,
	int *startIndex,
	const int selectionCriterium
) {
	if ( !improvement ) {
		printModel( "no improvement", selectionCriterium );
		*startIndex = 0;
		if ( !parameter->affection_status_phenotype ) {
			logger->info( "finalise is not implemented for continuous traits" );
			backwardModel=*this;
			improvement = saveguardbackwardstep( backwardModel, selectionCriterium);
			//not use makeMFFL in this case 
		} else {
		//diese Schritte sind eigentlich nur einmal nötig da, man sicher nichts
  //  gewinnt wenn man die ersten paar SNP oft und oft wiederholt
			logger->debug( "finalizeModelSelection last run with min of  parameter.PValueBorder or parameter.reset" );
			makeMFFL(
				min( PValueBorder, parameter->reset ),
				startIndex,
				selectionCriterium
			);
		}

  //  don't set to an explicit value because of memory leaks when the number 
  // of variables is very small!
		backwardModel = *this;
		improvement = saveguardbackwardstep( backwardModel, selectionCriterium );
		if ( improvement ) {
			logger->debug( "improvement after newstart forward step" );
			*this = backwardModel;	// REMARK<BB>: Erich had this line after the return
			return true;
			printModel( "final model", selectionCriterium );
		}
		return true;
	} else {
		logger->debug( "after forward step" );
		return false;
	}
}

bool Model::makeForwardStepLinear (
	Model *forwardModel,
	double* bestMSC,
	const size_t PValueBorder,
	int *startIndex,
	const int selectionCriterium
) {
	bool improvement = false;

	Model model3( *data_ );
 model3=*this;//when everthing fails this remains as result
 forwardModel=this;
 for(int ii=0;ii<20;ii++) //parameter 
        {	
		const int addedSNP = makeForwardStep( *forwardModel, PValueBorder, selectionCriterium );

	double locMSC=forwardModel->getMSC();
	   if(locMSC<*bestMSC)
	     {       	
			improvement = true;
                *bestMSC=locMSC;
		model3 =*forwardModel;
			logger->info(
				"better bigger model in iteration %u bestMSC=%f",
				ii,
				*bestMSC
			);
		*this=model3;
	     }
	   else
	     break;
        }
 *this =model3; //when something is updated in the loop then we will see it
	return improvement;
}
/*make ForwardStep for MData::selectModel for logistic Regression score version
 *  *this is forwardModel
 *
 */
bool Model::makeForwardStepLogistic (
	double* bestMSC,
	const size_t PValueBorder,
	int *startIndex,
	vector<size_t> score,
	const int selectionCriterium
) {
	bool improvement = false;
	if ( *startIndex < parameter->reset ) {		//ERICH wenn man 100 als Grenze nimmt sollte man das hier auch herabsetzen 
	 // 500 waren gut wenn man eine Grenze von 3000 hatte
	 // also docg eher PValueBorder/6
		*startIndex = 0;
		makeMFFL( PValueBorder, startIndex, score, selectionCriterium );
	} else {
	 //da auch 300 für 3000 
	 //daher 30 für 100
		*startIndex = max( 0, *startIndex - parameter->jump_back );
		makeMFFL( PValueBorder, startIndex, score, selectionCriterium );
	}
//REMOVED FOR SPEED  scoreTest();//more scoreTests

 double   locMSC=getMSC();
		  if(*bestMSC>locMSC)
		  { 
		improvement = true;
		   *bestMSC=locMSC;
		  }
	return improvement;
}
/*make ForwardStep for MData::selectModel for logistic Regression
 *  *this is forwardModel
 *
 */
bool Model::makeForwardStepLogistic (
	double* bestMSC,
	const size_t PValueBorder,
	int *startIndex,
	const int selectionCriterium
) {
	bool improvement = false;
	if( *startIndex < parameter->reset ) {		//500
		logger->info( "startIndex=%d, but 0 is used", *startIndex );
		*startIndex = 0;
		makeMFFL( PValueBorder, startIndex, selectionCriterium );
	} else {
		*startIndex = max( 0, *startIndex - parameter->jump_back );	//min because nothing prevents to be jump_back to be bigger than reset!
		logger->info(
			"startIndex=%d is used, with jump_back=%d",
			*startIndex,
			parameter->jump_back
		);
		makeMFFL( PValueBorder, startIndex, selectionCriterium );
	} 
//REMOVED FOR SPEED 
// scoreTest();//more scoreTests

 double   locMSC=getMSC();
		  if(*bestMSC>locMSC)
		  { 
		improvement = true;
		   *bestMSC=locMSC;
		  }
	return improvement;
}

/** makeMFFS make Fast Forward local fast forward 
 * */
bool Model::makeMFFS (
	const size_t PValueBorder,
	int* startIndex
) {
	int selectionCriterium = parameter->selectionCriterium;
         //FAST MultipleForward is needed locally!	
	 // TODO<BB>: remove unelegant parameter setting in code
        bool oldValue = parameter->ms_FastMultipleForwardStep;
        parameter->ms_FastMultipleForwardStep = true;

        const int altSize= getModelSize();
	makeMultiForwardStep( PValueBorder, selectionCriterium, startIndex );

	parameter->ms_FastMultipleForwardStep = oldValue;
	return true;
}

//and now with score test 
bool Model::makeMFFL(
	const size_t PValueBorder,
	int* startIndex,
	vector<size_t> score,
	const int selectionCriterium
) {       
//	int selectionCriterium=0;
         //MultipleForward is needed locally!	
   //     bool oldValue= parameter.ms_FastMultipleForwardStep;
        //parameter.ms_FastMultipleForwardStep=true;
        
       
        const int altSize= getModelSize();
makeMultiForwardStepScore ( PValueBorder,selectionCriterium,startIndex, score  );
	logger->info( "MFFL SCORE startIndex %d model size %u", *startIndex, altSize );
}


/**  makeMFFL make Fast Forward local but without changing to fast search*/

bool Model::makeMFFL(
	const size_t PValueBorder,
	int* startIndex,
	const int selectionCriterium
) {       
         //MultipleForward is needed locally!	
   //     bool oldValue= parameter.ms_FastMultipleForwardStep;
        //parameter.ms_FastMultipleForwardStep=true;
        
       
        const int altSize= getModelSize();
makeMultiForwardStep ( PValueBorder,selectionCriterium,startIndex  );
	logger->info( "MFFL startIndex %d model size %u", *startIndex, altSize );
}








//ScoreScoreScoreScoreScoreScoreScoreScoreScoreScoreScoreScoreScore
bool Model::makeMultiForwardStepScore (
	size_t PValueBorder,
	const int selectionCriterium,
	int* startIndex,
	vector<size_t> score
) {
	const size_t snps = data_->getSnpNo();
//if(score) 
	//getscores this means that we have logistic regression 
       // data_->getOrderedSNP should be something different for score
if(NULL==startIndex)
  {
    int dummy=0;	  
   startIndex=&dummy;
  }
  int returnIndex=*startIndex;
	logger->info( "Start Multiple-Forward-Step SCORE" );
  if (0==PValueBorder)
  {
		PValueBorder = snps;
		logger->info( "Default setting for PValuePorder: select all" );
  }
  //if modelsize is bigger than0 and selectionCriterium ~! 1
  //then new modus
int startSize=getModelSize();
int newSNP=0;

	logger->debug( "Model start size %u", startSize );
	if ( Parameter::selectionCriterium_BIC != selectionCriterium && 0 < startSize ) {
		logger->info( "new usage of Fastforward" );
		newSNP = parameter->ms_MaximalSNPsMultiForwardStep;
	} else {
		logger->info( "normal FastForward" );
		newSNP = parameter->ms_forward_step_max;
	}

 // }
  double oldBIC,newBIC;
  bool orig_affection_status_phenotype = true; //init just to init...

  // Fast multiple Forward Step for affection phenotypes: linear regression is used 
	if ( parameter->ms_FastMultipleForwardStep ) {
		logger->info( "Fast Multiple-Forward Score Used." );
		// TODO<BB>: repair ugly parameter setting in code
		orig_affection_status_phenotype = parameter->affection_status_phenotype;
		parameter->affection_status_phenotype = false;
		upToDateBetas_ = false;		// betas computed by a different regression
	}

	if ( parameter->affection_status_phenotype ) {
		Model NewModel( *data_ );
    oldBIC = computeMSC(selectionCriterium);
	
    for (
      int i = *startIndex ;
      getModelSize() </*=*/  startSize+/*parameter.ms_MaximalSNPsMultiForwardStep*/
	newSNP && i < PValueBorder && i < snps;
      ++i
    ) {
      // progress checking
		if ( 0 == i % 200 ) {
			logger->info(
				"Done %3.5f%%...",
				i / ( snps / 100.0 )
			);
		}

      NewModel = *this;
        // To quickly search whether SNP is already "in"
    const ModelIndex modelIndex( modelSnps_ ); 
    const size_t orderedSnpIndex = score[i]; //data_-> getOrderedSNP( i );
    if ( !modelIndex.contains( orderedSnpIndex ) )     {
	 NewModel.addSNPtoModel( orderedSnpIndex ); // Add SNPs according to p-Value
      
      // check if regression works properly, else try next SNP
         if (  NewModel.computeRegression() ) {
            newBIC = NewModel.computeMSC(selectionCriterium);
        
         // check if SNP is added to the Model
          if ( newBIC < oldBIC ) {
          *this = NewModel; //model is now updated 
			if ( parameter->detailed_selction ) {
				logger->info(
					"ADD SNP: %s ModelSize: %d MSC: %f",
					data_->getSNP( orderedSnpIndex ).getSnpId().c_str(),
	      				NewModel.getModelSize(),
					newBIC
				);
			}
          oldBIC = newBIC;	  
       }
	}//if snp is  in model	 
     returnIndex=i; 
     }
   }
  // Fast multiple Forward Step:  reset original parameter
  // reset parameter.affection_status_phenotype and compute 
		if ( parameter->ms_FastMultipleForwardStep ) {
			printModelNew();
			// TODO<BB>: It is unelegant and not thread-safe to manipulate the configuration object like this
			parameter->affection_status_phenotype = orig_affection_status_phenotype;
			if ( !computeRegression() ) {
				logger->error( "Fast Multiple-Forward: failed" );
				// TODO<BB>: Remove exit statement for library
				exit(12);
			}
		}
	}
	*startIndex = returnIndex++;
	if ( *startIndex >= snps ) {
		*startIndex = 0;
	}
	logger->info( "startIndex=%d", *startIndex );
	logger->info( "Finish Multiple-Forward-Step: %u SNPs have been added",  getModelSize() );
	return true;
}

/* this is the main forward step in the model creation 
 * */
bool Model::makeMultiForwardStep (
	size_t PValueBorder,
	const int selectionCriterium,
	int *startIndex,
        TBitset * exclusivedSNP,
        TBitset * goodSNPs
) {
	const size_t
		idvs = data_->getIdvNo(),
		snps = data_->getSnpNo();
//if(score) 
	//getscores this means that we have logistic regression 
       // data_->getOrderedSNP should be something different for score
	int dummy = 0;
if(NULL==startIndex)
  {
   startIndex=&dummy;
  }
  int returnIndex=*startIndex;
	logger->info( "Start Multiple-Forward-Step" );
  if (0==PValueBorder)
  {
		PValueBorder = snps;
		logger->info( "Default setting for PValuePorder: select all" );
  }
  //if modelsize is bigger than0 and selectionCriterium ~! 1
  //then new modus
int startSize=getModelSize();
int newSNP=0;

	logger->debug( "Model start size %u before MultiForwardStep", startSize );
	if ( Parameter::selectionCriterium_BIC != selectionCriterium && 0 < startSize) {
		logger->info( "new usage of Fastforward" );
		newSNP = parameter->ms_MaximalSNPsMultiForwardStep;
	} else {
		logger->info( "normal FastForward" );
		newSNP = parameter->ms_forward_step_max;
	}

  double oldBIC,newBIC;
  bool orig_affection_status_phenotype = true; //init just to init...

  // Fast multiple Forward Step for affection phenotypes: linear regression is used 
	// TODO<BB>: again (similar code above) remove the ugly parameter setting in code
	if ( parameter->ms_FastMultipleForwardStep ) {
		logger->info( "Fast Multiple-Forward Used." );
		orig_affection_status_phenotype = parameter->affection_status_phenotype;
		parameter->affection_status_phenotype = false;
		upToDateBetas_ = false;		// betas computed by a different regression
	}

	oldBIC = computeMSC(selectionCriterium);
	if ( ::isinf( oldBIC ) && oldBIC < 0.0 ) {
   		logger->info( "model fully explains observations" );
		return false;
	}
	if ( parameter->affection_status_phenotype ) {
		Model NewModel( *data_ ); //new Model to test SNPs
	
    for (
      int i = *startIndex ;
      getModelSize() </*=*/  startSize+/*parameter.ms_MaximalSNPsMultiForwardStep*/
      newSNP && i < PValueBorder /*data_->getSnpNo()*/;
      ++i
    ) {
//DEBUG cout<<endl<<"getModelSize() =  startSize+ newSNP"<<getModelSize()<<"="
//<< startSize+/*parameter.ms_MaximalSNPsMultiForwardStep*/ newSNP;

	if (goodSNPs != 0) {   //  GA
		if ( false == (*goodSNPs)[ data_->getOrderedSNP( i ) ] ) {
			continue;
		}
	}
      if (exclusivedSNP != 0) // GA
      {
	if ( true == (*exclusivedSNP)[ data_->getOrderedSNP(i)] )
        {
          continue;          
        }     
      }    
      // progress checking
		if ( 0 == i % 20 ) {
			logger->info(
				"Done %3.5f%%...",
				i / ( snps / 100.0 )
			);
		}

      NewModel = *this;
        // To quickly search whether SNP is already "in"
    const ModelIndex modelIndex( modelSnps_ ); 
	const size_t orderedSnpIndex = data_->getOrderedSNP( i );
    if ( !modelIndex.contains( orderedSnpIndex ) )     {
	 NewModel.addSNPtoModel( orderedSnpIndex ); // Add SNPs according to p-Value
      
      // check if regression works properly, else try next SNP
         if (  NewModel.computeRegression() ) {
            newBIC = NewModel.computeMSC(selectionCriterium);
         // check if SNP is added to the Model
          if ( newBIC < oldBIC ) {
				logger->debug(
					"relative improvement of new model is %f",
					(oldBIC-newBIC)/oldBIC
				);
          *this = NewModel; //model is now updated 
				if ( parameter->detailed_selction ) {
						logger->debug(
							"ADD SNP: %s ModelSize: %u MSC: %f",
							data_->getSNP( orderedSnpIndex ).getSnpId().c_str(),
							NewModel.getModelSize(),
							newBIC
						);
				}
          oldBIC = newBIC;
          if (exclusivedSNP != 0)    // GA
		(*exclusivedSNP)[ data_->getOrderedSNP( i ) ] = true;  
          }
	  
       }
	}//if snp is  in model	 
     returnIndex=i; 
     }
   } else { // linear models
	if ( idvs <= getModelSize() ) {
			logger->debug( "model already at maximum possible size before linear dependence" );
		return false;
	}
     // Compare: oraculateOptimalForwardStep() and TODO<BB>: redesign to reduce code duplication
     // Fast incremental linear regression calculator
     QRuncher qruncher( yVec );
     
     // Import the coefficient matrix into the QRuncher
     for ( size_t col = 0; col < xMat.countColumns(); ++col ) {
		qruncher.pushColumn( xMat.columnVector( col ) );
    }
    
    // To quickly search whether SNP is already "in"
    const ModelIndex modelIndex = getIndex();
    
    // Ignore first (fixed) columns not corresponding to SNPs
    // TODO<BB>: Handle this differently once xMat will be pre-transformed by the fixed columns' Householder vectors
    oldBIC = computeMSC(selectionCriterium , qruncher.calculateRSS() );//hopefully 1 is also BiC here
		logger->debug( "oldBIC=%f", oldBIC );
    for (
	size_t snpCol = *startIndex;
      getModelSize() <=  startSize + /*parameter.ms_MaximalSNPsMultiForwardStep*/ newSNP && snpCol < snps;
      ++snpCol
    ) {
			logger->info( "Done %3.5f%%...", snpCol / ( snps / 100.0 ) );
      
		const size_t orderedSnpIndex = data_->getOrderedSNP( snpCol );
      
      if (goodSNPs != 0)
      {
        if (orderedSnpIndex < 0 || orderedSnpIndex > goodSNPs->size())
        {
					logger->error(
						"orderedSnpIndex %d out of range",
						orderedSnpIndex
					);
					// TODO<BB>: avoid exit
          exit(-1);
        }
        if ((*goodSNPs)[orderedSnpIndex] == false)
          continue;
      }  

	if (exclusivedSNP != 0)
	{
		if ( true == (*exclusivedSNP)[orderedSnpIndex] )	// for GA
		{
			continue;          
		}
	}
      
      if ( modelIndex.contains( orderedSnpIndex ) ) continue;
      
      // Prepare new column
	AutoVector xVec( idvs );
	data_->getXcolumn( orderedSnpIndex, xVec );
      qruncher.pushColumn( xVec );
      
      newBIC = computeMSC( selectionCriterium, qruncher.calculateRSS() ); //1 instead of selection criterion
      // TODO<BB>: Here we should store the RSS in a ResultStore to avoid duplicate calculation
      
      if ( newBIC < oldBIC ) {
				logger->debug( "newBIC=%f", newBIC );
        addSNPtoModel( orderedSnpIndex );
      
	if (exclusivedSNP != 0)  // GA
	{
		(*exclusivedSNP)[orderedSnpIndex] = true;   
	}  
        oldBIC = newBIC;
	if ( ::isinf( oldBIC ) && oldBIC < 0.0 ) {
					logger->info( "model fully explains observations" );
		// a miracle has happened
		break;
	}
      } else {
        qruncher.popColumn();
      }

      returnIndex=snpCol;
    }
   }

  // Fast multiple Forward Step:  reset original parameter
  // reset parameter.affection_status_phenotype and compute 
	if ( parameter->ms_FastMultipleForwardStep ) {
		printModelNew();
		// TODO<BB>: It is unelegant and not thread-safe to manipulate the configuration object like this
		parameter->affection_status_phenotype = orig_affection_status_phenotype;
		if (! computeRegression() ) {
			logger->error( "Fast Multiple-Forward: failed" );
			// TODO<BB>: again remove exit
			exit(12);
		}
  }
   *startIndex=returnIndex++;
	if ( *startIndex >= snps ) {
		*startIndex = 0;
	}
	logger->info( "startIndex=", *startIndex );
	logger->info( "Finish Multiple-Forward-Step: %u SNPs have been added", getModelSize() );
	return true;
}

/** In such a step, one SNP is removed from the model.
* From the models with one SNP less, 
* the model with the lowest MJC is selected.
*/
int Model::makeBackwardStep ( Model &smallerModel ) {

	double 	compareMJC = DBL_MAX; // arbitrary large number
	Model NMin( *data_ );
	int 	removedSNP = -1;

	const bool useOracle = !parameter->affection_status_phenotype;	// linear model
	size_t optimalSNP;
	if ( useOracle ) {
		oraculateOptimalLinearBackwardStep( &optimalSNP );
	}
	// if oracle is used, for-loop degenerates to a single run-through with the optimal SNP
	for (
		int i = useOracle ? optimalSNP : 0;
		( !useOracle || optimalSNP == i ) && i < getModelSize();
		++i
	) {
		Model NTest( *data_ );	// construct new model on same MData
		NTest = *this; // copy current model
		NTest.removeSNPfromModel(i);
		
		// check if Regression works, else try next SNP
		if ( NTest.computeRegression() ) {
			//~ if (NTest.getMJC() < compareMJC)
			if (NTest.computeMSC() < compareMJC)
			{//here ERICH
				NMin = NTest;
				compareMJC = NTest.computeMSC();
				removedSNP = i;
			}
		}
	}
	
	// check if a valid model is returned, (0-model stays 0-model)

	if ( 0 < getModelSize() ) {
		 smallerModel = NMin;
	}
	return removedSNP;
}

/** In such a step, one SNP is added to the model.
* From the models with one SNP more, 
* the model with the lowest MJC is selected.
*/
int Model::makeForwardStep ( Model &biggerModel, const int boundSNP, const int selectionCriterium ) {
	double compareMSC =   getMSC(); // DBL_MAX;	// arbitrarily large number
	Model NMax( *data_ );
	int addedSNP=-1;
	// TODO<BB>: Avoid expensive copying of Model objects.
	// oracle bug
		const bool useOracle=false;

	size_t optimalSNP;
	if ( useOracle ) {
		oraculateOptimalLinearForwardStep( &optimalSNP, boundSNP );
	}
	// if oracle is used, for-loop degenerates to a single run-through with the optimal SNP
	for (
		int i = useOracle ? optimalSNP : 0;
		( !useOracle || optimalSNP == i ) && i < boundSNP;
		++i
	) {
		Model NTest( *data_ );	// construct new model on same MData
		NTest = *this;	// copy current model
		NTest.addSNPtoModel( data_->getOrderedSNP( i ) );
		
		// check if Regression works, else try next SNP
		if ( NTest.computeRegression() ) {
			if ( NTest.computeMSC( selectionCriterium ) < compareMSC ) {
				NMax = NTest;
				compareMSC = NTest.computeMSC( selectionCriterium );
				addedSNP = data_->getOrderedSNP( i );
			}
		}
	}

	// a bigger, meaningfull model is found
	if ( -1 != addedSNP ) { 
		biggerModel = NMax;
	}
	return addedSNP;
}


/** here a saveguarded backward step when 2 consecutive backwardsteps dosent improve the model stop the procedure this will save time by backward of big probems
  */
bool Model::saveguardbackwardstep (
	Model &smallerModel,
	const int selectionCriterium
) {
	Model backwardModel( *data_ );
	Model model3( *data_ );
model3=*this;

//this is for giving back the original model instead of the smallest one
double bestMSC=getMSC(); //from
	logger->info( "saveguardbackwardstep bestMSC=%f", bestMSC );
int breakfor=0;
bool improvment=false;
 // compute steps
 // improvment=false;
 for (int ii=getModelSize();ii>0;ii--)//0 has MSC 0 but this is for me easier
	 //than to implement the 0 case
	 //with an explicit 0 model
 { //reset
 	
		const int removedSNP = makeBackwardStepED( backwardModel, selectionCriterium );
 	*this=backwardModel;//wird ja auch in den Schritten verkleinert
 	//die nichts bringen werden
 double	locMSC=backwardModel.getMSC();
		logger->debug(
			"Modelsize=%u SNP(%u)=[] MSC=%f",
			backwardModel.getModelSize(),
			removedSNP+1,
			locMSC
		);
 
// if ther is an improvment then update the best model 	
 if( locMSC<=bestMSC)
 {    breakfor=0;
     	improvment=true;
         bestMSC=locMSC;
 			logger->debug( "better model found" );
  	model3 = backwardModel;
         *this=model3; 
// model3.printModel();
 
 }
 //else set the counter up and copy the backwardModel to *this
 else
 {//what if the criterion is positive?
			if ( bestMSC > 0 ) {
				// do nothing this could happen when you change from mBIC with 45 to mBIC2/,
				// and should only happen, when you start modelselection again with
				// a stronger criterion
	 		} else {
		 		++breakfor;
				if ( breakfor >= parameter->saveguardsteps ) {	//here one could set any value %the number here was 2
	 				break;
				}
			}
		}
	}
//after the fullbackward
      smallerModel=model3;
	if ( improvment ) {
		smallerModel.printModel(
			"best Model after FullBackward",
			selectionCriterium
		);
	}
*this=smallerModel;//restauration of original model;
		 	 return improvment;
}

/** In such a step, one SNP is removed from the model.
* From the models with one SNP less, 
* the model with the lowest MJC is selected.
*/
int Model::makeBackwardStepED ( Model &smallerModel, const int selectionCriterium ) {
	double 	compareMSC = DBL_MAX,returnVal=DBL_MAX; // arbitrary large number
	double   local=0;
	Model NMin( *data_ );
	int 	removedSNP = -1;
        //vector<double> MSC(getModelSize(),-1);//-1 is not a valid size
	//double MSC;;
	//const bool useOracle = !parameter.affection_status_phenotype;	// linear model
	const bool useOracle = false; //for an error
	size_t optimalSNP = 2000000000UL;
	
	if ( useOracle ) {
		returnVal =oraculateOptimalLinearBackwardStep( &optimalSNP );
		if ( DBL_MAX == returnVal ) {
			logger->debug( "oraculateOptimalLinearBackwardStep did not suggest anything" );
		}
	        if ( 2000000000UL == optimalSNP ) {
			logger->debug(
				"Regression failed by searching the optimal SNP, so optimalSNP=%u",
				optimalSNP
			);
		}
	}
	// if oracle is used, for-loop degenerates to a single run-through with the optimal SNP
	for (
		int i = useOracle ? optimalSNP : 0;
		( !useOracle || optimalSNP == i ) && i < getModelSize();
		++i
	) {
		Model NTest( *data_ );	// construct new model on same MData
		NTest = *this; // copy current model
		NTest.removeSNPfromModel(i);
//	  cerr<<"removed SNP="<<i<<endl;	
		// check if Regression works, else try next SNP
		if ( NTest.computeRegression() ) {
				local=NTest.computeMSC( selectionCriterium );
			if (local < compareMSC)
		{ 
				NMin = NTest;
				compareMSC= local;
				removedSNP = i;
			}
			//it cannot increase MSC globaly because ist starts with $\infty$
		} else {
			logger->error( "Regression failed by removing SNP %u", i );
		}
	}
	
	// check if a valid model is returned, (0-model stays 0-model)

	if ( 0 < getModelSize() ) {
		 smallerModel = NMin;
	}
	return removedSNP;
}



bool Model::makeMultiBackwardStep () {
	logger->info( "Start Multiple-Backward-Step" );
	
	Model BackwardModel( *data_ );
	int		removedSNP;// the removed SNP by a Backward-Step: !!! the integer is the position in the Model (rel) before the backward step !!!
	int		oldModelsize = getModelSize();
	bool	stop=false;
	// checks if we get lower MSCs 
	
	// from a given model, backward steps are performed until a further backwardstep does not lower the MSC
	// while there are SNPs left
	while ( !stop && 0 < getModelSize() ) {
		removedSNP = makeBackwardStepED(BackwardModel);	//compute Model with the lowest MSC with one SNP less
		if ( BackwardModel.computeMSC() <= computeMSC() ) // check if new Model has a lower MSC,
		{
			if ( parameter->detailed_selction ) {
				logger->debug(
					"REMOVE	SNP: %s ModelSize: %u MSC: %f",
					getSNPId( removedSNP ).c_str(),
					BackwardModel.getModelSize(),
					BackwardModel.computeMSC()
				);
			}

			// TODO<BB>: Avoid expensive copying of Model objects.
			*this = BackwardModel;			// the Model is updated
		}
		else
		{
			stop = true;
		}
	}

	logger->info(
		"Finished Multiple-Backward-Step: %u SNPs removed",
		oldModelsize-getModelSize()
	);
	return true;
}

bool Model::computeRegression () {
	if ( parameter->affection_status_phenotype ) {
		return computeLogRegression();
	} else {
		return computeLinRegression();
	}
}


bool Model::computeLinRegression () {
	const size_t cols = xMat.countColumns();
	QRuncher qruncher( yVec );
	bool solvable = true;
	for ( size_t col = 0; col < cols; ++col ) {
		solvable &= 0.0 != qruncher.pushColumn( xMat.columnVector( col ) );
	}
	
	if ( !solvable ) {
		upToDateBetas_ = false;
		return false;
	}

	qruncher.calculateCoefficients( beta );
	upToDateBetas_ = true;
	modelJudgingCriterion_ = qruncher.calculateRSS();
	return true;
}

void Model::scoreTest ( const string& extra ) {
	const size_t
		snps = data_->getSnpNo(),
		vars = getNoOfVariables();
	SortVec score;
	ScoreTestShortcut stsc( *data_ );
	stsc.ScoreTestShortcut::scoreTests( *this, 0u, snps, score );
	stringstream ss;
	ss << vars;
	ofstream S;
	S.exceptions( ofstream::eofbit | ofstream::failbit | ofstream::badbit );
	try{
		S.open( ( parameter->out_file_name + "_" + extra + ss.str() + ".score" ).c_str(), fstream::out );
		const size_t nscores = score.size();
		for ( size_t i = 0; i < nscores; ++i ) {
			S << score.getId( i ) << " " << score.getValue( i ) << endl;
		}
		S.close();
	}
	catch ( ofstream::failure e ) {
		logger->error( "Could not write score-file: %s", e.what() );
	}
}

void Model::scoreTestWithOneSNPless ( const size_t position, SortVec &score ) {
//foreach model SNP: remove them (make a  regression) and make scoretest
//take the result snps and add these snps, when the improve the model 
//then take the best of these snp in the new model.
	Model intermediateModel( *data_ );
	intermediateModel=*this;
	//SortVec score(size); should be created outside

    	intermediateModel.removeSNPfromModel( position );
	intermediateModel.computeLogRegression ();
        //DEBUG intermediateModel.printModel();
	ScoreTestShortcut stsc( *data_ );
	//this should only check some
	//snp around the removed snp
        //big windows make the prozess very very slow
	//HARDCODED for my computer this varies for different computer and probably
	//for different size of the XMat_ 
	size_t fenster = 25;
	if (35<getModelSize ())  //when you assume the models have to be in this size or when the computer is fast enough set it to somewhat higher
		fenster=10;
	if (45<getModelSize ())
		fenster=5;
	if (55<getModelSize ())
		fenster=1;
	const size_t
		snps = data_->getSnpNo(),
		start = modelSnps_[position] > fenster ? modelSnps_[position] - fenster : 0u,
		stop = min( modelSnps_[position] + fenster, snps );
	stsc.scoreTests( intermediateModel, start, stop - start, score );
}

bool Model::computeLogRegression () {
	AutoVector beta_array( beta );
	double 	loglik = 0.0;		        // the log-likelihood

	// logistic firth regression 
	if ( 
		logistffit(
			// Input
			xMat,
			yVec,
			// In+Output
			beta_array,
			// Output
			&loglik, 
			// Control Parameters
			parameter->logrC_maxit,
			parameter->logrC_maxhs,
			parameter->logrC_maxstep,
			parameter->logrC_lconv,
			parameter->logrC_gconv,
			parameter->logrC_xconv
		)
	) {
		modelJudgingCriterion_ = loglik;	//will be calculated in firth-fit
		beta.copy( beta_array );
		upToDateBetas_ = true;
		return true;	// logistic regression worked, betas are updated
	} else {
		logger->info( "Logistic regression cancelled: Fisher matrix not invertible" );
		return false;	// logistic regression failed
	}
}

double Model::computeSingleRegressorTest () {
	assert( 1 == getModelSize() );
	if ( parameter->affection_status_phenotype ) {
		return computeSingleLogRegressorTest();
	} else {
		return computeSingleLinRegressorTest();
	}
}

/** @see Model::computeLinRegression() for more detailed comments */
double Model::computeSingleLinRegressorTest () {
	// REMARK<BB>: The following code is similar to that in Model::computeLinRegression().
	// TODO: Unify both methods by splitting LinearModel and LogisticModel
	// and then permanently maintaining a QRuncher instead of XMat in the former.
	const size_t cols = xMat.countColumns();
	QRuncher qruncher( yVec );
	double absRii;	// |R_{col,col}| of current QR decomposition
	bool solvable = true;
	for ( size_t col = 0; col < cols; ++col ) {
		absRii = qruncher.pushColumn( xMat.columnVector( col ) );
		solvable &= 0.0 != absRii;
	}
	
	if ( !solvable ) {
		upToDateBetas_ = false;
		return 1.0;
		// if diagonal element 0: then all SNPs equal because zero variance in the cov matrix
		// and the MAF-requirement >0.05 is not satisfied
	}

	qruncher.calculateCoefficients( beta );
	upToDateBetas_ = true;
	modelJudgingCriterion_ = qruncher.calculateRSS();
	
	// compute test statistic for the regression coefficient of the SNP

	const double diff = data_->getIdvNo() - static_cast<double>( cols );
	if ( 0.0 >= diff ) {
		return 1.0;
	}
	const double div = sqrt( modelJudgingCriterion_ / diff ) / absRii;
	if ( 0 == div ) {
		return 0.0;
		// if RSS/MJC = 0, no error, extremely good fit, low P,
	}
	const double beta_i = beta.get( cols - 1 );
	const double test_stat = beta_i / div;
					
	// compute the p-value for the test-statistic. 	
	// to check: http://home.ubalt.edu/ntsbarsh/Business-stat/otherapplets/pvalues.htm#rtdist
	return 2 * gsl_cdf_tdist_P( -fabs(test_stat), diff );
	// TODO<BB>: it would be more efficient not to calculate all P values for cut-off,
	// but rather sort by critical value and use the inverse of the distribution function to cut off.
}

double Model::computeSingleLogRegressorTest () {
	const size_t snp = modelSnps_.at( 0 );	// mind the precondition model has size 1

	// compute the genotype frequency table of SNP snp
	GenotypeFreq freqOneSNP( *data_, snp );
	
	// choose test for single marker test for case-control data 
	switch ( parameter->singleMarkerTest ) {
		case Parameter::singleMarkerTest_CHI_SQUARE:
			return freqOneSNP.calculateChiSquare();
		case Parameter::singleMarkerTest_COCHRAN_ARMITAGE:
			return freqOneSNP.calculateCATT();
		default:
			throw Exception(
				"Implementation error:"
				" Unrecognised choice %d for single marker test.",
				parameter->singleMarkerTest
			);
	}
}

double Model::computeMSC ( const int selectionCriterium ) {
	return computeMSC( selectionCriterium, getMJC() );
}

double Model::computeMSC ( const int selectionCriterium, double mjc ) {
	if ( !upToDateBetas_ ) // check if betas_ and therefor also MJC is up-to-date, otherwise update
		// REMARK<BB>: pretty useless, since method parameter mjc is used below
	{
		computeRegression();
	}

	const size_t
		n = data_->getIdvNo(),
		p = parameter->nSNPKriterium,	// data_->getSnpNo(); das ist original
		q = getModelSize();

	// choose the likelihood part depending if the Data is quantitative or affection
	double LRT, d;
	if ( parameter->affection_status_phenotype ) {
		LRT = -2.0*( mjc - data_->getLL0M() );	// LRT = -2 log (likelihood(0-model)/likelihood(model)) the inverse is the right
	}
	else
	{
		LRT = n*log( mjc ) ;  // n * log(RSS)
	}
	
	logger->debug( "LRT=%f", LRT );

	switch ( selectionCriterium ) {
		case Parameter::selectionCriterium_BIC:
			msc = LRT + q * log(n);
			return msc;

		case Parameter::selectionCriterium_EBIC:
			if ( ::isnan( parameter->EBIC_gamma ) ) {
				// How the default value comes about, see Zhao Chen:
				// "A feature selection approach to case-control genome-wide association studies"
				parameter->EBIC_gamma = 1.0 - log( n ) / ( 2 * log( p ) );
			}
			msc = LRT + q * log(n) + 2 * parameter->EBIC_gamma * logFactorial.logChoose( p, q );
			return msc;

		case Parameter::selectionCriterium_mBIC_firstRound:
		case Parameter::selectionCriterium_mBIC:
			d = -2 * log(
				Parameter::selectionCriterium_mBIC_firstRound == selectionCriterium
				? parameter->mBIC_firstRound_expectedCausalSNPs
				: parameter->mBIC_expectedCausalSNPs
			);
			msc = LRT + q * ( log(n) + 2 * log(p) + d );
			return msc;

		case Parameter::selectionCriterium_mBIC2:
			d = -2 * log( 4 );
			msc = LRT + q * ( log(n) + 2 * log(p) + d ) - 2 * log_factorial(q);
			return msc;

		default:
			throw Exception(
				"Unimplemented selection criterium code %d.",
				selectionCriterium
			);
	}
}

double Model::getMSC() const {
	return msc;
}

/** Destructor */
Model::~Model () {
}

void Model::printModelNew() const {
	const size_t
		idvs = data_->getIdvNo(),
		covs = data_->getCovNo();
	ofstream	SNPL;
	
	SNPL.exceptions ( ofstream::eofbit | ofstream::failbit | ofstream::badbit ); // checks if Logfile can be writt
	try
	{
		SNPL.open( ( parameter->out_file_name + "_SNPListNew.txt" ).c_str(), fstream::out );

		AutoVector vec( idvs );

		for ( size_t i = 0; i < getModelSize(); ++i ) {
			const size_t snp = modelSnps_.at( i );
			data_->getXcolumn( snp, vec );
			SNPL << getSNPId( i ) <<" <- c(";
			for ( size_t idv = 0; idv < idvs; ++idv ) {
				SNPL << ( 0 < idv ? "," : "" );
				SNPL << vec.get( idv );
			}
			SNPL<< ")"<< endl;
		}

		for ( size_t i = 0; i < covs; ++i ) {
			SNPL << data_->getCovMatElementName( i ) << " <- c(";
			data_->getCovariateColumn( i, vec );
			for ( size_t idv = 0; idv < idvs; ++idv ) {
				if ( 0 < idv ) SNPL << ",";
				SNPL << vec.get( idv );
			}
			SNPL<< ")"<< endl;	
		}

		SNPL << "Y" <<" <- c(";
		for ( size_t idv = 0; idv < idvs; ++idv ) {
			if ( 0 < idv ) SNPL << ",";
			SNPL << yVec.get( idv );
		}
		SNPL<< ")"<< endl;
		SNPL <<endl;
		
		SNPL << "Intercept" <<" <- c(";
		for ( size_t idv = 0; idv < idvs; ++idv ) {
			if ( 0 < idv ) SNPL << ",";
			SNPL << 1;
		}
		SNPL<< ")"<< endl;
		
		SNPL << "X <- cbind( Intercept";
		for  ( size_t i = 0; i < getModelSize(); ++i ) {
			if ( 0 < i ) SNPL << ",";
			SNPL << getSNPId( i );
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

/**
 * @brief Print model information (snps and msc) on the screen. Just for testing GA
 */
ostream &operator << ( ostream &out, const Model &m ) {
  vector<size_t> v = m.modelSnps_;
  sort(v.begin(), v.end());
  out << "[";
  if (v.size() > 0)
  {
    out << v.at( 0 );
  }  
  for ( size_t i = 1; i < v.size(); ++i ) {
    out << ", " << v.at( i );
  }
  return out << "], msc: " << setprecision(8) << m.msc;
}  

/**
 * @brief gets snp at given position, 
 * @param pos relative position of snp at vector modelSnps_
 * @return snp
 */
size_t Model::getSNPat ( const size_t pos ) const {
  if (pos < modelSnps_.size())
    return modelSnps_[pos];
  else
  {
		logger->error(
			"getSNPat(%u) out of range [0,%u]",
			pos,
			modelSnps_.size()
		);
		// TODO<BB>: avoid exit in library code.
    exit(-1);
  }
}

vector<size_t> Model::getModelSnps() const {
	return modelSnps_;
}

/**
 * @brief Clears previus model and creates new model from given snps
 * @param snps - set of snps 
 * //WARNING I think there is faster way to create new model from given snps. I'm goint to do it. 
 */
void Model::createFromSNPs( const set<size_t>& snps ) { 
  clearModel();
  computeRegression();
  for ( set<size_t>::const_iterator it = snps.begin(); it != snps.end(); ++it ) {
    addSNPtoModel( *it );
  }
}

/**
 * @brief removes SNP form Model, 
 * @param oneSNP is SNP value at vector modelSnps_
 * @return true if SNP removed, and false if an error occors
 * ------------------------------------------------------------------------------
 */
bool Model::removeSNPValFromModel( const size_t oneSNP ) {
  vector<size_t>::iterator it = find(modelSnps_.begin(), modelSnps_.end(), oneSNP);
  return removeSNPfromModel(it - modelSnps_.begin());
}

void Model::clearModel () {
	modelSnps_.clear(); 
	initializeModel();
}

double Model::computeMSCfalseRegression ( const int selectionCriterium ) {
	vector<size_t> removedSnps;
	msc = computeMSCfalseRegression( selectionCriterium, removedSnps );
	return msc;
}
  
double Model::computeMSCfalseRegression (
	const int selectionCriterium,
	vector<size_t> &removedSnps
) {
  if (computeRegression() == false)
  {
    removedSnps.push_back(modelSnps_[modelSnps_.size() - 1]);
    if (removeSNPfromModel(modelSnps_.size() - 1) == false)
	logger->warning( "removeSNPfromModel failed!" );
    //cout << "try computeMSC for model (-1): " << *this << endl;
    //char cc; cout << "Press a key... "; cin >> cc;
		msc = computeMSCfalseRegression( selectionCriterium, removedSnps );
    
    size_t snp = removedSnps.back();
    
    removedSnps.pop_back();
    addSNPtoModel(snp);
    int n = data_->getIdvNo();
    int p = data_->getSnpNo();
    int q = getModelSize() + 1;   // + 1 for additional corelated snp
	double d = -2 * log( 4 );
    // !!!!!!1
    // czy tak obliczamy tylko dla mBIC2, czy dla wsszystkich?
    // czy d też ma bstyć w tym wzorze?
    double pen = q * (log(n) + 2 * log(p) + d ) - 2 * (log_factorial(q));     
    //cout << "msc = " << msc << ", add snp: " << snp << ", pen = " << pen << ", new msc = " << msc << endl;
    //char c; cout << "Press a key... "; cin >> c;
    msc += pen;
    return msc;
  }
  else
  {
		return computeMSC( selectionCriterium );
  }
}
//diese sollen jetzt auch noch das Kriterium manipulieren können manipulieren können
bool Model::selectModel (
	Model &startFromModel,
	size_t PValueBorder,
	const int maxModel,
	const int selectionCriterium
) {
	logger->info( "Score select model" );
	if ( getModelSize() > parameter->maximalModelSize ) {
		return false;
	}

	double best= getMSC();
	//if best is bigger than 0
	//than a smaller model is better or at least the 0 model which is empty
	 Model 	*backwardModel;
	 backwardModel=this;
//ORIGINAL value is now default   int PValueBorder =100,
	int *startIndex;
 	int dummy=0;
	startIndex=&dummy;
	PValueBorder = min( PValueBorder, 400u );	//override to high PValueBorders !!!
 	double bestMSC=getMSC(); //this one is the best up to now!

 	const size_t
		snps = data_->getSnpNo(),
		vars = getNoOfVariables();
	bool improvement = true, improvement2 = true, improvement3 = true;
	computeRegression();
	ScoreTestShortcut stsc( *data_ );
	SortVec score;
	stsc.scoreTests ( *this, 0u, snps, score );
	const size_t nscores = score.size();
	vector<size_t> Score( nscores );
	for ( size_t i = 0; i < nscores; ++i ) {
 		Score.at( i ) = score.getId( i );
	}

	while(
		( improvement = makeForwardStepLogistic( &bestMSC, PValueBorder, startIndex, Score, selectionCriterium ) )
		|| improvement2
		|| improvement3
	) {
		if( getModelSize() > maxModel ) break;
		startFromModel = *this;
		improvement2 = replaceModelSNPSCORE( selectionCriterium );
		improvement3 = saveguardbackwardstep( startFromModel, selectionCriterium );
		*this=startFromModel;
	}
	finalizeModelSelection(
		*backwardModel,
		improvement || improvement2,
		PValueBorder,
		startIndex,
		Score,
		selectionCriterium
	);
	return getMSC() < best;
}


bool Model::operator == ( const Model &m ) const {
  if (modelSnps_.size() != m.modelSnps_.size())
    return false;
  if (modelSnps_.size() == 0)
    return true;

  vector<size_t> v1 = modelSnps_;
  sort(v1.begin(), v1.end());

  vector<size_t> v2 = m.modelSnps_;
  sort(v2.begin(), v2.end());

  unsigned int i = 0;
  while (i < v1.size())
  {
    if (v1[i] != v2[i])
      return false;
    ++i;
  }
  return true;
}
