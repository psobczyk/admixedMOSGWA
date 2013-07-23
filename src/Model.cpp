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

#include "Model.hpp"
#include "GenotypeFreq.hpp"
#include "QRuncher.hpp"
#include "lookup/package.hpp"
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
bool DEBUG=false,DEBUG2=false;
////++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//  class Model
////++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



Model::Model ( const MData & mData ) : XMat_( NULL ), YVec_( NULL ), betas_( NULL ) {
  data_ = &mData;	// set reference to MData
  upToDateXMat_ = false;	// matrix not allocated
}


Model& Model::operator= ( const Model & orig ) {
	if ( &orig != this ) {
		// Destructor Part 
		if ( NULL != XMat_ ) {
			gsl_matrix_free( XMat_ );
			XMat_ = NULL;
		}
		if ( NULL != YVec_ ) {
			gsl_vector_free( YVec_ );
			YVec_ = NULL;
		}
		if ( NULL != betas_ ) {
			gsl_vector_free( betas_ );
			betas_ = NULL;
		}

		// Assigment Part
		data_ = orig.data_; 
		modelSnps_ = orig.modelSnps_; 
		msc = orig.msc;

		if ( NULL != orig.XMat_ ) {
			if(DEBUG==1)
			{cerr<<"Individuen="<<data_->getIdvNo()<<endl
			     <<"getNoOfVariables()="<<getNoOfVariables()<<endl
			     <<"orig.XMat_->size1 vom orig="<< orig.XMat_->size1<<endl
			     <<"orig.XMat_->size2 vom orig="<< orig.XMat_->size2<<endl;
			}
			XMat_ = gsl_matrix_alloc( data_->getIdvNo(), getNoOfVariables() );
			gsl_matrix_memcpy(XMat_, orig.XMat_);
		}
		if ( NULL != orig.YVec_ ) {
			YVec_ = gsl_vector_alloc( data_->getIdvNo() );
			gsl_vector_memcpy(YVec_, orig.YVec_);
		}
		if ( NULL != orig.betas_ ) {
			betas_ = gsl_vector_alloc( getNoOfVariables() );
			gsl_vector_memcpy(betas_, orig.betas_);
		}

		upToDateXMat_ = orig.upToDateXMat_;
		upToDateBetas_ = orig.upToDateBetas_;
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
	return 1 + parameter.covariables + getModelSize();
}

double Model::getMJC () const {
	return modelJudgingCriterion_;
}

double Model::getBeta ( const int i ) const {
	if ( i >= getNoOfVariables() ) {
		throw;
	} else {
		return gsl_vector_get( betas_, i );
	}	
}

string 	Model::getSNPId ( const snp_index_t i ) const {
		if ( 0 > i || getModelSize() <=i ) { cerr<<"The requested SNP of the model "<<i <<" is outside of [0, getModelSize()-1] Modelsize is actually "<<getModelSize()<<endl<<"you will see an SNP_false instead"<<endl;
	//	throw; //simple but bad because this string is only for reporting!
	return "ERROR SNP";
	} else {	
		return data_->getSNP( modelSnps_.at(i) )->getSnpId();
	}
}

void Model::sortSNPsAccordingBetas () {
	vector<size_t> SNP( getNoOfVariables(), 0 );
  SortVec sbetas(getNoOfVariables());
  vector<double> betas(getNoOfVariables(),0);
 for (int i=0;i<getNoOfVariables();i++) //this asumes that getNoOfVariables	
 { betas[i]=getBeta(i);
   SNP[i] = modelSnps_[i];
 }
sbetas.fillVec(getNoOfVariables(), &SNP[0],&betas[0],false);
}

void Model::addSNPtoModel ( const snp_index_t snp ) {

	// TODO<BB>: I do not see: How is adding a SNP duplicate avoided?
	// only when you use addmanySNP 
	modelSnps_.push_back( snp );

	// Model was previously empty and so it has to be initilized or the XMat is not up-to-date
	if ( 1 == getModelSize() or ! upToDateXMat_ ) {
		this->initializeModel();
	}
	else // an up-to-date XMat is expanded by one column 
	{	
		// Allocate Bigger Matrix 
		// TODO<BB>: Avoid too many allocations and copying, e.g. by use of modern linalg::Matrix
		gsl_matrix * NewXMat = gsl_matrix_alloc ( data_->getIdvNo(), getNoOfVariables() );
		gsl_vector * CopyV = gsl_vector_alloc ( data_->getIdvNo() ); //to copy columns
		gsl_vector *   NEWbetas_ =  gsl_vector_alloc( getNoOfVariables() );
		// copy columns of the old matrix;
		for ( int i = 0; i < getNoOfVariables() - 1; ++i ) {
			gsl_matrix_get_col (CopyV,XMat_,i);
			gsl_matrix_set_col (NewXMat,i,CopyV);
		}
	
		// add the new column at the end
		const Vector xVec = const_cast<MData*>( data_ )->getXcolumn( snp );
		for ( int i = 0; i < data_->getIdvNo(); ++i ) {
			gsl_matrix_set( NewXMat, i, getNoOfVariables() - 1, xVec.get( i ) );
		}
                for(int i=0;i<getNoOfVariables()-1;++i) //one less we have now a model with +1 variables
		     gsl_vector_set(NEWbetas_,i,gsl_vector_get(betas_,i));
gsl_vector_set(NEWbetas_, getNoOfVariables() - 1,0);
		// free old gsl-objects
		if ( NULL != XMat_ ) {
			gsl_matrix_free( XMat_ );
			XMat_ = NULL;
		}
		if ( NULL != betas_ ) {
			gsl_vector_free( betas_ );
			betas_ = NULL;
		}

		XMat_ = NewXMat; // use the new matrix
		upToDateXMat_= true;
                betas_ = NEWbetas_;
		upToDateBetas_= false; 

		gsl_vector_free (CopyV);
	}
}

/** replaces replaceSNPinModel 
 * set snp at position  
 * position is the place where from 0 to  getNoOfVariables()
 * that says that even covariates could be replaced!!!!
 *
 * */

bool Model::replaceSNPinModel ( const snp_index_t snp,  const snp_index_t  position ) {
if ( 0 == getModelSize())
{       cerr<<"replacing SNPs is not implemented for empty Models"<<endl;
       	return false;
}
if(0>position||getModelSize()<=position)
 {	cerr<< "replaceSNPinModel position is not in Model";
	return false;}
if(DEBUG2)
{
	cerr<<"XMat_->size1="<< XMat_->size1<<endl
        <<"XMat_->size2="<<XMat_->size2<<endl<<"getNoOfVariables="<<getNoOfVariables()<<endl;
	for (int j=0;j<100;j++)
{cerr<<"XMat"<<j<<" ";
	for ( int i = 1 + parameter.covariables; i < getNoOfVariables(); ++i ) {
	cerr<<gsl_matrix_get(XMat_,j,i)<<" ";
	}
cerr<<endl;
}
}

	const snp_index_t reset = 1 + parameter.covariables + position;		//position 0 is the first 
 //cout<<"reset="<<reset<<endl;
 upToDateXMat_= false;
 modelSnps_[reset-1]= snp;//
// cerr<<"reset="<<reset;
	// add the new column at the end
       // gsl_matrix_set_col(XMat_,position,data_->xMat.columnVector(snp));
		const Vector xVec = const_cast<MData*>( data_ )->getXcolumn( snp );
		for ( int i = 0; i < data_->getIdvNo(); ++i )//statt position reset 
			gsl_matrix_set( XMat_, i, reset, xVec.get( i ) );
			gsl_vector_set(betas_,reset,0); //0 is relativ good for an unknow variable.
		
if(DEBUG2)
{cerr<<"XMat_->size1="<< XMat_->size1<<endl
<<"XMat_->size2="<<XMat_->size2<<endl<<"getNoOfVariables("<<getNoOfVariables()<<endl;

for (int j=0;j<100;j++)
{cerr<<"XMat"<<j<<" ";
	for ( int i = 1 + parameter.covariables; i < getNoOfVariables(); ++i ) {
	cerr<<gsl_matrix_get(XMat_,j,i)<<" ";
	
}	cerr<<endl;
}
}
//for ( int i=0;i<getNoOfVariables();++i)
gsl_vector_set(betas_,reset,0);


upToDateXMat_= true;
upToDateBetas_= false; 
//Y not altered
return true;	
}

/** removes SNP from Model, SNP is relativ position at vector modelSnps_ this is covariate aware */
bool Model::removeSNPfromModel ( const snp_index_t snp ) {

	// check if SNP is a valid postition. 
	if ( 0 > snp || getModelSize() <= snp ) {
			return false;
	} else {
		modelSnps_.erase( modelSnps_.begin() + snp );

		if ( ! upToDateXMat_ ) {
			this->initializeModel();	
		} else {
			// TODO<BB>: Avoid too much allocation and copying
			gsl_matrix * NewXMat = gsl_matrix_alloc ( data_->getIdvNo(), getNoOfVariables() );	// Allocate Smaller Matrix 
			gsl_vector * CopyV = gsl_vector_alloc ( data_->getIdvNo() );	//to Copy columns
                        gsl_vector * NEWbetas_ = gsl_vector_alloc ( getNoOfVariables() );

			// Copy columns of the old Matrix
			// except 1 + parameter.covariables + snp == i
			for ( int i = 0; i <= getNoOfVariables(); ++i ) {	// compare <= because noOfVariables has been decremented
				if ( i < 1 + parameter.covariables + snp ) {
					gsl_matrix_get_col( CopyV, XMat_, i );
					gsl_matrix_set_col( NewXMat, i, CopyV );
					gsl_vector_set(NEWbetas_,i,gsl_vector_get(betas_,i));
				} else if ( i > 1 + parameter.covariables + snp ) {
					gsl_matrix_get_col( CopyV, XMat_, i );
					gsl_matrix_set_col( NewXMat, i-1, CopyV );
					gsl_vector_set(NEWbetas_,i-1,gsl_vector_get(betas_,i));//NEWbetas_(i)=betas(i+1)
				}
			}
				
			// deallocate old gsl-objects
			if ( NULL != XMat_ ) {
				gsl_matrix_free( XMat_ );
				XMat_ = NULL;
			}
			if ( NULL != betas_ ) {
				gsl_vector_free( betas_ );
				betas_ = NULL;
			}

			XMat_ = NewXMat;	// use the new matrix
			upToDateXMat_ = true;
			//betas_ = gsl_vector_alloc( getNoOfVariables() );	// allocate smaller Vector (less coefficients)
                         betas_=NEWbetas_;
			upToDateBetas_ = false;
		
			gsl_vector_free( CopyV );
		}
		return true;
	}
}


void Model::initializeModel () {
	// TODO<BB>: No mem leak any more; but shortcut free & alloc sequence.
	if ( NULL != XMat_ ) {
		gsl_matrix_free( XMat_ );
		XMat_ = NULL;
	} 
	
	XMat_ = gsl_matrix_alloc( data_->getIdvNo(), getNoOfVariables() ); // Allocate Matrix
	if ( NULL != YVec_ ) {
		gsl_vector_free( YVec_ );
		YVec_ = NULL;
	}
	YVec_ = gsl_vector_alloc( data_->getIdvNo() );
	if ( NULL != betas_ ) {
		gsl_vector_free( betas_ );
		betas_ = NULL;
	}
	betas_ = gsl_vector_alloc( getNoOfVariables() );
	 gsl_vector_set_zero (betas_);

	// Set XMat_
	// Columns of XMat_: intercept, covariables, SNP-Data
	Matrix xMat( *XMat_ );
	size_t col = 0;
	xMat.columnVector( col++ ).fill( 1.0 );		// intercept
	for ( size_t cov = 0; cov < parameter.covariables; ++cov ) {
		const Vector covVec = const_cast<MData*>( data_ )->getCovariateColumn( cov );
		Vector xVec = xMat.columnVector( col++ );
		xVec.copy( covVec );
	}
	for ( size_t modelSnp = 0; modelSnp < getModelSize(); ++modelSnp ) {
		const size_t snp = modelSnps_.at( modelSnp );
		const Vector genotypeVec = const_cast<MData*>( data_ )->getXcolumn( snp );
		Vector xVec = xMat.columnVector( col++ );
		xVec.copy( genotypeVec );
	}
	assert( getNoOfVariables() == col );

	// Set YVec_
	Vector yVec = Vector( *YVec_ );
	yVec.copy( const_cast<MData*>( data_ )->getY() );

	upToDateXMat_ = true;	// XMat_ is uptodate 
	upToDateBetas_ = false; // the betas_ are not uptodate, they are just initialiesed
}

//**adds many SNP to the Model

void Model::addManySNP ( vector<snp_index_t> selected ) {
	//tests
	if (
		0 == selected.size()
		||
		data_->getSnpNo () < selected.size()
	) {
		printLOG("addManySNP selected SNP vector is of size 0 or is bigger than the number of SNP");exit(1);
	}
	//add further testing here
	// dubletts in selected
	// using < as comparison 
	
	//instead using sort use SortVec
	sort (selected.begin(), selected.end()); // sorts all
	if ( 0 > selected[0] ) {
		printLOG( "addmanySNP   negative index in selected" ); return;
	}
	
	if ( data_->getSnpNo () < selected.back() ) {
		printLOG( "addmanySNP index in selected is bigger as Number of SNP" );return;
	}
	vector<snp_index_t> finalselect;
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
		//printLOG( int2str( finalselect[i] ) );
		addSNPtoModel(finalselect[i]);
	}
}

	/** set \beta could be used for generating Y
            the \betas should be  gsl_vector* Beta      
 */

void Model::setBeta ( gsl_vector* & beta)
{
 if (0==beta->size)
 {printLOG("setBeta vector Beta is empty maybe an error");return;}
 if (getNoOfVariables()!=beta->size)
{//Achtung der intercept term nicht vergessen!
 printLOG("setBeta vector Beta has not the same size as NoOfVariables "+ int2str(getNoOfVariables()));return;
}
else
   gsl_vector_memcpy (betas_,beta);
   upToDateBetas_= true;
}
/***
 * Ybinary() calls expXbeta
 *
 */
void Model::Ybinary()
{expXbeta();  //this should be the binary case
}


/***
 *Ycontinous() 
 * the beta is set in advance
 * set an ran_gaussian this should be the continous!
 *
 */
void Model::Ycontinous(){
  gsl_vector* TVec  = gsl_vector_alloc ( (data_-> getIdvNo ()));//TVec=zeros(1,getIdvNo)
  int	     p  = (data_-> getIdvNo ()); //p=zeros(1,getIdvNo)
  gsl_blas_dgemv( CblasNoTrans, 1, XMat_, betas_, 0, TVec );//TVec=XMat_*betas
 for (int i=0; i<data_-> getIdvNo ();++i)
  cout << gsl_vector_get(TVec,i) <<" ";
  cout << endl;
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
	try {Y.open( ( parameter.out_file_name + ".yvm" ).c_str(), fstream::out );}
	catch  ( ofstream::failure e/*xception*/ ) {
		cerr << "Could not open additional pheno file for octave" << ( parameter.out_file_name + ".yvm" ).c_str()<< endl;}
// the header line is not that good when to use as  an input to octave
// but heres she is:
Y<<"FID IID ";
for(int j=1;j<=parameter.replications; j++)//0 is not good 
Y<<j<<" ";
Y<<endl; 
//this was the first line of yvm

    // the header line is not that good when to use as  an input to octave
    // but heres she is:
for (int i=0;i<p;++i)
   { FID=data_->getFID(i);
     ID=data_->getID(i);
     Y <<FID << " "<<ID<<" "; //hier werden die Variablen aufgerufen, das geschieht jede Zeile
	 for(int j=0;j<parameter.replications;++j)
	         
		 Y<<(gsl_vector_get(TVec,i)+gsl_ran_gaussian(r,1))<<" ";
            Y<<endl; //after all a newline
   }
 try{Y.close();}
    catch ( ofstream::failure e/*xception*/ ) {
		cerr << "Could not  additional pheno file for octave" << ( parameter.out_file_name + ".yvm" ).c_str()<< endl;}
 gsl_rng_free (r); 
//YVec_(i)=0+gausian(r,1)
 /*
  for(int i=0;i<  p->size;++i)
	  gsl_vector_set (YVec_,i,gsl_vector_get(TVec,i)+gsl_ran_gaussian(r,1));
	   
  gsl_rng_free (r);
  */
 gsl_vector_free( TVec ); 

}


/** expXbeta calculates the Y values used in simulations 
 */
void Model::expXbeta () {
//	printLOG("my first BLAS call!");
  gsl_vector* TVec  = gsl_vector_alloc ( (data_-> getIdvNo ()));
/** 
    p will be the result of:
    $$ p=\frac{e^{\betaX}}{1-\betaX}$$
    or equivalent  $$\frac{1}{1+e^{-\betaX}:u}$$
*/
  gsl_vector* p  = gsl_vector_alloc ( (data_-> getIdvNo ()));
//These functions compute the matrix-vector product and sum y = \alpha op(A) x + \beta y, where op(A) = A, A^T, A^H for TransA = CblasNoTrans, CblasTrans, CblasConjTrans.
//DEBUGCODE 
for(int i=0; i< betas_->size;++i)
		cout<<"betas["<<i<<"]="<<gsl_vector_get( betas_,i)<<endl;

 
//  for(int i=0; i<data_->getIdvNo();++i)
//{for(int j=0; j< getNoOfVariables();++j)
//	cout <<"["<<i<<","<<j<<"]="<<gsl_matrix_get(XMat_,i,j)<<" ";
//cout<<endl;}
/*DEBUGCODE */
 //double general matrix v vector:w multiplication
	gsl_blas_dgemv( CblasNoTrans, -1, XMat_, betas_, 0, TVec ); //-X*beta
/*
8.3.8 Vector operations
-----------------------
 -- Function: int gsl_vector_add (gsl_vector * A, const gsl_vector * B)
     This function adds the elements of vector B to the elements of
     vector A, a'_i = a_i + b_i. The two vectors must have the same
     length.
-- Function: int gsl_vector_div (gsl_vector * A, const gsl_vector * B)
     This function divides the elements of vector A by the elements of
     vector B, a'_i = a_i / b_i. The two vectors must have the same
     length.
-- Function: int gsl_vector_scale (gsl_vector * A, const double X)
     This function multiplies the elements of vector A by the constant
     factor X, a'_i = x a_i.
 -- Function: int gsl_vector_add_constant (gsl_vector * A, const double
          X)
     This function adds the constant value X to the elements of the
     vector A, a'_i = a_i + x.
*/
	double lokalValue=0;
	double mean=0;
for (int i=0;i<TVec->size;++i)	
      { lokalValue=1/(1+exp(gsl_vector_get(TVec,  i)));
        //DEBUG	      cout<< lokalValue<<","<<endl;
	gsl_vector_set(TVec,i,lokalValue);
        mean+=lokalValue;
            //  printLOG( double2str(gsl_vector_get (TVec,  i)));
      }
     mean/=TVec->size;
     //mean=0.5-mean; //this should do the trick, no longer the mean but the difference to the mean
printLOG("the mean should be near to 0.5 mean="+ double2str(mean));
//   for (int i=0;i<TVec->size;++i)
//{gsl_vector_set(TVec,i,gsl_vector_get(TVec,  i)+0);//mean);
//	printLOG( double2str(gsl_vector_get (TVec,  i)));
//}
//gsl_vector_set(betas_,0,mean); //beta reset 
/*gsl_blas_dgemv( CblasNoTrans, -1, XMat_, betas_, 0, TVec );
for (int i=0;i<TVec->size;++i)
{gsl_vector_set(TVec,i,1/(1+exp(gsl_vector_get(TVec,  i))));
		printLOG( double2str(gsl_vector_get (TVec,  i)));
}
*/	//  TVec=0.5-mean(TVec)
//
/** -- Function: int gsl_vector_memcpy (gsl_vector * DEST, const
          gsl_vector * SRC)
*/
// gsl_vector* out  = gsl_vector_alloc ( (data_-> getIdvNo ()));
 // if parameter.replication == 0 set him to 1
 if (0==parameter.replications)
{ parameter.replications=1;
     printLOG("replication has to be set, it is now 1");}
 gsl_matrix* out  = gsl_matrix_alloc ( (data_-> getIdvNo ()),parameter.replications); //generate 1000 replications
/* only for the random number generator */
 const gsl_rng_type * T; 
 gsl_rng * r;//unknown
 gsl_rng_env_setup();
 T = gsl_rng_default;
 //
 r = gsl_rng_alloc (T);

 //long int seed = time (NULL) * getpid();
 gsl_rng_set (r, time (NULL));
/* only for the random number generator */

 for(int i=0;i<  p->size;++i)
	 for(int j=0;j<parameter.replications;++j)
 //gsl_vector_set (out,i,(gsl_ran_flat (r,0.0,1.0 )>gsl_vector_get(TVec,i))?1:0);
  gsl_matrix_set (out,i,j,(gsl_ran_flat (r,0.0,1.0 )>gsl_vector_get(TVec,i))?1:0);
	 gsl_rng_free (r);//    

//printLOG("now the resulting p:");
//*PRINTING when needed
//Now the printing of 
string FID,ID;
ofstream Y;
	Y.exceptions ( ofstream::eofbit | ofstream::failbit | ofstream::badbit );
	try {Y.open( ( parameter.out_file_name + ".yvm" ).c_str(), fstream::out );}
	catch  ( ofstream::failure e/*xception*/ )
	{
		cerr << "Could not open additional pheno file for octave" << ( parameter.out_file_name + ".yvm" ).c_str()<< endl;}
               // the header line is not that good when to use as  an input to octave
               // but heres she is:
               Y<<"FID IID ";
               for(int j=1;j<=parameter.replications; j++)
                   Y<<j<<" ";
               Y<<endl; 
               //this was the first line of yvm
               for (int i=0;i<p->size;++i)
                   { FID=data_->getFID(i);
                     ID=data_->getID(i);	    
	             Y <<FID << " "<<ID<<" "; //hier werden die Variablen aufgerufen, das geschieht jede Zeile
                     for(int j=0;j<parameter.replications; j++)
                         Y<<gsl_matrix_get (out,  i, j)<<" ";
                     Y<<endl;
		   }      
        try{Y.close();}
        catch ( ofstream::failure e/*xception*/ )
	{cerr << "Could not close additional pheno file for octave" << ( parameter.out_file_name + ".txt" ).c_str()<< endl;}

gsl_matrix_get_col (YVec_,out, 0); // the first column
// gsl_vector_memcpy (YVec_,out); //copy out to YVec_ 
 gsl_matrix_free( out ); 
 gsl_vector_free( p );

}
void Model::getYvec ( vector<bool> &out ) {
	if ( NULL == YVec_ ) {
		cerr	<< "YVec_ is NULL I take the Y from MDATA getYvalue"
			<< endl
			<< "and YVec_-1"
			<<endl;
		YVec_ = gsl_vector_alloc( data_->getIdvNo() );	//some initialisation for the YVec_
		Vector yVec = Vector( *YVec_ );
		yVec.copy( const_cast<MData*>( data_ )->getY() );
		// BB: previously it was gsl_vector_set(YVec_, i, data_->getYvalue(i)-1 );}
	}
	if ( out.size() < YVec_->size ) {
		out.resize(YVec_->size);//out
	}
	for ( size_t idv = 0; idv<YVec_->size; ++idv ) {
		out[idv] = gsl_vector_get( YVec_, idv );
	}
	cout << "getYvec";
}

//printYvec 
// check = false print data_->getYvalue(i)
// else 
// check =true print gsl_vector_get(YVec_,i)
void Model::printYvec ( const bool check ) {
	ofstream Y;
	Y.exceptions ( ofstream::eofbit | ofstream::failbit | ofstream::badbit );
	try {
		Y.open( ( parameter.out_file_name + "Yvecout" ).c_str(), fstream::out );
		if ( check ) {
			for ( size_t idv = 0; idv < YVec_->size; ++idv ) {
				Y << gsl_vector_get( YVec_, idv ) << "\n";
			}
		} else {
			const Vector yVec = const_cast<MData*>( data_ )->getY();
			for ( size_t idv = 0; idv < data_->getIdvNo(); ++idv ) {
	   			Y << yVec.get( idv ) <<endl;
			}
		}
		Y.close();
	}
	catch ( ofstream::failure e ) {
		cerr	<< "Could not write Yvecout file"
			<< ( parameter.out_file_name + "Yvecout" ).c_str()
			<< endl;
	}
}

//REM 
void Model::printModel ( const string& out, const string& filemodifier ) {
	stringstream ss; // to save output
	ofstream OUT; // output model to file

	// generate output
	ss << out << endl;
	ss << "ModelSize: " << getModelSize() << "\tMSC: "<< computeMSC() <<"\tMJC: "<< getMJC() << endl;
	ss << "Nr \tSNPId       \tChr \tPos   \tbeta    \tp-Value Single Marker Test" <<endl;

	for ( int i = 0; i < getModelSize(); ++i ) {
		ss << modelSnps_.at(i)<< "\t" 
				<< getSNPId(i) << "\t"
				<< (data_->getSNP(modelSnps_.at(i)))->getChromosome()<<"\t"
				<< setw(10)<<(data_->getSNP(modelSnps_.at(i)))->getBasePairPosition()<<"\t"   //the setw for formating 
				<< getBeta( 1 + parameter.covariables + i ) << "\t"
				<< (data_->getSNP(modelSnps_.at(i)))->getSingleMarkerTest()
				<< endl;
	}
        if(0==getModelSize())
		ss<<"\tIntercept \t \t \t is not available empty model"<<endl;
	else
	ss << "\tIntercept \t \t \t" << getBeta(0)<< endl;

	for ( int i = 0; i < parameter.covariables ; ++i ) {
		ss << "\t" << data_->getCovMatElementName(i) << "\t\t\t"<< getBeta( i + 1 ) << endl;
	}

	// output to screen
	cout << ss.str();
	
	// output to file
	try {
		OUT.open( ( parameter.out_file_name +filemodifier+ ".mod" ).c_str(), ios::out );
		OUT << ss.str();
		OUT.close();
	} catch( ofstream::failure ) {
		printLOG( "Could not write Modelfile \"" + parameter.out_file_name +filemodifier+ ".mod\"." );
	}
}



void Model::printModelInR () const {
	vector<string> output;
	for ( int i=0; i < getModelSize(); ++i ) {
		output.push_back(getSNPId(i));
	}
	data_->printSelectedSNPsInR(output);
}

/*does not print the elements from the model but something from xMat instead*/
void Model::printModelInMatlab (const string& dummy  ) const {
	vector<string> output;
	for ( int i=0; i < getModelSize(); ++i ) {
		output.push_back(getSNPId(i));
	}
//	string dummy("");
	data_->printSelectedSNPsInMatlab(output,dummy);//debug this it messes something up
}

/*bool Model::replaceModel(const int typeNr)
{while(replaceModelSNPbyNearCorrelated(const int typeNr))
 improvment=saveguardbackwardstep( backwardModel);

}
*/
/* the old 1 SNP at a time routine for Correlated */
bool  Model::replaceModelSNPbyNearCorrelated1(const int typeNr)

{  //Model NewModel( *data_ );
    printLOG("replaceModelSNPbyNearCorrelated 1 dim variant");
    Model model0(*data_);
    bool improve=false;
    computeMSC(typeNr);
    cout <<typeNr<<endl;
    double bestMSC=getMSC();
   // cerr<<"bestMSC"<<bestMSC<<endl;//DEBUG
    int fenster=6;
    double  threshold=0.1; //should select every model
    double val =DBL_MAX; 
    vector <snp_index_t> zwischen ;
        for(int i=getModelSize()-1; i >=1; --i)//es schaut so aus als ob nur die letzte Variable also die die neu hinzugekommen ist überprüft werdfen soll
 //int i=0;
 // 
        {               model0=*this;
		if(i==getModelSize()-1) fenster=10;//25 is mybe to large
             printStronglyCorrelatedSnps2(i,threshold,zwischen,fenster);
             for(int j=0;j<zwischen.size();++j)
             {  
                    model0.replaceSNPinModel ( zwischen[j],  i ); //geht oder nicht!
                    val=DBL_MAX;
                    model0.computeRegression(); //regression should be calculated!
                    double val= model0.computeMSC(typeNr);

                if(val<1.0002*bestMSC) //saveguard against rouning errors in logistic regression
                      { double  alt =bestMSC;
			 bestMSC=val;
                       //  cerr<<"!!!!!!!!!!!!!!!!!!!Better Model for "<<i<<" with SNP "<<zwischen[j]<<" neues MSC"<<bestMSC<<" "<<val<<endl;
			 printLOG("!!!!!!!!!!Better Model at position " + int2str(i+1) +" SNP= "
			 +int2str( modelSnps_.at(i)) + " is replaced with " + int2str(zwischen[j])
			 +" oldMSC="+ double2str(alt)+ " newMSC ="+ double2str(bestMSC));
                         model0.printModel("Replaced_inter");
                        *this=model0;   model0.computeMSC(typeNr);
                         improve=true; 
                      }
             }
        }
	return improve;
}
/** replaceModelSNPbyNearCorrelate 
 * type Nr isa the Number of the sel criteria 
 * near is set by fenster
 * and correlated is set by threshold 
 * everything other should be fixed
 */
bool Model::replaceModelSNPbyNearCorrelated(const int typeNr)

{  //Model NewModel( *data_ ); 
    printLOG("replaceModelSNPbyNearCorrelated fenster=10 threshold=0.2");
    bool improve=false;
    Model model0(*data_);
    
    computeMSC(typeNr);
    cout <<typeNr<<endl;
    double bestMSC=getMSC();
   // cerr<<"bestMSC"<<bestMSC<<endl;//DEBUG
    int fenster=6;//smaller is better? 10 original  
    double	threshold=0.1; //should select ever model
    double val =DBL_MAX ;
    vector <snp_index_t> zwischen ;
    //jetzt sollen 2 benachbarte SNP zur gleichen Zeit variert werden
    vector <snp_index_t> zwischen2 ;
       	for(int i=getModelSize()-1; i >=1; i-=2) //model size but index runs fro 0 to modelsize-1
 //int i=0;	
	{  		model0=*this;
	     printStronglyCorrelatedSnps2(i,threshold,zwischen,fenster);
             printStronglyCorrelatedSnps2(i-1,threshold,zwischen2,fenster,true);//all =true every SNP allso the SNP for which the variable is selected
	     for(int j=0;j<zwischen.size();++j)
	     {	
                    model0.replaceSNPinModel ( zwischen[j],  i );
                    for(int k=0;k<zwischen2.size();++k)
		    { model0.replaceSNPinModel ( zwischen2[k],  i-1 ); //der mit dem 1 verkleinerten index soll getestet werden
		    val=DBL_MAX;
		    model0.computeRegression(); //regression should be calculated!
		     val= model0.computeMSC(typeNr);
                       //cerr<<"val="<<val<<" "<< zwischen[j]<<","<<zwischen2[k]<<endl;
	             if(val<1.0001*bestMSC)
	              { bestMSC=val;			 
		         printLOG("!!!!!!!!!!!!!!!!!!!Better Model for "+int2str(i)+" whith SNP "+int2str(zwischen[j])+" and SNP "+int2str( zwischen2[k]) +" new MSC"+double2str(bestMSC)+" "+double2str(val));
	                         model0.printModel("Replaced_inter");
		        *this=model0;
			
	
		       cout<<"MSC"<<getMSC();
                       improve=true;

		      }
		      }
	     } 
	}
         computeMSC(typeNr);
//	 printModelInMatlab("after");
//	printModel("Replaced");
return improve;
}

/**
 printStronglyCorrelatedSnps2 should  use later a  faster implementation of corr
 and should be more flexible than printStronglyCorrelatedSnps 
 */

void Model::printStronglyCorrelatedSnps2 (const int which_snp, const double threshold,vector<snp_index_t> &zwischen, int fenster, bool all) const {
	//fenster=40;
	zwischen.clear();
	//zwischen.reserve(2*fenster+1);
		double	 abscor=0;
	unsigned int nSNP=data_->getSnpNo();
	if (which_snp > getModelSize() || which_snp<0)

		cerr<<"printStronglyCorrelatedSnps2 called for snp index outside of the model"<<endl;
//	for ( int i = 0; i < getModelSize(); ++i )
int i=  which_snp;
	{
		int	 low=max(0,(int)modelSnps_.at(i)-fenster);//unsigned is wrong  !
		int	 high=min(nSNP,modelSnps_.at(i)+fenster);
		int line=0;
                 for (int j =low  ;j<high;j++)
		    {   abscor = fabs( data_->computeCorrelation( modelSnps_.at(i), j ) );
		         if ( abscor  >= threshold )

			 {    if(!all){ if( j!=modelSnps_.at(i))
				 zwischen.push_back(j);}

                                 zwischen.push_back(j);//when not
				// cout<<j<<";";
				++line;
			 }
		    }	 

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
        const int fenster=400; //left and right side of the fenster
        const int reserve=2*fenster+1; // if fenster=1 then there is the left and the right and the middle
	double                  zwischen[2*reserve];  //only here for hdf5
	const char                    *S[reserve]={NULL};
	int 					count=0;			// count strongly correlated SNPs for a given SNP
	double					abscor;			// |Correlation| for two SNPs
         //ED
	 //
	unsigned int nSNP=data_->getSnpNo();//für max
//cerr<<"nSNP"<<nSNP<<endl;
//HDF5 
stringstream name;
hid_t file,fid,dataset,space,dset,/*memtype,*/ filetype;
herr_t status;
hsize_t dim[]={0,2};//dim für Vector der korrelierte
hsize_t di[]={0};
//if (0<getModelSize())
//{	
file=H5Fcreate((parameter.out_file_name + extra + "Corr.h5").c_str(),H5F_ACC_TRUNC,H5P_DEFAULT,H5P_DEFAULT); //standart 
//HDF5
for ( int i = 0; i < getModelSize(); ++i ) {
		count =0;
	
		// search for strongly correlated SNPs
		// PARALLELISIERBAR kann GEWALTIG GRO
		//for (int j = 0; j < data_->getSnpNo(); ++j )
		//ternary operator :make
		//(logical expression)?if_true : if false
	
	//for (int j = 0; j<nSNP;j++)
		size_t low = max( 0, (int)modelSnps_.at(i)-fenster );//unsigned is wrong  !
		size_t high = min( nSNP, modelSnps_.at(i) + fenster );
		int      line=0;
		
		for ( size_t j = low; j < high; ++j ) {	//für 0U Literale füur unsigned int	
				if(false)
				abscor =  data_->computeCorrelation( modelSnps_.at(i), j ) ;
		         	else//sometime one wants the signed version 
                                abscor = fabs( data_->computeCorrelation( modelSnps_.at(i), j ) ); // compute correlation between model SNP and SNP j
			// add if correlation is big enough
			 if ( abscor  >= threshold ) {
			 // if (fabs(abscor)>=threshold) {//for the signed version one is interesstet in the high correlated snps
				StrongCor.insert(pair<double, int>(abscor, j));
				S[line]=(data_->getSNP(j))->getSnpId().c_str();
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
	//ED funktioniert nicht; dataset = H5Dcreate( file, name.str().c_str(), H5T_NATIVE_DOUBLE, fid, H5P_DEFAULT );
        /** write the data to the dataset
	 *  and then close nothing because it works 
	 */
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


			ss 	<< (*it).second  //das 2te von pair<double,int>
				<< "\t" << (data_->getSNP((*it).second))->getSnpId() 
				<< "\t" << (data_->getSNP((*it).second))->getChromosome() 
				<< "\t" << (data_->getSNP((*it).second))->getBasePairPosition() 
				<< "\t" << (*it).first << "\n";
		}
		StrongCor.clear();
	
	}
	
	// output to screen
	cout << ss.str();
	
	// output to file
	try {
        status = H5Fclose(file);//the file will closed!
		printLOG("close "+parameter.out_file_name +"_"+extra +"_Corr.txt with status:"+ int2str(status));
		OUT.open( ( parameter.out_file_name + "_" + extra + "_Corr.txt" ).c_str(), ios::out );
		OUT << ss.str();
		OUT.close();

	} catch( ofstream::failure ) {
		printLOG( "Could not write the correlation file \"" + parameter.out_file_name + extra  + "_Corr.txt\"." );
}
//}//if model is 0;	
}

/**
 * replaceModelSNPbyNearFromCAT 
 * if an SNP s enters the Model in the sequence of the accending p-values of the CAT test,
 * occasional the first variable wich enters improves the model  but there is an variable t situated near to s 
 *  which is even better than the selected variable.
 *  currentPosition could be set to 0 search all or the psoition from the forwardstep
 * */
bool Model::replaceModelSNPbyNearFromCAT( int currrentPosition, int PValueBorder,const int typeNr)
{ //current position is the position returned from fast forward,
  // you have to use this befor ressetting
	//find the position of the SNP which has entered the Model;
	//the last SNP in the model will be checked here!
	if (0==getModelSize())
		return false;
	bool improve=false;
	int fenster=50; //be bold
	const int grace=5000;//war 500
    const   int nSNP=data_->getSnpNo();
	printLOG("replaceModelSNPbyNearFromCAT  currentPosition" + int2str(currrentPosition) +" fenster="+ int2str(fenster)+" grace="+int2str(grace)+ " letzter Model SNP ="+int2str( modelSnps_.at(getModelSize()-1)));
//	printLOG(int2str( modelSnps_.at(getModelSize()-1)));
	double bestMSC=getMSC(); //the msc in the current model is the best one up to now
	Model model0(*data_);
	const unsigned int ref=data_-> getOrderedSNP( currrentPosition );
	//cerr<<"ref"<<ref<<endl;
//	cerr<<"getModelSize(="<<getModelSize()<<endl;
/*random permutation
 */
	const size_t N=getModelSize();
	//DEBUG cerr<<"N="<<getModelSize()<<endl; //DEBUG
	const gsl_rng_type *T;
	gsl_rng *r;
	gsl_permutation *a =gsl_permutation_alloc (N);
	gsl_permutation_init(a);
	//DEBUG gsl_permutation_fprintf(stdout,a," %u"); //DEBUG
	cerr<<endl; //DEBUG
    //gsl_permutation *q =gsl_permutation_alloc (N);
	gsl_rng_env_setup();
 T = gsl_rng_default;
 r = gsl_rng_alloc (T);
// gsl_permutation_init(a);
gsl_ran_shuffle (r, a->data, N, sizeof (size_t));//reference to data!!!

//gsl_permutation_fprintf (stdout, a, " %u");
//	for(int jo=getModelSize()-1; jo >=0; jo--)
for(int jo=0; jo<getModelSize(); jo++)
{ int j=gsl_permutation_get(a,jo);
	//DEBUG cerr<<"a["<<jo<<"]="<<j<<endl;
	for(snp_index_t i=0/*cu/rrrentPosition*/;i<min(max(PValueBorder+grace,1000),nSNP);++i) //saveguard against overrun of nSNP and only 500SNP in large Problems, with 1000 in most cases all relevant SNP will found
	        { //cerr<<(abs(ref-data_-> getOrderedSNP(i))<50)<<",";
		       	if (
				abs( (int)modelSnps_.at(j) - (int) data_->getOrderedSNP(i) )<fenster
				&& data_->getOrderedSNP(i) != (int)modelSnps_.at(j)
			) //the original model should not be considered!

			{
		          model0=*this;
			  model0.replaceSNPinModel (data_-> getOrderedSNP(i) ,  j ); //counting from 0
       	                  double val=DBL_MAX;
                          model0.computeRegression(); //regression should be calculated!
                          val= model0.computeMSC(typeNr);
                    //  printLOG("ModelSNP="+ int2str(j) + "POS="+ int2str(data_-> getOrderedSNP(i))+"val=" +double2str(val));

                  if(bestMSC>0?val<0.9999*bestMSC:val<1.0002*bestMSC)
		       	//saveguard against rouning errors in logistic regression
		       //2 versions for < and > 0 
			// REMARK<BB>: To me something like ( val < bestMSC - 0.001 ) would make more sense.
                      { double  alt =bestMSC;
			 bestMSC=val;
			 printLOG("!!!!!!!!!!!!!!!!!!!Better Model at position " + int2str(j) +" SNP= "
			 +int2str( modelSnps_.at(j)) + " is replaced with " + int2str(data_-> getOrderedSNP(i))
			 +" oldMSC="+ double2str(alt)+ " newMSC ="+ double2str(bestMSC));
                         model0.printModel("Replaced_inter");
                        *this=model0;   model0.computeMSC(typeNr);
                         improve=true; 
		      }
                      }
		}
}//rand order

gsl_permutation_free (a);
gsl_rng_free (r);//
return improve;


}//replaceModelSNPbyNearFromCAT ends
bool Model::replaceModelSNPbyNearFromSCORE( int currrentPosition, int PValueBorder,vector<int> SCORE,const int typeNr)
{ //current position is the position returned from fast forward,
  // you have to use this befor ressetting
	//find the position of the SNP which has entered the Model;
	//the last SNP in the model will be checked here!
	bool improve=false;
	int fenster=50; //be bold
	const int grace=5000; //war 500
        const int nSNP=data_->getSnpNo();
	printLOG("replaceModelSNPbyNearFromSCORE  currentPosition" + int2str(currrentPosition) +" fenster="+ int2str(fenster)+" grace="+int2str(grace)+ " letzter Model SNP ="+int2str( modelSnps_.at(getModelSize()-1)));
	printLOG(int2str( modelSnps_.at(getModelSize()-1)));
	double bestMSC=getMSC(); //the msc in the current model is the best one up to now
	Model model0(*data_);
	const unsigned int ref=data_-> getOrderedSNP( currrentPosition ); 	//cerr<<"ref"<<ref<<endl;
//	model0=*this;
//	model0.printModel("bevor SCORE swaps");

	for(int j=getModelSize()-1; j >=0; j--)
	for(snp_index_t i=0/*currrentPosition*/;i<min(max(PValueBorder+grace,1000),nSNP);++i) //saveguard against overrun of nSNP and only 500SNP in large Problems, with 1000 in most cases all relevant SNP will found
	        {//////ERICH was bedeutet  getOrderedSNP(SCORE[i]))??????
				if (
					abs( (int) modelSnps_.at(j) - (int) data_->getOrderedSNP(SCORE[i]) )<fenster
					&&
					data_->getOrderedSNP(SCORE[i])!=(int)modelSnps_.at(j))
					//alle die im Fenster sind werden 
			{// cerr<<endl<<(data_-> getOrderedSNP(i))<<",";
				model0=*this;
//				cerr<<data_-> getOrderedSNP(SCORE[i])<<endl;
//					cerr<<SCORE[i]<<endl;
			  model0.replaceSNPinModel (data_-> getOrderedSNP(SCORE[i]) ,  j ); //counting from 0
       	                  double val=DBL_MAX;
                          model0.computeRegression(); //regression should be calculated!
                          val= model0.computeMSC(typeNr);

                  if(bestMSC>0?val>0.9999*bestMSC:val<1.0002*bestMSC)
		       	//saveguard against rouning errors in logistic regression
		       //2 versions for < and > 0 
                      { double  alt =bestMSC;
			 bestMSC=val;
			 printLOG("!!!!!!!!!!!!!!!!!!!Better Model at position " + int2str(j) +" SNP= "
			 +int2str( modelSnps_.at(j)) + " is replaced with " + int2str(data_-> getOrderedSNP(SCORE[i]))
			 +" oldMSC="+ double2str(alt)+ " newMSC ="+ double2str(bestMSC));
                         model0.printModel("Replaced_inter");
                        *this=model0;   model0.computeMSC(typeNr);
                         improve=true; 
		      }
                      }
		}
return improve;


}//replaceModelSNPbyNearFromCAT ends
bool Model::replaceModelSNPSCORE( ) //vector<int> SCORE,const int typeNr)
{cerr<<"NEW"<<"NEW"<<"NEW"<<"NEW"<<"NEW"<<"NEW"<<"NEW"<<"NEW"<<endl;
 bool improve=false; //we don't know weather we could improve
 int fenster=50; //be bold
 const int grace=5000;
 const int nSNP=data_->getSnpNo();
 double bestMSC=getMSC(); //the msc in the current model is the best one up to now
 //cerr<<"bestMSC="<<bestMSC<<endl;
 Model model0(*data_);
 SortVec score(200); //Warning 50+50+1 that should be variable
 	for(int j=getModelSize()-1; j >=0; j--)
	{//replace
//		  cerr<<modelSnps_[j]<<endl;
	int nscores=	scoreTestWithOneSNPless(j, score);
//	cerr<<"num_of_scores="<<nscores<<endl;

		for(snp_index_t i=0;i<nscores;++i) 
		{       model0=*this;
//       		cerr<<"SNP="<<score.getId(i);
                        model0.replaceSNPinModel (score.getId(i) ,  j );

	                double val=DBL_MAX;
                          model0.computeRegression(); //regression should be calculated!
                          val= model0.computeMSC(0);  //should be mBIC2
       //                 cerr<<" mBIC="<<val<<endl;
                
			if(bestMSC>0?val>0.9999*bestMSC:val<1.0002*bestMSC)
	                      { double  alt =bestMSC;
				 bestMSC=val;
				 printLOG("!!!!!!!!!!!!!!!!!!!Better Model at Modelposition " + int2str(j) + " SCOREPostition=" +int2str(i) +" SNP= "
				 +int2str( modelSnps_.at(j)) + " is replaced with " + int2str(score.getId(i))
				 +" oldMSC="+ double2str(alt)+ " newMSC ="+ double2str(bestMSC));
	                         model0.printModel("Replaced_inter");
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


double Model::oraculateOptimalLinearForwardStep( snp_index_t *snp, size_t bound ) const {
	// First, intermediate step towards using linalg instead plain GSL
	const Matrix xMat( *XMat_ );
	const Vector yVec( *YVec_ );

	// Fast incremental linear regression calculator
	QRuncher qruncher( yVec );

	// Import the coefficient matrix into the QRuncher
	for ( size_t col = 0; col < xMat.countColumns(); ++col ) {
		// TODO: Take care of return value against adding linearly dependent columns
		qruncher.pushColumn( const_cast<Matrix&>( xMat ).columnVector( col ) );
	}

	// To quickly search whether SNP is already "in"
	const ModelIndex modelIndex = getIndex();

	// Search best to add
	double bestRSS = DBL_MAX;
	snp_index_t bestSNP;
	for ( snp_index_t snpCol = 0; snpCol < bound; ++snpCol ) {
		const snp_index_t orderedSnpIndex = data_-> getOrderedSNP( snpCol );
		if ( modelIndex.contains( orderedSnpIndex ) ) continue;

		// Prepare new column
		const Vector xVec = const_cast<MData*>( data_ )->getXcolumn( orderedSnpIndex );
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
* @see oraculateOptimalLinearForwardStep( size_t, snp_index_t* )
*/
double Model::oraculateOptimalLinearBackwardStep( snp_index_t *snp ) const {
	// First, intermediate step towards using linalg instead plain GSL
	const Matrix xMat( *XMat_ );
	const Vector yVec( *YVec_ );

	// Fast linear regression calculator with Bernhard's backward shortcut algorithm
	QRuncher qruncher( yVec );

	// Import the coefficient matrix into the QRuncher
	for ( size_t col = 0; col < xMat.countColumns(); ++col ) {
		// TODO: Take care of return value against adding linearly dependent columns
		qruncher.pushColumn( const_cast<Matrix&>( xMat ).columnVector( col ) );
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
		cerr<<rss<<";";
		if ( rss < bestRSS ) {
			bestRSS = rss;
			bestCOL = col;
		}
	}
	if ( bestRSS < DBL_MAX ) {
		*snp = modelSnps_.at( bestCOL );
		cerr <<"bestCol"<<bestCOL<<"modelSnps_.at( bestCOL )"<<*snp<<endl;
	}
	return bestRSS;
}
/** finalizeModelSelection()
 *   finalization work call some function and quit
 */
bool Model::finalizeModelSelection(Model &backwardModel,int JJ, bool improvment, int PValueBorder, int *startIndex, vector<int> score)
{
  if (false==improvment)      
  { printModel( "no improvment");
    int dummy=0; *startIndex=dummy;
    if(!parameter.affection_status_phenotype)//quantitative
      { cout<<"not implemented"<<endl;backwardModel=*this;
 improvment=saveguardbackwardstep( backwardModel);//not use makeMFFL in this case 
	} else {
		makeMFFL(min(500, PValueBorder),startIndex,score);
	}
    //don't set to an explicit value because of memory leaks when the number of variables is very small!
    JJ++;
	   //REMOVED FOR SPEED    printModelInMatlab(int2str(JJ)); 
	   backwardModel=*this;
	      improvment= saveguardbackwardstep( backwardModel);
	      //what when improvement? set return stop =false
	      if(improvment==true)
	      {	printLOG("finalizeModelSelection");
		return true  ; //if only 1 snp was added and the backward step was succesfull than this is anerror
		//but if let
		*this=backwardModel;
	      }

	       
               printModel( "final model");
	      //REMOVED FOR SPEED printStronglyCorrelatedSnps( 0.99, int2str(parameter.in_values_int) + "the_result" );
	      // printYvec(true);
	     //REMOVED FOR SPEED if(parameter.affection_status_phenotype)
	     //  scoreTest(); //only if case control structure
	      return true;
  }
else

{ 	printLOG("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!after forward step");
//	cerr<<addedSNP<<endl<<"MSC"<<forwardModel->computeMSC();
	return false;
	}
}// finalizeModelSelection()


bool Model::finalizeModelSelection(Model &backwardModel,int JJ, bool improvment, int PValueBorder, int *startIndex)

{
  if (false==improvment)      
  { printModel( "no improvment");
    int dummy=0; *startIndex=dummy;
    if(!parameter.affection_status_phenotype)//quantitative
      { cout<<"not implemented"<<endl;backwardModel=*this;
 improvment=saveguardbackwardstep( backwardModel);//not use makeMFFL in this case 
	} else {
		//diese Schritte sind eigentlich nur einmal nötig da, man sicher nichts
  //  gewinnt wenn man die ersten paar SNP oft und oft wiederholt
   printLOG("finalizeModelSelection last run with min of  parameter.PValueBorder or parameter.reset ");
	   makeMFFL(min(PValueBorder, parameter.reset),startIndex);
	}

  //  don't set to an explicit value because of memory leaks when the number 
  // of variables is very small!
  //  JJ;
 
	   //REMOVED FOR SPEED    printModelInMatlab(int2str(JJ)); 
	   backwardModel=*this;
	      improvment= saveguardbackwardstep( backwardModel);
	      //what when improvement? set return stop =false
	      if(improvment==true)
	      {	printLOG("improvment after newstart forward step");
		return true  ; //if only 1 snp was added and the backward step was succesfull than this is anerror
		//but if let
		*this=backwardModel;
		printModel( "final model");
	      }

	       
              
	      //REMOVED FOR SPEED printStronglyCorrelatedSnps( 0.99, int2str(parameter.in_values_int) + "the_result" );
	      // printYvec(true);
	     //REMOVED FOR SPEED if(parameter.affection_status_phenotype)
	     //  scoreTest(); //only if case control structure
	      return true;
  }
else

{    // cerr<<improvment<<endl;
	printLOG("after forward step");
//	cerr<<addedSNP<<endl<<"MSC"<<forwardModel->computeMSC();
	return false;
	}
}// finalizeModelSelection()
//makeForwardStepLinear the currentModle ist this
//
bool  Model::makeForwardStepLinear(Model *forwardModel, int JJ, double* bestMSC, int PValueBorder,int *startIndex)
{bool improvment=false; //we don't know from an succes up to now.
	  cerr<<"linear regression______________________";

 Model model3(*data_);
 model3=*this;//when everthing fails this remains as result
 forwardModel=this;
 for(int ii=0;ii<20;ii++) //parameter 
        {	
int	 addedSNP = makeForwardStep( *forwardModel, PValueBorder );

	double locMSC=forwardModel->getMSC();
//	 swapHelper=forwardModel;
	 //*this =swapHelper; //currentModel in the caller gets the local information of swapHelper ?
	   if(locMSC<*bestMSC)
	     {       	
		improvment=true;
                *bestMSC=locMSC;
		model3 =*forwardModel;
	        printLOG("better bigger  Model in Iter"+int2str(ii)+" bestMSC="+double2str(*bestMSC));
		*this=model3;
	     }
	   else
	     break;
	//currentModel=&model3;
        }
 *this =model3; //when something is updated in the loop then we will see it
return improvment;
}
/*make ForwardStep for MData::selectModel for logistic Regression score version
 *  *this is forwardModel
 *
 */
bool  Model::makeForwardStepLogistic(int JJ, double* bestMSC, int PValueBorder,int *startIndex, vector<int> score)
{ bool improvment=false;
  printLOG("par--------------------------------------");  //DEBUG
int dummy= *startIndex;
//cout<<"dummy="<<dummy<<endl;
 if( dummy<parameter.reset)//ERICH wenn man 100 als Grenze nimmt sollte man das hier auch herabsetzen 
	 // 500 waren gut wenn man eine Grenze von 3000 hatte
	 // also docg eher PValueBorder/6
 { //cerr<<"startIndex="<<*startIndex<< "but 0 is used"<<endl;
   dummy=0;*startIndex=dummy;
   makeMFFL( PValueBorder,startIndex,score);JJ++;
   //REMOVED FOR SPEED printModelInMatlab(int2str(JJ));
 }
 else if (dummy>=parameter.reset) //these values are only guesses
	 //da auch 300 für 3000 
	 //daher 30 für 100
 {  dummy=max(0,dummy-parameter.jump_back);*startIndex=dummy; //better 500 (?)
     //    cerr<<"StartIndex= "<<*startIndex<<" is used"<<endl;
    makeMFFL( PValueBorder,startIndex,score);JJ++;
//REMOVED FOR SPEED    printModelInMatlab(int2str(JJ));
 } 
//REMOVED FOR SPEED  scoreTest();//more scoreTests

 double   locMSC=getMSC();
		  if(*bestMSC>locMSC)
		  { 
 		   improvment=true;
		   *bestMSC=locMSC;
		  }
 return improvment;
}
/*make ForwardStep for MData::selectModel for logistic Regression
 *  *this is forwardModel
 *
 */
bool  Model::makeForwardStepLogistic(int JJ, double* bestMSC, int PValueBorder,int *startIndex)
{ bool improvment=false;
  printLOG("pValue ordered--------------------------------------");  //DEBUG
int dummy= *startIndex;
//cout<<"dummy="<<dummy<<endl;
 if( dummy< parameter.reset) //500 
 { printLOG("startIndex="+int2str(*startIndex)+ "but 0 is used");
   dummy=0;*startIndex=dummy;
   makeMFFL( PValueBorder,startIndex);JJ++;
  //DEBUG  printModelInMatlab(int2str(JJ));
 }
 else if (dummy>= parameter.reset) //these values are only guesses

 {  dummy=max(0,dummy-parameter.jump_back);*startIndex=dummy; //min because nothing prevents to be jump_back to be bigger than reset!
         printLOG("StartIndex= "+int2str(*startIndex)+" is used, with jump_back="+int2str(parameter.jump_back));
    makeMFFL( PValueBorder,startIndex);JJ++;
//REMOVED FOR SPEED    printModelInMatlab(int2str(JJ));
 } 
//REMOVED FOR SPEED 
// scoreTest();//more scoreTests

 double   locMSC=getMSC();
		  if(*bestMSC>locMSC)
		  { 
 		   improvment=true;
		   *bestMSC=locMSC;
		  }
 return improvment;
}
/** makeMFFS make Fast Forward local fast forward 
 * */

bool Model::makeMFFS(int PValueBorder, int* startIndex)
{       //if (startIndex==NULL)
	    //startIndex
	int selectionCriterium=0;
         //FAST MultipleForward is needed locally!	
        bool oldValue= parameter.ms_FastMultipleForwardStep;
        parameter.ms_FastMultipleForwardStep=true;
        
       
        const int altSize= getModelSize();
makeMultiForwardStep ( PValueBorder,selectionCriterium, startIndex);

 parameter.ms_FastMultipleForwardStep=oldValue;
 return true;
}
//and now with score test 
bool Model::makeMFFL(int PValueBorder, int* startIndex, vector<int> score)
{       
	int selectionCriterium=0;
         //MultipleForward is needed locally!	
   //     bool oldValue= parameter.ms_FastMultipleForwardStep;
        //parameter.ms_FastMultipleForwardStep=true;
        
       
        const int altSize= getModelSize();
makeMultiForwardStepScore ( PValueBorder,selectionCriterium,startIndex, score  );
cout<<"MFFL startIndex SCORE"<<*startIndex<<"Model Size="<<altSize<<endl;

 //parameter.ms_FastMultipleForwardStep=oldValue;
}


/**  makeMFFL make Fast Forward local but without changing to fast search*/

bool Model::makeMFFL(int PValueBorder, int* startIndex)
{       
	int selectionCriterium=0;
         //MultipleForward is needed locally!	
   //     bool oldValue= parameter.ms_FastMultipleForwardStep;
        //parameter.ms_FastMultipleForwardStep=true;
        
       
        const int altSize= getModelSize();
makeMultiForwardStep ( PValueBorder,selectionCriterium,startIndex  );
cout<<"MFFL startIndex"<<*startIndex<<"Model Size="<<altSize<<endl;

 //parameter.ms_FastMultipleForwardStep=oldValue;
}








//ScoreScoreScoreScoreScoreScoreScoreScoreScoreScoreScoreScoreScore
bool Model::makeMultiForwardStepScore ( int PValueBorder, int selectionCriterium, int* startIndex, vector<int> score)
	{
//if(score) 
	//getscores this means that we have logistic regression 
       // data_->getOrderedSNP should be something different for score
if(NULL==startIndex)
  {
    int dummy=0;	  
   startIndex=&dummy;
  }
  int returnIndex=*startIndex;
  printLOG( "Start Multiple-Forward-Step" );
  if (0==PValueBorder)
  {
   PValueBorder=data_->getSnpNo(); 
   printLOG( "Default setting for PValuePorder: select all");
  }
  //if modelsize is bigger than0 and selectionCriterium ~! 1
  //then new modus
int startSize=getModelSize();
int newSNP=0;

  cout<<"ModelSize before MultiForwardStep="<<startSize<<endl;

  if(1!=selectionCriterium && 0<startSize)
		  {cout<<"new usage of Fastforward"<<endl;
		  newSNP= parameter.ms_MaximalSNPsMultiForwardStep;
		  }
  else
                  { printLOG("normal FastForward");
                  newSNP=parameter.ms_forward_step_max;
		  }

 // }
  double oldBIC,newBIC;
  bool orig_affection_status_phenotype = true; //init just to init...

  // Fast multiple Forward Step for affection phenotypes: linear regression is used 
  if ( parameter.ms_FastMultipleForwardStep ) {
    printLOG("Fast Multiple-Forward Score Used.");
    orig_affection_status_phenotype = parameter.affection_status_phenotype;
    parameter.affection_status_phenotype = false;
    upToDateBetas_= false; // betas computed by a different regression
  }

  if ( parameter.affection_status_phenotype ) {
    Model NewModel( *data_ ); //new Model to test SNPs
    oldBIC = computeMSC(selectionCriterium);
	
    for (
      int i = *startIndex ;
      getModelSize() </*=*/  startSize+/*parameter.ms_MaximalSNPsMultiForwardStep*/
      newSNP && i <  PValueBorder /*data_->getSnpNo()*/;
      ++i
    ) {
      // progress checking
      if (0==i%20 )
      {printf ("\rDone %3.5f%%...", i / (data_->getSnpNo()/100.0));
      fflush (stdout);}
	
      NewModel = *this;
        // To quickly search whether SNP is already "in"
    const ModelIndex modelIndex( modelSnps_ ); 
    const snp_index_t orderedSnpIndex = score[i]; //data_-> getOrderedSNP( i );
    if ( !modelIndex.contains( orderedSnpIndex ) )     {
	 NewModel.addSNPtoModel( orderedSnpIndex ); // Add SNPs according to p-Value
      
      // check if regression works properly, else try next SNP
         if (  NewModel.computeRegression() ) {
            newBIC = NewModel.computeMSC(selectionCriterium);
        
         // check if SNP is added to the Model
          if ( newBIC < oldBIC ) {
          *this = NewModel; //model is now updated 
          if ( parameter.detailed_selction ) {
            printLOG(
            //~ "Step " + int2str(i) + " of " + int2str(data_->getSnpNo()) +": "+
              "ADD    SNP: " + data_->getSNP(orderedSnpIndex )->getSnpId() + " ModelSize: "
	      + int2str( NewModel.getModelSize() ) + " MSC: " + double2str(newBIC)
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
  if ( parameter.ms_FastMultipleForwardStep ) {
    printModelNew();
    // TODO<BB>: It is unelegant and not thread-safe to manipulate the configuration object like this
    parameter.affection_status_phenotype = orig_affection_status_phenotype;
    
    if (! computeRegression() ) {
      printLOG( "Fast Multiple-Forward: failed" );
      exit(12);
    }
    //upToDateBetas_=false;
  }
 }  
   *startIndex=returnIndex++;
   if (*startIndex>=data_->getSnpNo())
	*startIndex=0;
   cout<<"startIndex="<< *startIndex<<endl<<endl;
  printLOG( "Finish Multiple-Forward-Step: " + int2str( getModelSize() ) + " SNPs have been added" );
  return true;
}















/* this is the main forward step in the model creation 
 * */

bool Model::makeMultiForwardStep ( int PValueBorder, int selectionCriterium, int* startIndex, std::set<snp_index_t> *exclusivedSNP) {
//if(score) 
	//getscores this means that we have logistic regression 
       // data_->getOrderedSNP should be something different for score
	int dummy = 0;
if(NULL==startIndex)
  {
   startIndex=&dummy;
  }
  int returnIndex=*startIndex;
  printLOG( "Start Multiple-Forward-Step" );
  if (0==PValueBorder)
  {
   PValueBorder=data_->getSnpNo(); 
   printLOG( "Default setting for PValuePorder: select all");
  }
  //if modelsize is bigger than0 and selectionCriterium ~! 1
  //then new modus
int startSize=getModelSize();
int newSNP=0;

  cout<<"ModelSize before MultiForwardStep="<<startSize<<endl;

  if(1!=selectionCriterium && 0<startSize)
		  {cout<<"new usage of Fastforward"<<endl;
		  newSNP= parameter.ms_MaximalSNPsMultiForwardStep;
		  }
  else
                  { printLOG("normal FastForward");
                  newSNP=parameter.ms_forward_step_max;
		  }

  double oldBIC,newBIC;
  bool orig_affection_status_phenotype = true; //init just to init...

  // Fast multiple Forward Step for affection phenotypes: linear regression is used 
  if ( parameter.ms_FastMultipleForwardStep ) {
    printLOG("Fast Multiple-Forward Used.");
    orig_affection_status_phenotype = parameter.affection_status_phenotype;
    parameter.affection_status_phenotype = false;
    upToDateBetas_= false; // betas computed by a different regression
  }

	oldBIC = computeMSC(selectionCriterium);
	if ( ::isinf( oldBIC ) && oldBIC < 0.0 ) {
   		printLOG( "model fully explains observations" );
		return false;
	}
  if ( parameter.affection_status_phenotype ) {
    Model NewModel( *data_ ); //new Model to test SNPs
	
    for (
      int i = *startIndex ;
      getModelSize() </*=*/  startSize+/*parameter.ms_MaximalSNPsMultiForwardStep*/
      newSNP && i <  PValueBorder /*data_->getSnpNo()*/;
      ++i
    ) {
//DEBUG cout<<endl<<"getModelSize() =  startSize+ newSNP"<<getModelSize()<<"="
//<< startSize+/*parameter.ms_MaximalSNPsMultiForwardStep*/ newSNP;

      if (exclusivedSNP != 0)
      {
        if (exclusivedSNP->find( data_->getOrderedSNP(i)) != exclusivedSNP->end())   // for GA
        {
          continue;          
        }     
      }    
      // progress checking
      if (0==i%20) 
      {printf ("\rDone %3.5f%%...", i / (data_->getSnpNo()/100.0));
      fflush (stdout);}
	
      NewModel = *this;
        // To quickly search whether SNP is already "in"
    const ModelIndex modelIndex( modelSnps_ ); 
    const snp_index_t orderedSnpIndex = data_-> getOrderedSNP( i );
    if ( !modelIndex.contains( orderedSnpIndex ) )     {
	 NewModel.addSNPtoModel( orderedSnpIndex ); // Add SNPs according to p-Value
      
      // check if regression works properly, else try next SNP
         if (  NewModel.computeRegression() ) {
            newBIC = NewModel.computeMSC(selectionCriterium);
        
         // check if SNP is added to the Model
          if ( newBIC < oldBIC ) {
		  printLOG("improvement of new model is :" +double2str((oldBIC-newBIC)/oldBIC));

		   

          *this = NewModel; //model is now updated 
          if ( parameter.detailed_selction ) {
            printLOG(
            //~ "Step " + int2str(i) + " of " + int2str(data_->getSnpNo()) +": "+
              "ADD    SNP: " + data_->getSNP(orderedSnpIndex )->getSnpId() + " ModelSize: "
	      + int2str( NewModel.getModelSize() ) + " MSC: " + double2str(newBIC)
            );
          }
          oldBIC = newBIC;
          if (exclusivedSNP != 0)
            exclusivedSNP->insert(orderedSnpIndex);   // for GA
          }
	  
       }
	}//if snp is  in model	 
     returnIndex=i; 
     }
   } else { // linear models
	if ( data_->getIdvNo() <= getModelSize() ) {
   		printLOG( "model already at maximum possible size before linear dependence" );
		return false;
	}
     // Compare: oraculateOptimalForwardStep() and TODO<BB>: redesign to reduce code duplication
     
     // First, intermediate step towards using linalg instead plain GSL
     const Matrix xMat( *XMat_ );
     const Vector yVec( *YVec_ );
     
     // Fast incremental linear regression calculator
     QRuncher qruncher( yVec );
     
     // Import the coefficient matrix into the QRuncher
     for ( size_t col = 0; col < xMat.countColumns(); ++col ) {
	qruncher.pushColumn( const_cast<Matrix&>( xMat ).columnVector( col ) );
    }
    
    // To quickly search whether SNP is already "in"
    const ModelIndex modelIndex = getIndex();
    
    // Ignore first (fixed) columns not corresponding to SNPs
    // TODO<BB>: Handle this differently once xMat will be pre-transformed by the fixed columns' Householder vectors
    oldBIC = computeMSC(selectionCriterium , qruncher.calculateRSS() );//hopefully 1 is also BiC here
   printLOG( "oldBIC________________"+double2str(oldBIC));
    for (
      snp_index_t snpCol = *startIndex;
      getModelSize() <=  startSize + /*parameter.ms_MaximalSNPsMultiForwardStep*/ newSNP && snpCol < data_->getSnpNo();
      ++snpCol
    ) {
      printf ("\rDone %3.5f%%...", snpCol / (data_->getSnpNo()/100.0));
      
      const snp_index_t orderedSnpIndex = data_-> getOrderedSNP( snpCol );
      
      if (exclusivedSNP != 0)
      {
        if ((exclusivedSNP->find( orderedSnpIndex)!= exclusivedSNP->end()))   // for GA
        {
          continue;          
        }     
      }
      
      if ( modelIndex.contains( orderedSnpIndex ) ) continue;
      
      // Prepare new column
	const Vector xVec = const_cast<MData*>( data_ )->getXcolumn( orderedSnpIndex );
      qruncher.pushColumn( xVec );
      
      newBIC = computeMSC( selectionCriterium, qruncher.calculateRSS() ); //1 instead of selection criterion
      // TODO<BB>: Here we should store the RSS in a ResultStore to avoid duplicate calculation
      
      if ( newBIC < oldBIC ) {
 printLOG("newBIC________________"+double2str(newBIC));
        addSNPtoModel( orderedSnpIndex );
      
        if (exclusivedSNP != 0)
        {
          exclusivedSNP->insert(orderedSnpIndex);   // for GA        
        }  
        oldBIC = newBIC;
	if ( ::isinf( oldBIC ) && oldBIC < 0.0 ) {
   		printLOG( "model fully explains observations" );
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
  if ( parameter.ms_FastMultipleForwardStep ) {
    printModelNew();
    // TODO<BB>: It is unelegant and not thread-safe to manipulate the configuration object like this
    parameter.affection_status_phenotype = orig_affection_status_phenotype;
    
    if (! computeRegression() ) {
      printLOG( "Fast Multiple-Forward: failed" );
      exit(12);
    }
    //upToDateBetas_=false;
  }
   *startIndex=returnIndex++;
   if (*startIndex>=data_->getSnpNo())
	*startIndex=0;
   printLOG("startIndex="+ int2str(*startIndex));;
  printLOG( "Finish Multiple-Forward-Step: " + int2str( getModelSize() ) + " SNPs have been added" );
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

	const bool useOracle = !parameter.affection_status_phenotype;	// linear model
	snp_index_t optimalSNP;
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
				compareMJC= NTest.computeMSC() ;//.computeMSC(2);
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
int Model::makeForwardStep ( Model &biggerModel, const int boundSNP ) {
	double compareMSC =   getMSC(); // DBL_MAX;	// arbitrarily large number
	Model NMax( *data_ );
	int addedSNP=-1;
	// TODO<BB>: Avoid expensive copying of Model objects.
	// oracle bug
		const bool useOracle=false;
	if (!parameter.affection_status_phenotype)
	const bool useOracle = true;// !parameter.affection_status_phenotype;	// linear model
  
	
	snp_index_t optimalSNP;
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
		NTest.addSNPtoModel( data_-> getOrderedSNP(i) );
		
		// check if Regression works, else try next SNP
		if ( NTest.computeRegression() ) {
			if ( NTest.computeMSC(0) < compareMSC ) {
				NMax = NTest;
				compareMSC= NTest.computeMSC(0);//there was a 2 mBic instead of mBIC2
double compareMSC1= NTest.computeMSC(1);
double compareMSC2= NTest.computeMSC(2);
			        cerr<<"makeForwardStep compareMSC="<<	compareMSC<<","<<compareMSC1<<compareMSC1<<compareMSC2<<","<<compareMSC2<<endl;
				addedSNP = data_-> getOrderedSNP(i);
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
bool Model::saveguardbackwardstep(Model &smallerModel)//model3, Model currentModel, Model backwardModel, Model swapHelper)
{//cerr<<"saveguardbackwardstep"<<endl;
  printModel("saveguardbackwardstep");
//reset parameter.expected_causal_snps1
//bool reset=false;
//int reset_causal_snps_=parameter.expected_causal_snps;
//  if(2==parameter.expected_causal_snps)
//{//everthing is ok
//
//}
//else
//{reset=true;
//	 parameter.expected_causal_snps=2;
//	computeMSC();//calulate the new MSC value for the model
//}
//init of model
Model backwardModel(*data_);
Model model1(*data_);	
Model model3(*data_);
model3=*this;

Model *swapHelper=this;
//this is for giving back the original model instead of the smallest one
Model currentModel(*data_);
currentModel=*this;
double bestMSC=getMSC(); //from
printLOG(" bestMSC="+double2str(getMSC()));
int breakfor=0;
bool improvment=false;
 // compute steps
 // improvment=false;
 for (int ii=getModelSize();ii>1;ii--)//currentModel to sm
 { //reset
 	
 int	removedSNP = makeBackwardStepED( backwardModel );
 	*this=backwardModel;//wird ja auch in den Schritten verkleinert
 	//die nichts bringen werden
 double	locMSC=backwardModel.getMSC();
         swapHelper = &backwardModel;
        // *this =*swapHelper;
 	printLOG("Modelsize="+int2str(backwardModel.getModelSize())
 	     +" SNP("+ int2str(removedSNP+1) +")=[] "
     	     +"MSC()=" +double2str(locMSC));
 
// if ther is an improvment then update the best model 	
 if( locMSC<=bestMSC)
 {    breakfor=0;
     	improvment=true;
         bestMSC=locMSC;
 		printLOG("Better Model");
  	model3 = backwardModel;
         *this=model3; 
        // best =	swapHelper;
// model3.printModel();
 
 }
 //else set the counter up and copy the backwardModel to *this
 else
 {breakfor++;
  if (breakfor>=2)
 	{break; }
       
 }
 }
//after the fullbackward
      smallerModel=model3;
 if(improvment)     smallerModel.printModel("best Model after FullBackward"); //print only if improved
		       
//
//
*this=smallerModel;//restauration of original model;
		 	 return improvment;
}
/** In such a step, one SNP is removed from the model.
* From the models with one SNP less, 
* the model with the lowest MJC is selected.
*/
int Model::makeBackwardStepED ( Model &smallerModel ) {
//       printLOG("BackwardStep ED");
	double 	compareMSC = DBL_MAX,returnVal=DBL_MAX; // arbitrary large number
	double   local=0;
	Model NMin( *data_ );
	int 	removedSNP = -1;
        //vector<double> MSC(getModelSize(),-1);//-1 is not a valid size
	//double MSC;;
	//const bool useOracle = !parameter.affection_status_phenotype;	// linear model
	const bool useOracle = false; //for an error
	 snp_index_t optimalSNP=2000000000UL;
	
	if ( useOracle ) {
		returnVal =oraculateOptimalLinearBackwardStep( &optimalSNP );
		if (DBL_MAX== returnVal)
                    cerr<<"oraculateOptimalLinearBackwardStep fails"<<endl;
	        if (optimalSNP==2000000000UL)
		cerr<<"Regression failed by searching the optimal SNP -1 SNP"<<endl;
		cerr<<"optimalSNP="<<optimalSNP<<endl; 
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
                                 local=NTest.computeMSC();
	//ED DEBUG 	    cerr<<" SNP("<< modelSnps_[i]<<")="<<local<<endl;
			if (local < compareMSC)
		{ 
				NMin = NTest;
				compareMSC= local;
				removedSNP = i;
			}
			//it cannot increase MSC globaly because ist starts with $\infty$
		}
		else
		 {//false
		  cerr<<"Regression failed by removing the "<<i<<"SNP"<<endl;
		 }
	}
	
	// check if a valid model is returned, (0-model stays 0-model)

	if ( 0 < getModelSize() ) {
		 smallerModel = NMin;
	}
	return removedSNP;
}



bool Model::makeMultiBackwardStep () {

	printLOG( "Start Multiple-Backward-Step" );
	
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
			if ( parameter.detailed_selction ) {
				printLOG(
						//to state the removed SNP, we need the old Model
						
						"REMOVE	SNP: " + getSNPId(removedSNP) + " ModelSize: " + int2str( BackwardModel.getModelSize() ) + " MSC: " + double2str(BackwardModel.computeMSC())
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

	printLOG( "Finished Multiple-Backward-Step: " + int2str(oldModelsize-getModelSize()) + " SNPs removed" );
	return true;
	
}

bool Model::computeRegression () {
	
	if (!upToDateXMat_) // check if the current XMat is correct
	{
		this->initializeModel();
	}
	
	// perform regression according to data
	if ( parameter.affection_status_phenotype ) {	
		return computeLogRegression();
	}
	else
	{
		return computeLinRegression();
	}
	
}


bool Model::computeLinRegression () {
	// allocate auxiliary gsl objects 
	gsl_vector* beta  = gsl_vector_alloc( getNoOfVariables() );
	gsl_matrix* TMat  = gsl_matrix_alloc( getNoOfVariables(), getNoOfVariables() );
	gsl_vector* TVec  = gsl_vector_alloc( getNoOfVariables() );
	gsl_vector* YCopy = gsl_vector_alloc( data_->getIdvNo() );
	
	// copy YVec_, since some gsl methodes "destroys" (do not seperate between in/output) input data
	gsl_vector_memcpy (YCopy, YVec_);
	
	//Function: int gsl_blas_dgemm (CBLAS_TRANSPOSE_t TransA, CBLAS_TRANSPOSE_t TransB, double alpha, const gsl_matrix * A, const gsl_matrix * B, double beta, gsl_matrix * C)
	//These functions compute the matrix-matrix product and sum C = \alpha op(A) op(B) + \beta C where op(A) = A, A^T, A^H for TransA = CblasNoTrans, CblasTrans, CblasConjTrans and similarly for the parameter TransB. 
	gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1, XMat_, XMat_, 0, TMat); // TMat= X^T * X
	
	//Function: int gsl_blas_dgemv (CBLAS_TRANSPOSE_t TransA, double alpha, const gsl_matrix * A, const gsl_vector * x, double beta, gsl_vector * y)
	//These functions compute the matrix-vector product and sum y = \alpha op(A) x + \beta y, where op(A) = A, A^T, A^H for TransA = CblasNoTrans, CblasTrans, CblasConjTrans.
	gsl_blas_dgemv(CblasTrans,1,XMat_,YVec_,0,TVec);	// TVec = X^T * Y
	
	//+++++++
	// solve lineare system:
	
	//Function: int gsl_linalg_HH_solve (gsl_matrix * A, const gsl_vector * b, gsl_vector * x)
    //This function solves the system A x = b directly using House/bacholder transformations. On output the solution is stored in x and b is not modified. The matrix A is destroyed by the Householder transformations. 
	
	// problem: is system solveable?
	
	gsl_set_error_handler_off();
	int status = gsl_linalg_HH_solve(TMat,TVec, beta);  // slove (X^T * X) * beta = X^T * Y
	
	// Check if gsl_linalg_HH_solve returns a solution, (status != 0), otherwise clear memory and return. 
	if (status)
	{
		gsl_matrix_free (TMat);
		gsl_vector_free (TVec);
		gsl_vector_free (beta);
		gsl_vector_free (YCopy);
		upToDateBetas_=false;
		return false;
	}
	
	gsl_set_error_handler(NULL);
	
	// lineare system solved ! 
	//+++++++
	
	// return betas_
	gsl_vector_memcpy (betas_,beta);
	upToDateBetas_=true; // betas_ computed, therefore uptodate
	
	//+++++++
	// compute residual sum of squares
	
	// compute residuals
	gsl_blas_dgemv( CblasNoTrans,-1,XMat_,beta,1,YCopy); // YCopy is now the error (-1)*X*beta+Y = e, would destroy YVec_
	
	
	//Function: int gsl_blas_ddot (const gsl_vector * x, const gsl_vector * y, double * result)
	//These functions compute the scalar product x^T y for the vectors x and y, returning the result in result. 
	gsl_blas_ddot (YCopy, YCopy, &modelJudgingCriterion_); // e^T e = RSS, residual sum of squares, saved in MJC
	
	// deallocate auxiliary gsl objects 
	gsl_matrix_free (TMat);
	gsl_vector_free (TVec);
	gsl_vector_free (beta);
	gsl_vector_free (YCopy);
	return true;
	
}
/** basic functionality*/
bool Model::scoreTest(string extra){ 
         
	ScoreTestShortcut stsc( *data_);
	int size=data_->getSnpNo();
	SortVec score(size);
	     stsc.ScoreTestShortcut::scoreTests ( *this, score );
ofstream S;
	stringstream ss;//create a stringstream
   ss << getNoOfVariables();//add number to the stream
   //return a string with the contents of the stream
	S.exceptions ( ofstream::eofbit | ofstream::failbit | ofstream::badbit ); // checks if Logfile can be written
	try{ S.open( ( parameter.out_file_name + "_"+ extra +  ss.str()  + ".score" ).c_str(), fstream::out );
		for(int i=0;i<=size-getNoOfVariables();i++)
	S<<score.getId(i)<<" "<< score.getValue(i) <<endl; //see sortvec.cpp
		S.close();
	}
	catch
		(ofstream::failure e)
	{
		cerr << "Could not write score-File" <<endl;
	}	


	
	//     for(int i=0;i<max(1000,size-getNoOfVariables());i++)
	//S<<i<<"="<<score.getId(i)<<endl;
return 	true;
}

size_t Model::scoreTestWithOneSNPless ( size_t position, SortVec &score )
{
//foreach model SNP: remove them (make a  regression) and make scoretest
//take the result snps and add these snps, when the improve the model 
//then take the best of these snp in the new model.
Model intermediateModel(*data_);
	intermediateModel=*this;
	//SortVec score(size); should be created outside

    	intermediateModel.removeSNPfromModel( position );
	intermediateModel.computeLogRegression ();
        //DEBUG intermediateModel.printModel();
        ScoreTestShortcut stsc( *data_);
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
		size=data_->getSnpNo(),
		start = modelSnps_[position] > fenster ? modelSnps_[position] - fenster : 0u,
		stop = min( modelSnps_[position] + fenster, size ),
		nscore = stsc.scoreTests( intermediateModel, score, start, stop );
        //then try this snp for replacing
        return nscore; //does this make sense? 	

//ERICH
} 

bool Model::computeLogRegression () {	
	// the plain C code for logistic firth-regression "logistf package" for R from Georg Heinze is used
	// therefore, the C++ objects have to be converted to C arrays
	// and all parameters for the C function have to be set  
//~ cout << "in log reg"<< endl;
	int 	i, j;
	int intboolfalse = 0;		// plain C has no type "boolen", 0 is "false", every other int is "true"
	int intbooltrue	= 1;
	
	int 	n = data_->getIdvNo();	// # of rows
	int 	k = getNoOfVariables();	// # of columns
	//double 	x[n*k]; 		// an array representing XMat_
	int		y[n];		// an array representing YVec_ // implicit cast! YVec_ is double, but should only contain 0, 1
	double 	beta_array[k];		// the coefficicents
	int 	col_fit[k];				
	double 	pi[n];
	double 	weight[n];
	double 	offset[n];
	double 	var[k*k];
	double 	H[n];
	double 	loglik=0;		        // the log-likelihood
	int 	iter, evals;
	double 	lchange, ret_max_U_star, ret_max_delta; 
	
	for(i=0; i < n; i++) 
	{
	//	for (j=0; j < k; j++)
	//	{
	//		x[i*k + j] = gsl_matrix_get(XMat_, i, j);// copy XMat_ to array
	//	}
		y[i] = int(gsl_vector_get(YVec_, i)); 		// copy YVec_ to array	 
		offset[i] = 0;			                // set offset to 0
		weight[i]=1;		    	                // set weights to 1
		
	}
	
	for (j=0; j < k; j++)
	{
		beta_array[j] = gsl_vector_get(betas_,j); 	//  initialize according to old model //intialize beta as 0
	//	cout<<beta_array[j]<<",";  //DEBUG
		col_fit[j] = 1; 	// use every column
	}
               //cout<<endl; DEBUG
	// logistic firth regression 
	if ( 
		logistffit(
					// Input
					&k, // # of rows of X (# of variables) BB: I'd rather call that "columns"
					&n, // # of columns of X (# of samples) BB: I'd rather call that "rows"
					XMat_,//x, // 
					y,  //
					// Output
					beta_array, 
					var, // 
					pi, //
					H, // 
					&loglik, 
					&iter, // # of main iterations;
					&evals, 
					&lchange, // the change in the log-likelihood between steps
					&ret_max_U_star, //
					&ret_max_delta,			
					//// optional Input
					weight, // the weights
					offset, 	// offset
					&intbooltrue,	//  use firth-regression
					col_fit, // a "boolean" vector indication which colums to use
					beta_array,	// initial values for beta
					&intboolfalse,	// only evaluate likelihood
					//// Control Parameters
					&parameter.logrC_maxit,	
					&parameter.logrC_maxhs,	
					&parameter.logrC_maxstep,
					&parameter.logrC_lconv,
					&parameter.logrC_gconv,
					&parameter.logrC_xconv	
				)
		)
	{
		modelJudgingCriterion_ = loglik; //will be calculated in firth-fit	
	        // return betas
		for (j=0; j < k; j++)
		{
			gsl_vector_set(betas_, j, beta_array[j]);
		}
		
		upToDateBetas_= true;
		return true; // logistic regression worked, betas are updated
	}
	else
	{
		//~ cout << "iter: "<< iter<<endl;
		return false;	// logistic regression failed
	}	

}



double Model::computeSingleRegressorTest ( const snp_index_t snp ) {
	
	// Intitialisation
	
	// check if we have an 1-SNP model
	if ( 1 != getModelSize() ) {
		modelSnps_.clear();
		modelSnps_.push_back( snp );
		initializeModel();
	}
	else // replace in 1-SNP model the old SNP with new SNP snp
	{
		modelSnps_.at(0) = snp;
		const Vector xVec = const_cast<MData*>( data_ )->getXcolumn( snp );
		// TODO<BB>: Use vectorial replacement instead of loop
		for ( int i = 0; i < data_->getIdvNo(); ++i ) {
			// the SNP is in the last column noOfVariables - 1		
			gsl_matrix_set( XMat_, i, getNoOfVariables() - 1, xVec.get( i ) );
		}
	}
	
	// XMat_ intialised and uptodate, betas_ are now false
	upToDateXMat_ = true;
	upToDateBetas_ = false;
	
	// perform single marker test according to data
	if( parameter.affection_status_phenotype ) {
		return computeSingleLogRegressorTest( snp );
	}
	else
	{
		return computeSingleLinRegressorTest( snp );
	}
}



/** @see Model::computeLinRegression() for more detailed comments */
double Model::computeSingleLinRegressorTest ( const snp_index_t snp ) {
	
	// allocate auxiliary gsl objects
	gsl_vector * beta = gsl_vector_alloc( getNoOfVariables() );
	gsl_matrix * TMat =gsl_matrix_alloc( getNoOfVariables(), getNoOfVariables() );
	gsl_matrix * InvMat = gsl_matrix_alloc( getNoOfVariables(), getNoOfVariables() );
	gsl_vector * TVec = gsl_vector_alloc( getNoOfVariables() );
	gsl_vector * YCopy = gsl_vector_alloc( data_->getIdvNo() );
	gsl_vector_memcpy( YCopy, YVec_ );
	
	gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1, XMat_, XMat_, 0, TMat); // TMat= X^T * X
	gsl_blas_dgemv(CblasTrans, 1, XMat_, YVec_, 0, TVec);	// TVec = X^T * Y
	
	//++++++++++ 
	// invert matrix
	gsl_permutation * p = gsl_permutation_alloc ( getNoOfVariables() );	// needed for inversion of matrix
	
	// Need to invert TMat, we obtain InvMat = (X^T * X)^(-1);
	int s=0;
	gsl_linalg_LU_decomp(TMat, p, &s);//changes TMat
	
	// Check if Matix is singular: 
	// 		-> one diagonal-element is 0
	// needed to run gsl_linalg_LU_invert(TMat, p, InvMat) without error
	
	for ( int i = 0; i < getNoOfVariables(); ++i ) {
		if ( 0  == gsl_matrix_get (TMat, i, i) )
		{
			// the matrix is not invertabele. 
			printLOG("Error single linear regressor test for SNP "+ (data_->getSNP( snp ))->getSnpId()+ ": Regressionmatrix not invertable!" );
		
			gsl_matrix_free (TMat);
			gsl_matrix_free (InvMat);
			gsl_vector_free (TVec);
			gsl_vector_free (beta);
			gsl_vector_free (YCopy);
			gsl_permutation_free (p);
			
			return 1; // no meaningfull result possible, so we set the p-value to 1
		}
	}
	
	gsl_linalg_LU_invert(TMat, p, InvMat);
		
	// matrix inverted 
	//++++++++++
	 
	gsl_blas_dgemv(CblasNoTrans,1,InvMat,TVec,0,beta); 	// beta = (X^T * X)^(-1)* X^T * Y
	
	gsl_blas_dgemv(CblasNoTrans,-1,XMat_,beta,1,YCopy); // YCopy is now the error (-1)*X*beta+Y = e, would destroy YVec_
	gsl_blas_ddot (YCopy, YCopy, &modelJudgingCriterion_); // e^T e = RSS
	
	gsl_vector_memcpy (betas_,beta);
	upToDateBetas_ = true;
	
	// compute test statistic for the regressionscoefficient of the SNP
	double diff, div, test_stat;
	double x_ii = gsl_matrix_get( InvMat, getNoOfVariables() - 1, getNoOfVariables() - 1 );
	double beta_i = gsl_vector_get( beta, getNoOfVariables() - 1 );

	// Free Memory befor return
	gsl_matrix_free (TMat);
	gsl_matrix_free (InvMat);
	gsl_vector_free (TVec);
	gsl_vector_free (beta);
	gsl_vector_free (YCopy);
	gsl_permutation_free (p);
	
	
	
	diff = data_->getIdvNo() - getNoOfVariables();	//Intercept not counted 
	if ( diff <= 0 ) 
		{return 1;} 
	else 
	{ 
		div = sqrt( (modelJudgingCriterion_ * x_ii)/diff);
	}
	if ( div == 0 ) 
		{return 0;} // falls RSS/MJC = 0, kein Fehler, extrem hohes T, niedriges P, falls diagonalelement 0 :dann alle SNPs gleich (weil keine Varianz (aus der Kovarianzmatrix)
		// und die MAF-Bedigung >0.05 ist nicht erfüllt 	
	else
	{
		test_stat = (beta_i)/div;
		
		// compute the p-value for the test-statistic. 	
		// to checkt: http://home.ubalt.edu/ntsbarsh/Business-stat/otherapplets/pvalues.htm#rtdist
		return ( gsl_cdf_tdist_P ( -fabs(test_stat), diff ) + gsl_cdf_tdist_Q ( fabs(test_stat), diff )); // F(-|x|) +  (1-F(|x|)) 
	}

}



double Model::computeSingleLogRegressorTest ( const snp_index_t snp ) {
	// compute the genotype frequency table of SNP snp
	GenotypeFreq freqOneSNP( *data_, snp );
	
	// choose test for single marker test for case-control data 
	if ( parameter.cc_SingleMarkerTest_ChiSq ) {
		return freqOneSNP.calculateChiSquare();
	} else if ( parameter.cc_SingleMarkerTest_CATT ) {
		return freqOneSNP.calculateCATT();		
	} else {	// default
		return freqOneSNP.calculateChiSquare();
	}	

}

double Model::computeMSC ( const int typeNr  ) {
if(DEBUG)	cerr<<"MSC variant="<< typeNr<<endl;
	return computeMSC( typeNr, getMJC() );
}

double Model::computeMSC ( const int typeNr, double mjc ) {
//	cout<<"typeNr"<<typeNr<<endl;
	if ( !upToDateBetas_ ) // check if betas_ and therefor also MJC is up-to-date, otherwise update
	{
		computeRegression();
	}

	int n = data_->getIdvNo();
	int p = data_->getSnpNo();
	//int p = 780675; // p as magic number
	int q = getModelSize();
	//double d = -2 * log( parameter.ms_ExpectedCausalSNPs );
	
	// choose the likelihood part depending if the Data is quantitative or affection
	double LRT;
	if ( parameter.affection_status_phenotype ) {
		LRT = (-2.0)*( mjc - data_->getLL0M()); // LRT = -2 log (likelihood(0-model)/likelihood(model))
	}
	else
	{
		LRT = n*log( mjc ) ;  // n * log(RSS)
	}
	

	// compute the modelselection criterion
	switch (typeNr)
	{
		case 1: /*BIC*/ 
		        msc = LRT + q*log(n);
			return  msc;
		break;
		
		case 2: /*just LRT, for internal comparisons (Backward Step)*/
		        msc = LRT;
			return msc;
		break;

		default: /*mBIC2*/
		        double d = -2 * log( parameter.ms_ExpectedCausalSNPs );
		//	cout<<"q"<<q<<"p="<<p<<"d="<<d<<"LRT="<<LRT<<endl;

		        msc = LRT + q*(log(n) + 2* log(p) + d ) - 2*(log_factorial(q));	
			return  msc;
		break;
	}	
	
}

double Model::getMSC() const {
	return msc;
}

/** Destructor */
Model::~Model () {
  clearModel();     // I put lines below to clearModel()
}

void Model::printModelNew() const {
	ofstream	SNPL;
	
	SNPL.exceptions ( ofstream::eofbit | ofstream::failbit | ofstream::badbit ); // checks if Logfile can be writt
	try
	{
		SNPL.open( ( parameter.out_file_name + "_SNPListNew.txt" ).c_str(),  fstream::out );
		
		
		for ( size_t i = 0; i < getModelSize(); ++i ) {
			const size_t snp = modelSnps_.at(i);
			const Vector xVec = const_cast<MData*>( data_ )->getXcolumn( snp );
			SNPL << getSNPId(i) <<" <- c(";
			for ( size_t idv = 1; idv  < data_->getIdvNo(); ++idv ) {
				SNPL << ( 0 < idv ? "," : "" );
				SNPL << xVec.get( idv );
			}
			SNPL<< ")"<< endl;
			
		}

		for ( int i = 0; i < parameter.covariables; ++i ) {
			SNPL << data_->getCovMatElementName(i)<<" <- c(";
			const Vector covVec = const_cast<MData*>( data_ )->getCovariateColumn( i );
			SNPL <<  covVec.get( 0 );
			for (int j= 1; j < data_->getIdvNo(); j++)
			{
				SNPL << "," << covVec.get( j );
			}
			SNPL<< ")"<< endl;	
		
		}

			SNPL << "Y" <<" <- c(";
			SNPL<< gsl_vector_get(YVec_,0);
			for (int j= 1; j < data_->getIdvNo(); j++)
			{
				SNPL <<","<< gsl_vector_get(YVec_,j);
			}
			SNPL<< ")"<< endl;
		SNPL <<endl;
		
			SNPL << "Intercept" <<" <- c(1";
			for (int j= 1; j < data_->getIdvNo(); j++)
			{
				SNPL <<","<<1;
			}
			SNPL<< ")"<< endl;
		
		//~ SNPL << "summary(lm(Y~Dummy1 + Dummy2 + Dummy3 +";
		SNPL <<"X <- cbind( Intercept";
		for  (int i = 0; i< getModelSize(); i++)
		{
			//~ if ( it_SNPList != SNPList.begin() )
			SNPL << ",";;
			SNPL <<  getSNPId(i);
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

/**
 * @brief Print model information (snps and msc) on the screen. Just for testing GA
 */
ostream &operator << ( ostream &out, const Model &m ) {
	vector<snp_index_t>::const_iterator it = m.modelSnps_.begin();
  out << "[";
  if (m.modelSnps_.size() > 0)
  {
    out << *it;
    ++it;
  }  
  for (; it != m.modelSnps_.end(); ++it)
    out << ", " << *it;
  return out << "], msc: " << m.msc;
}  

/**
 * @brief gets snp at given position, 
 * @param pos relative position of snp at vector modelSnps_
 * @return snp
 */
snp_index_t Model::getSNPat ( const snp_index_t pos ) const {
  if (pos < modelSnps_.size())
    return modelSnps_[pos];
  else
  {
    cerr << "getSNPat(...): " << pos << "out of range [0, " << modelSnps_.size() << "]" << endl;
    exit(-1);
  }
}

vector<snp_index_t> Model::getModelSnps() const {
	return modelSnps_;
}

/**
 * @brief Clears previus model and creates new model from given snps
 * @param snps - set of snps 
 * //WARNING I think there is faster way to create new model from given snps. I'm goint to do it. 
 */
void Model::createFromSNPs( const set<snp_index_t>& snps ) { 
  clearModel();
  computeRegression();
  for ( set<snp_index_t>::const_iterator it = snps.begin(); it != snps.end(); ++it ) {
    addSNPtoModel( *it );
  }
}

/**
 * @brief removes SNP form Model, 
 * @param oneSNP is SNP value at vector modelSnps_
 * @return true if SNP removed, and false if an error occors
 * ------------------------------------------------------------------------------
 */
bool Model::removeSNPValFromModel( const snp_index_t oneSNP ) {
  vector<snp_index_t>::iterator it = find(modelSnps_.begin(), modelSnps_.end(), oneSNP);
  return removeSNPfromModel(it - modelSnps_.begin());
}

/**------------------------------------------------------------------------------
 * @brief Deallocates memory and makes model ready for create new model with createFromSNPs(...)
 * WARNING Please check
 * -----------------------------------------------------------------------------
 */ 
void Model::clearModel()
{
  upToDateXMat_ = false; // matrix not allocated
  modelSnps_.clear(); 
  upToDateBetas_ = false;    // WARNING Is this line OK? If it's not, then delete, please.

  if ( NULL != XMat_ ) {
    gsl_matrix_free (XMat_);
    XMat_ = NULL;
  }
  if ( NULL != YVec_ ) {
    gsl_vector_free (YVec_);
    YVec_ = NULL;
  }
  if ( NULL != betas_ ) {
    gsl_vector_free (betas_);
    betas_ = NULL;
  }
}

double Model::computeMSCfalseRegression(const int typeNr)
{
  vector<snp_index_t> removedSnps;
  msc = computeMSCfalseRegression(typeNr, removedSnps);
  return msc;
}
  
double Model::computeMSCfalseRegression(const int typeNr, vector<snp_index_t> &removedSnps)
{
  if (computeRegression() == false)
  {
    //cout << "try computeMSC for model: " << *this << "modelel size: " << modelSnps_.size() << endl;
    removedSnps.push_back(modelSnps_[modelSnps_.size() - 1]);
    //cout << "remove snp: " << modelSnps_.size() - 1 << endl;;
    if (removeSNPfromModel(modelSnps_.size() - 1) == false)
      cout << "removeSNPfromModel failed!" << endl;
    //cout << "try computeMSC for model (-1): " << *this << endl;
    //char cc; cout << "Press a key... "; cin >> cc;
    msc = computeMSCfalseRegression(typeNr, removedSnps);
    
    snp_index_t snp = removedSnps.back();
    
    removedSnps.pop_back();
    addSNPtoModel(snp);
    int n = data_->getIdvNo();
    int p = data_->getSnpNo();
    int q = getModelSize() + 1;   // + 1 for additional corelated snp
    double d = -2 * log( parameter.ms_ExpectedCausalSNPs );
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
    return computeMSC(typeNr);
  }
}
bool Model::selectModel(Model &startFromModel,int  PValueBorder,int maxModel)
{
printLOG("SCORE SCORE SELECT MODEL");
if(getModelSize()>parameter.maximalModelSize)
	return false; //nothing will be done, when ModelSize is >49
cerr<<"S1"<<endl;
	double best= getMSC();
	 Model 	*backwardModel;
	 backwardModel=this;
//ORIGINAL value is now default   int PValueBorder =100,
	int *startIndex;
 	int dummy=0;
	    startIndex=&dummy;
	    PValueBorder =2000; //override to high PValueBorders !!!
 	double bestMSC=getMSC(); //this one is the best up to now!

 ScoreTestShortcut stsc( *data_);
 	int size=data_->getSnpNo(),JJ=0;
	bool improvment=true,improvment2=true, stop=false;
	SortVec score(size);
 this->computeRegression();//for SCORE maybe not calculated?
        stsc.ScoreTestShortcut::scoreTests ( *this, score );
 vector<int> Score(size);
        for(int i=0;i<=size-getNoOfVariables();i++)
 		Score[i]= score.getId(i);	
cerr<<"S2"<<endl;
	while (improvment&&!(getModelSize()>min(parameter.maximalModelSize,maxModel))) //||improvment
 {cerr<<getModelSize()<<endl;
   while(improvment=makeForwardStepLogistic(JJ, &bestMSC,  PValueBorder, startIndex,Score))
   { if(getModelSize()>min(parameter.maximalModelSize,maxModel)) break;
	   startFromModel=*this; //see //startFromModel is the inputvariable
	   	 if (getModelSize()>min(parameter.maximalModelSize,maxModel))
		 break;
	  //replaceModelSNPbyNearFromSCORE(*startIndex , PValueBorder,Score);
          replaceModelSNPSCORE();
   improvment2=saveguardbackwardstep( startFromModel );
   *this=startFromModel;//if improvment2
   }

 } 
  stop=finalizeModelSelection( *backwardModel, JJ,  improvment||improvment2,  PValueBorder,  startIndex,Score);
  if( getMSC()<best)
	 return true;
  else 
	 return false; 
 
}

void Model::checkallSNPS ()
{       cerr<<"WARNUNG NUR EIN SPEZIALFALL FUNKTIONIERT MIT DIESER FUNKTION"<<endl;
	Model model0(*data_);
	model0=*this;
	int nSNPs=	data_->getSnpNo();
	vector<int> snps( nSNPs );//SCORE
        vector<double> scores( nSNPs);//SCORE

        snp_index_t orderedSnpIndex =50505;// data_-> getOrderedSNP( 0 );
//snps[0]=orderedSnpIndex;
model0.addSNPtoModel(orderedSnpIndex);
model0.addSNPtoModel(53580);
//Erich maybe score need a Regression ??,
if (  model0.computeRegression() ) 
          model0.computeMSC(0);
	ScoreTestShortcut stsc( *data_);//SCORE

	SortVec score(nSNPs);//SCORE
stsc.ScoreTestShortcut::scoreTests ( model0, score );//SCORE
	 vector<int> Score(nSNPs);//SCORE
        for(int i=0;i<=nSNPs-model0.getNoOfVariables();i++)//SCORE/
 		Score[i]= score.getId(i);//SCORE
	int sModel=model0.getModelSize()-1;
orderedSnpIndex = data_-> getOrderedSNP( 1 );
model0.addSNPtoModel(orderedSnpIndex);
sModel=model0.getModelSize()-1;
//dies soll die Auswirkung der Richtigen SNPS auf position 1 des Model
//relativ zum an und für sich gewählten den data_-> getOrderedSNP( 0)
//sichtbar machen
snps[0]=orderedSnpIndex;
 // check if regression works properly, else try next SNP
         if (  model0.computeRegression() ) 
            scores[0] = model0.computeMSC(0);

for (int i=2;i<1000/*nSNPs*/;i++)
{
	orderedSnpIndex = data_-> getOrderedSNP( i );
	snps[i-1]=orderedSnpIndex;
	model0.replaceSNPinModel(orderedSnpIndex,sModel);
    if (  model0.computeRegression() ) 
            scores[i-1] = model0.computeMSC(0);
//cerr<<snps[i]<<" "<<scores[i]<<endl;
}
//SortVec(sModel,&snps[0],&scores[0],false);
 ofstream Y;
cerr<< parameter.out_file_name + "1Models"<<endl;
	Y.exceptions ( ofstream::eofbit | ofstream::failbit | ofstream::badbit );
	try {Y.open( ( parameter.out_file_name + "1Models" ).c_str(), fstream::out );}
		catch  ( ofstream::failure e/*xception*/ ){
		cerr <<"Could not  write 1 Model file file " << (parameter.out_file_name + "1Models"  ).c_str()<< endl; }
 
			for(int i=0;i<2000/*nSNPs*/;i++)
			{//cerr<<i;
			Y<<snps[i]<<"\t"<<scores[i]<<"\t"<<Score[i]<<endl;
			}
	

 try{Y.close();}
    catch ( ofstream::failure e/*xception*/ ) {
		cerr << "Could not close 1 Model file file  " << (parameter.out_file_name + "1Models"  ).c_str()<< endl;}
		}

