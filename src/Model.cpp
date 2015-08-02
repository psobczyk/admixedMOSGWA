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
#include "Exception.hpp"
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
bool DEBUG=false,DEBUG2=false,DEBUG3=false;
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
	return 1 + data_->getCovNo() + getModelSize();
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
		const size_t idvs = data_->getIdvNo();
		AutoVector xVec( idvs );
		data_->getXcolumn( snp, xVec );
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

bool Model::printXmat(const string& extra){}
//  	char mat [data_->getIdvNo()][getModelSize()];
// 	// genotype-data
// 	for(int i=0;i<data_->getIdvNo();++i){
//		const Vector genotypes = data_->getX().rowVector( i );
//		for ( int j = 0; j < getModelSize(); ++j ) {
//			mat[i][j+ parameter.dummy_covariables + parameter.covariables + 1]=
//			   	genotypes.get( modelSnps_.at(j) );}
//    //init of genotype matrix
//}	cerr<<endl;
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
//    hsize_t  cdims[2];
// 
//    int      idx;
//    int      i,j, numfilt;
//   // int      buf[DIM0][DIM1];
//  // int      rbuf [DIM0][DIM1];
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
//	   data->getIdvNo()][getModelSize()]*/
//    dims[0] = data_->getIdvNo();
//    dims[1] = getModelSize();
//    dataspace_id = H5Screate_simple (2, dims, NULL); //2 ist der Rang
//
//    plist_id  = H5Pcreate (H5P_DATASET_CREATE);
//
//    /* Dataset must be chunked for compression */
//    cdims[0] = 20;
//    cdims[1] = getModelSize();
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
//    status = H5Pclose (plist_id);
//    status = H5Fclose (file_id);
//
//
//}

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
	for (int j=0;j<30;j++)
{cerr<<"XMat"<<j<<" ";
	for ( int i = 1 + data_->getCovNo(); i < getNoOfVariables(); ++i ) {
	cerr<<gsl_matrix_get(XMat_,j,i)<<" ";
	}
cerr<<endl;
}
}

	const snp_index_t reset = 1 + data_->getCovNo() + position;		//position 0 is the first 
 //cout<<"reset="<<reset<<endl;
 upToDateXMat_= false;
	modelSnps_.at( position ) = snp;
// cerr<<"reset="<<reset;
	// add the new column at the end
       // gsl_matrix_set_col(XMat_,position,data_->xMat.columnVector(snp));
		const size_t idvs = data_->getIdvNo();
		AutoVector xVec( idvs );
		data_->getXcolumn( snp, xVec );
		// TODO: Refactor to use Vector.copy
		for ( int i = 0; i < idvs; ++i ) {	//instead position reset 
			gsl_matrix_set( XMat_, i, reset, xVec.get( i ) );
		}
			gsl_vector_set(betas_,reset,0); //0 is relativ good for an unknow variable.
		
if(DEBUG2)
{cerr<<"XMat_->size1="<< XMat_->size1<<endl
<<"XMat_->size2="<<XMat_->size2<<endl<<"getNoOfVariables("<<getNoOfVariables()<<endl;

for (int j=0;j<30;j++)
{cerr<<"XMat"<<j<<" ";
	for ( int i = 1 + data_->getCovNo(); i < getNoOfVariables(); ++i ) {
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
				if ( i < 1 + data_->getCovNo() + snp ) {
					gsl_matrix_get_col( CopyV, XMat_, i );
					gsl_matrix_set_col( NewXMat, i, CopyV );
					gsl_vector_set(NEWbetas_,i,gsl_vector_get(betas_,i));
				} else if ( i > 1 + data_->getCovNo() + snp ) {
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
	for ( size_t cov = 0; cov < data_->getCovNo(); ++cov ) {
		Vector xVec = xMat.columnVector( col++ );
		data_->getCovariateColumn( cov, xVec );
	}
	for ( size_t modelSnp = 0; modelSnp < getModelSize(); ++modelSnp ) {
		const size_t snp = modelSnps_.at( modelSnp );
		Vector xVec = xMat.columnVector( col++ );
		data_->getXcolumn( snp, xVec );
	}
	assert( getNoOfVariables() == col );

	// Set YVec_
	Vector yVec = Vector( *YVec_ );
	data_->getY( yVec );

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
     //mean=0.5-mean; //this should do the trick, no longer the mean but the difference to the meansetting            setting
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
		data_->getY( yVec );
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
			const size_t idvs = data_->getIdvNo();
			AutoVector yVec( idvs );
			data_->getY( yVec );
			for ( size_t idv = 0; idv < idvs; ++idv ) {
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
void Model::printModel (
	const string& out,
	const int selectionCriterium,
	const string& filemodifier
) {
	stringstream ss; // to save output
	ofstream OUT; // output model to file

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
				<< getBeta( 1 + data_->getCovNo() + i ) << "\t"
				<< data_->getSingleMarkerTestAt( modelSnps_.at(i) )
				<< endl;
	}
        if(0==getModelSize())
		ss<<"\tIntercept \t \t \t is not available empty model"<<endl;
	else
	ss << "\tIntercept \t \t \t" << getBeta(0)<< endl;

	for ( int i = 0; i < data_->getCovNo() ; ++i ) {
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
        const int fenster=400; //keine Fenster mehr left and right side of the fenster
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
// REMARK<BB> 0 should not be a reserved value, it is a valid SNP index.
bool Model::replaceModelSNPbyNearFromCAT (
	int currentPosition,
	int PValueBorder,
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
	const int nSNP = data_->getSnpNo();
	printLOG(
		"replaceModelSNPbyNearFromCAT currentPosition=" + int2str(currentPosition)
		+ " window=" + int2str(fenster)
		+ " grace=" + int2str(grace)
		+ " last model SNP =" + int2str( modelSnps_.at(getModelSize()-1))
	);
	double bestMSC=getMSC(); //the msc in the current model is the best one up to now
	Model model0(*data_);
	const unsigned int ref=data_-> getOrderedSNP( currentPosition );
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
gsl_ran_shuffle (r, a->data, N, sizeof (size_t));//reference to data!!!

//gsl_permutation_fprintf (stdout, a, " %u");
//	for(int jo=getModelSize()-1; jo >=0; jo--)
for(int jo=0; jo<getModelSize(); jo++)
{ int j=gsl_permutation_get(a,jo);
	//DEBUG cerr<<"a["<<jo<<"]="<<j<<endl;
	for(snp_index_t i=0/*curentPosition*/;i<min(max(PValueBorder+grace,1000),nSNP);++i) //saveguard against overrun of nSNP and only 500SNP in large Problems, with 1000 in most cases all relevant SNP will found
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
				val = model0.computeMSC( selectionCriterium );
                    //  printLOG("ModelSNP="+ int2str(j) + "POS="+ int2str(data_-> getOrderedSNP(i))+"val=" +double2str(val));

                  if(val<bestMSC-0.0001)
		       	//saveguard against rouning errors in logistic regression
                      { double  alt =bestMSC;
			 bestMSC=val;
				printLOG(
					"Better Model at position " + int2str(j)
					+ ": SNP= " +int2str( modelSnps_.at(j) )
					+ " is replaced with " + int2str( data_-> getOrderedSNP(i) )
					+" oldMSC=" + double2str(alt)
					+ " newMSC="+ double2str(bestMSC)
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
 const int nSNP=data_->getSnpNo();
 double bestMSC=getMSC(); //the msc in the current model is the best one up to now
 Model model0(*data_);
 SortVec score(200); //Warning 50+50+1 that should be variable
 	for(int j=getModelSize()-1; j >=0; j--)
	{//replace
	int nscores=	scoreTestWithOneSNPless(j, score);
		for(snp_index_t i=0;i<nscores;++i) 
		{       model0=*this;
                        model0.replaceSNPinModel (score.getId(i) ,  j );

	                double val=DBL_MAX;
                          model0.computeRegression(); //regression should be calculated!
			val= model0.computeMSC( selectionCriterium );
                
			if(bestMSC>0?val<0.9999*bestMSC:val<1.0002*bestMSC)
	                      { double  alt =bestMSC;
				 bestMSC=val;
				 printLOG("!!!!!!!!!!!!!!!!!!!Better Model at Modelposition " + int2str(j) + " SCOREPostition=" +int2str(i) +" SNP= "
				 +int2str( modelSnps_.at(j)) + " is replaced with " + int2str(score.getId(i))
				 +" oldMSC="+ double2str(alt)+ " newMSC ="+ double2str(bestMSC));
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


double Model::oraculateOptimalLinearForwardStep( snp_index_t *snp, size_t bound ) const {
	// First, intermediate step towards using linalg instead plain GSL
	const Matrix xMat( *XMat_ );
	const Vector yVec( *YVec_ );
	AutoVector xVec( data_->getIdvNo() );

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
bool Model::finalizeModelSelection (
	Model &backwardModel,
	bool improvement,
	int PValueBorder,
	int *startIndex,
	vector<int> score,
	const int selectionCriterium
) {
	if ( !improvement ) {
		printModel( "no improvement", selectionCriterium );
		*startIndex = 0;
		if( !parameter.affection_status_phenotype ) {
			cerr << "finalise with score not implemented for continuous traits" << endl;
			backwardModel=*this;
			improvement = saveguardbackwardstep( backwardModel, selectionCriterium);
			// not use makeMFFL in this case 
		} else {
			makeMFFL(
				max( parameter.reset, PValueBorder ),
				startIndex,
				score,
				selectionCriterium
			);
		}
		//don't set to an explicit value because of memory leaks when the number of variables is very small!
		backwardModel = *this;
		improvement = saveguardbackwardstep( backwardModel, selectionCriterium );
		if ( improvement ) {
			printLOG( "finalizeModelSelection" );
			*this=backwardModel;	// REMARK<BB>: Erich had this line after the return … ?
			return true;
		}

		printModel( "final model", selectionCriterium );
		return true;
  }
else

{ 	printLOG("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!after forward step");
	return false;
	}
}// finalizeModelSelection()


bool Model::finalizeModelSelection (
	Model &backwardModel,
	bool improvement,
	int PValueBorder,
	int *startIndex,
	const int selectionCriterium
) {
	if ( !improvement ) {
		printModel( "no improvement", selectionCriterium );
		*startIndex = 0;
		if ( !parameter.affection_status_phenotype ) {
			cerr << "finalise not implemented for continuous traits" << endl;
			backwardModel=*this;
			improvement = saveguardbackwardstep( backwardModel, selectionCriterium);
			//not use makeMFFL in this case 
		} else {
		//diese Schritte sind eigentlich nur einmal nötig da, man sicher nichts
  //  gewinnt wenn man die ersten paar SNP oft und oft wiederholt
   printLOG("finalizeModelSelection last run with min of  parameter.PValueBorder or parameter.reset ");
			makeMFFL(
				min( PValueBorder, parameter.reset ),
				startIndex,
				selectionCriterium
			);
		}

  //  don't set to an explicit value because of memory leaks when the number 
  // of variables is very small!
		backwardModel = *this;
		improvement = saveguardbackwardstep( backwardModel, selectionCriterium );
		if ( improvement ) {
			printLOG( "improvement after newstart forward step" );
			*this = backwardModel;	// REMARK<BB>: Erich had this line after the return
			return true;
			printModel( "final model", selectionCriterium );
		}

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
bool Model::makeForwardStepLinear (
	Model *forwardModel,
	double* bestMSC,
	int PValueBorder,
	int *startIndex,
	const int selectionCriterium
) {
	bool improvement = false;

 Model model3(*data_);
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
	        printLOG("better bigger  Model in Iter"+int2str(ii)+" bestMSC="+double2str(*bestMSC));
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
	int PValueBorder,
	int *startIndex,
	vector<int> score,
	const int selectionCriterium
) {
	bool improvement = false;
int dummy= *startIndex;
//cout<<"dummy="<<dummy<<endl;
 if( dummy<parameter.reset)//ERICH wenn man 100 als Grenze nimmt sollte man das hier auch herabsetzen 
	 // 500 waren gut wenn man eine Grenze von 3000 hatte
	 // also docg eher PValueBorder/6
 { //cerr<<"startIndex="<<*startIndex<< "but 0 is used"<<endl;
   dummy=0;*startIndex=dummy;
		makeMFFL( PValueBorder, startIndex, score, selectionCriterium );
 }
 else if (dummy>=parameter.reset) //these values are only guesses
	 //da auch 300 für 3000 
	 //daher 30 für 100
 {  dummy=max(0,dummy-parameter.jump_back);*startIndex=dummy; //better 500 (?)
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
	int PValueBorder,
	int *startIndex,
	const int selectionCriterium
) {
	bool improvement = false;
int dummy= *startIndex;
//cout<<"dummy="<<dummy<<endl;
 if( dummy< parameter.reset) //500 
 { printLOG("startIndex="+int2str(*startIndex)+ "but 0 is used");
   dummy=0;*startIndex=dummy;
		makeMFFL( PValueBorder, startIndex, selectionCriterium );
 }
 else if (dummy>= parameter.reset) //these values are only guesses

 {  dummy=max(0,dummy-parameter.jump_back);*startIndex=dummy; //min because nothing prevents to be jump_back to be bigger than reset!
         printLOG("StartIndex= "+int2str(*startIndex)+" is used, with jump_back="+int2str(parameter.jump_back));
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

bool Model::makeMFFS(int PValueBorder, int* startIndex)
{       //if (startIndex==NULL)
	    //startIndex
	int selectionCriterium = parameter.selectionCriterium;
         //FAST MultipleForward is needed locally!	
        bool oldValue= parameter.ms_FastMultipleForwardStep;
        parameter.ms_FastMultipleForwardStep=true;
        
       
        const int altSize= getModelSize();
makeMultiForwardStep ( PValueBorder,selectionCriterium, startIndex);

 parameter.ms_FastMultipleForwardStep=oldValue;
 return true;
}
//and now with score test 
bool Model::makeMFFL(
	int PValueBorder,
	int* startIndex,
	vector<int> score,
	const int selectionCriterium
) {       
//	int selectionCriterium=0;
         //MultipleForward is needed locally!	
   //     bool oldValue= parameter.ms_FastMultipleForwardStep;
        //parameter.ms_FastMultipleForwardStep=true;
        
       
        const int altSize= getModelSize();
makeMultiForwardStepScore ( PValueBorder,selectionCriterium,startIndex, score  );
cout<<"MFFL startIndex SCORE"<<*startIndex<<"Model Size="<<altSize<<endl;

 //parameter.ms_FastMultipleForwardStep=oldValue;
}


/**  makeMFFL make Fast Forward local but without changing to fast search*/

bool Model::makeMFFL(
	int PValueBorder,
	int* startIndex,
	const int selectionCriterium
) {       
         //MultipleForward is needed locally!	
   //     bool oldValue= parameter.ms_FastMultipleForwardStep;
        //parameter.ms_FastMultipleForwardStep=true;
        
       
        const int altSize= getModelSize();
makeMultiForwardStep ( PValueBorder,selectionCriterium,startIndex  );
cout<<"MFFL startIndex"<<*startIndex<<"Model Size="<<altSize<<endl;

 //parameter.ms_FastMultipleForwardStep=oldValue;
}








//ScoreScoreScoreScoreScoreScoreScoreScoreScoreScoreScoreScoreScore
bool Model::makeMultiForwardStepScore (
	int PValueBorder,
	const int selectionCriterium,
	int* startIndex,
	vector<int> score
) {
//if(score) 
	//getscores this means that we have logistic regression 
       // data_->getOrderedSNP should be something different for score
if(NULL==startIndex)
  {
    int dummy=0;	  
   startIndex=&dummy;
  }
  int returnIndex=*startIndex;
  printLOG( "Start Multiple-Forward-Step SCORE" );
  if (0==PValueBorder)
  {
   PValueBorder=data_->getSnpNo(); 
   printLOG( "Default setting for PValuePorder: select all");
  }
  //if modelsize is bigger than0 and selectionCriterium ~! 1
  //then new modus
int startSize=getModelSize();
int newSNP=0;

  cout<<"ModelSize ="<<startSize<<endl;

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
	newSNP && i < PValueBorder && i < data_->getSnpNo();
      ++i
    ) {
      // progress checking
      if (0==i%200 )
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
              "ADD    SNP: " + data_->getSNP(orderedSnpIndex ).getSnpId() + " ModelSize: "
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

bool Model::makeMultiForwardStep (
	int PValueBorder,
	const int selectionCriterium,
	int *startIndex,
        TBitset * exclusivedSNP,
        TBitset * goodSNPs
) {
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

	if (goodSNPs != 0) {   //  GA
		if ( false == (*goodSNPs)[data_->getOrderedSNP(i)] ) {
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
              "ADD    SNP: " + data_->getSNP(orderedSnpIndex ).getSnpId() + " ModelSize: "
	      + int2str( NewModel.getModelSize() ) + " MSC: " + double2str(newBIC)
            );
          }
          oldBIC = newBIC;
          if (exclusivedSNP != 0)    // GA
            (*exclusivedSNP)[data_->getOrderedSNP(i)] = true;  
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
      
      if (goodSNPs != 0)
      {
        if (orderedSnpIndex < 0 || orderedSnpIndex > goodSNPs->size())
        {
          cerr << "ERROR: orderedSnpIndex: " << orderedSnpIndex << endl;
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
	AutoVector xVec( data_->getIdvNo() );
	data_->getXcolumn( orderedSnpIndex, xVec );
      qruncher.pushColumn( xVec );
      
      newBIC = computeMSC( selectionCriterium, qruncher.calculateRSS() ); //1 instead of selection criterion
      // TODO<BB>: Here we should store the RSS in a ResultStore to avoid duplicate calculation
      
      if ( newBIC < oldBIC ) {
 printLOG("newBIC________________"+double2str(newBIC));
        addSNPtoModel( orderedSnpIndex );
      
	if (exclusivedSNP != 0)  // GA
	{
		(*exclusivedSNP)[orderedSnpIndex] = true;   
	}  
        oldBIC = newBIC;
	if ( ::isinf( oldBIC ) && oldBIC < 0.0 ) {
   		printLOG( "model fully explains observations" );
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
			if ( NTest.computeMSC( selectionCriterium ) < compareMSC ) {
				NMax = NTest;
				compareMSC = NTest.computeMSC( selectionCriterium );
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
bool Model::saveguardbackwardstep (
	Model &smallerModel,
	const int selectionCriterium
) {
Model backwardModel(*data_);
Model model3(*data_);
model3=*this;

//this is for giving back the original model instead of the smallest one
double bestMSC=getMSC(); //from
printLOG(" bestMSC="+double2str(bestMSC));
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
// model3.printModel();
 
 }
 //else set the counter up and copy the backwardModel to *this
 else
 {//what if the criterion is positive?
	 if (bestMSC>0)
		; //do nothing this could happen when you change from mBIC with 45 to mBIC2/, and should only happen, when you start modelselection again with 
		 //a stronger criterion
	 else	 
	 {breakfor++;
  if (breakfor>=parameter.saveguardsteps) //here one could set any value %the number here was 2
 	{break; }
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
				local=NTest.computeMSC( selectionCriterium );
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
		size=data_->getSnpNo()-1,//otherwise 1 to big 
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
	int 	n = data_->getIdvNo();	// # of rows
	int 	k = getNoOfVariables();	// # of columns
	//double 	x[n*k]; 		// an array representing XMat_
	int		y[n];		// an array representing YVec_ // implicit cast! YVec_ is double, but should only contain 0, 1
	double 	beta_array[k];		// the coefficicents
	double 	pi[n];
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
	}
	
	for (j=0; j < k; j++)
	{
		beta_array[j] = gsl_vector_get(betas_,j); 	//  initialize according to old model //intialize beta as 0
	//	cout<<beta_array[j]<<",";  //DEBUG
	}
               //cout<<endl; DEBUG
	// logistic firth regression 
	if ( 
		logistffit(
					// Input
					k, // # of rows of X (# of variables) BB: I'd rather call that "columns"
					n, // # of columns of X (# of samples) BB: I'd rather call that "rows"
					XMat_,//x, // 
					y,  //
					// Output
					beta_array, 
					pi, //
					H, // 
					&loglik, 
					&iter, // # of main iterations;
					&evals, 
					&lchange, // the change in the log-likelihood between steps
					&ret_max_U_star, //
					&ret_max_delta,			
					//// optional Input
					true,	//  use firth-regression
					beta_array,	// initial values for beta
					false,	// only evaluate likelihood
					//// Control Parameters
					parameter.logrC_maxit,	
					parameter.logrC_maxhs,	
					parameter.logrC_maxstep,
					parameter.logrC_lconv,
					parameter.logrC_gconv,
					parameter.logrC_xconv	
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
	// Initialisation

	// check if we have an 1-SNP model
	if ( 1 != getModelSize() ) {
		modelSnps_.clear();
		modelSnps_.push_back( snp );
		initializeModel();
	}
	else // replace in 1-SNP model the old SNP with new SNP snp
	{
		const size_t idvs = data_->getIdvNo();
		AutoVector xVec( idvs );
		modelSnps_.at(0) = snp;
		data_->getXcolumn( snp, xVec );
		// TODO<BB>: Use vectorial replacement instead of loop
		for ( size_t i = 0; i < idvs; ++i ) {
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
	// REMARK<BB>: The following code is similar to that in Model::computeLinRegression(),
	// but uses LU and the explicit Inverse.
	// TODO: Unify both approaches and use QRuncher.

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
			printLOG("Error single linear regressor test for SNP "+ (data_->getSNP( snp )).getSnpId()+ ": Regressionmatrix not invertable!" );
		
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
		// TODO<BB>: it would be more efficient not to calculate all P values for cut-off,
		// but rather sort by critical value and use the inverse of the distribution function to cut off.
	}

}



double Model::computeSingleLogRegressorTest ( const snp_index_t snp ) {
	// compute the genotype frequency table of SNP snp
	GenotypeFreq freqOneSNP( *data_, snp );
	
	// choose test for single marker test for case-control data 
	switch ( parameter.singleMarkerTest ) {
		case Parameter::singleMarkerTest_CHI_SQUARE:
			return freqOneSNP.calculateChiSquare();
		case Parameter::singleMarkerTest_COCHRAN_ARMITAGE:
			return freqOneSNP.calculateCATT();
		default:
			throw Exception(
				"Implementation error: Unrecognised choice %d for single marker test."
			);
	}
}

double Model::computeMSC ( const int selectionCriterium ) {
	return computeMSC( selectionCriterium, getMJC() );
}

double Model::computeMSC ( const int selectionCriterium, double mjc ) {
	if ( !upToDateBetas_ ) // check if betas_ and therefor also MJC is up-to-date, otherwise update
	{
		computeRegression();
	}

	const size_t
		n = data_->getIdvNo(),
		p = parameter.nSNPKriterium,	// data_->getSnpNo(); das ist original
		q = getModelSize();

	// choose the likelihood part depending if the Data is quantitative or affection
	double LRT, d;
	if ( parameter.affection_status_phenotype ) {
		LRT = (-2.0)*( mjc - data_->getLL0M()); // LRT = -2 log (likelihood(0-model)/likelihood(model)) the inverse is the right 
	}
	else
	{
		LRT = n*log( mjc ) ;  // n * log(RSS)
	}
	
if(DEBUG) cerr<<"LRT="<<LRT<<endl;

	switch ( selectionCriterium ) {
		case Parameter::selectionCriterium_BIC:
			msc = LRT + q * log(n);
			return msc;

		case Parameter::selectionCriterium_EBIC:
			if ( ::isnan( parameter.EBIC_gamma ) ) {
				// How the default value comes about, see Zhao Chen:
				// "A feature selection approach to case-control genome-wide association studies"
				parameter.EBIC_gamma = 1.0 - log( n ) / ( 2 * log( p ) );
			}
			msc = LRT + q * log(n) + 2 * parameter.EBIC_gamma * logFactorial.logChoose( p, q );
			return msc;

		case Parameter::selectionCriterium_mBIC_firstRound:
		case Parameter::selectionCriterium_mBIC:
			d = -2 * log(
				Parameter::selectionCriterium_mBIC_firstRound == selectionCriterium
				? parameter.mBIC_firstRound_expectedCausalSNPs
				: parameter.mBIC_expectedCausalSNPs
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
  clearModel();     // I put lines below to clearModel()
}

void Model::printModelNew() const {
	const size_t idvs = data_->getIdvNo();
	ofstream	SNPL;
	
	SNPL.exceptions ( ofstream::eofbit | ofstream::failbit | ofstream::badbit ); // checks if Logfile can be writt
	try
	{
		SNPL.open( ( parameter.out_file_name + "_SNPListNew.txt" ).c_str(),  fstream::out );

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

		for ( size_t i = 0; i < data_->getCovNo(); ++i ) {
			SNPL << data_->getCovMatElementName(i) << " <- c(";
			data_->getCovariateColumn( i, vec );
			for ( size_t idv = 0; idv < idvs; ++idv ) {
				if ( 0 < idv ) SNPL << ",";
				SNPL << vec.get( idv );
			}
			SNPL<< ")"<< endl;	
		}

		const Vector yVec( *YVec_ );
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
  vector<snp_index_t> v = m.modelSnps_;
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

double Model::computeMSCfalseRegression ( const int selectionCriterium ) {
	vector<snp_index_t> removedSnps;
	msc = computeMSCfalseRegression( selectionCriterium, removedSnps );
	return msc;
}
  
double Model::computeMSCfalseRegression (
	const int selectionCriterium,
	vector<snp_index_t> &removedSnps
) {
  if (computeRegression() == false)
  {
    removedSnps.push_back(modelSnps_[modelSnps_.size() - 1]);
    if (removeSNPfromModel(modelSnps_.size() - 1) == false)
      cout << "removeSNPfromModel failed!" << endl;
    //cout << "try computeMSC for model (-1): " << *this << endl;
    //char cc; cout << "Press a key... "; cin >> cc;
		msc = computeMSCfalseRegression( selectionCriterium, removedSnps );
    
    snp_index_t snp = removedSnps.back();
    
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
	int PValueBorder,
	const int maxModel,
	const int selectionCriterium
) {
	printLOG( "Score select model" );
	if ( getModelSize() > parameter.maximalModelSize ) return false;

	double best= getMSC();
	//if best is bigger than 0
	//than a smaller model is better or at least the 0 model which is empty
	 Model 	*backwardModel;
	 backwardModel=this;
//ORIGINAL value is now default   int PValueBorder =100,
	int *startIndex;
 	int dummy=0;
	startIndex=&dummy;
	PValueBorder =min(PValueBorder,400); //override to high PValueBorders !!!
 	double bestMSC=getMSC(); //this one is the best up to now!

 ScoreTestShortcut stsc( *data_);
 	int size = data_->getSnpNo();
	bool improvement = true, improvement2 = true, improvement3 = true;
	SortVec score(size);
 this->computeRegression();//for SCORE maybe not calculated?
        stsc.ScoreTestShortcut::scoreTests ( *this, score );
 vector<int> Score(size);
        for(int i=0;i<=size-getNoOfVariables();i++)
 		Score[i]= score.getId(i);	

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

  vector<snp_index_t> v1 = modelSnps_;
  sort(v1.begin(), v1.end());

  vector<snp_index_t> v2 = m.modelSnps_;
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
