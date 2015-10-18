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

#include "AModel.hpp"
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

#include <algorithm> //min und max


using namespace std;
using namespace linalg;
using namespace lookup;
using namespace logging;

////++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//  class Model
////++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



AModel::AModel ( const AMdata& mData )
:
	data_( &mData ),
	xMat( mData.getIdvNo(), 1 + mData.getCovNo() ),
	cMat( mData.getIdvNo(), 1 + mData.getCnpNo() ),
	yVec( mData.getIdvNo() ),
	beta( 1 + mData.getCovNo() ),
	gamma( 1 + mData.getCnpNo() )
{
	const size_t covs = data_->getCovNo();
	size_t col = 0;
	xMat.columnVector( col++ ).fill( 1.0 );		// intercept is a column of 1
	for ( size_t cov = 0; cov < covs; ++cov ) {
		Vector xVec = xMat.columnVector( col++ );
		data_->getCovariateColumn( cov, xVec );
	}
	for ( size_t cov = 0; cov < covs; ++cov ) {
		Vector zVec = cMat.columnVector( cov );
		data_->getCcolumn( cov, zVec );
	}
	data_->getY( yVec );
	initializeModel();
}


AModel::~AModel () {}

void AModel::initializeModel () {
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
	gamma.upSize( cols );
	gamma.fill( 0.0 );
	upToDateBetas_ = false;
	upToDateGammas_ = false;
}

// getters
ModelIndex AModel::getIndex () const {
	return ModelIndex( modelSnps_ );
}

int AModel::getModelSize () const {
	return modelSnps_.size();
}

int AModel::getNoOfVariables () const {
	return 1 + data_->getCovNo() + getModelSize();
}

double AModel::getMJC () const {
	return modelJudgingCriterion_;
}

double AModel::getBeta ( const int i ) const {
	if ( i >= getNoOfVariables() ) {
		throw;
	} else {
		return beta.get( i );
	}	
}

double AModel::getGamma ( const int i ) const {
	if ( i >= getNoOfVariables() ) {
		throw;
	} else {
		return gamma.get( i );
	}	
}

string 	AModel::getSNPId ( const size_t i ) const {
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

void AModel::sortSNPsAccordingBetas () {
  vector<size_t> SNP( getNoOfVariables(), 0 );
  SortVec sbetas(getNoOfVariables());
  vector<double> betas(getNoOfVariables(),0);
 for (int i=0;i<getNoOfVariables();i++) //this assumes that getNoOfVariables	
 { betas[i]=getBeta(i);
   SNP[i] = modelSnps_[i];
 }
sbetas.fillVec(getNoOfVariables(), &SNP[0],&betas[0],false);
}

void AModel::setBeta ( const Vector& newBeta ) {
	beta.copy( newBeta );
	upToDateBetas_ = true;
}



bool AModel::operator == ( const AModel &m ) const {
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
