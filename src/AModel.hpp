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

#ifndef AMODEL_HPP
#define AMODEL_HPP

#include <string>
#include <iostream>
#include "Parameter.hpp"
#include "Helpfull.hpp"
#include "AMdata.hpp"
#include "MData.hpp"
#include "lookup/ModelIndex.hpp"
#include "linalg/AutoVector.hpp"
#include "linalg/AutoMatrix.hpp"
#include "util/LogFactorial.hpp"

#include <set>      
#include <iterator> 

class MData;        

/** For a given MData set, Model implements a statistic model on MData,
* that is a subset of SNPs of MData.
* To handle this set of SNPs further information is stored.
* SNPs can be added and removed from the model, and properties of the model can be computed.
*
* TODO<BB>: Split class Model into subclasses for linear and logistic model
* so that the user of a Model will always get a subclass optimised for the specific algorithm.
* E.g. linear models could fully exploit {@link QRuncher} capabilities
* without a need to carry and copy matrices needed for logistic regression.
*/
class AModel {

	friend std::ostream &operator << ( std::ostream &out, const AModel &m );
  
private:
	/** a pointer to the AMData whereupon the model is based */
	const AMdata* data_;

	/** the position-no of the SNPs in the Model (for faster access), size is modelSize_ */
	std::vector<size_t> modelSnps_;

	/** the design matrix X for the model (data.getIdvNo() rows, noOfVariables columns) */
	linalg::AutoMatrix xMat;

	/** the ancestry matrix C for the model (data.getIdvNo() rows, noOfVariables columns) */
	linalg::AutoMatrix cMat;

	/** the target vector Y (size is data.getIdvNo()) */
	linalg::AutoVector yVec;

	/** the regression coefficients (size is noOfVariables_) */
	linalg::AutoVector beta;

	/** the regression coefficients for ancestry matrix (size is noOfVariables_) */
	linalg::AutoVector gamma;

	/** states if betas_ is up-to-date */
	bool upToDateBetas_;

        /** states if betas_ is up-to-date */
	bool upToDateGammas_;

	/** Residual sum of squares (RSS) for quantitive traits, log-likelihood for case-control */
	double modelJudgingCriterion_;

	/** For initializing GA population. Set of used SNPs. SNPs may belong to one model */

	/** Helps calculating factorials and combinations.
	* TODO: do it in a thread-safe way.
	*/
	util::LogFactorial logFactorial;

	/** msc of the model. calculateMSC() set up msc value */
	double msc;                         //  FOR GA. 

// public methods:	
public:
  /** Constructor: Model must be base on some MData */
  AModel ( const AMdata & mData );
  
  /** Destructor */
  ~AModel ();
  

  void initializeModel ();

  void sortSNPsAccordingBetas ();
	/** Suggest an optimal SNP to be added.
	* Precondition: The Model does not yet contain all SNPs below <code>bound</code>.
	* @param snp points to a variable to store the selected SNP index.
	* @returns the residual sum of squares for the best of the considered models.
	*/
	double oraculateOptimalLinearForwardStep( size_t *snp, size_t bound ) const;


	/** Suggest an optimal SNP to be removed from the model.
	* Precondition: The Model contains at least one SNP.
	* @param snp points to a variable to store the selected SNP index.
	* @returns the residual sum of squares for the best of the considered models.
	*/
	double oraculateOptimalLinearBackwardStep( size_t *snp ) const;

	// regression for model depending on the type
	bool computeLinRegression ();

	
	/** Calculate P value of a 1 SNP model.
	* Precondition: Model holds exactly one SNP.
	*/
	double computeSingleLinRegressorTest ();

	/** Calculate P value of a 1 ZNP model.
	* Precondition: Model holds exactly one ZNP.
	*/
	double computeSingleLinAncestryRegressorTest ();

	// Getters

	/** Retrieve the multi-index identifying the model.
	* It contains exactly the indices (in {@link MData}) of the chosen SNPs, i.e. regression variables.
	*/
	lookup::ModelIndex getIndex () const;

	int getModelSize () const;

	int getNoOfVariables () const;

	/** Get the model judging criterion. */
	double getMJC () const;

	/** Get the beta of the model */
	double getBeta ( const int i ) const;

	/** Get the gamma of the model */
	double getGamma ( const int i ) const;

	/** get the name of a SNP by its relativ position in Model */
	std::string getSNPId ( const size_t i ) const;

        // Setters
        
        /** set \beta could be used for generating Y
            the \betas should be  gsl_vector* Beta     
          */

	void setBeta ( const linalg::Vector& newBeta );

	// Output
	/** outputs the Model information, string out is a title */
	void printModel (
		const std::string& out = "",
		const int selectionCriterium = Parameter::selectionCriterium_mBIC2,
		const std::string& filemodifier = ""
	);


	bool selectModel (
		AModel &startFromModel,
		size_t PValueBorder = 100,
		const int maxMode
l=parameter->maximalModelSize,
		const int selectionCriterium = Parameter::selectionCriterium_mBIC2
	);

	bool operator == ( const AModel &m ) const;

	bool operator != ( const AModel &m ) const { return !(*this == m); }
};

#endif	/* AMODEL_HPP */
