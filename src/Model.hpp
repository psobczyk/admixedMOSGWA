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

#ifndef MODEL_HPP
#define MODEL_HPP

#include <string>
#include <iostream>
#include "Parameter.hpp"
#include "Helpfull.hpp"
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
class Model {

	friend std::ostream &operator << ( std::ostream &out, const Model &m );
  
private:

	/** a pointer to the MData whereupon the model is based */
	const MData* data_;

	/** the position-no of the SNPs in the Model (for faster access), size is modelSize_ */
	vector<size_t> modelSnps_;

	/** the design matrix X for the model (data.getIdvNo() rows, noOfVariables columns) */
	linalg::AutoMatrix xMat;
public:
	/** the target vector Y (size is data.getIdvNo()) */
	linalg::AutoVector yVec;
private:
	/** the regression coefficients (size is noOfVariables_) */
	linalg::AutoVector beta;
public:
	/** states if betas_ is up-to-date */
	bool upToDateBetas_;

	/** Residual sum of squares (RSS) for quantitive traits, log-likelihood for case-control */
	double modelJudgingCriterion_;

	/** For initializing GA population. Set of used SNPs. SNPs may belong to one model */

	/** Helps calculating factorials and combinations.
	* TODO: do it in a thread-safe way.
	*/
	util::LogFactorial logFactorial;

	/** msc of the model. calculateMSC() set up msc value */
	double msc;                         //  FOR GA. 
	
	/** load data from MData (plink-format) to Model (gsl-format) */
	void initializeModel ();

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

	bool computeLogRegression ();

	/** Calculate P value of a 1 SNP model.
	* Precondition: Model holds exactly one SNP.
	*/
	double computeSingleLinRegressorTest ();

	/** Calculate selection value of a 1 SNP model.
	* Precondition: Model holds exactly one SNP.
	*/
       	double computeSingleLogRegressorTest ();

//++++++++++++
// public methods:
public:
	/** Constructor: Model must be base on some MData */
	Model ( const MData & mData );

	/** Assignment operator */
	DEPRECATED( Model& operator= ( const Model & orig ); );

	/** Destructor */
	~Model ();

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

	void printModelNew () const;

	bool replaceModelSNPbyNearFromCAT (
		int currentPosition,
		const size_t PValueBorder,
		const int selectionCriterium = Parameter::selectionCriterium_mBIC2
	);

	/** Try to replace using conditional score test for each model SNP removed */
	bool replaceModelSNPSCORE ( const int selectionCriterium );

	/** outputs for every SNP in the model SNPs with a correlation above const double threshold */
	void printStronglyCorrelatedSnps ( const double threshold, string extra ="" ) const;
	/* sortSNPsAccordingBetas does exactly what is says*/
	void sortSNPsAccordingBetas();
	// change model Size
	/** add SNP to Model, snp is absolute postion in MData::snps_ */
	void addSNPtoModel ( const size_t snp );
      bool replaceSNPinModel ( const size_t snp,  const size_t position );
		void addManySNP ( std::vector<size_t> selected );
  
	/** Continous Y
	* @deprecated test data generation will be removed from class Model
	*/
	DEPRECATED( void Ycontinous(); );
  
	/** remove SNP from Model, snp is relativ position at vector modelSnps_
	* @returns true if SNP removed, and false if an error occors */
	bool removeSNPfromModel ( const size_t snp );

	/** adds the SNPs minimzing the MSC.
	* Model &biggerModel is the new model, boundSNP is a border on the SNPs to add
	* @returns the absolut positon of the added SNP, or -1 on error */
	int makeForwardStep ( Model &biggerModel, const int boundSNP, const int selectionCriterium );
 
	/** computes the Best Model with a SNP smaller,
	* Model &smallerModel is the new model
	* @returns the relativ position of the removed SNP, or -1 on error */
	int makeBackwardStep ( Model &smallerModel );
	/*special Version of the multiforward step*/
	bool finalizeModelSelection (
		Model &backwardModel,
		bool improvement,
		const size_t PValueBorder,
		int *startIndex,
		vector<int> score,
		const int selectionCriterium = Parameter::selectionCriterium_mBIC2
	);
	bool finalizeModelSelection (
		Model &backwardModel,
		bool improvment,
		const size_t PValueBorder,
		int *startIndex,
		const int selectionCriterium = Parameter::selectionCriterium_mBIC2
	);
	bool makeForwardStepLinear (
		Model *forwardModel,
		double* bestMSC,
		const size_t PValueBorder,
		int *startIndex,
		const int selectionCriterium = Parameter::selectionCriterium_mBIC2
	);
	/**  makeForwardStepLogistic replaces the code in selectModel*/
	bool makeForwardStepLogistic (
		double *bestMSC,
		const size_t PValueBorder,
		int *startIndex,
		const int selectionCriterium = Parameter::selectionCriterium_mBIC2
	);
/**makeForwardStepLogistic score version */
	bool makeForwardStepLogistic (
		double *bestMSC,
		const size_t PValueBorder,
		int *startIndex,
		vector<int> score,
		int selectionCriterium = Parameter::selectionCriterium_mBIC2
	);

        bool makeMFFS (
		const size_t PValueBorder,
		int *startIndex
	); //with linear

	bool makeMFFL (
		const size_t PValueBorder,
		int *startIndex,
		const int selectionCriterium = Parameter::selectionCriterium_mBIC2
	);

	bool makeMFFL (
		const size_t PValueBorder,
		int *startIndex,
		vector<int> score,
		int criterium = 0
	); //with logistic+score

        bool makeMultiForwardStepScore (
		size_t PValueBorder,
		int selectionCriterium,
		int* startIndex,
		vector<int> scores
	);

        /** makeMultiForwardStep take the PValueBorder, an selection Criterium, the default Value is 1 for BIC in the initial ForwardStep and an  an exclusivedSNP set,*/ 
	bool makeMultiForwardStep (
		size_t PValueBorder = 0,
		int selectionCriterium = 1,
		int *startIndex = NULL,
		TBitset * exclusivedSNP = 0,
		TBitset * goodSNPs = 0
	);

	bool makeMultiBackwardStep ();

	bool saveguardbackwardstep (
		Model &smallerModel,
		const int selectionCriterium = Parameter::selectionCriterium_mBIC2
	);
	/** makeBackwardStepED  a variation of makeBackwardStep */
        int makeBackwardStepED ( Model &smallerModel, const int criterium=0 );

	bool selectModel (
		Model &startFromModel,
		size_t PValueBorder = 100,
		int maxModel=parameter->maximalModelSize,
		const int selectionCriterium = Parameter::selectionCriterium_mBIC2
	);
	/** Compute Regression for Model.
	* @returns false on error */
	bool computeRegression ();

        /** scoreTest runs the score tests from Bernhard */
	void scoreTest ( const string& extra = "" );

	size_t scoreTestWithOneSNPless ( size_t SNPindex, SortVec &score );
        
	/** computes single marker test for a one-SNP model.
	* Precondition: <code>this</code> holds exactly one SNP.
	* @returns the p-value of the according teststatistic */
	double computeSingleRegressorTest ();

	/** print the yVec from Model to file. */
	void printYvec ();

	/** Calculate the model selection criterion. */
	double computeMSC ( const int selectionCriterium = Parameter::selectionCriterium_mBIC2 );

	/** Calculate the model selection criterion from the model judging criterion.
	* @see computeMSC ( int )
	*/
	double computeMSC ( const int selectionCriterium, double mjc );

	/** @brief Returns msc of the model. Remember to use calculateMSC() before getMSC() (for GA)
      @returns msc of the model
  */
	double getMSC() const;

	/** @brief Gets the snp at a given position (for GA)
      @param pos - snp position
  */
	size_t getSNPat( const size_t pos ) const;

	/** @brief Returns snps vector of the model (for GA). */
	std::vector<size_t> getModelSnps() const;

	/** @brief Creates new model from given snps, snps are absolut postion in MData::snps_ (for GA) */
	void createFromSNPs ( const std::set<size_t>& snps );

	/** @brief Dealocates memory and makes the object ready to create new model (for GA)*/
	void clearModel();

	/** @brief Removes oneSNP form Model, oneSNP is SNP value (for GA)
  * @returns value = true if SNP removed, and false if an error occors 
  */
	bool removeSNPValFromModel( const size_t oneSNP );
  
  
  /** @brief Computes MSC for model with correlated snps. Use if computeRegression() returns false
   *TODO This function is not optimised
   */
	double computeMSCfalseRegression ( const int selectionCriterium = Parameter::selectionCriterium_mBIC2 );
  
	double computeMSCfalseRegression ( const int selectionCriterium, std::vector<size_t> &removedSnps );

	bool operator == ( const Model &m ) const;

	bool operator != ( const Model &m ) const { return !(*this == m); }
};

#endif	/* MODEL_HPP */
