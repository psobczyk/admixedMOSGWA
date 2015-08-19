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

#ifndef MDATA_HPP
#define MDATA_HPP

#include <string>
#include <vector>
#include <memory>

#include "Parameter.hpp"
#include "Helpfull.hpp"
#include "SortVec.hpp"
#include "Individual.hpp"
#include "SNP.hpp"
#include "linalg/AutoVector.hpp"
#include "linalg/AutoMatrix.hpp"
#include "io/Input.hpp"
#include "io/CachedInput.hpp"
#include "deprecate.h"

#include "Model.hpp" // ag
class Model;

using namespace std;	// TODO<BB>: remove namespace import in *.hpp

/** Stores all the information from the input plink files.
* Data is read and access methodes are provided.
* Single marker test and model selection can be performed.
*/
class MData {

	/** Where to fetch input data from. E.g. Genotype data is read only on demand. */
	io::Input *input;

	/** Input cache; only used for genotype vectors. */
	std::auto_ptr<io::CachedInput> inputCache;

	/** Whether the input was given as <code>NULL</code> and thus has been allocated. */
	const bool allocateInput;

	/** the name of the target trait */
	string Y_name_;
	
	/** the names of the Covariates */
	vector <string> covNames;

	/** Matrix where the covariables are stored.
	* It has parameter.covariables columns and getIdvNo() rows.
	*/
	linalg::AutoMatrix covMat;

	/** the number of cases (for affection phenotype) */
	size_t caseNo_;

	/** the number of controls (for affection phenotype) */
	size_t contNo_;

	/** Holds the single marker tests for all SNPs */
	vector <double> singleMarkerTestResult;

	/** the SNPs, sorted according to ascending single marker tests */
	SortVec snp_order_;

	/** the log-likelihood of the zero-model (without any SNPs) */
	double loglikelihood0Model_;

	/** Check input data */
	void checkData ();

	/** Remove the given individual from the data. */
	void removeIndividual ( const size_t idv );

public:
	// not private anymore because of setting Y values on the fly
	/** for affection-type targets, stores the log-likelihood of the zero-model */
	void setLL0M ( const double ll );

	/** Default Constructor: reads the input-files in, sets parameters, deals with missing phenotypes.
	* @param input provides access to input data for the transition of the MOSGWA architecture.
	* If given, the pointer must be valid as long as the <code>MData</code> is alive.
	* Otherwise, input is read according to the preference settings in {@link Parameter}.
	* @see Parameter
	*/
	MData ( io::Input * const externalInput = NULL );

	/** Destructor: clean up */
	~MData();

	/** get the /number of SNPs in the data (the number of variables) */
	size_t getSnpNo () const;

	/** get the number of individuals in the data (the sample size) */
	size_t getIdvNo () const;
        std::string getFID ( const size_t index ) const;
        std::string getID ( const size_t index ) const;
	/** Access a column of the regression matrix.
	* It has {@link MData::getIdvNo()} rows and {@link MData::getSnpNo()} columns.
	* The entries are -1 (negative homozygote), 0 (heterozygote), +1 (positive homozygote)
	* for the respective individual and SNP, or <code>NaN</code> for missing data.
	* The returned matrix column provides a snapshot of currently stored data.
	* Its behaviour becomes undefined
	* if the content of MData is changed after the call of <code>getXcolumn()</code>.
	* WARNING: Modifying the X matrix will most likely invalidate already calculated model selection criteria.
	*/
	void getXcolumn ( const size_t snp, linalg::Vector& vector ) const;

	/** Get the regression target vector.
	* The returned Vector provides a snapshot.
	* Its behaviour becomes undefined if the content of MData is changed after the call of getY().
	* WARNING: Modifying the Y vector will most likely invalidate already calculated model selection criteria.
	*/
	void getY ( linalg::Vector& vector ) const;

	/** for an integer snp a pointer to the snp-th SNP is returned */
	const SNP & getSNP ( const size_t snp ) const;

	/** For an integer i, return the (absolute) postition of the SNP with the snp-th lowest p-value. */
	size_t getOrderedSNP ( const size_t snp ) const;

	/** Get the number of covariate vectors in the data. */
	size_t getCovNo () const;

	/** Access a column of the covariate matrix.
	* The returned <code>Vector</code> provides a snapshot.
	* Its behaviour becomes undefined if the content of MData is changed after the call.
	* WARNING: Modifying the matrix will most likely invalidate already calculated model selection criteria.
	*/
	void getCovariateColumn ( const size_t cov, linalg::Vector& vector ) const;

	/** return the name of the cov-th covaribale */
	const string& getCovMatElementName( const size_t cov ) const;

	// for affection-type targets
	/** # of cases in for affection-type target (should not be used for quantitive traits, intialised with 0) */
	size_t getCaseNo () const;

	/** # of controls in for affection-type target (should not be used for quantitive traits, intialised with 0) */
	size_t getContNo () const;

	/** returns the log-likelihood for affection-type targets.
	* @deprecated A {@link search::Search} algorithm should calculate model regressions itself.
	* This is not the role of MData.
	* Also mind using {@link lookup::ResultStore}.
	*/
	DEPRECATED( double getLL0M () const );

//+++++++++++++

	// general methods
	/** Computes the correlation between two SNPs (loci) */
	double computeCorrelation ( const size_t locus1, const size_t locus2 ) const;

	// model selection
	/** single marker test are performed, and saved (in object SNP).
		The SNPs are sorted w.r.t ascending p-values in SortVec snp_order_.
		The ordered SNPs are output into file "*outname*_IT.txt"
	* @deprecated A {@link search::Search} algorithm should calculate model regressions itself.
	* Dealing with {@link Model}s is not the role of MData.
	* Also mind using {@link lookup::ResultStore}.
	*/
	DEPRECATED( void calculateIndividualTests() );
/**calculate PValuePorder needs parameter.ms_MaximalPValueForwardStep but set in the conf-file
 *and of course the MData variable.
 determines the SNPs with p-Value < parameter.ms_MaximalPValueForwardStep
 these are tested in the Forward Step
 */

	size_t calculatePValueBorder() const;

	/** model selection is performed
	* @deprecated Formulate this as a {@link search::Search} algorithm.
	* Dealing with {@link Model}s is not the role of MData.
	* PValueBorder could now set directly it should be 0<PValueBorder<=nSNPs
	* maxModel is usefull to limit the mximal modelsize, respect to an older call of selecModel with different
	*  ExpectedCausalSNPs (this is a variable used with mBIC2
	*/
	DEPRECATED( bool selectModel() );
	bool selectModel (
		Model * inputModel,
		size_t PValueBorder,
		int maxModel = parameter->maximalModelSize,
		const int selectionCriterium = Parameter::selectionCriterium_mBIC2
	);
	//+++++++++++++
	// methods for testing CAUTION using print_____ on big datasets
	/** TESTING // screen-output of the informations stored for the SNPs */
	void printSNPs ();

	/** a special version for Artur */
	/** TESTING // given a list of SNP-Names, the Genotype of the SNP
	 is extracted, the SNPs, Covariables, Intercept are written to
	 an R-linear model. Was used to check correctness of linear regression
	 but may be helpful in general to extract all the date of a model... */
	void printSelectedSNPsInR ( vector<string> SNPList ) const;

	/** The same in Matlab extra is the number in the yvm file */
	void printSelectedSNPsInMatlab ( vector<string> SNPList, string extra ) const;

  /** Artur's new code:
   * @brief Gets snp_order value. It's needed for multiForwardStep
   * @param i - position 
   */
  double getSnp_order_Value ( const size_t i ) const;
  
	/** Access the Single Marker Test result for the given SNP. */
	double getSingleMarkerTestAt ( const size_t index ) const;

  /** Artur's new code:
   * @brief sets single markert test for i-th snp with value
   * @param i - position of the snps_
   * @param value - value of single market test
   */
  void setSingleMarkerTestAt ( const size_t index, const double value );
  
  /** Artur's new code:
   * @brief fills vector
   */
  void fillSnp_order_Vec ( const size_t snpNo, size_t* SNPList, double* TestStat );
  
  
  /** Artur's new code:
   * @brief writes ordered snps to file. File name is set up in parameters class.
   */
  void printSnpOrder();
  //that should work
// template<class Vals> void  sortingPermutation(const Vals& values, std::vector<int>& v);
  /** findSNPIndex is here to add many SNP from a list of SNPNames instead of a vector of int's */
  void findSNPIndex ( vector<string>& SNPNames, vector<size_t>& index ) const;
};
#endif
