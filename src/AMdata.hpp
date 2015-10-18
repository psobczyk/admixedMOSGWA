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

#ifndef AMDATA_HPP
#define AMDATA_HPP

#include <string>
#include "Individual.hpp"
#include "SNP.hpp"
#include "linalg/Vector.hpp"

/** Stores all the information required for an admixture GWAS.
* Data is read and access methodes are provided.
*/
class AMdata {

public:
  /** Default Constructor: reads the input-files in, sets parameters, deals with missing phenotypes.
   * @param input provides access to input data for the transition of the MOSGWA architecture.
   * If given, the pointer must be valid as long as the <code>MData</code> is alive.
   * Otherwise, input is read according to the preference settings in {@link Parameter}.
   * @see Parameter
   */
  AMdata ();
  
  /** Destructor: clean up */
  ~AMdata ();
  
  /** get the number of individuals in the data (the sample size) */
  size_t getIdvNo () const;
  
  /** get the index-th individual */
	Individual getIndividual ( const size_t index ) const;

	/** get the /number of SNPs in the data (the number of variables) */
	size_t getSnpNo () const;

	/** for an integer snp a pointer to the snp-th SNP is returned */
	const SNP & getSNP ( const size_t snp ) const;

	/** Access a column of the regression matrix.
	* It has {@link MData::getIdvNo()} rows and {@link MData::getSnpNo()} columns.
	* The entries are -1 (negative homozygote), 0 (heterozygote), +1 (positive homozygote)
	* for the respective individual and SNP, or <code>NaN</code> for missing data.
	* The returned matrix column provides a snapshot of currently stored data.
	* Its behaviour becomes undefined
	* if the content of MData is changed after the call of <code>getXcolumn()</code>.
	*/
	void getXcolumn ( const size_t snp, linalg::Vector& vector ) const;

	/** get the /number of CNPs in the data (the number of C matrix variables) */
	size_t getCnpNo () const;

	/** for an integer snp a pointer to the snp-th CNP is returned */
	const SNP & getCNP ( const size_t snp ) const;

	void getCcolumn ( const size_t snp, linalg::Vector& vector ) const;

	/** Get the number of covariate vectors in the data. */
	size_t getCovNo () const;

	/** return the name of the cov-th covaribale */
	const std::string& getCovMatElementName( const size_t cov ) const;

	/** Access a column of the covariate matrix.
	* The returned <code>Vector</code> provides a snapshot.
	* Its behaviour becomes undefined if the content of MData is changed after the call.
	*/
	void getCovariateColumn ( const size_t cov, linalg::Vector& vector ) const;

	/** Get the regression target vector.
	* The returned Vector provides a snapshot.
	* Its behaviour becomes undefined if the content of MData is changed after the call of getY().
	*/
	void getY ( linalg::Vector& vector ) const;
};
#endif
