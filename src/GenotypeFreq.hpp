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

#ifndef GENOTYPEFREQ_HPP
#define GENOTYPEFREQ_HPP

#include "MData.hpp"

/** The genotype frequency table for case-control phenotypes. */
class GenotypeFreq {

private:
	size_t
		/** case with genotype "0" */
		r0_,
		/** case with genotype "1" */
		r1_,
		/** case with genotype "2" */
		r2_,
		/** # of cases */
		r_,
		/** control with genotype "0" */
		s0_,
		/** control with genotype "1" */
		s1_,
		/** control with genotype "2" */
		s2_,
		/** # of controls */
		s_,
		/** total # of genotype "0" */
		n0_,
		/** total # of genotype "1" */
		n1_,
		/** total # of genotype "2" */
		n2_;
		/** # of individuals */

	const size_t n_;

public:
	/** Number of cases with genotype 0 */
	size_t get_r0 () const { return r0_; }

	/** Number of cases with genotype 1 */
	size_t get_r1 () const { return r1_; }

	/** Number of cases with genotype 2 */
	size_t get_r2 () const { return r2_; }

	/** Number of cases total */
	size_t get_r () const { return r_; }

	/** Number of controls with genotype 0 */
	size_t get_s0 () const { return s0_; }

	/** Number of controls with genotype 1 */
	size_t get_s1 () const { return s1_; }

	/** Number of controls with genotype 2 */
	size_t get_s2 () const { return s2_; }

	/** Number of controls total */
	size_t get_s () const { return s_; }

	/** Total number of genotype 0 */
	size_t get_n0 () const { return n0_; }

	/** Total number of genotype 1 */
	size_t get_n1 () const { return n1_; }

	/** Total number of genotype 2 */
	size_t get_n2 () const { return n2_; }

	/** Total number of individuals */
	size_t get_n () const { return n_; }

	/** computes for a given SNP (as int of its position) the genotype frequency table */
	GenotypeFreq ( const MData & mData, const size_t snp );

	/** returns p-Value for Pearson-Chi-Square Test */
	double calculateChiSquare () const;

	/** a auxiliary for double GenotypeFreq::calculateChiSquare() */
	double chiSquareSummand ( const double obs, const double totalRow, const double totalCol, const double total ) const;

	/** returns p-Value for Cochran-Armitage trend Test */
	double calculateCATT () const;
};

#endif
