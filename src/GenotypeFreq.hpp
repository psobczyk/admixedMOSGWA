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

#ifndef GENOTYPEFREQ_HPP
#define GENOTYPEFREQ_HPP

#include "MData.hpp"

/** The genotype frequency table for case-control phenotypes. */
class GenotypeFreq {

private:
	int r0_; // Case with genotype "0"
	int r1_; // Case with genotype "1"
	int r2_; // Case with genotype "2"
	int r_;	 // # of Cases
	int s0_; // Control with genotype "0"
	int s1_; // Control with genotype "1"
	int s2_; // Control with genotype "2"
	int s_;	 // # of Controls
	int n0_; // total # of genotype "0"
	int n1_; // total # of genotype "1"
	int n2_; // total # of genotype "2"
	int n_;	 // # of individuals

public:
	/** Number of cases with genotype 0 */
	int get_r0 () const { return r0_; }

	/** Number of cases with genotype 1 */
	int get_r1 () const { return r1_; }

	/** Number of cases with genotype 2 */
	int get_r2 () const { return r2_; }

	/** Number of cases total */
	int get_r () const { return r_; }

	/** Number of controls with genotype 0 */
	int get_s0 () const { return s0_; }

	/** Number of controls with genotype 1 */
	int get_s1 () const { return s1_; }

	/** Number of controls with genotype 2 */
	int get_s2 () const { return s2_; }

	/** Number of controls total */
	int get_s () const { return s_; }

	/** Total number of genotype 0 */
	int get_n0 () const { return n0_; }

	/** Total number of genotype 1 */
	int get_n1 () const { return n1_; }

	/** Total number of genotype 2 */
	int get_n2 () const { return n2_; }

	/** Total number of individuals */
	int get_n () const { return n_; }

	/** computes for a given SNP (as int of its position) the genotype frequency table */
	GenotypeFreq ( const MData & mData, const int snp );

	/** returns p-Value for Pearson-Chi-Square Test */
	double calculateChiSquare () const;

	/** a auxiliary for double GenotypeFreq::calculateChiSquare() */
	double chiSquareSummand ( const double obs, const double totalRow, const double totalCol, const double total ) const;

	/** returns p-Value for Cochran-Armitage trend Test */
	double calculateCATT () const;
};
#endif
