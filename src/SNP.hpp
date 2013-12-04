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

#ifndef SNP_HPP
#define SNP_HPP

#include <string>

/** Single Nucleotide Polymorphism. */
class SNP {

	/** from 0 to 26: where
	* 0 = unplaced,
	* 1-22 autosomal,
	* 23 = X chromosome,
	* 24 = Y chromosome,
	* 25 = pseudo-autosomal region of X,
	* 26 = MT Mitochondrial
	*/
	unsigned int chromosomeId;

	/** name of the SNP */
	std::string snpId;

	//in morgans
	double geneticDistance;

	//in bp units
	int basePairPosition;

	// the alleles represented by 0 and 1
	char allele1;
	char allele2;

	// to save the outcome of the accoring single marker test
	double singleMarkerTest;

public:

	/** Construct a SNP. */
	SNP (
		const unsigned int chromosomeId,
		const std::string &snpId,
		const double geneticDistance,
		const int basePairPosition,
		const char allele1,
		const char allele2
	);

	/** Construct a SNP. */
	SNP (
		const std::string &chromosomeId,
		const std::string &snpId,
		const double geneticDistance,
		const int basePairPosition,
		const char allele1,
		const char allele2
	);

	// Setters
	void setSingleMarkerTest ( const double singleMarkerTest );

	// Getters
	int getChromosome () const;
	std::string getSnpId () const;
	double getGeneticDistance () const;
	int getBasePairPosition () const;
	char getAllele1 () const;
	char getAllele2 () const;
	double getSingleMarkerTest () const;
};

/** Output a {@link SNP}. */
std::ostream& operator<< ( std::ostream& s, const SNP& snp );

#endif
