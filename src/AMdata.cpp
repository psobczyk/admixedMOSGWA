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

#include "AMdata.hpp"
#define	IDV_NO 3
#define	SNP_NO 2
#define	CNP_NO 1
#define	COV_NO 1

using namespace std;

const Individual individuals[] = {
	Individual( "Bierwetter", "Bernhard", "Camillo", "Praxedis", Individual::MALE ),
	Individual( "Wetter", "Praxedis", "Gianluca", "Caecilia", Individual::FEMALE ),
	Individual( "Bier", "Camillo", "Gonzago", "Pandora", Individual::MALE )
};

const SNP snps[] = {
	SNP( 26, string( "SNP26" ), 26.26, 2626, 'A', 'C' ),
	SNP( 17, string( "SNP17" ), 17.17, 1717, 'C', 'T' )
};

const SNP cnps[] = {
	SNP( 9, string( "CNP9" ), 9.9, 99, 'G', 'A' )
};

const string covariates[2] = { "Age", "Sex" };

AMdata::AMdata () {}

AMdata::~AMdata () {}

size_t AMdata::getIdvNo () const {
	return IDV_NO;
}

Individual AMdata::getIndividual ( const size_t index ) const {
	return individuals[index];
}

size_t AMdata::getSnpNo () const {
	return SNP_NO;
}

const SNP & AMdata::getSNP ( const size_t snp ) const {
	return snps[snp];
}

void AMdata::getXcolumn ( const size_t snp, linalg::Vector& vector ) const {
	const double X[IDV_NO][SNP_NO] = {
		1, 0,
		-1, 1,
		0, 1
	};
	for ( size_t i = 0; i < IDV_NO; ++i ) {
		vector.set( i, X[i][snp] );
	}
}

size_t AMdata::getCnpNo () const {
	return CNP_NO;
}

const SNP & AMdata::getCNP ( const size_t cnp ) const {
	return cnps[cnp];
}

void AMdata::getCcolumn ( const size_t cnp, linalg::Vector& vector ) const {
	const double C[IDV_NO][CNP_NO] = {
		-1,
		1,
		1,
	};
	for ( size_t i = 0; i < IDV_NO; ++i ) {
		vector.set( i, C[i][cnp] );
	}
}

size_t AMdata::getCovNo () const {
	return COV_NO;
}

const std::string& AMdata::getCovMatElementName( const size_t cov ) const {
	return covariates[cov];
}

void AMdata::getCovariateColumn ( const size_t cov, linalg::Vector& vector ) const {
	const double Cov[IDV_NO][COV_NO] = {
		-1,
		1,
		-1,
	};
	for ( size_t i = 0; i < IDV_NO; ++i ) {
		vector.set( i, Cov[i][cov] );
	}
}

void AMdata::getY ( linalg::Vector& vector ) const {
	const double Y[IDV_NO] = {
		0, 0, 1
	};
	vector.fill( Y );
}
