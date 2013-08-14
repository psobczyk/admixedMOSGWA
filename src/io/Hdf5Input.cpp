/********************************************************************************
 *	This file is part of the MOSGWA program code.				*
 *	Copyright ©2012–2013, Bernhard Bodenstorfer.				*
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

#include "Hdf5Input.hpp"
#include "../Exception.hpp"
#include <cassert>

using namespace std;
using namespace hdf5;
using namespace linalg;

namespace io {

	const char
		* const Hdf5Input::snpListPath = "/single_nucleotide_polymorphisms",
		* const Hdf5Input::individualListPath = "/individuals",
		* const Hdf5Input::genotypeMatrixPath = "/genome_matrix",
		* const Hdf5Input::covariateListPath = "/covariates",
		* const Hdf5Input::covariateMatrixPath = "/covariate_matrix",
		* const Hdf5Input::phenotypeVectorPath = "/phenotypes";

	Hdf5Input::Hdf5Input ( const char * const filename, const bool useCovariates )
		:
		file( filename ),
		phenotypes( file, phenotypeVectorPath ),
		genotypesTransposed( file, genotypeMatrixPath )
	{
		// Read about individuals
		{
			StringList individualList( file, individualListPath );
			const size_t individualCount = individualList.countDimensions();
			vector<string> individualNames( individualCount );
			individualList.readAll( individualNames.data() );
			individuals.reserve( individualCount );		// hint space requirements
			for ( size_t individualIndex = 0; individualIndex < individualCount; ++individualIndex ) {
				const string individualId = individualNames.at( individualIndex );
				const Individual individual(
					"",
					individualId,
					"",
					"",
					Individual::MISSING
				);
				individuals.push_back( individual );
			}
		}

		// Read about SNPs
		{
			StringList snpList( file, snpListPath );
			const size_t snpCount = snpList.countDimensions();
			vector<string> snpNames( snpCount );
			snpList.readAll( snpNames.data() );
			snps.reserve( snpCount );	// hint space requirements
			for ( size_t snpIndex = 0; snpIndex < snpCount; ++snpIndex ) {
				const string snpId = snpNames.at( snpIndex );
				const size_t idLength = snpId.length();
				size_t
					i,
					chromosomeIdLength = 0,
					positionStringStart = 0;
				for ( i = 0; i < idLength; ++i ) {
					const char c = snpId[i];
					if ( '_' == c && 0 == positionStringStart ) {
						chromosomeIdLength = i;
						positionStringStart = i + 1;
					} else if ( c < '0' || '9' < c ) {
						throw Exception(
							"HDF5 input file \"%s\" dataset \"%s\""
							" SNP[%d] has bad character in name \"%s\" at position %l;"
							" expecting decimals indicating chromosome and position,"
							" separated by a single underscore.",
							file.getName().c_str(),
							snpListPath,
							snpIndex,
							snpId.c_str(),
							i
						);
					}
				}
				const string chromosomeId( snpId, 0, chromosomeIdLength );
				const unsigned long position = strtoul( snpId.c_str() + positionStringStart, NULL, 10 );
				const SNP snp( chromosomeId, snpId, 0.0, position, 0, 0 );
				snps.push_back( snp );
			}
		}

		if ( useCovariates ) {
			StringList covariateList( file, covariateListPath );
			const size_t covariateCount = covariateList.countDimensions();
			covariates.resize( covariateCount );
			covariateList.readAll( covariates.data() );
			covariatesTransposedPtr.reset( new DoubleTable( file, covariateMatrixPath ) );
		} else {
			// else no info is read thus leaving covariate count 0
		}

		// Tests for input data consistency
		if ( countSnps() != genotypesTransposed.countRows() ) {
			throw Exception(
				"HDF5 input file \"%s\""
				" dataset \"%s\" has %u dimensions"
				" and dataset \"%s\" has major %u dimensions"
				" but both numbers should be equal.",
				filename,
				snpListPath,
				countSnps(),
				genotypeMatrixPath,
				genotypesTransposed.countRows()
			);
		}
		if ( countIndividuals() != genotypesTransposed.countColumns() ) {
			throw Exception(
				"HDF5 input file \"%s\""
				" dataset \"%s\" has %u dimensions"
				" and dataset \"%s\" has minor %u dimensions"
				" but both numbers should be equal.",
				filename,
				individualListPath,
				countIndividuals(),
				genotypeMatrixPath,
				genotypesTransposed.countColumns()
			);
		}
		if ( countIndividuals() != phenotypes.countDimensions() ) {
			throw Exception(
				"HDF5 input file \"%s\""
				" dataset \"%s\" has %u dimensions"
				" and dataset \"%s\" has %u dimensions"
				" but both numbers should be equal.",
				filename,
				individualListPath,
				countIndividuals(),
				phenotypeVectorPath,
				phenotypes.countDimensions()
			);
		}

		if ( useCovariates ) {
			if ( countIndividuals() != covariatesTransposedPtr->countColumns() ) {
				throw Exception(
					"HDF5 input file \"%s\""
					" dataset has %u rows"
					" and transposed covariates dataset \"%s\" has minor %u dimensions"
					" but both numbers should be equal.",
					filename,
					countIndividuals(),
					covariateMatrixPath,
					covariatesTransposedPtr->countColumns()
				);
			}
			if ( countCovariates() != covariatesTransposedPtr->countRows() ) {
				throw Exception(
					"HDF5 input file \"%s\""
					" dataset \"%s\" has %u dimensions"
					" and dataset \"%s\" has %u dimensions"
					" but both numbers should be equal.",
					filename,
					covariateListPath,
					countCovariates(),
					covariateMatrixPath,
					covariatesTransposedPtr->countRows()
				);
			}
		}
	}

	void Hdf5Input::retrieveGenotypeVector ( const size_t snpIndex, Vector& v ) {
		const size_t
			rows = genotypesTransposed.countColumns(),	// equals rows
			cols = genotypesTransposed.countRows();		// equals cols
		assert( snpIndex < cols );
		assert( rows == v.countDimensions() );
		vector<double> array( rows );
		genotypesTransposed.readRow( snpIndex, array.data() );
		v.fill( array.data() );
	}

	void Hdf5Input::retrievePhenotypeVector ( Vector& v ) {
		const size_t dims = phenotypes.countDimensions();	// equals rows
		assert( dims == v.countDimensions() );
		vector<double> array( dims );
		phenotypes.readAll( array.data() );
		v.fill( array.data() );
	}

	void Hdf5Input::retrieveCovariateVector ( const size_t covIndex, linalg::Vector& v ) {
		const size_t
			rows = covariatesTransposedPtr->countColumns(), // equals rows
			covs = covariatesTransposedPtr->countRows();	// equals covs
		assert( covIndex < covs );
		assert( rows == v.countDimensions() );
		vector<double> array( rows );
		covariatesTransposedPtr->readRow( covIndex, array.data() );
		v.fill( array.data() );
	}

	Hdf5Input::~Hdf5Input () {}

}
