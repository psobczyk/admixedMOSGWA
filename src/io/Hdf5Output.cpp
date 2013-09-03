/********************************************************************************
 *	This file is part of the MOSGWA program code.				*
 *	Copyright ©2013, Bernhard Bodenstorfer.					*
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

#include "Hdf5Output.hpp"
#include "Hdf5Constants.hpp"
#include "../Exception.hpp"
#include <vector>
#include <cassert>
#include <cctype>

using namespace std;
using namespace linalg;
using namespace hdf5;
using namespace io::Hdf5Constants;

namespace io {

	Hdf5Output::Hdf5Output (
		const char * const filename,
		const size_t individualCount,
		const size_t snpCount,
		const size_t covariateCount
	)
		:
		individualCount( individualCount ),
		snpCount( snpCount ),
		covariateCount( covariateCount ),
		file( filename, true ),
		individualList( file, individualListPath, individualCount ),
		snpList( file, snpListPath, snpCount ),
		covariateList( file, covariateListPath, covariateCount ),
		phenotypes( file, phenotypeVectorPath, individualCount ),
		genotypesTransposed( file, genotypeMatrixPath, snpCount, individualCount ),
		covariatesTransposed( file, covariateMatrixPath, covariateCount, individualCount )
	{}

	void Hdf5Output::setIndividuals ( const Individual * individuals ) {
		vector<string> individualNames( individualCount );
		for (
			size_t individualIndex = 0;
			individualIndex < individualCount;
			++individualIndex
		) {
			// TODO: save also the family ID using HDF5 compound types
			string individualName = individuals[ individualIndex ].getIndividualID().c_str();
			individualNames.at( individualIndex ) = individualName;
		}
		individualList.writeAll( individualNames.data() );
	}

	void Hdf5Output::storePhenotypeVector ( const Vector& v ) {
		assert( v.countDimensions() == individualCount );
		vector<double> array( individualCount );
		for (
			size_t individualIndex = 0;
			individualIndex < individualCount;
			++individualIndex
		) {
			// eliminate stride (could perhaps be done directly in HDF5)
			array.at( individualIndex ) = v.get( individualIndex );
		}
		phenotypes.writeAll( array.data() );
	}

	void Hdf5Output::setSnps ( const SNP * snps ) {
		vector<string> snpIds( snpCount );
		for (
			size_t snpIndex = 0;
			snpIndex < snpCount;
			++snpIndex
		) {
			string snpId = snps[ snpIndex ].getSnpId();
			enum { NO_SEPARATOR, SEPARATOR, INVALID } status = NO_SEPARATOR;
			for (
				const char* cursor = snpId.c_str();
				0 != *cursor;
				++cursor
			) {
				const char c = *cursor;
				if ( NO_SEPARATOR == status && '_' == c ) {
					status = SEPARATOR;
				} else if ( isdigit( c ) ) {
					// a digit is o.k.
				} else {
					throw Exception(
						"HDF5 output file \"%s\" dataset \"%s\""
						" SNP[%d] has identifier \"%s\"."
						" However,"
						" the current MOSGWA implementation requires"
						" that the identifier are"
						" two decimal numbers"
						" indicating chromosome and position,"
						" separated by a single underscore.",
						file.getName().c_str(),
						snpListPath,
						snpIndex,
						snpId.c_str()
					);
				}
			}
			snpIds.at( snpIndex ) = snpId;
		}
		snpList.writeAll( snpIds.data() );
	}

	void Hdf5Output::storeGenotypeVector ( const size_t snpIndex, const Vector& v ) {
		assert( snpIndex < snpCount );
		assert( v.countDimensions() == individualCount );
		vector<double> array( individualCount );
		for (
			size_t individualIndex = 0;
			individualIndex < individualCount;
			++individualIndex
		) {
			// eliminate stride (could perhaps be done directly in HDF5)
			array.at( individualIndex ) = v.get( individualIndex );
		}
		genotypesTransposed.writeRow( snpIndex, array.data() );
	}

	void Hdf5Output::setCovariates ( const string * covariates ) {
		vector<string> covariateNames( covariateCount );
		for (
			size_t covariateIndex = 0;
			covariateIndex < covariateCount;
			++covariateIndex
		) {
			string covariateName = covariates[ covariateIndex ];
			covariateNames.at( covariateIndex ) = covariateName;
		}
		covariateList.writeAll( covariateNames.data() );
	}

	void Hdf5Output::storeCovariateVector ( const size_t covIndex, const Vector& v ) {
		assert( covIndex < covariateCount );
		assert( v.countDimensions() == individualCount );
		vector<double> array( individualCount );
		for (
			size_t individualIndex = 0;
			individualIndex < individualCount;
			++individualIndex
		) {
			// eliminate stride (could perhaps be done directly in HDF5)
			array.at( individualIndex ) = v.get( individualIndex );
		}
		covariatesTransposed.writeRow( covIndex, array.data() );
	}

	Hdf5Output::~Hdf5Output () {}

}
