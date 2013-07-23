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
#include <vector>
#include <hdf5.h>

using namespace std;
using namespace linalg;

namespace io {

	const char
		* const Hdf5Input::snpListPath = "/single_nucleotide_polymorphisms",
		* const Hdf5Input::individualListPath = "/individuals",
		* const Hdf5Input::genotypeMatrixPath = "/genome_matrix",
		* const Hdf5Input::phenotypeVectorPath = "/phenotypes";

	Hdf5Input::Hdf5Input ( const char * const filename )
		:
		fileId( filename ),
		snpNames( fileId, snpListPath ),
		individualNames( fileId, individualListPath ),
		genotypesTransposed( fileId, genotypeMatrixPath ),
		phenotypes( fileId, phenotypeVectorPath )
	{
		// Tests for input data consistency
		if ( snpNames.countDimensions() != genotypesTransposed.countRows() ) {
			throw Exception(
				"HDF5 input file \"%s\""
				" dataset \"%s\" has %u dimensions"
				" and dataset \"%s\" has major %u dimensions"
				" but both numbers should be equal.",
				filename,
				snpListPath,
				snpNames.countDimensions(),
				genotypeMatrixPath,
				genotypesTransposed.countRows()
			);
		}
		if ( individualNames.countDimensions() != genotypesTransposed.countColumns() ) {
			throw Exception(
				"HDF5 input file \"%s\""
				" dataset \"%s\" has %u dimensions"
				" and dataset \"%s\" has minor %u dimensions"
				" but both numbers should be equal.",
				filename,
				individualListPath,
				individualNames.countDimensions(),
				genotypeMatrixPath,
				genotypesTransposed.countColumns()
			);
		}
		if ( individualNames.countDimensions() != phenotypes.countDimensions() ) {
			throw Exception(
				"HDF5 input file \"%s\""
				" dataset \"%s\" has %u dimensions"
				" and dataset \"%s\" has %u dimensions"
				" but both numbers should be equal.",
				filename,
				individualListPath,
				individualNames.countDimensions(),
				phenotypeVectorPath,
				phenotypes.countDimensions()
			);
		}
	}

	size_t Hdf5Input::countSnps () {
		return snpNames.countDimensions();
	}

	size_t Hdf5Input::countIndividuals () {
		return individualNames.countDimensions();
	}

	SNP Hdf5Input::getSnp ( const size_t snpIndex ) {
		const string snpId = snpNames.readOne( snpIndex );
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
					fileId.getName(),
					snpListPath,
					snpIndex,
					snpId.c_str(),
					i
				);
			}
		}
		const string chromosomeId( snpId, 0, chromosomeIdLength );
		const unsigned long position = strtoul( snpId.c_str() + positionStringStart, NULL, 10 );
		return SNP( chromosomeId, snpId, 0.0, position, 0, 0 );
	}

	Individual Hdf5Input::getIndividual ( const size_t individualIndex ) {
		const string individualId = individualNames.readOne( individualIndex );

		const Individual individual(
			"",
			individualId,
			"",
			"",
			Individual::MISSING,
			phenotypes.readOne( individualIndex )
		);

		return individual;
	}

	void Hdf5Input::retrieveGenotypesIntoVector ( const size_t snpIndex, Vector& v ) {
		const size_t
			rows = genotypesTransposed.countColumns(),	// equals rows
			cols = genotypesTransposed.countRows();		// equals cols
		assert( snpIndex < cols );
		assert( rows == v.countDimensions() );
		vector<double> array( rows );
		genotypesTransposed.readRow( snpIndex, array.data() );
		v.fill( array.data() );
	}

	void Hdf5Input::retrievePhenotypesIntoVector ( Vector& v ) {
		const size_t dims = phenotypes.countDimensions();	// equals rows
		assert( dims == v.countDimensions() );
		vector<double> array( dims );
		phenotypes.readAll( array.data() );
		v.fill( array.data() );
	}

	Hdf5Input::~Hdf5Input () {}

}
