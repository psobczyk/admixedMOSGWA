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

#include "DoubleTable.hpp"
#include "../../Exception.hpp"
#include <cassert>

using namespace std;

namespace hdf5 {

	DoubleTable::DoubleTable ( File& file, const std::string& objectPath )
		: Object<2>( file, objectPath )
	{}

	size_t DoubleTable::countRows () const {
		return countDimensions( 0 );		// major in C: array[major][minor]
	}

	size_t DoubleTable::countColumns () const {
		return countDimensions( 1 );		// minor in C: array[major][minor]
	}

	void DoubleTable::readRow ( const size_t row, double* array ) {
		assert( row < countRows() );

		const hsize_t
			start[2] = { row, 0 },
			count[2] = { 1, countColumns() };
		// mark selection on a copy only to keep it thread-safe
		Dataspace<2> dataspaceCopy(
			H5Scopy( dataspace.getId() ),
			dataspace.getName() + "[clone]"
		);
		const herr_t status = H5Sselect_hyperslab(
			dataspaceCopy.getId(),
			H5S_SELECT_SET,
			start, NULL,
			count, NULL
		);

		Dataspace<1> memSpace(
			H5Screate_simple( 1, &count[1], NULL ),
			"creature of simple row space"
		);
		const herr_t status2 = H5Dread(
			dataset.getId(),
			H5T_NATIVE_DOUBLE,
			memSpace.getId(),
			dataspaceCopy.getId(),
			H5P_DEFAULT,
			array
		);
		if ( 0 > status2 ) {
			throw Exception(
				"HDF5 \"%s\" dataset read row[%u] of length %u failed.",
				dataset.getName().c_str(),
				row,
				countColumns()
			);
		}
	}
}
