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

#ifndef _IO_HDF5OBJECT_TPL_
#define _IO_HDF5OBJECT_TPL_

#include "Hdf5Object.hpp"
#include "../Exception.hpp"
#include <vector>
#include <cassert>

using namespace std;

namespace io {

	template<size_t D> Hdf5Object<D>::Hdf5Object ( Hdf5FileId& fileId, const std::string& objectPath )
		:
		datasetId( fileId, objectPath ),
		datatypeId( datasetId ),
		dataspaceId( datasetId )
	{
		const int dims = H5Sget_simple_extent_ndims( dataspaceId.id );
		if ( 0 > dims ) {
			throw Exception(
				"HDF5 \"%s\" dataspace failed to determine dimension.",
				dataspaceId.getName()
			);
		} else if ( D != dims ) {
			throw Exception(
				"HDF5 \"%s\" dataspace has %u dimensions, but should have %u.",
				dataspaceId.getName(),
				dims,
				D
			);
		}
		hsize_t h5size[D];
		if ( 0 > H5Sget_simple_extent_dims( dataspaceId.id, h5size, NULL ) ) {
			throw Exception(
				"HDF5 \"%s\" dataspace failed to determine dimension sizes.",
				dataspaceId.getName()
			);
		}

		volume = 1;
		for ( size_t dim = 0; dim < D; ++dim ) {
			volume *= ( size[dim] = h5size[dim] );
			assert( h5size[dim] == size[dim] );	// guard against overflow from hsize_t to size_t
		}
		// guard against overflow for the product
		bool nullVolume = false;
		size_t volumeTest = volume;
		for ( size_t dim = 0; dim < D; ++dim ) {
			nullVolume |= 0 == size[dim];
			volumeTest /= size[dim];
		}
		assert( nullVolume || 1 == volumeTest );
	}

	template<size_t D> size_t Hdf5Object<D>::countItems () const {
		return volume;
	}

	template<size_t D> void Hdf5Object<D>::readAll ( double* array ) {
		if ( 0 > H5Dread( datasetId.id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, array ) ) {
			throw Exception(
				"HDF5 \"%s\" dataset \"%s\" read failed.",
				datasetId.getName()
			);
		}
	}
}

#endif	/* _IO_HDF5OBJECT_TPL_ */
