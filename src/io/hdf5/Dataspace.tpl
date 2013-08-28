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

#ifndef _IO_HDF5_DATASPACE_TPL_
#define _IO_HDF5_DATASPACE_TPL_

#include "Dataspace.hpp"
#include "../../Exception.hpp"
#include <cassert>

namespace hdf5 {

	template<size_t D> Dataspace<D>::Dataspace ( const std::string& name, const size_t size[D] )
		:
		Id(
			H5Screate( H5S_SIMPLE ),
			name
		)
	{
		if ( 0 > getId() ) {
			throw Exception(
				"HDF5 \"%s\" simple dataspace creation failed.",
				getName().c_str()
			);
		}
		hsize_t hSize[D];
		for ( size_t i = 0; i < D; ++i ) {
			hSize[i] = size[i];
			assert( size[i] == hSize[i] );	// guard against type overflow
		}
		if ( 0 > H5Sset_extent_simple( getId(), D, hSize, NULL ) ) {
			throw Exception(
				"HDF5 \"%s\" dataspace set %u sizes failed.",
				getName().c_str(),
				D
			);
		}
		// initialise the remainder back from HDF5 object
		initFromId();
		for ( size_t i = 0; i < D; ++i ) {
			assert( size[i] == this->size[i] );
		}
	}

	template<size_t D> Dataspace<D>::Dataspace ( const hid_t id, const std::string& name )
		:
		Id( id, name )
	{
		initFromId();
	}

	template<size_t D> Dataspace<D>::Dataspace ( Dataset& dataset )
		:
		Id(
			H5Dget_space( dataset.getId() ),
			dataset.getName()
		)
	{
		initFromId();
	}

	template<size_t D> void Dataspace<D>::initFromId () {
		if ( 0 > getId() ) {
			throw Exception(
				"HDF5 \"%s\" get dataspace failed.",
				getName().c_str()
			);
		}
		const int dims = H5Sget_simple_extent_ndims( getId() );
		if ( 0 > dims ) {
			throw Exception(
				"HDF5 \"%s\" dataspace failed to determine dimension.",
				getName().c_str()
			);
		} else if ( D != dims ) {
			throw Exception(
				"HDF5 \"%s\" dataspace has %u dimensions, but should have %u.",
				getName().c_str(),
				dims,
				D
			);
		}
		hsize_t h5size[D];
		if ( 0 > H5Sget_simple_extent_dims( getId(), h5size, NULL ) ) {
			throw Exception(
				"HDF5 \"%s\" dataspace failed to determine dimension sizes.",
				getName().c_str()
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

	template<size_t D> size_t Dataspace<D>::countDimensions ( const size_t dim ) const {
		assert( dim < D );
		return size[dim];
	}

	template<size_t D> size_t Dataspace<D>::countItems () const {
		return volume;
	}

	template<size_t D> Dataspace<D>::~Dataspace () {
		if ( 0 > H5Sclose( getId() ) ) {
			throw Exception(
				"HDF5 \"%s\" dataspace close failed.",
				getName().c_str()
			);
		}
	}
}

#endif	/* _IO_HDF5_DATASPACE_TPL_ */
