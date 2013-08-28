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

#ifndef _IO_HDF5_OBJECT_TPL_
#define _IO_HDF5_OBJECT_TPL_

#include "Object.hpp"
#include "../../Exception.hpp"
#include <vector>
#include <cassert>

namespace hdf5 {

	template<size_t D> Object<D>::Object (
		File& file,
		const std::string& objectPath,
		const bool isString,
		const size_t size[D]
	)
		:
		dataset(
			H5Dcreate2(
				file.getId(),
				objectPath.c_str(),
				isString ? Datatype::varString.getId() : H5T_NATIVE_DOUBLE,
				Dataspace<D>( objectPath, size ).getId(),
				H5P_DEFAULT,
				H5P_DEFAULT,
				H5P_DEFAULT
			),
			objectPath
		),
		datatype( dataset ),
		dataspace( dataset )
	{
		if ( 0 > dataset.getId() ) {
			throw Exception(
				"HDF5 \"%s\" %u-dimensional %s dataset creation failed.",
				dataset.getName().c_str(),
				D,
				isString ? "var-length string" : "double"
			);
		}
	}

	template<size_t D> Object<D>::Object ( File& file, const std::string& objectPath )
		:
		dataset( file, objectPath ),
		datatype( dataset ),
		dataspace( dataset )
	{
	}

	template<size_t D> size_t Object<D>::countDimensions ( const size_t dim ) const {
		return dataspace.countDimensions( dim );
	}

	template<size_t D> size_t Object<D>::countItems () const {
		return dataspace.countItems();
	}

	template<size_t D> void Object<D>::readAll ( double* array ) {
		if ( 0 > H5Dread( dataset.getId(), H5T_NATIVE_DOUBLE, dataspace.getId(), dataspace.getId(), H5P_DEFAULT, array ) ) {
			throw Exception(
				"HDF5 \"%s\" dataset \"%s\" read all as double precision numbers failed.",
				dataset.getName().c_str()
			);
		}
	}

	template<size_t D> void Object<D>::writeAll ( const double* array ) {
		if ( 0 > H5Dwrite(
			dataset.getId(),
			H5T_NATIVE_DOUBLE,
			dataspace.getId(),
			dataspace.getId(),
			H5P_DEFAULT,
			array
		) ) {
			throw Exception(
				"HDF5 \"%s\" dataset \"%s\" write all as double precision numbers failed.",
				dataset.getName().c_str()
			);
		}
	}

	template<size_t D> void Object<D>::readAll ( std::string* array ) {

		// prepare type
		const size_t datatypeSize = datatype.size() + 1;	// + 1 for trailing \000
		const bool isVarString = datatype.isVariableString();
		Datatype memType( H5Tcopy( H5T_C_S1 ), "copy of H5T_C_S1" );
		if ( 0 > H5Tset_size( memType.getId(), isVarString ? H5T_VARIABLE : datatypeSize ) ) {
			throw isVarString
				? Exception( "HDF5 set string size %u failed.", datatypeSize )
				: Exception( "HDF5 set string size variable failed." );
		}

		const size_t items = countItems();
		if ( isVarString ) {
			std::vector<char*> buffer( items );
			if ( 0 > H5Dread( dataset.getId(), memType.getId(), dataspace.getId(), dataspace.getId(), H5P_DEFAULT, buffer.data() ) ) {
				throw Exception(
					"HDF5 \"%s\" dataset read all as variable length strings failed.",
					dataset.getName().c_str()
				);
			}
			for ( size_t i = 0; i < items; ++i ) {
				const char* value = buffer.at( i );
				array[i] = value;
			}
			if ( 0 > H5Dvlen_reclaim( memType.getId(), dataspace.getId(), H5P_DEFAULT, buffer.data() ) ) {
				throw Exception(
					"HDF5 \"%s\" dataset reclaim buffer failed for %u variable length strings.",
					dataset.getName().c_str(),
					items
				);
			}
		} else {
			const size_t totalSize = datatypeSize * items;
			std::vector<char> buffer( totalSize );
			if ( 0 > H5Dread( dataset.getId(), memType.getId(), H5S_ALL, H5S_ALL, H5P_DEFAULT, buffer.data() ) ) {
				throw Exception(
					"HDF5 \"%s\" dataset read all as fixed length strings failed.",
					dataset.getName().c_str()
				);
			}
			for ( size_t i = 0; i < items; ++i ) {
				const char* value = & buffer.at( datatypeSize * i );
				array[i] = value;
			}
		}
	}

	template<size_t D> void Object<D>::writeAll ( const std::string* array ) {
		const size_t items = countItems();
		std::vector<const char*> buffer( items );
		for ( size_t i = 0; i < items; ++i ) {
			const char* value = array[i].c_str();
			buffer.at( i ) = value;
		}

		if ( 0 > H5Dwrite(
			dataset.getId(),
			Datatype::varString.getId(),
			dataspace.getId(),
			dataspace.getId(),
			H5P_DEFAULT,
			buffer.data()
		) ) {
			throw Exception(
				"HDF5 \"%s\" dataset write all as variable length strings failed.",
				dataset.getName().c_str()
			);
		}
	}
}

#endif	/* _IO_HDF5_OBJECT_TPL_ */
