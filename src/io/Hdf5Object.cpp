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

#include "Hdf5Object.hpp"
#include "../Exception.hpp"
#include <vector>
#include <cassert>

using namespace std;

namespace io {

	Hdf5Id::Hdf5Id ( const hid_t id, const string& name ) : id( id ), name( name ) {}

	const char * Hdf5Id::getName () const {
		return name.c_str();
	}

	/** If the abstract destructor has no implementation, g++ yields:
	* "undefined reference to `io::Hdf5Id::~Hdf5Id()'".
	*/
	Hdf5Id::~Hdf5Id () {}

	Hdf5FileId::Hdf5FileId ( const string& path ) : Hdf5Id(
		H5Fopen( path.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT ),
		path
	) {
		if ( 0 > id ) {
			throw Exception(
				"HDF5 \"%s\" open input file failed.",
				getName()
			);
		}
	}

	Hdf5FileId::~Hdf5FileId () {
		if ( 0 > H5Fclose( id ) ) {
			throw Exception(
				"HDF5 \"%s\" file close failed.",
				getName()
			);
		}
	}

	Hdf5DatasetId::Hdf5DatasetId ( Hdf5FileId& fileId, const string& objectPath ) : Hdf5Id(
		H5Dopen2( fileId.id, objectPath.c_str(), H5P_DEFAULT ),
		fileId.name + objectPath
	) {
		if ( 0 > id ) {
			throw Exception(
				"HDF5 \"%s\" file does not contain dataset \"%s\" (use command \"h5ls %s\" to list the data sets).",
				fileId.getName(),
				objectPath.c_str(),
				fileId.getName()
			);
		}
		const hid_t creationProperties = H5Dget_create_plist( id );
		const htri_t filtersAvailable = H5Pall_filters_avail( creationProperties );
		H5Pclose( creationProperties );
		if ( 0 > filtersAvailable ) {
			throw Exception(
				"HDF5 \"%s\" dataset use of filters cannot be determined.",
				getName()
			);
		} else if ( ! filtersAvailable ) {
			throw Exception(
				"HDF5 \"%s\" dataset uses filters which are not available on your platform (use command \"h5ls -v %s\" to find out which filters are used).",
				getName(),
				fileId.getName()
			);
		}
	}

	Hdf5DatasetId::~Hdf5DatasetId () {
		if ( 0 > H5Dclose( id ) ) {
			throw Exception(
				"HDF5 \"%s\" dataset close failed.",
				getName()
			);
		}
	}

	Hdf5DatatypeId::Hdf5DatatypeId ( Hdf5DatasetId& datasetId ) : Hdf5Id(
		H5Dget_type( datasetId.id ),
		datasetId.name
	) {
		if ( 0 > id ) {
			throw Exception(
				"HDF5 \"%s\" dataset failed to determine datatype.",
				getName()
			);
		}
	}

	Hdf5DatatypeId::Hdf5DatatypeId ( const hid_t id, const string& name ) : Hdf5Id(
		id,
		name
	) {
		if ( 0 > id ) {
			throw Exception(
				"HDF5 \"%s\" datatype creation failed.",
				getName()
			);
		}
	}

	size_t Hdf5DatatypeId::size () {
		const size_t datatypeSize = H5Tget_size( id );
		if ( 0 >= datatypeSize ) {
			throw Exception(
				"HDF5 \"%s\" datatype failed to get its size.",
				getName()
			);
		}
		return datatypeSize;
	}

	bool Hdf5DatatypeId::isVariableString () {
		const htri_t isVarString = H5Tis_variable_str( id );
		if ( 0 > isVarString ) {
			throw Exception(
				"HDF5 \"%s\" datatype failed to determine whether it is a variable length string.",
				getName()
			);
		}
		return isVarString;
	}

	Hdf5DatatypeId::~Hdf5DatatypeId () {
		if ( 0 > H5Tclose( id ) ) {
			throw Exception(
				"HDF5 \"%s\" datatype close failed.",
				getName()
			);
		}
	}

	Hdf5DataspaceId::Hdf5DataspaceId ( Hdf5DatasetId& datasetId ) : Hdf5Id(
		H5Dget_space( datasetId.id ),
		datasetId.name
	) {
		if ( 0 > id ) {
			throw Exception(
				"HDF5 \"%s\" dataset failed to determine dataspace.",
				getName()
			);
		}
	}

	Hdf5DataspaceId::Hdf5DataspaceId ( const hid_t id, const string& name ) : Hdf5Id(
		id,
		name
	) {
		if ( 0 > id ) {
			throw Exception(
				"HDF5 \"%s\" dataspace creation failed.",
				getName()
			);
		}
	}

	Hdf5DataspaceId::~Hdf5DataspaceId () {
		if ( 0 > H5Sclose( id ) ) {
			throw Exception(
				"HDF5 \"%s\" dataspace close failed.",
				getName()
			);
		}
	}

	Hdf5StringList::Hdf5StringList ( Hdf5FileId& fileId, const std::string& objectPath )
		: Hdf5Object<1>( fileId, objectPath )
	{}

	size_t Hdf5StringList::countDimensions () const {
		return size[0];
	}

	string Hdf5StringList::readOne ( const size_t index ) {
		assert( index < countDimensions() );

		// prepare type
		const size_t datatypeSize = datatypeId.size() + 1;	// + 1 for trailing \000
		const bool isVarString = datatypeId.isVariableString();
		Hdf5DatatypeId memType( H5Tcopy( H5T_C_S1 ), "copy of H5T_C_S1" );
		if ( 0 > H5Tset_size( memType.id, isVarString ? H5T_VARIABLE : datatypeSize ) ) {
			throw isVarString
				? Exception( "HDF5 set string size %u failed.", datatypeSize )
				: Exception( "HDF5 set string size variable failed." );
		}

		// prepare space
		const hsize_t coordinates[1] = { index };
		assert( coordinates[0] == index );	// guard against assignment overflow
		// REMARK<BB>: select breaks thread-safety.
		if ( 0 > H5Sselect_elements( dataspaceId.id, H5S_SELECT_SET, 1, coordinates ) ) {
			throw Exception(
				"HDF5 \"%s\" dataset select element %u failed.",
				datasetId.getName(),
				index
			);
		}
		Hdf5DataspaceId memSpace( H5Screate( H5S_SCALAR ), "creature of H5S_SCALAR" );

		// read
		string value;
		if ( isVarString ) {
			char* buffer[1];
			if ( 0 > H5Dread( datasetId.id, memType.id, memSpace.id, dataspaceId.id, H5P_DEFAULT, buffer ) ) {
				throw Exception(
					"HDF5 \"%s\" dataset read variable length string[%u] failed.",
					datasetId.getName(),
					index
				);
			}
			value = string( buffer[0] );
			if ( 0 > H5Dvlen_reclaim( memType.id, memSpace.id, H5P_DEFAULT, buffer ) ) {
				throw Exception(
					"HDF5 \"%s\" dataset reclaim buffer for variable length string[%u] failed.",
					datasetId.getName(),
					index
				);
			}
		} else {
			vector<char> buffer( datatypeSize );
			if ( 0 > H5Dread( datasetId.id, memType.id, memSpace.id, dataspaceId.id, H5P_DEFAULT, buffer.data() ) ) {
				throw Exception(
					"HDF5 \"%s\" dataset read fixed length string[%u] failed.",
					datasetId.getName(),
					index
				);
			}
			value = string( buffer.data(), 0, datatypeSize );
		}

		return value;
	}

	Hdf5DoubleList::Hdf5DoubleList ( Hdf5FileId& fileId, const std::string& objectPath )
		: Hdf5Object<1>( fileId, objectPath )
	{}

	size_t Hdf5DoubleList::countDimensions () const {
		return size[0];
	}

	double Hdf5DoubleList::readOne ( const size_t index ) {
		assert( index < countDimensions() );

		// prepare space
		const hsize_t coordinates[1] = { index };
		assert( coordinates[0] == index );	// guard against assignment overflow
		// REMARK<BB>: select breaks thread-safety.
		if ( 0 > H5Sselect_elements( dataspaceId.id, H5S_SELECT_SET, 1, coordinates ) ) {
			throw Exception(
				"HDF5 \"%s\" dataset select element %u failed.",
				datasetId.getName(),
				index
			);
		}
		Hdf5DataspaceId memSpace( H5Screate( H5S_SCALAR ), "creature of H5S_SCALAR" );

		// read
		double value;
		if ( 0 > H5Dread( datasetId.id, H5T_NATIVE_DOUBLE, memSpace.id, dataspaceId.id, H5P_DEFAULT, &value ) ) {
			throw Exception(
				"HDF5 \"%s\" dataset read double[%u] failed.",
				datasetId.getName(),
				index
			);
		}
		return value;
	}

	Hdf5DoubleTable::Hdf5DoubleTable ( Hdf5FileId& fileId, const std::string& objectPath )
		: Hdf5Object<2>( fileId, objectPath )
	{}

	size_t Hdf5DoubleTable::countRows () const {
		return size[0];		// major in C: array[major][minor]
	}

	size_t Hdf5DoubleTable::countColumns () const {
		return size[1];		// minor in C: array[major][minor]
	}

	void Hdf5DoubleTable::readRow ( const size_t row, double* array ) {
		assert( row < countRows() );

		const hsize_t
			start[2] = { row, 0 },
			count[2] = { 1, countColumns() };
		// REMARK<BB>: select breaks thread-safety.
		const herr_t status = H5Sselect_hyperslab( dataspaceId.id, H5S_SELECT_SET, start, NULL, count, NULL );

		Hdf5DataspaceId memSpace( H5Screate_simple( 1, &count[1], NULL ), "creature of simple row space" );
		const herr_t status2 = H5Dread( datasetId.id, H5T_NATIVE_DOUBLE, memSpace.id, dataspaceId.id, H5P_DEFAULT, array );
		if ( 0 > status2 ) {
			throw Exception(
				"HDF5 \"%s\" dataset read row[%u] of length %u failed.",
				datasetId.getName(),
				row,
				countColumns()
			);
		}
	}
}
