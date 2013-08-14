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

#include "Datatype.hpp"
#include "../../Exception.hpp"

using namespace std;

namespace hdf5 {

	const Datatype Datatype::varString ( "hdf5::Datatype::varString", true );

	Datatype::Datatype ( const string& name, const bool isString )
		: Id(
			isString ? H5Tcopy( H5T_C_S1 ) : H5Tcopy( H5T_NATIVE_DOUBLE ),
			name
		)
	{
		if ( 0 > getId() ) {
			throw Exception(
				"HDF5 \"%s\" datatype %s copy failed.",
				getName().c_str(),
				isString ? "string" : "double"
			);
		}
		if ( isString && 0 > H5Tset_size( getId(), H5T_VARIABLE ) ) {
			throw Exception(
				"HDF5 \"%s\" string type set variable-length failed.",
				getName().c_str()
			);
		}
	}

	Datatype::Datatype ( const hid_t id, const string& name ) : Id(
		id,
		name
	) {
		if ( 0 > getId() ) {
			throw Exception(
				"HDF5 \"%s\" datatype creation failed.",
				getName().c_str()
			);
		}
	}

	Datatype::Datatype ( Dataset& dataset ) : Id(
		H5Dget_type( dataset.getId() ),
		dataset.getName()
	) {
		if ( 0 > getId() ) {
			throw Exception(
				"HDF5 \"%s\" dataset failed to determine datatype.",
				getName().c_str()
			);
		}
	}

	size_t Datatype::size () {
		const size_t datatypeSize = H5Tget_size( getId() );
		if ( 0 >= datatypeSize ) {
			throw Exception(
				"HDF5 \"%s\" datatype failed to get its size.",
				getName().c_str()
			);
		}
		return datatypeSize;
	}

	bool Datatype::isVariableString () {
		const htri_t isVarString = H5Tis_variable_str( getId() );
		if ( 0 > isVarString ) {
			throw Exception(
				"HDF5 \"%s\" datatype failed to determine whether it is a variable length string.",
				getName().c_str()
			);
		}
		return isVarString;
	}

	Datatype::~Datatype () {
		if ( 0 > H5Tclose( getId() ) ) {
			throw Exception(
				"HDF5 \"%s\" datatype close failed.",
				getName().c_str()
			);
		}
	}

}
