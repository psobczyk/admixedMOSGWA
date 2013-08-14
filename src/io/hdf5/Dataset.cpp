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

#include "Dataset.hpp"
#include "../../Exception.hpp"

using namespace std;

namespace hdf5 {

	Dataset::Dataset ( const hid_t id, const string& name )
	: Id( id, name ) {
		if ( 0 > getId() ) {
			throw Exception(
				"HDF5 \"%s\" dataset creation failed.",
				getName().c_str()
			);
		}
	}

	Dataset::Dataset ( File& file, const string& objectPath ) : Id(
		H5Dopen2( file.getId(), objectPath.c_str(), H5P_DEFAULT ),
		file.getName() + objectPath
	) {
		if ( 0 > getId() ) {
			throw Exception(
				"HDF5 \"%s\" file does not contain dataset \"%s\" (use command \"h5ls %s\" to list the data sets).",
				file.getName().c_str(),
				objectPath.c_str(),
				file.getName().c_str()
			);
		}
		const hid_t creationProperties = H5Dget_create_plist( getId() );
		const htri_t filtersAvailable = H5Pall_filters_avail( creationProperties );
		H5Pclose( creationProperties );
		if ( 0 > filtersAvailable ) {
			throw Exception(
				"HDF5 \"%s\" dataset use of filters cannot be determined.",
				getName().c_str()
			);
		} else if ( ! filtersAvailable ) {
			throw Exception(
				"HDF5 \"%s\" dataset uses filters which are not available on your platform (use command \"h5ls -v %s\" to find out which filters are used).",
				getName().c_str(),
				file.getName().c_str()
			);
		}
	}

	Dataset::~Dataset () {
		if ( 0 > H5Dclose( getId() ) ) {
			throw Exception(
				"HDF5 \"%s\" dataset close failed.",
				getName().c_str()
			);
		}
	}

}
