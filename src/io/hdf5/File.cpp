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

#include "File.hpp"
#include "../../Exception.hpp"

using namespace std;

namespace hdf5 {

	File::File ( const string& path, const bool create ) : Id(
		create
		? H5Fcreate( path.c_str(), H5F_ACC_EXCL, H5P_DEFAULT, H5P_DEFAULT )
		: H5Fopen( path.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT ),
		path
	) {
		if ( 0 > getId() ) {
			throw Exception(
				"HDF5 \"%s\" %s file failed.",
				getName().c_str(),
				create ? "create output" : "open input"
			);
		}
	}

	File::~File () {
		if ( 0 > H5Fclose( getId() ) ) {
			throw Exception(
				"HDF5 \"%s\" file close failed.",
				getName().c_str()
			);
		}
	}

}
