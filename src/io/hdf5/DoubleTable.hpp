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

#ifndef _IO_HDF5_DOUBLETABLE_HPP_
#define _IO_HDF5_DOUBLETABLE_HPP_

#include "Object.hpp"
#include <string>
#include <hdf5.h>

namespace hdf5 {

	/** Two-dimensional table of numbers in HDF5. */
	class DoubleTable : public Object<2> {

		void transferRow (
			const size_t row,
			double* array,
			herr_t function ( hid_t, hid_t, hid_t, hid_t, hid_t, void* ),
			const char * description
		);

		public:

		/** Create table in HDF5 file. */
		DoubleTable ( File& file, const std::string& objectPath, const size_t rows, const size_t columns );

		/** Open existing table in HDF5 file. */
		DoubleTable ( File& fileId, const std::string& objectPath );

		/** Retrieve the number of rows. */
		size_t countRows () const;

		/** Retrieve the number of columns. */
		size_t countColumns () const;

		/** Retrieve one row (= majour index).
		* @param row specifies which row to read.
		*/
		void readRow ( const size_t row, double* array );

		/** Store one row (= majour index).
		* @param row specifies which row to write.
		*/
		void writeRow ( const size_t row, const double* array );
	};

}

#endif	/* _IO_HDF5_DOUBLETABLE_HPP_ */
