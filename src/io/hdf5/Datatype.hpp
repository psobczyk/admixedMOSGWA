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

#ifndef _IO_HDF5_DATATYPE_HPP_
#define _IO_HDF5_DATATYPE_HPP_

#include "Dataset.hpp"

namespace hdf5 {

	/** RAII wrapper for datatype ids. */
	class Datatype : public Id {

		public:

		static const Datatype varString;

		/** Create a double or var-length string data type. */
		Datatype ( const std::string& name, const bool isString );

		/** Wrap an already existing data type <code>hid_t</code> for certain destruction. */
		Datatype ( const hid_t id, const std::string& name );

		/** Wrap the data type of the given data set. */
		Datatype ( Dataset& dataset );

		/** Retrieve the data type size. */
		size_t size ();

		/** Tell whether the data type is variable sized string. */
		bool isVariableString ();

		virtual ~Datatype ();
	};

}

#endif	/* _IO_HDF5_DATATYPE_HPP_ */
