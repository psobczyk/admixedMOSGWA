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

#ifndef _IO_HDF5_OBJECT_HPP_
#define _IO_HDF5_OBJECT_HPP_

#include "Datatype.hpp"
#include "Dataspace.hpp"

namespace hdf5 {

	/** Practical access adapter for HDF5 data. */
	template<size_t D> class Object {

		protected:

		Dataset dataset;
		Datatype datatype;
		Dataspace<D> dataspace;

		public:

		/** Wraps the object in a HDF5 file identified by the given path. */
		Object ( File& file, const std::string& objectPath );

		/** Creates in a HDF5 file the object identified by the given path.
		* @param isString specifies whether the object should hold variable-length strings
		* or otherwise double precision numbers.
		* @param size specifies the array dimensions from majour to minour
		*/
		Object (
			File& file,
			const std::string& objectPath,
			const bool isString,
			const size_t size[D]
		);

		/** Get the size of the given dimension. */
		size_t countDimensions ( const size_t dim ) const;

		/** Get the number of data items. */
		size_t countItems () const;

		/** Retrieve numerical data in one big lump.
		* The array size must be at least {@link #countItems}
		* and the type <code>double<code> must match the stored HDF5-type.
		*/
		void readAll ( double* array );

		/** Retrieve textual data in one big lump.
		* The array size must be at least {@link #countItems}
		* and the type <code>std::string<code> must match the stored HDF5-type.
		*/
		void readAll ( std::string* array );

		/** Store numerical data in one big lump.
		* The array size must be at least {@link #countItems}
		* and the type <code>double<code> must match the stored HDF5-type.
		*/
		void writeAll ( const double* array );

		/** Store textual data in one big lump.
		* The array size must be at least {@link #countItems}
		* and the type <code>std::string<code> must match the stored HDF5-type.
		*/
		void writeAll ( const std::string* array );
	};

}

// C++ compilers want template functions source code available when the template functions are used
#include "Object.tpl"

#endif	/* _IO_HDF5_OBJECT_HPP_ */
