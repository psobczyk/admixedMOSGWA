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

#ifndef _IO_HDF5OBJECT_HPP_
#define _IO_HDF5OBJECT_HPP_

#include <string>
#include <hdf5.h>

namespace io {

	/** Resource-acquisition is initialisation (RAII) wrapper for <code>hid_t</code>. */
	struct Hdf5Id {

		/** The wrapped identifier. */
		const hid_t id;

		/** Name from opening and for exception messages. */
		const std::string name;

		/** Take ownership of and wrap the given identifier. */
		Hdf5Id ( const hid_t id, const std::string& name );

		/** Convenience access to the name. */
		const char * getName () const;

		/** Release the wrapped object using the appropriate close. */
		virtual ~Hdf5Id () = 0;
	};

	/** RAII wrapper for file ids. */
	struct Hdf5FileId : public Hdf5Id {
		Hdf5FileId ( const std::string& path );
		virtual ~Hdf5FileId ();
	};

	/** RAII wrapper for dataset ids. */
	struct Hdf5DatasetId : public Hdf5Id {
		Hdf5DatasetId ( Hdf5FileId& fileId, const std::string& name );
		virtual ~Hdf5DatasetId ();
	};

	/** RAII wrapper for datatype ids. */
	struct Hdf5DatatypeId : public Hdf5Id {

		/** Wrap an already existing data type <code>hid_t</code> for certain destruction. */
		Hdf5DatatypeId ( const hid_t id, const std::string& name );

		/** Wrap the data type of the given data set. */
		Hdf5DatatypeId ( Hdf5DatasetId& datasetId );

		/** Retrieve the data type size. */
		size_t size ();

		/** Tell whether the data type is variable sized string. */
		bool isVariableString ();

		virtual ~Hdf5DatatypeId ();
	};

	/** RAII wrapper for dataspace ids. */
	struct Hdf5DataspaceId : public Hdf5Id {

		/** Wrap an already existing data space <code>hid_t</code> for certain destruction. */
		Hdf5DataspaceId ( const hid_t id, const std::string& name );

		/** Wrap the data space of the given data set. */
		Hdf5DataspaceId ( Hdf5DatasetId& datasetId );

		virtual ~Hdf5DataspaceId ();
	};

	/** Practical access adapter for HDF5 data. */
	template<size_t D> class Hdf5Object {

		protected:

		Hdf5DatasetId datasetId;
		Hdf5DatatypeId datatypeId;
		Hdf5DataspaceId dataspaceId;

		/** Volume is total array size, size holds sizes of individual dimensions. */
		size_t
			volume,
			size[D];

		public:

		/** Wraps the object in a HDF5 file identified by the given path. */
		Hdf5Object ( Hdf5FileId& fileId, const std::string& objectPath );

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
	};

	struct Hdf5StringList : public Hdf5Object<1> {
		Hdf5StringList ( Hdf5FileId& fileId, const std::string& objectPath );

		/** Retrieve the number of strings. */
		size_t countDimensions () const;
	};

	struct Hdf5DoubleList : public Hdf5Object<1> {
		Hdf5DoubleList ( Hdf5FileId& fileId, const std::string& objectPath );

		/** Retrieve the number of strings. */
		size_t countDimensions () const;
	};

	struct Hdf5DoubleTable : public Hdf5Object<2> {
		Hdf5DoubleTable ( Hdf5FileId& fileId, const std::string& objectPath );

		/** Retrieve the number of rows. */
		size_t countRows () const;

		/** Retrieve the number of columns. */
		size_t countColumns () const;

		/** Retrieve one row. */
		void readRow ( const size_t row, double* array );
	};

}

// C++ compilers want template functions source code available when the template functions are used
#include "Hdf5Object.tpl"

#endif	/* _IO_HDF5OBJECT_HPP_ */
