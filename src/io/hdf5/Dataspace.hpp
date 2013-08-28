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

#ifndef _IO_HDF5_DATASPACE_HPP_
#define _IO_HDF5_DATASPACE_HPP_

#include "Dataset.hpp"

namespace hdf5 {

	/** RAII wrapper for dataspace ids. */
	template<size_t D> class Dataspace : public Id {

		size_t
			size[D],
			volume;

		/** Construction helper for creation from existing data space. */
		void initFromId ();

		public:

		/** Create a <code>D</code>-dimensional rectangular data space. */
		Dataspace ( const std::string& name, const size_t size[D] );

		/** Wrap an already existing data space <code>hid_t</code> for certain destruction. */
		Dataspace ( const hid_t id, const std::string& name );

		/** Wrap the data space of the given data set. */
		Dataspace ( Dataset& dataset );

		/** Get the size of the given dimension. */
		size_t countDimensions ( const size_t dim ) const;

		/** Get number of cells in the data space. */
		size_t countItems () const;

		virtual ~Dataspace ();
	};

}

// C++ compilers want template functions source code available when the template functions are used
#include "Dataspace.tpl"

#endif	/* _IO_HDF5_DATASPACE_HPP_ */
