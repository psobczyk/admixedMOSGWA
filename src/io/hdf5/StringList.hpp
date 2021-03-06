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

#ifndef _IO_HDF5_STRINGLIST_HPP_
#define _IO_HDF5_STRINGLIST_HPP_

#include "Object.hpp"
#include <string>
#include <hdf5.h>

namespace hdf5 {

	/** One-dimensional table of strings in HDF5. */
	class StringList : public Object<1> {

		public:

		/** Create list in HDF5 file. */
		StringList ( File& file, const std::string& objectPath, const size_t size );

		/** Open existing list in HDF5 file. */
		StringList ( File& fileId, const std::string& objectPath );

		/** Retrieve the number of strings. */
		size_t countDimensions () const;
	};

}

#endif	/* _IO_HDF5_STRINGLIST_HPP_ */
