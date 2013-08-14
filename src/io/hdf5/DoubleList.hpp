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

#ifndef _IO_HDF5_DOUBLELIST_HPP_
#define _IO_HDF5_DOUBLELIST_HPP_

#include "Object.hpp"
#include <string>
#include <hdf5.h>

namespace hdf5 {

	/** One-dimensional table of numbers in HDF5. */
	class DoubleList : public Object<1> {

		public:

		DoubleList ( File& fileId, const std::string& objectPath );

		/** Retrieve the number of strings. */
		size_t countDimensions () const;
	};

}

#endif	/* _IO_HDF5_DOUBLELIST_HPP_ */
