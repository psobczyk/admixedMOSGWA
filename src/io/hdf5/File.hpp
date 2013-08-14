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

#ifndef _IO_HDF5_FILE_HPP_
#define _IO_HDF5_FILE_HPP_

#include "Id.hpp"

namespace hdf5 {

	/** RAII wrapper for file ids. */
	class File : public Id {

		public:

		File ( const std::string& path, const bool create = false );

		virtual ~File ();
	};

}

#endif	/* _IO_HDF5_FILE_HPP_ */
