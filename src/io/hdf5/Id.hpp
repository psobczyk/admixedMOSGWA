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

#ifndef _IO_HDF5_ID_HPP_
#define _IO_HDF5_ID_HPP_

#include <string>
#include <hdf5.h>

/** Resource-acquisition is initialisation (RAII) wrappers for accessing HDF5 files. */
namespace hdf5 {

	/** Base class RAII wrapper for <code>hid_t</code>. */
	class Id {

		/** The wrapped identifier. */
		const hid_t id;

		/** Name from opening and for exception messages. */
		const std::string name;

		/** Forbid copy constructor. */
		Id ( const Id& id );

		/** Forbid assigmnent operator. */
		Id& operator= ( const Id& id );

		public:

		/** Take ownership of and wrap the given identifier. */
		Id ( const hid_t id, const std::string& name );

		/** Access the stored HDF5 identifier. */
		hid_t getId () const;

		/** Access the name (used for exception messages). */
		const std::string& getName () const;

		/** Release the wrapped object using the appropriate close. */
		virtual ~Id () = 0;
	};

}

#endif	/* _IO_HDF5_ID_HPP_ */
