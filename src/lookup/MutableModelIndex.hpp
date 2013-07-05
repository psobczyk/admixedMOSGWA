/********************************************************************************
 *	This file is part of the MOSGWA program code.				*
 *	Copyright ©2012–2013, Bernhard Bodenstorfer.				*
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

#ifndef MUTABLE_MODEL_INDEX_HPP
#define MUTABLE_MODEL_INDEX_HPP

#include "ModelIndex.hpp"

namespace lookup {

	/** A mutable {@link ModelIndex}. */
	class MutableModelIndex : public ModelIndex {

		public:

		/** Copy constructor, to avoid the default. */
		MutableModelIndex ( const MutableModelIndex& original );

		/** Copy from {@link ModelIndex} constructor. */
		MutableModelIndex ( const ModelIndex& original );

		/** Assignment operator, to avoid the default. */
		MutableModelIndex& operator= ( const MutableModelIndex& original );

		/** Assignment from {@link ModelIndex} operator. */
		MutableModelIndex& operator= ( const ModelIndex& original );
	};

}

#endif	/* MUTABLE_MODEL_INDEX_HPP */
