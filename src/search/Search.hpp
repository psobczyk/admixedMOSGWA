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

#ifndef SEARCH_SEARCH_HPP
#define SEARCH_SEARCH_HPP

#include "../logging/Logger.hpp"
#include "../Parameter.hpp"
#include "../Model.hpp"

/** Search algorithms. These try to find optimal {@link Model}s for given {@link io::Input} data.
* @author Bernhard Bodenstorfer
*/
namespace search {

	/** Abstract base class for search algorithms. */
	class Search {

		protected:

		/** Construct search machinery. */
		Search ();

		public:

		/** Perform the search. */
		virtual void run () = 0;

		/** Retrieve a pointer to the winning model for copying and then further processing.
		* This method should be called only after {@link run}.
		* The returned pointer will most likely point to an internally stored {@link Model}.
		* It may not be valid after <code>this</code> {@link Search} is destructed.
		* Since <code>this</code> is responsible,
		* clients calling this method should not try to deallocate the memory pointed to.
		* @returns a pointer to the best found model as yet.
		*/
		virtual const Model* result () = 0;

		/** Destruct and potentially free resources. */
		virtual ~Search ();
	};

}

#endif	/* SEARCH_SEARCH_HPP */
