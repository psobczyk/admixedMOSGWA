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

#include "../io/Input.hpp"
#include "../lookup/ResultStore.hpp"

/** Search algorithms. These try to find optimal {@link Model}s for given {@link Input} data.
* @author Bernhard Bodenstorfer
*/
namespace search {

	/** Abstract base class for search algorithms. */
	class Search {

		protected:

		/** The source of data with regard to which to search. */
		io::Input& input;

		/** Already calculated model selection criteria. */
		lookup::ResultStore& resultStore;

		/** Construct with reference to data and result store. */
		Search ( io::Input& input, lookup::ResultStore& resultStore );

		public:

		/** Perform the search. */
		virtual void run () = 0;

		/** Destruct and potentially free resources. */
		virtual ~Search ();
	};

}

#endif	/* SEARCH_SEARCH_HPP */
