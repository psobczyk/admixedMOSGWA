/********************************************************************************
 *	This file is part of the MOSGWA program code.				*
 *	Copyright ©2015, Bernhard Bodenstorfer.					*
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

#ifndef SEARCH_EGREEDY_PACKAGE_HPP
#define SEARCH_EGREEDY_PACKAGE_HPP

#include "../Search.hpp"
#include "../../MData.hpp"
#include "../../Model.hpp"

/** Elaborated greedy model search.
* @author Erich Dolejsi, Bernhard Bodenstorfer
*/
namespace egreedy {

	/** Implements an elaborate greedy model search algorithm.
	*/
	class ElaboratedGreedy : public search::Search {

		/** Data needed by the model below. */
		MData data;

		/** The best model yet found. */
		Model model;

		public:

		/** Set up the search environment. */
		ElaboratedGreedy ();

		/** Run the search. */
		virtual void run ();

		/** Retrieve the winning model.
		* This should be called only after {@link run}.
		*/
		virtual const Model* result ();
	};

}

#endif
