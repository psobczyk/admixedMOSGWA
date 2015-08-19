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

#ifndef SEARCH_GENETICA_PACKAGE_HPP
#define SEARCH_GENETICA_PACKAGE_HPP

#include "../Search.hpp"

/** Genetic algorithm linear model search.
* @author Artur Gola, Bernhard Bodenstorfer
*/
namespace genetica {

	/** Implements a linear model search genetic algorithm.
	* Parameters used from <code>parameter</code>:
	* @param outNo - it is number of output.
	* It is used to run GA in a loop and save all the results in the different files.
	* @param modelsFileName - it is the file name which contains information about an initial population.
	* If a modelsFileName is empty, GA creates a new initial population
	*/
	class GeneticAlgorithm : public search::Search {

		public:

		/** Set up the search environment. */
		GeneticAlgorithm ();

		/** @brief Runs genetic algorithm as a model selection method. */
		virtual void run ();

	};

}

#endif
