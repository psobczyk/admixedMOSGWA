
/********************************************************************************
 *	This file is part of the MOSGWA program code.				*
 *	Copyright ©2015, Bernhard Bodenstorfer, Piotr Sobczyk.  		*
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

#include "AdmixtureSearch.hpp"
#include "MData.hpp"
#include "Model.hpp"

using namespace std;
using namespace logging;

AdmixtureSearch::AdmixtureSearch () 
  : data(), tests( data ), model( data ) 
{
  logger->info( "Prepared setting for elaborated greedy search" );
}
  
/** Run the search. */
virtual void AdmixtureSearch::run () {
  logger->info( "Start greedy search for admixture model" );
  
  // calculate single marker test (need for modelselection)
  data.calculateIndividualTests();

  model.computeRegression();
  data.setLL0M( model.getMJC() );
  data.selectModel(
		   &model,
		   &tests,
		   parameter->PValueBorder,
		   parameter->maximalModelSize,
		   Parameter::selectionCriterium_mBIC_firstRound
		   );
  model.printModel( "First round result" );
  data.selectModel(
		   &model,
		   &tests,
		   5000,
		   parameter->maximalModelSize,
		   parameter->selectionCriterium
		   );
  logger->info( "Finish greedy search for adherancy model" );
}
  
/** Retrieve the winning model.
 * This should be called only after {@link run}.
 */
virtual const AdmixtureSearch::AModel* result () {
	return &model;
}
