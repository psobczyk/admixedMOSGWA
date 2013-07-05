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

#ifndef LOOKUP_RESULT_STORE_HPP
#define LOOKUP_RESULT_STORE_HPP

#include <map>
#include "MutableModelIndex.hpp"

namespace lookup {

	/** Stores already calculated model selection criteria to find the best and avoid duplicate calculation.
	* TODO<BB>: Take it to the next step: Prepare this class for multi-threading.
	* Then the get-MSC-methods may need to return a synchronising handle instead of a <code>double</code>
	* to access the result or otherwise coordinate its calculation.
	*/
	class ResultStore {

		/** Stores results. */
		std::map<ModelIndex,double> linearMSC, logisticMSC;

		/** Stores optimal model indices. */
		MutableModelIndex optimalLinearModel, optimalLogisticModel;

		/** Stores optimal model selection criteria. */
		double optimalLinearMSC, optimalLogisticMSC;

		/** Forbid assignment operator. */
		ResultStore& operator= ( const ResultStore& );

		/** Forbid copy constructor. */
		ResultStore ( const ResultStore& );

		public:

		/** Set the model selection criterion for linear regression.
		* @see getLinearMSC( ModelIndex& );
		*/
		void setLinearMSC ( const ModelIndex& index, const double msc );

		/** Retrieve the model selection criterion for linear regression.
		* @returns <code>NaN</code> if the value is not (yet) available.
		* TODO<BB>: Future architecture outlook:
		* For multi-threading, an object should be returned
		* which either gives the <code>double</code> value
		* or can provide a lock object which is either locked (i.e. somebody is calculating)
		* or free (an offer to calculate and let potential other interested threads wait for the result).
		*/
		double getLinearMSC ( const ModelIndex& index ) const;

		/** Retrieve so far optimal set of variables with regard to the linear selection criterion.
		* Precondition: {@link ResultStore::setLinearMSC} must have been called before.
		*/
		const ModelIndex getOptimalLinearModel () const;

		/** Set the model selection criterion for logistic regression.
		* @see getLogisticMSC( ModelIndex& );
		*/
		void setLogisticMSC ( const ModelIndex& index, const double msc );

		/** Retrieve the model selection criterion for logistic regression.
		* @see getLinearMSC( ModelIndex& );
		*/
		double getLogisticMSC ( const ModelIndex& index ) const;

		/** Retrieve so far optimal set of variables with regard to the logistic selection criterion.
		* Precondition: {@link ResultStore::setLogisticMSC} must have been called before.
		*/
		const ModelIndex getOptimalLogisticModel () const;
	};

}

#endif	/* LOOKUP_RESULT_STORE_HPP */
