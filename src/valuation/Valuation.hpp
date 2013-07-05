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

#ifndef VALUATION_VALUATION_HPP
#define VALUATION_VALUATION_HPP

#include "../lookup/ModelIndex.hpp"

/** Valuation algorithms. These judge how good a {@link lookup::ModelIndex} matches given {@link io::Input} data.
* @author Bernhard Bodenstorfer
*/
namespace valuation {

	/** Abstract base class for valuation algorithms. */
	class Valuation {

		public:

		/** Judge a given model. */
		virtual double valuate ( const lookup::ModelIndex& modelIndex ) = 0;

		/** Destruct and potentially free resources. */
		virtual ~Valuation ();
	};

}

#endif	/* VALUATION_VALUATION_HPP */
