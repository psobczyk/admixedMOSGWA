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
