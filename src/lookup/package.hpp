#ifndef PACKAGE_LOOKUP_HPP
#define PACKAGE_LOOKUP_HPP

#include "ModelIndex.hpp"

/** Facilities to store already calculated {@link Model} selection criteria for later reference
* to avoid unnecessary recalculation.
* @author Bernhard Bodenstorfer
*/
namespace lookup {

	/** Output to <code>cout</code> for testing and debugging. */
	void printModelIndex ( const ModelIndex& mi );

}

#endif	/* PACKAGE_LOOKUP_HPP */
