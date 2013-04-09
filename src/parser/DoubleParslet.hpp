#ifndef _DOUBLEPARSLET_HPP_
#define _DOUBLEPARSLET_HPP_

#include "ParsletTemplate.hpp"

namespace parser {

	/** Configures a <code>double</code> variable. */
	class DoubleParslet : public ParsletTemplate<double> {

		public:

		/** Construct an instance for a given <code>double</code> variable.
		* @param variable to be configured.
		*/
		DoubleParslet ( double &variable );

		virtual bool parse ( const char* &text );

		virtual const char* type ();
	};

}

#endif /* _DOUBLEPARSLET_HPP_ */
