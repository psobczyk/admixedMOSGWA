#ifndef _INTPARSLET_HPP_
#define _INTPARSLET_HPP_

#include "ParsletTemplate.hpp"

namespace parser {

	/** Configures an <code>int</code> variable. */
	class IntParslet : public ParsletTemplate<int> {

		public:

		/** Construct an instance for a given <code>int</code> variable.
		* @param variable to be configured.
		*/
		IntParslet ( int &variable );

		virtual bool parse ( const char* &text );

		virtual const char* type ();
	};

}

#endif /* _INTPARSLET_HPP_ */
