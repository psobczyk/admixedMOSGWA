#ifndef _BOOLPARSLET_HPP_
#define _BOOLPARSLET_HPP_

#include "ParsletTemplate.hpp"

namespace parser {

	/** Configures a boolean variable with syntactic keys <code>false</code> and <code>true</code>. */
	class BoolParslet : public ParsletTemplate<bool> {

		public:

		/** Construct an instance for a given <code>bool</code> variable.
		* @param variable to be configured.
		*/
		BoolParslet ( bool &variable );

		virtual bool parse ( const char* &text );

		virtual const char* type ();
	};

}

#endif /* _BOOLPARSLET_HPP_ */
