#ifndef _VARIABLEPARSLET_HPP_
#define _VARIABLEPARSLET_HPP_

#include "NameParslet.hpp"

namespace parser {

	/** Configures a <code>string</code> variable to hold a variable name. */
	class VariableParslet : public NameParslet {

		public:

		/** Construct an instance for a given variable name holder variable.
		* @param variable to be configured.
		*/
		VariableParslet ( string &variable );

		virtual bool parse ( const char* &text );

		virtual const char* type ();
	};

}

#endif /* _VARIABLEPARSLET_HPP_ */
