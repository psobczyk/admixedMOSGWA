#ifndef _STRINGPARSLET_HPP_
#define _STRINGPARSLET_HPP_

#include "ParsletTemplate.hpp"

#include <string.h>
#include <string>

using namespace std;

namespace parser {

	/** Configures an <code>string</code> variable. */
	class StringParslet : public ParsletTemplate<string> {

		public:

		/** Construct an instance for a given <code>string</code> variable.
		* @param variable to be configured.
		*/
		StringParslet ( string &variable );

		virtual bool parse ( const char* &text );

		virtual const char* type ();
	};

}

#endif /* _STRINGPARSLET_HPP_ */
