#ifndef _NAMEPARSLET_HPP_
#define _NAMEPARSLET_HPP_

#include "ParsletTemplate.hpp"
#include <string>

using namespace std;

namespace parser {

	/** Abstract base class for parsing variable and section names */
	class NameParslet : public ParsletTemplate<string> {

		protected:

		/** Parse a variable name into a <code>string</code>. */
		string parseName ( const char* &text );

		public:

		/** Configure a <code>string</code> name of a variable.
		* @param variable to be configured.
		*/
		NameParslet ( string &variable );
	};

}

#endif /* _NAMEPARSLET_HPP_ */
