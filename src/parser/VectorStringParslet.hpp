#ifndef _VECTORSTRINGPARSLET_HPP_
#define _VECTORSTRINGPARSLET_HPP_

#include "ParsletTemplate.hpp"
#include <vector>
#include <string.h>
#include <string>

using namespace std;

namespace parser {

	/** Configures an <code>string</code> variable. */
	class VectorStringParslet : public ParsletTemplate<vector<string> > {

		public:

		/** Construct an instance for a given <code>string</code> variable.
		* @param variable to be configured.
		*/
		VectorStringParslet (vector< string> &variable );

		virtual bool parse ( const char* &text );

		virtual const char* type ();

		/** Advance the text pointer to the next non-whitespace character. */
		void skipWhitespace ( const char* &text );
	};

}

#endif /* _VECTORSTRINGPARSLET_HPP_ */
