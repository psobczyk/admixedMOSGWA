#ifndef _SECTIONPARSLET_HPP_
#define _SECTIONPARSLET_HPP_

#include "NameParslet.hpp"

namespace parser {

	/** Configures a <code>string</code> variable to hold a section name. */
	class SectionParslet : public NameParslet {

		public:

		/** Construct an instance for a given section name holder variable.
		* @param variable to be configured.
		*/
		SectionParslet ( string &variable );

		virtual bool parse ( const char* &text );

		virtual const char* type ();
	};

}

#endif /* _SECTIONPARSLET_HPP_ */
