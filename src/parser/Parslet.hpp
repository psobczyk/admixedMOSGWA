#ifndef _PARSLET_HPP_
#define _PARSLET_HPP_

namespace parser {

	/** Abstract class to parse a value according to its expected type. */
	class Parslet {

		public:

		/** Parse from a string.
		* @param text is parsed
		* and upon success advanced to after the last character parsed
		* @returns whether anything has been be parsed,
		* i.e. whether <code>text</code> has been moved
		*/
		virtual bool parse ( const char* &text ) = 0;

		/** Informs about the expected variable type for log messages. */
		virtual const char* type () = 0;

		/** Empty placeholder for valuable subclass destructors. */
		virtual ~Parslet ();
	};

}

#endif /* _PARSLET_HPP_ */
