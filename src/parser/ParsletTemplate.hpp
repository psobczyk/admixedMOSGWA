#ifndef _PARSLETTEMPLATE_HPP_
#define _PARSLETTEMPLATE_HPP_

#include "Parslet.hpp"

namespace parser {

	/** Abstract class template to help implementing {@link Parslet}. */
	template <class T> class ParsletTemplate : public Parslet {

		private:

		/** The variable to be configured. */
		T &variable;

		public:

		/** Construct an instance for a given variable.
		* @param variable to be configured.
		*/
		ParsletTemplate ( T &variable );

		/** Get the current variable value. */
		T get () const;

		/** Set the value of the configured variable. */
		void set ( const T value );
	};

}

// C++ compilers want template functions source code available when the template functions are used
#include "ParsletTemplate.cpp"

#endif /* _PARSLETTEMPLATE_HPP_ */
