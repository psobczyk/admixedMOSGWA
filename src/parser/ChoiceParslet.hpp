#ifndef _CHOICEPARSLET_HPP_
#define _CHOICEPARSLET_HPP_

#include "ParsletTemplate.hpp"
#include <string>
#include <map>

using namespace std;

namespace parser {

	/** Configures a variable via syntactic choices taken from a {@link std::map}. */
	template <class T> class ChoiceParslet : public ParsletTemplate<T> {

		private:

		/** Maps possible syntactic choices to variable values */
		const map < const string, T > choices;

		/** Describes the kind of {@link Parslet} this is.
		* @see #type()
		*/
		string typeString;

		public:

		/** Construct an instance for a given variable and given map from choices to values.
		* @param variable to be configured.
		* @param choices maps syntactic keys to variable values
		*/
		ChoiceParslet ( T &variable, const map < const string, T > &choices );

		virtual bool parse ( const char* &text );

		virtual const char* type ();
	};

}

// C++ compilers want template functions source code available when the template functions are used
#include "ChoiceParslet.cpp"

#endif /* _CHOICEPARSLET_HPP_ */
