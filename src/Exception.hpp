#ifndef EXCEPTION_HPP
#define EXCEPTION_HPP

#include <string>
#include <exception>

/** Error reporting mechanism.
* @author Bernhard Bodenstorfer
* @see std::exception
*/
class Exception : public std::exception {

	std::string message;

	public:

	/** Initialise with a given message.
	* The arguments work like for <code>printf</code>.
	*/
	Exception ( const char * format, ... );

	/** Retrieve the message. */
	virtual const char * what () const throw ();

	virtual ~Exception () throw ();

};

#endif	/* EXCEPTION_HPP */
