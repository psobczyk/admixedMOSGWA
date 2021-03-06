/********************************************************************************
 *	This file is part of the MOSGWA program code.				*
 *	Copyright ©2012–2013, Bernhard Bodenstorfer.				*
 *										*
 *	This program is free software; you can redistribute it and/or modify	*
 *	it under the terms of the GNU General Public License as published by	*
 *	the Free Software Foundation; either version 3 of the License, or	*
 *	(at your option) any later version.					*
 *										*
 *	This program is distributed in the hope that it will be useful,		*
 *	but WITHOUT ANY WARRANTY; without even the implied warranty of		*
 *	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.			*
 *	See the GNU General Public License for more details.			*
 ********************************************************************************/

#ifndef EXCEPTION_HPP
#define EXCEPTION_HPP

#include <string>
#include <exception>

extern const size_t MAX_MSG_SIZE;

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
