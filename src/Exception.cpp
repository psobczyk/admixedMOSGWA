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

#include "Exception.hpp"
#include <cstdio>
#include <cstdarg>

using namespace std;

const size_t MAX_MSG_SIZE = 1024;

Exception::Exception ( const char * format, ... ) {
	char buffer[MAX_MSG_SIZE];
	va_list arguments;
	va_start( arguments, format );
	const size_t needed = vsnprintf( buffer, MAX_MSG_SIZE, format, arguments );
	message = buffer;
	va_end( arguments );
}

const char * Exception::what () const throw () {
	return message.c_str();
}

Exception::~Exception () throw () {}
