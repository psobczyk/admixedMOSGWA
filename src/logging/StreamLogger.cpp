/********************************************************************************
 *	This file is part of the MOSGWA program code.				*
 *	Copyright ©2015, Bernhard Bodenstorfer.					*
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

#include "StreamLogger.hpp"

using namespace std;

namespace logging {

	StreamLogger::StreamLogger ( ostream& stream ) {
		setStream( stream );
	}

	void StreamLogger::setStream ( ostream& stream ) {
		this->stream = &stream;
	}

	void StreamLogger::write ( const char * const text ) {
		( *stream ) << text << endl;
	}

	StreamLogger::~StreamLogger () {}

}
