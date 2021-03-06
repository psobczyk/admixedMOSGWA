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

#include "FileLogger.hpp"

using namespace std;

namespace logging {

	FileLogger::FileLogger ( const char * const filename )
	:
		file( filename, ios_base::app )
	{
		file.exceptions ( ofstream::eofbit | ofstream::failbit | ofstream::badbit );
		setStream( file );
	}

	FileLogger::~FileLogger () {
		file.close();
	}

}
