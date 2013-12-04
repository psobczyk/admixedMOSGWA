/********************************************************************************
 *	This file is part of the MOSGWA program code.				*
 *	Copyright ©2011–2013, Erich Dolejsi, Bernhard Bodenstorfer.		*
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

#ifndef LOG_HPP
#define LOG_HPP

#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>

#include "Parameter.hpp"
#include "Helpfull.hpp"

using namespace std;

/** defines the File-stream (out) Log */
extern ofstream LOG;

/** print string s to the Log-file */
void printLOG( const string s );

/** print the logo at the start of the programm */
void printStartScreen();

#endif
