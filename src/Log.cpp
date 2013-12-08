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

#include "Log.hpp"
#include "buildinfo.hpp"

/** get and convert time */
string timestamp () {
	const time_t ltime=time(NULL);	//get current calendar time
	// TODO<BB>: inspect localtime memory allocation
	const tm* now = localtime(&ltime);	// struct for the day, year....

	// format time information
	const string outTime = int2strPadWith( now->tm_mday, 2, '0' ) + "."
		+ int2strPadWith( now->tm_mon+1, 2, '0' ) + "."
		+ int2str( now->tm_year+1900 ) + " "
		+ int2strPadWith( now->tm_hour, 2, '0' ) + ":"
		+ int2strPadWith( now->tm_min, 2, '0' ) + ":"
		+ int2strPadWith( now->tm_sec, 2, '0' );

	return outTime;
}

/** Print string to the Log-file */
void printLOG ( const string s ) {
	// print to logfile
	const string time = timestamp();
	LOG << time << ": " << s << endl;
	LOG.flush();

	// print to screen
	if ( ! parameter.silent ) {
		cout << time << ": " << s << endl;
		cout.flush();
	}
}

/** Print ASCII-Art-Logo for MOSGWA. */
void printStartScreen () {
	cout << endl;
	cout << " _____ _____ _____ _____ _ _ _ _____ " << endl;
	cout << "|     |     |   __|   __| | | |  _  |\tModel Selection" << endl;
	cout << "| | | |  |  |__   |  |  | | | |     |\tfor Genome-wide" << endl;
	cout << "|_|_|_|_____|_____|_____|_____|__|__|\tAssociations" << endl;
	//cout << "Build time: " << buildinfo::timestamp << endl;
	cout << endl;
}
