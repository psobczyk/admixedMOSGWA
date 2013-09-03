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
	cout << "Build time: " << buildinfo::timestamp << endl;
	cout << endl;
}
