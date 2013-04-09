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
