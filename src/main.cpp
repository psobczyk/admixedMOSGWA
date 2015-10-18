/********************************************************************************
 *	This file is part of the MOSGWA program code.				*
 *	Copyright ©2011–2015, Erich Dolejsi, Bernhard Bodenstorfer.		*
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


// Piotr Sobczyk

#include "AModel.hpp"
#include <iostream>
#include "AMdata.hpp"
#include "linalg/AutoVector.hpp"

using namespace std;

/** Application entry point */
int main ( const int argn, const char *argv[] ) {
  std::cout << " _____ _____ _____ _____ _ _ _ _____ \n";
  std::cout << "|     |     |   __|   __| | | |  _  |\n";
  std::cout << "| | | |  |  |__   |  |  | | | |     |\n";
  std::cout << "|_|_|_|_____|_____|_____|_____|__|__|\n";
  std::cout << "Model Selection for Genom-Wide Associations\n";

  const  AMdata data;

  std::cout << data.getSnpNo ( ) << endl;
  
  AModel model ( data );

  std::cout <<  model.getBeta( 1 ) << endl;
  std::cout <<  model.getGamma( 1 ) << endl;

linalg::Vector wektor;
wektor[0] = 1;
wektor[1] = 2;

  return 0;
}

