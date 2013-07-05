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

#include <string>
#include <map>
#include <iostream>

#include "Parslet.hpp"
#include "ParsletTemplate.hpp"
#include "NameParslet.hpp"
#include "SectionParslet.hpp"
#include "VariableParslet.hpp"
#include "BoolParslet.hpp"
#include "IntParslet.hpp"
#include "DoubleParslet.hpp"
#include "StringParslet.hpp"
#include "ChoiceParslet.hpp"
#include "ConfigParser.hpp"

using namespace std;
using namespace parser;

/** Test */
int main () {
	const char* text;

	string c( "Section" );
	SectionParslet cp(c);
	cout << cp.type() << endl;
	cout << c << endl;
	cout << cp.get() << endl;
	cp.set( "Wo?" );
	cout << c << endl;
	cout << cp.get() << endl;
	text = "[Terminal] 1";
	cp.parse( text );
	cout << c << endl;
	cout << cp.get() << "\t" << text << endl;

	string v( "Variable" );
	VariableParslet vp(v);
	cout << vp.type() << endl;
	cout << v << endl;
	cout << vp.get() << endl;
	vp.set( "Wo?" );
	cout << v << endl;
	cout << vp.get() << endl;
	text = "Binal Kuubah";
	vp.parse( text );
	cout << v << endl;
	cout << vp.get() << "\t" << text << endl;

	bool b(0);
	BoolParslet bp( b );
	cout << bp.type() << endl;
	cout << b << endl;
	cout << bp.get() << endl;
	bp.set( 8 );
	cout << b << endl;
	cout << bp.get() << endl;
	text = "trueth";
	bp.parse( text );
	cout << b << endl;
	cout << bp.get() << "\t" << text << endl;

	int i(0);
	IntParslet ip( i );
	cout << ip.type() << endl;
	cout << i << endl;
	cout << ip.get() << endl;
	ip.set( 8 );
	cout << i << endl;
	cout << ip.get() << endl;
	text = "-1234-3.6";
	ip.parse( text );
	cout << i << endl;
	cout << ip.get() << "\t" << text << endl;

	double d(0);
	DoubleParslet dp( d );
	cout << dp.type() << endl;
	cout << d << endl;
	cout << dp.get() << endl;
	dp.set( 8 );
	cout << d << endl;
	cout << dp.get() << endl;
	text = "-1234.6e+5-7";
	dp.parse( text );
	cout << d << endl;
	cout << dp.get() << "\t" << text << endl;

	string s( "Radix" );
	StringParslet sp(s);
	cout << sp.type() << endl;
	cout << s << endl;
	cout << sp.get() << endl;
	sp.set( "Wo?" );
	cout << s << endl;
	cout << sp.get() << endl;
	text = "\"Abrakadabr\"a.";
	sp.parse( text );
	cout << s << endl;
	cout << sp.get() << "\t" << text << endl;

	StringParslet sp2( *new string( "Nathan" ) );	// Leads to mem leak
	cout << sp2.get() << endl;
	text = "\"Wombit\"0";
	sp2.parse( text );
	cout << sp2.get() << "\t" << text << endl;

	map< const string, int > choice;
	choice["one"] = 1;
	choice["two"] = 2;
	choice["three"] = 3;
	ChoiceParslet<int> icp( i, choice );
	cout << icp.type() << endl;
	cout << i << endl;
	cout << icp.get() << endl;
	icp.set( 2 );
	cout << i << endl;
	cout << icp.get() << endl;
	text = "threethousand";
	icp.parse( text );
	cout << i << endl;
	cout << icp.get() << "\t" << text << endl;

	ConfigParser cfp;
	cfp.declare( "section", "section", c );
	cfp.declare( "section", "variable", v );
	cfp.declare( "section", "boolean", b );
	cfp.declare( "section", "integer", i );
	cfp.declare( "section", "double", d );
	cfp.declare( "section", "string", s );
	cfp.parse( cin, "stdin" );
	cout << "section" << c << endl
		<< "variable" << v << endl
		<< "boolean" << b << endl
		<< "integer" << i << endl
		<< "double" << d << endl
		<< "string" << s << endl;
}
