#include "VectorStringParslet.hpp"
#include  <iostream>
using namespace std;
namespace parser {

	VectorStringParslet::VectorStringParslet ( vector<string> &variable ) : ParsletTemplate<vector<string> >( variable ) {}

	bool VectorStringParslet::parse ( const char* &text ) {
		vector<string> strv;
		cerr<<"1:";
		if ( '{' != *text ) return false;
		text++;
			skipWhitespace(text); //erase withespace
                 cerr<<"0:"<<text<<endl;
		 bool debug1='}'!= *text;
cerr<< debug1 <<endl;
		while ('}'!= *text )
		{//skipWhitespace(text); 
		if ( '"' != *text ) return false;
		cerr<<"2:";
		text++;
		cerr<<"3:"<<text<<endl;
		const char* endQuote = strchr( text, '"' );
		if ( endQuote == text ) return false;
		 strv.push_back(string( text , endQuote - text /* - 1*/ ));
		 cerr<<"5:"<<string( text , endQuote - text /*- 1*/ )<<endl;
		text = endQuote + 1;
cerr<<"3:"<<text<<endl;
                skipWhitespace(text);
                if ( ',' == *text)
		         text++; //skip , 
		if ('}' == *text)
		;	//to nothing
                skipWhitespace(text); //maybe withespaces here!
		}
      text++;
		set(strv); 	cerr<<"2:";
		return true;
	}

	const char* VectorStringParslet::type () {
		return "vectorstring";
	}

	void VectorStringParslet::skipWhitespace ( const char* &text ) {
		while ( isspace( *text ) ) ++text;
	}


}	
