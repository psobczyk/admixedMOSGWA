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
