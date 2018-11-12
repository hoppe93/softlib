
#include <cstdarg>
#include <cstdio>
#include <cstdlib>
#include <exception>
#include <string>
#include <softlib/SOFTLibException.h>

using namespace std;

/**
 * Constructor.
 */
SOFTLibException::SOFTLibException() { }
SOFTLibException::SOFTLibException(const string& msg) {
	this->errormsg = msg;
}
const char *SOFTLibException::what() const throw() {
	return this->errormsg.c_str();
}
const string& SOFTLibException::whats() const throw() {
	return this->errormsg;
}
