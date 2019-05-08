/**
 * Implementation of the configuration script interpreter.
 */

#include <cstdlib>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include <softlib/Configuration.h>
#include <softlib/Configuration/ConfigurationScript.h>

using namespace std;

/**
 * Load configuration from a file.
 *
 * filename: Name of file to load configuration from.
 */
void ConfigurationScript::FromFile(const string& filename) {
	if (!this->filename.empty())
		throw ConfigurationException("A configuration file has already been opened with this ConfigurationScript object.");

	ifstream cfile(filename);
	stringstream buffer;

	if (!cfile.is_open())
		throw ConfigurationException("Unable to open configuration file: " + filename);
	else this->filename = filename;

	buffer << cfile.rdbuf();
	cfile.close();

	string contents = buffer.str();
	FromString(contents, filename);
}

/**
 * Load configuration from stdin.
 */
void ConfigurationScript::FromStdin() {
#	define _FROMSTDIN_BUFSIZE 1000
	char *buf = NULL;
	size_t bufl = 0;

	while (!cin.eof()) {
		buf = (char*)realloc(buf, sizeof(char)*(bufl+_FROMSTDIN_BUFSIZE+1));
		cin.read(buf+bufl, _FROMSTDIN_BUFSIZE);
		bufl += cin.gcount();
	}

	buf[bufl] = 0;

	string b = string(buf);
	FromString(b, "stdin");
}

/**
 * Load configuration from a string.
 *
 * config:         String containing configuration text
 * (optional) src: Name of source of this string (i.e. filename).
 */
ConfigBlock& ConfigurationScript::FromString(const string& config, const string& src) {
    this->file = src;

	// Convert the string to a stream of tokens
	vector<ConfigToken*> tkns = Lex(config);

	// Interpret and generate root block
	this->root = Interpret(tkns);

    // Delete config tokens
    for (vector<ConfigToken*>::iterator it = tkns.begin(); it != tkns.end(); it++)
        delete *(it);

	return this->root;
}


