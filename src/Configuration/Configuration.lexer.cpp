/**
 * Implementation of private methods for
 * the 'Configuration' class.
 */

#include <fstream>
#include <iostream>
#include <string>
#include <sstream>
#include <vector>

#include <softlib/Configuration.h>
#include <softlib/SOFTLibException.h>

using namespace std;

void Configuration::initbuffer(const string& buf) {
	this->cnfbuffer = buf;
	this->cnfbufferpos = 0;
	this->cnfbufferlength = buf.length();
    this->linecounter = 1;
    this->charpos = 1;
	this->instring = false;
    this->ininclude = false;
    this->prevws = false;
}
/**
 * Obtain the next non-whitespace character.
 * If 'instring == true', all characters are
 * returned, even whitespace.
 */
char Configuration::nextc() {
	char c;
	if (this->cnfbufferpos >= this->cnfbufferlength) return 0;

	do {
		c = this->cnfbuffer[this->cnfbufferpos++];
		this->charpos++;

		if (this->instring) return c;
		// Make sure whitespace and comments are ignored
		switch (c) {
			case '#':
				// Ignore all characters til end-of-line (or file)
				do {
					if (this->cnfbufferpos >= this->cnfbufferlength)
						break;

					c = this->cnfbuffer[this->cnfbufferpos++];
				} while (c != 0 && c != '\n');

				// Handle the newline or EOF in a new iteration
				this->cnfbufferpos--;

				break;
			case '\n':
				this->linecounter++;
				this->charpos = 1;
                // Intended fall-through
                [[fallthrough]];
			case ' ':
			case '\t':
				if (this->prevws) break;
				else {
					this->prevws = true;
					return ' ';	/* All white-space is turned into regular space */
				}
			default:
				this->prevws = false;
				return c;
		}

	} while (c != 0 && this->cnfbufferpos < this->cnfbufferlength);

	return 0;
}
ConfigToken *Configuration::gettkn(string& s, ConfigToken *prev) {
	// First we check for simple tokens
	if (s == "=") {
		return new ConfigToken(CTKN_EQUALS, s, prev, this->linecounter, this->charpos, this->file);
	} else if (s == "{") {
		return new ConfigToken(CTKN_BLOCK_START, s, prev, this->linecounter, this->charpos, this->file);
	} else if (s == "}") {
		return new ConfigToken(CTKN_BLOCK_END, s, prev, this->linecounter, this->charpos, this->file);
	} else if (s == ";") {
		return new ConfigToken(CTKN_END_STATEMENT, s, prev, this->linecounter, this->charpos, this->file);
    }

	// Next, check if previous token was an equal sign or block type
	// (then we interpret this token as a value or block name)
	if (prev != NULL) {
		if (prev->GetType() == CTKN_EQUALS || prev->GetType() == CTKN_VALUE) {
			return new ConfigToken(CTKN_VALUE, s, prev, this->linecounter, this->charpos, this->file);
		} else if (prev->GetType() == CTKN_BLOCKTYPE) {
			return new ConfigToken(CTKN_BLOCKNAME, s, prev, this->linecounter, this->charpos, this->file);
		}
	}

	// Finally, we check if this is a block type
	if (IsBlockType(s)) {
		return new ConfigToken(CTKN_BLOCKTYPE, s, prev, this->linecounter, this->charpos, this->file);
	}
	
	// Else, it must be a name.
	return new ConfigToken(CTKN_NAME, s, prev, this->linecounter, this->charpos, this->file);
}
void Configuration::pushbuffer(string& buffer, vector<ConfigToken*>& tkns) {
	if (buffer.empty()) return;

    ConfigToken *ct;
	if (tkns.empty())
		//tkns.push_back(gettkn(buffer, NULL));
        ct = gettkn(buffer, nullptr);
	else
		//tkns.push_back(gettkn(buffer, &(tkns.back())));
        ct = gettkn(buffer, tkns.back());

    tkns.push_back(ct);
}

/**
 * Continue reading from the file with name 'name'
 * as if it was part of the current stream. When
 * reading is finished, returns the state to the
 * original file.
 *
 * name: Name of configuration file to read.
 */
vector<ConfigToken*> Configuration::include_file(const string& name) {
    ifstream cfile(name);
    stringstream buf;
    string oldname, oldbuf;
    size_t bufpos, buflen, linec, charc;

    if (!cfile.is_open())
        throw ConfigurationException(
            this->linecounter, this->charpos, this->filename,
            "Unable to open included configuration file '%s'.",
            name.c_str()
        );

    // Save current state
    oldname = this->filename;
    this->filename = name;
    
    oldbuf = this->cnfbuffer;
    bufpos = this->cnfbufferpos;
    buflen = this->cnfbufferlength;
    linec  = this->linecounter;
    charc  = this->charpos;

    // Load file
    buf << cfile.rdbuf();
    cfile.close();

    // Lex the file contents
    string contents = buf.str();
    vector<ConfigToken*> tkns = Lex(contents);

    // Remove end-of-file token
    tkns.pop_back();

    // Restore old state
    this->charpos = charc;
    this->linecounter = linec;
    this->cnfbufferlength = buflen;
    this->cnfbufferpos = bufpos;
    this->cnfbuffer = oldbuf;
    this->filename = oldname;

    return tkns;
}

/**
 * Perform lexical analysis on the configuration script.
 */
vector<ConfigToken*> Configuration::Lex(const string& str) {
	vector<ConfigToken*> tkns;
	ConfigToken tkn;
	char c;
	string buffer;

	initbuffer(str);

	while ((c=nextc())) {
		// If we're in string, append all characters...
		if (this->instring && c != '"') {
			buffer += c;
			continue;
		} else if (this->ininclude && c != '>') {
            buffer += c;
            continue;
        }

		switch (c) {
			case '"': {
				pushbuffer(buffer, tkns);
				buffer.clear();
				this->instring = !this->instring;
			} break;
            case '<': {
                pushbuffer(buffer, tkns);
                buffer.clear();
                this->ininclude = true;
            } break;
            case '>': {
                this->ininclude = false;
                vector<ConfigToken*> t = include_file(buffer);
                tkns.insert(tkns.end(), t.begin(), t.end());
                buffer.clear();
            } break;
			case '=':
			case '{':
			case '}':
			case ';': {
				pushbuffer(buffer, tkns);
				// Clear buffer and replace with current character
				buffer.assign(1, c);
				pushbuffer(buffer, tkns);
				buffer.clear();
			} break;
			case ',':		// Treat commas as spaces for endless syntactic possibilities
			case ' ': {
				pushbuffer(buffer, tkns);
				buffer.clear();
            } break;
			default: {
				buffer += c;
            } break;
		}
	}

	// Push remaining buffer
    if (!buffer.empty())
        pushbuffer(buffer, tkns);

	// Append an End-of-File token
	if (tkns.empty())
		tkns.push_back(new ConfigToken(CTKN_EOF, "<EOF>", NULL, this->linecounter, this->charpos, this->file));
	else
		tkns.push_back(new ConfigToken(CTKN_EOF, "<EOF>", tkns.back(), this->linecounter, this->charpos, this->file));
	
	return tkns;
}

