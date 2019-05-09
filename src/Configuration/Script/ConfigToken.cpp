/**
 * Implementation of the 'ConfigToken' class,
 * representing a lexical object in the
 * configuration file.
 */

#include <string>
#include <softlib/Configuration.h>
#include <softlib/Configuration/ConfigurationScript.h>

using namespace std;

ConfigToken::ConfigToken() { }
ConfigToken::ConfigToken(const ConfigToken &ct) {
    this->type = ct.GetType();
    this->value = ct.GetValue();
    this->parent = ct.Previous();
    this->line = ct.GetLine();
    this->charpos = ct.GetCharpos();
    this->filename = ct.GetFilename();

    if (parent != nullptr)
        this->parent->SetNext(this);
}
ConfigToken::ConfigToken(conftoken_t type, const string& value, ConfigToken *parent, size_t line, size_t charpos, const string& filename) {
	this->type = type;
	this->value = value;
	this->parent = parent;
    this->line = line;
    this->charpos = charpos;
    this->filename = filename;

	if (parent != nullptr)
		this->parent->SetNext(this);
}
ConfigToken::~ConfigToken() { }

