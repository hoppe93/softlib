/**
 * Implementation of the 'Configuration' class, for
 * loading settings from a configuration file.
 */

#include <cstdlib>
#include <string>
#include <vector>

#include <softlib/Configuration.h>

using namespace std;

/****************************************************************
 ****************************************************************
 **
 ** PUBLIC METHODS
 **
 ****************************************************************
 ****************************************************************
 */
/**
 * Constructors & Destructors
 */
Configuration::Configuration() {
	confblock_t t = RegisterBlockType("/root");
	this->root = ConfigBlock(t, "<root>");
}
Configuration::Configuration(const Configuration& conf) {
    this->filename = conf.GetFilename();
    this->root = ConfigBlock(conf.GetRootBlock());
    this->blocktypes = conf.GetBlockTypes();
    this->errorflag = conf.HasError();
}
Configuration::~Configuration() { }

/**
 * Check if the given string is the name of
 * any block type.
 *
 * name: Name of potential block type.
 */
bool Configuration::IsBlockType(const string& name) {
	for (vector<string>::iterator it = blocktypes.begin(); it != blocktypes.end(); it++) {
		if (name == *it)
			return true;
	}

	return false;
}

/**
 * Merge this Configuration object with the given
 * Configuration object. This is equivalent to
 * mergeing the root ConfigBlocks of each
 * Configuration objects.
 *
 * conf:     Configuration object to merge into this object.
 * allowNew: If true, allows completely new settings
 *           to be defined by 'conf'. Otherwise, returns
 *           a list of strings containing the names of
 *           all settings that were exclusively defined
 *           by 'conf'.
 *
 * RETURNS a null pointer if allowNew = true, or if
 * allowNew = false and all settings defined in 'conf'
 * were also defined in this Configuration. If
 * allowNew = false and some settings were exclusively
 * defined in 'conf', a list containing the names of
 * all settings exclusively defined 'conf' is returned.
 */
vector<string> *Configuration::Merge(Configuration& conf, bool allowNew) {
    ConfigBlock cb = conf.GetRootBlock();
    return this->root.Merge(cb, allowNew);
}

/**
 * Register a ConfigBlock type.
 * Returns the ID of the newly registered block.
 *
 * name: Name of the block type to register.
 */
confblock_t Configuration::RegisterBlockType(const string& name) {
	// Make sure type has not already been registered
	if (IsBlockType(name))
		throw SOFTLibException("The block-type '" + name + "' has already been registered.");

	blocktypes.push_back(name);
	return blocktypes.size()-1;
}

/**
 * Set the error flag, indicating that one or
 * more errors occured.
 */
void Configuration::SetErrorFlag() { this->errorflag = true; }

/**
 * Given the name of a block type, returns
 * the corresponding 'confblock_t' identifier.
 *
 * name: Name of the block type
 */
confblock_t Configuration::StringToBlockType(const string& name) {
	confblock_t i=0;
	for (vector<string>::iterator it = blocktypes.begin(); it != blocktypes.end(); it++, i++) {
		if (*it == name) return i;
	}

	return 0;		// TODO maybe throw an exception?
}

