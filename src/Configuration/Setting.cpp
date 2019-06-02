
#include <string>
#include <vector>
#include <iostream>

#include <softlib/Configuration.h>
#include <softlib/SOFTLibException.h>

using namespace std;

Setting::~Setting() {}

/**
 * Get the name of this setting.
 */
const string& Setting::GetName() const { return this->name; }

/**
 * Check if this setting has the given name
 *
 * name: Name to compare to
 */
bool Setting::HasName(const string& name) const { return (this->name == name); }

/**
 * Converts a string to a boolean value, if possible.
 * Otherwise, throws an exception.
 *
 * val: String to convert.
 */
bool Setting::StringToBool(const string& value) {
	if (value == "yes" || value == "Yes" || value == "YES" ||
		value == "true" || value == "True" || value == "TRUE")
		return true;
	else if (value == "no" || value == "No" || value == "NO" ||
			 value == "false" || value == "False" || value == "FALSE")
		return false;
	else {
		throw SOFTLibException("Setting value is not a boolean.");
	}
}

/**
 * Mark this Setting object as "touched",
 * i.e. that it's value has been used.
 */
void Setting::Touch() {
    this->touched = true;
}

/**
 * Check if this setting has been touched.
 * If not, this indicates that the setting
 * is unrecognized and invalid.
 */
bool Setting::Touched() const {
    return this->touched;
}

