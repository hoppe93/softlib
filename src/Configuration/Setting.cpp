
#include <string>
#include <vector>
#include <iostream>

#include <softlib/Configuration.h>
#include <softlib/SOFTLibException.h>

using namespace std;

/**
 * Constructors & Destructors
 */
Setting::Setting(const Setting &s) : Setting(s.GetName(), s.GetTextVector()) { }
Setting::Setting(const string& name, const string& value) {
	this->name = name;
	this->values.push_back(value);
}
Setting::Setting(const string& name, const vector<string> &values) {
    this->name = name;
    this->values = values;
}
Setting::~Setting() { }

/**
 * Append a value to the list of values
 * of this Setting.
 *
 * val: String containing the raw value
 */
Setting& Setting::AppendValue(const string& val) {
	this->values.push_back(val);
	return *this;
}

/**
 * Get the value as a boolean
 */
bool Setting::GetBool(unsigned int index) {
	string& value = this->values.at(index);

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
 * Get the value as a 64-bit signed integer.
 */
int64_t Setting::GetInteger(unsigned int index) { return GetInteger64(index); }

/**
 * Get the value as a 32-bit signed/unsigned integer.
 */
int32_t Setting::GetInteger32(unsigned int index) {
    return stoi(this->values.at(index));
}
uint32_t Setting::GetUnsignedInteger32(unsigned int index) {
    return stoul(this->values.at(index));
}

/**
 * Get the value as a 64-bit signed integer.
 */
int64_t Setting::GetInteger64(unsigned int index) {
    return stoll(this->values.at(index));
}
uint64_t Setting::GetUnsignedInteger64(unsigned int index) {
    return stoull(this->values.at(index));
}

/**
 * Get the name of this setting.
 */
const string& Setting::GetName() const { return this->name; }

/**
 * Get the value as a scalar
 */
slibreal_t Setting::GetScalar(unsigned int index) {
    if (is_same<slibreal_t, float>::value)
        return stof(this->values.at(index), NULL);
    else if (is_same<slibreal_t, double>::value)
        return stod(this->values.at(index), NULL);
    else if (is_same<slibreal_t, long double>::value)
        return stold(this->values.at(index), NULL);
    else
        throw SOFTLibException("Invalid configuration of 'slibreal_t'.");
}

/**
 * Get the value as a string
 */
string& Setting::GetString(unsigned int index) { return this->values.at(index); }

/**
 * Get the value as a numeric vector
 */
vector<double> Setting::GetNumericVector() {
	size_t l = this->values.size();
	vector<double> ret = vector<double>(l, 0.0);

	for (size_t i = 0; i < l; i++) {
		ret[i] = stod(this->values[i], NULL);
	}

	return ret;
}

/**
 * Get the value as a string vector
 */
const vector<string> Setting::GetTextVector() const {
	return this->values;
}

/**
 * Get the number of elements in 'values'
 */
size_t Setting::GetNumberOfValues() {
	return this->values.size();
}

/**
 * Check if this setting has the given name
 *
 * name: Name to compare to
 */
bool Setting::HasName(const string& name) const { return (this->name == name); }

/**
 * Overwrite the values of this setting.
 * 
 * values: New values to assign to this setting.
 */
void Setting::OverwriteValues(const vector<string>& values) {
	this->values = values;
}

/**
 * Check if this setting is a valid boolean.
 * This function will return false unless
 * 'values' contains exactly one element.
 */
bool Setting::IsBool() const {
	if (this->values.size() != 1) return false;
    return IsBool(0);
}
bool Setting::IsBool(unsigned int i) const {
	const string& value = this->values.at(i);

	if (value == "yes"   || value == "Yes"   || value == "YES" ||
		value == "no"    || value == "No"    || value == "NO" ||
		value == "true"  || value == "True"  || value == "TRUE" ||
		value == "false" || value == "False" || value == "FALSE")
		return true;
	else return false;
}

/**
 * Check if this setting is a valid integer.
 */
bool Setting::IsInteger() const { return IsInteger64(); }
bool Setting::IsInteger(unsigned int i) const { return IsInteger64(i); }
bool Setting::IsInteger32() const {
    if (this->values.size() != 1) return false;
    return IsInteger32(0);
}
bool Setting::IsInteger32(unsigned int i) const {
    const string& value = this->values.at(i);
    size_t end = 0;

    try {
        stoi(value, &end);
    } catch (invalid_argument& e) {
        return false;
    } catch (out_of_range& e) {
        return false;
    }

    return (end == value.length());
}
bool Setting::IsUnsignedInteger32() const {
    if (this->values.size() != 1) return false;
    return IsUnsignedInteger32(0);
}
bool Setting::IsUnsignedInteger32(unsigned int i) const {
    const string& value = this->values.at(i);
    size_t end = 0;

    try {
        stoul(value, &end);
    } catch (invalid_argument& e) {
        return false;
    } catch (out_of_range& e) {
        return false;
    }

    return (end == value.length());
}
bool Setting::IsInteger64() const {
    if (this->values.size() != 1) return false;
    return IsInteger64(0);
}
bool Setting::IsInteger64(unsigned int i) const {
    const string& value = this->values.at(i);
    size_t end = 0;

    try {
        stoll(value, &end);
    } catch (invalid_argument& e) {
        return false;
    } catch (out_of_range& e) {
        return false;
    }

    return (end == value.length());
}
bool Setting::IsUnsignedInteger64() const {
    if (this->values.size() != 1) return false;
    return IsUnsignedInteger64(0);
}
bool Setting::IsUnsignedInteger64(unsigned int i) const {
    const string& value = this->values.at(i);
    size_t end = 0;

    try {
        stoull(value, &end);
    } catch (invalid_argument& e) {
        return false;
    } catch (out_of_range& e) {
        return false;
    }

    return (end == value.length());
}

/**
 * Check if this setting is a valid scalar.
 * This function returns true if the number of
 * elements of the value field of this setting
 * is exactly 1 and if that value is a valid
 * numeric string.
 */
bool Setting::IsScalar() const {
	if (this->values.size() != 1) return false;
    return IsScalar(0);
}
bool Setting::IsScalar(unsigned int i) const {
	const string& value = this->values.at(i);

	try {
		stod(value, NULL);
	} catch (invalid_argument& e) {
		return false;
	}

	return true;
}

/**
 * Check if this setting is a valid numeric vector.
 * 
 * (optional) n: Number of elements expected.
 */
bool Setting::IsNumericVector() const {
	size_t l = this->values.size();
	for (size_t i = 0; i < l; i++) {
		try {
			stod(this->values[i], NULL);
		} catch (invalid_argument& e) {
			return false;
		}
	}

	return true;
}
bool Setting::IsNumericVector(unsigned int n) const {
    if (this->values.size() != n)
        return false;
    else
        return IsNumericVector();
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

