
#include <string>
#include <vector>
#include <iostream>

#include <softlib/Configuration.h>
#include <softlib/Configuration/ConfigurationScript.h>
#include <softlib/SOFTLibException.h>

using namespace std;

/**
 * Constructors & Destructors
 */
SettingScript::SettingScript(const Setting *s) : SettingScript(s->GetName(), s->GetTextVector()) { }
SettingScript::SettingScript(const string& name, const string& value) {
	this->name = name;
	this->values.push_back(value);
}
SettingScript::SettingScript(const string& name, const vector<string> &values) {
    this->name = name;
    this->values = values;
}
SettingScript::~SettingScript() { }

/**
 * Append a value to the list of values
 * of this SettingScript.
 *
 * val: String containing the raw value
 */
SettingScript *SettingScript::AppendValue(const string& val) {
	this->values.push_back(val);
	return this;
}

/**
 * Create a copy of this setting.
 */
Setting *SettingScript::Copy() {
    return new SettingScript(this);
}

/**
 * Get the value as a boolean
 */
bool SettingScript::GetBool(unsigned int index) {
	string& value = this->values.at(index);
    return StringToBool(value);
}

/**
 * Get the value as a 64-bit signed integer.
 */
int64_t SettingScript::GetInteger(unsigned int index) { return GetInteger64(index); }

/**
 * Get the value as a 32-bit signed/unsigned integer.
 */
int32_t SettingScript::GetInteger32(unsigned int index) {
    return stoi(this->values.at(index));
}
uint32_t SettingScript::GetUnsignedInteger32(unsigned int index) {
    return stoul(this->values.at(index));
}

/**
 * Get the value as a 64-bit signed integer.
 */
int64_t SettingScript::GetInteger64(unsigned int index) {
    return stoll(this->values.at(index));
}
uint64_t SettingScript::GetUnsignedInteger64(unsigned int index) {
    return stoull(this->values.at(index));
}

/**
 * Get the value as a scalar
 */
slibreal_t SettingScript::GetScalar(unsigned int index) {
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
string& SettingScript::GetString(unsigned int index) { return this->values.at(index); }

/**
 * Get the value as a numeric vector
 */
vector<double> SettingScript::GetNumericVector() {
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
const vector<string> SettingScript::GetTextVector() const {
	return this->values;
}

/**
 * Get the number of elements in 'values'
 */
size_t SettingScript::GetNumberOfValues() {
	return this->values.size();
}

/**
 * Overwrite the values of this setting.
 * 
 * s: Setting, the value of which should override
 *    this setting's value.
 */
void SettingScript::OverwriteValues(const Setting *s) {
	this->values = s->GetTextVector();
}

/**
 * Check if this setting is a valid boolean.
 * This function will return false unless
 * 'values' contains exactly one element.
 */
bool SettingScript::IsBool() const {
	if (this->values.size() != 1) return false;
    return IsBool(0);
}
bool SettingScript::IsBool(unsigned int i) const {
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
bool SettingScript::IsInteger() const { return IsInteger64(); }
bool SettingScript::IsInteger(unsigned int i) const { return IsInteger64(i); }
bool SettingScript::IsInteger32() const {
    if (this->values.size() != 1) return false;
    return IsInteger32(0);
}
bool SettingScript::IsInteger32(unsigned int i) const {
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
bool SettingScript::IsUnsignedInteger32() const {
    if (this->values.size() != 1) return false;
    return IsUnsignedInteger32(0);
}
bool SettingScript::IsUnsignedInteger32(unsigned int i) const {
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
bool SettingScript::IsInteger64() const {
    if (this->values.size() != 1) return false;
    return IsInteger64(0);
}
bool SettingScript::IsInteger64(unsigned int i) const {
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
bool SettingScript::IsUnsignedInteger64() const {
    if (this->values.size() != 1) return false;
    return IsUnsignedInteger64(0);
}
bool SettingScript::IsUnsignedInteger64(unsigned int i) const {
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
bool SettingScript::IsScalar() const {
	if (this->values.size() != 1) return false;
    return IsScalar(0);
}
bool SettingScript::IsScalar(unsigned int i) const {
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
bool SettingScript::IsNumericVector() const {
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
bool SettingScript::IsNumericVector(unsigned int n) const {
    if (this->values.size() != n)
        return false;
    else
        return IsNumericVector();
}

