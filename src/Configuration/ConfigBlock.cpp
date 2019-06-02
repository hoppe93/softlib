
#include <string>
#include <vector>

#include <softlib/Configuration.h>
#include <softlib/SOFTLibException.h>

using namespace std;

/**
 * Constructor & Destructor
 */
ConfigBlock::ConfigBlock() { }
ConfigBlock::ConfigBlock(const ConfigBlock *cb) {
    this->type = cb->GetType();
    this->name = cb->GetName();
    this->secondary_type = cb->GetSecondaryType();
    this->parent = cb->GetParent();

    // Copy sub-blocks
    vector<ConfigBlock*> sub = cb->GetAllSubBlocks();
    for (unsigned int i = 0; i < sub.size(); i++)
        subblocks.push_back(new ConfigBlock(sub.at(i)));

    // Copy settings
    this->settings = cb->GetAllSettings();
    for (Setting *s : this->settings)
        settings.push_back(s->Copy());
}
ConfigBlock::ConfigBlock(confblock_t type, const string& name, const string& secondary_type, ConfigBlock *parent) {
	this->type = type;
	this->name = name;
    this->secondary_type = secondary_type;
	this->parent = parent;
}
ConfigBlock::~ConfigBlock() {
    // Delete settings
    for (Setting *s : this->settings) {
        delete s;
    }

    for (ConfigBlock *cb : this->subblocks) {
        delete cb;
    }
}

/**
 * Add setting to the list of settings
 * in this block from an existing setting.
 *
 * name: Name of the setting.
 * values: List of values.
 */
Setting *ConfigBlock::AddSetting(Setting *set) {
    if (HasSetting(set->GetName())) {
        // Replace setting
        Setting *s = GetSetting(set->GetName());
        s->OverwriteValues(set);

        delete set;

        return s;
    } else {
        this->settings.push_back(set);
        return set;
    }
}

/**
 * Add a sub-block to this ConfigBlock and
 * return the newly added object.
 */
ConfigBlock *ConfigBlock::AddSubBlock(confblock_t type, const string& name, const string& secondary_block) {
    return AddSubBlock(new ConfigBlock(type, name, secondary_block, this));
}

/**
 * Add an existing sub-block to this ConfigBlock and
 * return the newly added object.
 */
ConfigBlock *ConfigBlock::AddSubBlock(ConfigBlock *cb) {
    if (HasSubBlock(cb->GetType(), cb->GetName())) {
        return GetConfigBlock(cb->GetType(), cb->GetName());
    } else {
        this->subblocks.push_back(cb);
        return cb;
    }
}

/**
 * Get a list of all settings of this ConfigBlock.
 */
vector<Setting*> ConfigBlock::GetAllSettings() const {
	return this->settings;
}

/**
 * Get a list of all sub-blocks of this ConfigBlock.
 */
vector<ConfigBlock*> ConfigBlock::GetAllSubBlocks() const {
	return this->subblocks;
}

/**
 * Get a list of all sub-blocks of a specific type.
 *
 * type: Type of blocks to retrieve.
 */
vector<ConfigBlock*> ConfigBlock::GetAllSubBlocksOfType(confblock_t type) const {
    vector<ConfigBlock*> blocks;

    for (ConfigBlock *cb : this->subblocks) {
        if (cb->HasType(type))
            blocks.push_back(cb);
    }

    return blocks;
}

/**
 * Get the sub-configuration block of type 'type'
 * with name 'name'. Returns nullptr if not found.
 *
 * type: Type ID of ConfigBlock to get.
 * name: Name of ConfigBlock to get.
 */
ConfigBlock *ConfigBlock::GetConfigBlock(confblock_t type, const string& name) {
    for (ConfigBlock *cb : this->subblocks) {
        if (cb->HasTypeAndName(type, name))
            return cb;
    }

	return nullptr;
}

/**
 * Get the name of this ConfigBlock.
 */
const string& ConfigBlock::GetName() const { return name; }

/**
 * Get the secondary type of the ConfigBlock.
 */
const string& ConfigBlock::GetSecondaryType() const {
    if (secondary_type.empty())
        return name;
    else
        return secondary_type;
}

/**
 * Return the parent of this ConfigBlock.
 */
ConfigBlock *ConfigBlock::GetParent() const {
	return this->parent;
}

/**
 * Get the setting with name 'name'.
 * Returns nullptr if not found.
 *
 * name: Of setting to retrieve.
 */
Setting *ConfigBlock::GetSetting(const string& name) {
    for (Setting *s : this->settings) {
        if (s->HasName(name)) {
            s->Touch();
            return s;
        }
    }

	return nullptr;
}
Setting *ConfigBlock::GetSetting(const char *name) {
	string s = name;
	return GetSetting(s);
}

/**
 * Get the type of this ConfigBlock.
 */
confblock_t ConfigBlock::GetType() const { return type; }

/**
 * Returns a list of pointers to all settings
 * that were untouched.
 */
vector<Setting*> ConfigBlock::GetUntouchedSettings() {
    vector<Setting*> untouched;

    for (Setting *s : this->settings) {
        if (!s->Touched())
            untouched.push_back(s);
    }

    return untouched;
}

/**
 * Checks whether this ConfigBlock contains
 * the setting with name 'name'.
 *
 * name: Name of setting to check for.
 */
bool ConfigBlock::HasSetting(const string& name) {
    for (Setting *s : this->settings) {
        if (s->HasName(name))
            return true;
    }

    return false;
}

/**
 * Check if this ConfigBlock has a sub-block of
 * the given type and name.
 *
 * type: Type of ConfigBlock to check for.
 * name: Name of ConfigBlock go check for.
 */
bool ConfigBlock::HasSubBlock(confblock_t type, const string& name) {
    for (ConfigBlock *cb : this->subblocks) {
        if (cb->HasTypeAndName(type, name))
            return true;
    }

	return false;
}

/**
 * Check if this ConfigBlock has any sub-blocks.
 *
 * Optional parameters:
 *   type: Checks if ConfigBlock contains sub-blocks of type 'type'.
 */
bool ConfigBlock::HasSubBlocks() { return (subblocks.size() > 0); }
bool ConfigBlock::HasSubBlocks(confblock_t type) {
    for (ConfigBlock *cb : this->subblocks) {
        if (cb->HasType(type))
            return true;
    }

	return false;
}

/**
 * Check if this ConfigBlock has the given type.
 *
 * type: Type to compare to.
 */
bool ConfigBlock::HasType(confblock_t type) const { return (this->type == type); }

/**
 * Check if this block has the given type and name
 *
 * type: Type to compare to
 * name: Name to compare to
 */
bool ConfigBlock::HasTypeAndName(confblock_t type, const string& name) const {
	return (this->type == type && this->name == name);
}

/**
 * Check if this block has settings which remain
 * untouched, indicating that invalid settings have
 * been provided.
 *
 * RETURNS true if there were untouched settings in this
 * ConfigBlock. Returns false if all settings were touched.
 */
bool ConfigBlock::HasUntouchedSettings() const {
    for (Setting *s : this->settings) {
        if (!s->Touched())
            return true;
    }

    return false;
}

/**
 * Check if this ConfigBlock has parent.
 */
bool ConfigBlock::HasParent() const {
	return (this->parent != nullptr);
}

/**
 * Merge this ConfigBlock with the given ConfigBlock.
 * Any settings not already existing in this object
 * will be created (only if allowNew = true), while any
 * existing settings will be overwritten in favour of
 * the values assigned in the given object (always).
 *
 * cb:       ConfigBlock to copy settings from.
 * allowNew: Allow new settings to be defined by cb.
 *    If true, any setting defined by cb that does not
 *    already exist in this object will be created.
 *    If false, this method will return a list of
 *    strings containing the names of all settings
 *    that were in cb but not in this ConfigBlock.
 *
 * RETURNS 'nullptr' if allowNew = true, or if allowNew = false
 *    and all settings defined in cb were also defined
 *    in this ConfigBlock object. If allowNew = false
 *    and some settings were defined only in cb,
 *    returns a list of strings representing the
 *    names of all settings exclusively defined in cb.
 *
 * NOTE: The returned vector should be deleted after use.
 */
vector<string> *ConfigBlock::Merge(ConfigBlock *cb, bool allowNew) {
	ConfigBlock *tb;
	vector<Setting*> set = cb->GetAllSettings();
	unsigned int l = set.size(), i, j;
	vector<string> *unknown = nullptr;

	// Merge settings
	for (i = 0; i < l; i++) {
		if (this->HasSetting(set[i]->GetName()))
			this->MergeSetting(set[i]);
		else if (allowNew) {
			this->AddSetting(set[i]);
		} else {
			if (unknown == nullptr) unknown = new vector<string>();
			unknown->push_back(set[i]->GetName());
		}
	}

	// Merge sub-blocks
	vector<ConfigBlock*> sb = cb->GetAllSubBlocks();
	l = sb.size();
	for (i = 0; i < l; i++) {
		if (this->HasSubBlock(sb[i]->GetType(), sb[i]->GetName())) {
			tb = this->GetConfigBlock(sb[i]->GetType(), sb[i]->GetName());
			vector<string> *u = tb->Merge(sb[i], allowNew);

			if (u != nullptr) {
				if (unknown == nullptr) unknown = u;
				else {
					for (j = 0; j < u->size(); j++)
						unknown->push_back(u->at(j));

					delete u;
				}
			}
		} else if (allowNew) {
			this->AddSubBlock(new ConfigBlock(sb[i]));
		} else {
			if (unknown == nullptr) unknown = new vector<string>();
			unknown->push_back(sb[i]->GetName());
		}
	}

	return unknown;
}
void ConfigBlock::MergeSetting(Setting *s) {
    for (Setting *set : this->settings) {
        if (set->HasName(s->GetName())) {
            set->OverwriteValues(s);
            break;
        }
    }
}

/**
 * Set parent block of this ConfigBlock.
 *
 * parent: The ConfigBlock object to use as parent.
 */
void ConfigBlock::SetParent(ConfigBlock *parent) {
	this->parent = parent;
}

/**
 * Overload for comparing ConfigBlocks with '==' operator.
 * Equality is defined as both objects having the same name and type.
 *
 * rhs: The ConfigBlock object to compare to.
 */
bool ConfigBlock::operator==(const ConfigBlock& rhs) const {
	return this->HasTypeAndName(rhs.GetType(), rhs.GetName());
}

/**
 * Overload for accessing settings using the [] operator.
 * If the setting does not exist, a SOFTLibException is
 * thrown.
 *
 * setting: The name of the setting to return.
 */
string ConfigBlock::operator[](const string& setting) {
	Setting *set = this->GetSetting(setting);

	if (set == nullptr)
		throw SOFTLibException("No setting with name '%s' in block '%s'.", setting.c_str(), this->name.c_str());
	else return set->GetString();
}

