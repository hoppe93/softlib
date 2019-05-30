#ifndef _CONFIGURATION_H
#define _CONFIGURATION_H

#include <initializer_list>
#include <string>
#include <vector>

#include <softlib/config.h>
#include <softlib/ConfigurationException.h>

typedef size_t confblock_t;

/************************************************
 * PUBLIC OBJECTS                               *
 ************************************************/
/* Class representing an individual setting */
class Setting {
protected:
	std::string name;
    bool touched=false;

public:
    Setting() {}
    Setting(const Setting*);
	/*Setting(const std::string&, const std::string&);
    Setting(const std::string&, const std::vector<std::string>&);*/
	virtual ~Setting();

	bool HasName(const std::string&) const;
	const std::string& GetName() const;

    virtual Setting *Copy() = 0;
    virtual void OverwriteValues(const Setting*) = 0;
    virtual size_t GetNumberOfValues() = 0;

    // Get value(s)
	virtual bool GetBool(unsigned int index=0) = 0;
    virtual int64_t GetInteger(unsigned int index=0) = 0;
    virtual int32_t GetInteger32(unsigned int index=0) = 0;
    virtual uint32_t GetUnsignedInteger32(unsigned int index=0) = 0;
    virtual int64_t GetInteger64(unsigned int index=0) = 0;
    virtual uint64_t GetUnsignedInteger64(unsigned int index=0) = 0;
	virtual slibreal_t GetScalar(unsigned int index=0) = 0;
	virtual std::string& GetString(unsigned int index=0) = 0;
	virtual std::vector<slibreal_t> GetNumericVector() = 0;
	virtual const std::vector<std::string> GetTextVector() const = 0;

    // Check value type
	virtual bool IsBool() const = 0;
	virtual bool IsBool(unsigned int) const = 0;
    virtual bool IsInteger() const = 0;
    virtual bool IsInteger(unsigned int) const = 0;
    virtual bool IsInteger32() const = 0;
    virtual bool IsInteger32(unsigned int) const = 0;
    virtual bool IsInteger64() const = 0;
    virtual bool IsInteger64(unsigned int) const = 0;
    virtual bool IsUnsignedInteger32() const = 0;
    virtual bool IsUnsignedInteger32(unsigned int) const = 0;
    virtual bool IsUnsignedInteger64() const = 0;
    virtual bool IsUnsignedInteger64(unsigned int) const = 0;
	virtual bool IsScalar() const = 0;
	virtual bool IsScalar(unsigned int) const = 0;
	virtual bool IsNumericVector() const = 0;
    virtual bool IsNumericVector(unsigned int) const = 0;

    void Touch();
    bool Touched() const;

    static bool StringToBool(const std::string&);
};

/* Class representing a settings block */
class ConfigBlock {
	confblock_t type;
	std::string name;
    std::string secondary_type;
	ConfigBlock *parent;
	std::vector<ConfigBlock*> subblocks;
	std::vector<Setting*> settings;

public:
	ConfigBlock();
    ConfigBlock(const ConfigBlock&);
	ConfigBlock(confblock_t, const std::string&, const std::string& stype="", ConfigBlock *parent=nullptr);
	~ConfigBlock();

	/*Setting& AddSetting(const std::string&, const std::string&);
    Setting& AddSetting(const std::string&, const std::vector<std::string>&);*/
	Setting *AddSetting(Setting*);
	ConfigBlock *AddSubBlock(confblock_t, const std::string&, const std::string& stype="");
	ConfigBlock *AddSubBlock(ConfigBlock&);
	std::vector<Setting*> GetAllSettings() const;
	std::vector<ConfigBlock> GetAllSubBlocks() const;
	std::vector<ConfigBlock> GetAllSubBlocksOfType(confblock_t) const;
	ConfigBlock *GetConfigBlock(confblock_t, const std::string&);
	Setting *GetSetting(const std::string&);
	Setting *GetSetting(const char*);
	const std::string& GetName() const;
	confblock_t GetType() const;
    const std::string& GetSecondaryType() const;
    std::vector<Setting*> GetUntouchedSettings();
	bool HasSetting(const std::string&);		// Check if setting with specified name exists
	bool HasSubBlock(confblock_t, const std::string&);
	bool HasSubBlocks();			// Check if Configuration has any subblocks
	bool HasSubBlocks(confblock_t);// Check if sub-block with specified blocktype exists
	bool HasType(confblock_t) const;	// Check if this block has the specified type
	bool HasTypeAndName(confblock_t, const std::string&) const;// Check if this block has the given type and name
    bool HasUntouchedSettings() const;
	ConfigBlock *GetParent() const;
	bool HasParent() const;
	std::vector<std::string> *Merge(ConfigBlock&, bool allowNew=false);
	void MergeSetting(Setting*);
	void SetParent(ConfigBlock&);

	bool operator==(const ConfigBlock& rhs) const;
	std::string& operator[](const std::string&);
};

/* Root class, representing the configuration file */
class Configuration {
protected:
	std::string filename;
	ConfigBlock root;
	std::vector<std::string> blocktypes;
    bool errorflag=false;

public: 
	Configuration();
    Configuration(const Configuration&);
	virtual ~Configuration();

	virtual void FromFile(const std::string&) = 0;

    bool HasError() const { return this->errorflag; }
	ConfigBlock GetRootBlock() const { return this->root; }
	bool IsBlockType(const std::string&);
    std::vector<std::string> *Merge(Configuration&, bool allowNew=false);
	confblock_t RegisterBlockType(const std::string&);
    void SetErrorFlag();
	confblock_t StringToBlockType(const std::string&);

    std::vector<std::string> GetBlockTypes() const { return this->blocktypes; }
    std::string GetFilename() const { return this->filename; }
};

#endif/*_CONFIGURATION_H*/
