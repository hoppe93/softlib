#ifndef _CONFIGURATION_SFILE_H
#define _CONFIGURATION_SFILE_H

#include <vector>
#include <string>

#include <softlib/Configuration.h>

class SettingSFile : public Setting {
public:
    enum DataType {
        UNSPECIFIED,
        BOOL,
        INT32,
        INT64,
        UINT32,
        UINT64,
        REAL,
        STRING
    };
    enum DataType dtype = DataType::UNSPECIFIED;

    union {
        std::vector<bool> *b;
        std::vector<int32_t> *i32;
        std::vector<int64_t> *i64;
        std::vector<uint32_t> *u32;
        std::vector<uint64_t> *u64;
        std::vector<slibreal_t> *r;
        std::vector<std::string> *s;
    } data;

    SettingSFile(std::vector<bool> *b) { SetValueBool(b); }
    SettingSFile(std::vector<int32_t> *i) { SetValueInt32(i); }
    SettingSFile(std::vector<int64_t> *i) { SetValueInt64(i); }
    SettingSFile(std::vector<uint32_t> *u) { SetValueUInt32(u); }
    SettingSFile(std::vector<uint64_t> *u) { SetValueUInt64(u); }
    SettingSFile(std::vector<slibreal_t> *s) { SetValueReal(s); }
    SettingSFile(std::vector<std::string> *s) { SetValueString(s); }
	~SettingSFile();

    void DeleteCurrent();
	virtual void OverwriteValues(const Setting*) override;

    virtual Setting *Copy() override;
	virtual size_t GetNumberOfValues() const override;

    // Get value(s)
	virtual bool GetBool(unsigned int index=0) const override;
    virtual int64_t GetInteger(unsigned int index=0) const override;
    virtual int32_t GetInteger32(unsigned int index=0) const override;
    virtual uint32_t GetUnsignedInteger32(unsigned int index=0) const override;
    virtual int64_t GetInteger64(unsigned int index=0) const override;
    virtual uint64_t GetUnsignedInteger64(unsigned int index=0) const override;
	virtual slibreal_t GetScalar(unsigned int index=0) const override;
	virtual std::string GetString(unsigned int index=0) const override;
	virtual std::vector<slibreal_t> GetNumericVector() const override;
	virtual const std::vector<std::string> GetTextVector() const override;

    // Check value type
	virtual bool IsBool() const override;
	virtual bool IsBool(unsigned int) const override;
    virtual bool IsInteger() const override;
    virtual bool IsInteger(unsigned int) const override;
    virtual bool IsInteger32() const override;
    virtual bool IsInteger32(unsigned int) const override;
    virtual bool IsInteger64() const override;
    virtual bool IsInteger64(unsigned int) const override;
    virtual bool IsUnsignedInteger32() const override;
    virtual bool IsUnsignedInteger32(unsigned int) const override;
    virtual bool IsUnsignedInteger64() const override;
    virtual bool IsUnsignedInteger64(unsigned int) const override;
	virtual bool IsScalar() const override;
	virtual bool IsScalar(unsigned int) const override;
	virtual bool IsNumericVector() const override;
    virtual bool IsNumericVector(unsigned int) const override;

    // Setters
    void SetValueBool(std::vector<bool>*);
    void SetValueInt32(std::vector<int32_t>*);
    void SetValueInt64(std::vector<int64_t>*);
    void SetValueUInt32(std::vector<uint32_t>*);
    void SetValueUInt64(std::vector<uint64_t>*);
    void SetValueReal(std::vector<slibreal_t>*);
    void SetValueString(std::vector<std::string>*);
};
class ConfigurationSFile : public Configuration {
    ConfigurationSFile() : Configuration() {}
    ConfigurationSFile(const Configuration& c) : Configuration(c) {}
    virtual ~ConfigurationSFile() {}

    virtual void FromFile(const std::string&) override;
};

#endif/*_CONFIGURATION_SFILE_H*/
