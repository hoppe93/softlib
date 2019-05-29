#ifndef _CONFIGURATION_SFILE_H
#define _CONFIGURATION_SFILE_H

#include <vector>
#include <string>

#include <softlib/Configuration.h>

class ConfigurationSFile : public Configuration {
    ConfigurationSFile() : Configuration() {}
    ConfigurationSFile(const Configuration& c) : Configuration(c) {}
    virtual ~ConfigurationSFile() {}

    virtual void FromFile(const std::string&) override;
};

#endif/*_CONFIGURATION_SFILE_H*/
