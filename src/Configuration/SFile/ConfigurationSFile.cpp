/**
 * Implementation of the Configuration SFile interface.
 */

#include <vector>
#include <string>

#include <softlib/Configuration.h>
#include <softlib/ConfigurationException.h>
#include <softlib/Configuration/ConfigurationSFile.h>
#include <softlib/SFile.h>

/**
 * Load a configuration from the given input file
 * (using the SFile interface).
 *
 * filename: Name of input file to load.
 */
void ConfigurationSFile::FromFile(const string& filename) {
    if (!this->filename.empty())
        throw ConfigurationException("A configuration file has already been opened with this ConfigurationSFile object.");

    SFile *sf = SFile::Create(filename, SFILE_MODE_READ);
    this->filename = filename;

    FromSFile(sf);
}

void ConfigurationSFile::FromSFile(const SFile *sf) {
}

