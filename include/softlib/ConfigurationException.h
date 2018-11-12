#ifndef _CONFIGURATION_EXCEPTION_H
#define _CONFIGURATION_EXCEPTION_H

#include <string>
#include <softlib/SOFTLibException.h>

class ConfigurationException : public SOFTLibException {
    private:
        std::string filename;
        size_t linenumber=0, charpos=0;
    public:
        ConfigurationException(const std::string& msg) : SOFTLibException(msg) {}
        template<typename ... Args>
        ConfigurationException(const std::string& msg, Args&& ... args) : SOFTLibException(msg, std::forward<Args>(args) ...) {}
        template<typename ... Args>
        ConfigurationException(
            size_t line, size_t charpos, const std::string& file, 
            const std::string& msg, Args&& ... args
        ) {
            this->filename = file;
            this->linenumber = line;
            this->charpos = charpos;

            ConstructErrorMessage("%s:%zu:%zu: "+msg, file.c_str(), line, charpos, std::forward<Args>(args) ...);
        }

        size_t GetLine() const { return this->linenumber; }
        size_t GetCharpos() const { return this->charpos; }
};

#endif/*_CONFIGURATION_EXCEPTION_H*/
