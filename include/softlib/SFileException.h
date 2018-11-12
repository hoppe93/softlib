#ifndef _SFILE_EXCEPTION_H
#define _SFILE_EXCEPTION_H

#include <exception>
#include <string>
#include <softlib/SOFTLibException.h>

class SFileException : public SOFTLibException {
	public:
        template<typename ... Args>
		SFileException(const std::string& msg, Args&& ... args) : SOFTLibException(msg, std::forward<Args>(args) ...) {}
};

#endif/*_SFILE_EXCEPTION_H*/
