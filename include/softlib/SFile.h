#ifndef _SFILE_H
#define _SFILE_H

#include "config.h"
#include <string>

typedef unsigned long long int sfilesize_t;

enum sfile_mode {
	SFILE_MODE_READ,
	SFILE_MODE_UPDATE,
	SFILE_MODE_WRITE
};

enum sfile_type {
	SFILE_TYPE_UNKNOWN,
	SFILE_TYPE_HDF5,
	SFILE_TYPE_MATLAB,
	SFILE_TYPE_SDT
};

class SFile {
	public:
		std::string filename;
		enum sfile_mode mode;
		enum sfile_type filetype;

		static SFile *Create(const std::string&, enum sfile_mode);
		static SFile *Create(const std::string&, enum sfile_mode, enum sfile_type);
		static enum sfile_type GetFileType(const std::string&);
		static enum sfile_type TypeOfFile(const std::string&);

		virtual void Close() = 0;
		virtual double GetAttributeScalar(const std::string&, const std::string&) = 0;
		virtual std::string *GetAttributeString(const std::string&, const std::string&) = 0;
		virtual double **GetDoubles(const std::string&, sfilesize_t*) = 0;
        virtual double *GetDoubles1D(const std::string&, sfilesize_t*) = 0;
		virtual std::string *GetString(const std::string&) = 0;
		virtual void Open(const std::string&, enum sfile_mode) = 0;
		virtual void WriteArray(const std::string&, double**, sfilesize_t, sfilesize_t) = 0;
		virtual void WriteAttribute_scalar(const std::string&, const std::string&, double) = 0;
		virtual void WriteAttribute_string(const std::string&, const std::string&, const std::string&) = 0;
		virtual void WriteImage(const std::string&, double**, sfilesize_t) = 0;
		virtual void WriteList(const std::string&, double*, sfilesize_t) = 0;
		virtual void WriteString(const std::string&, const std::string&) = 0;

		// Functions implemented in the base class
		double *GetList(const std::string&, sfilesize_t*);
		double GetScalar(const std::string&);
        void WriteScalar(const std::string&, double);
};

#endif/*_SFILE_H*/
