#ifndef _SFILE_H
#define _SFILE_H

#include "config.h"
#include <functional>
#include <map>
#include <string>
#include <vector>

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
    private:
        template<typename T>
        T __GetSingle(const std::string&, T**(SFile::*)(const std::string&, sfilesize_t*), SFile*);
	public:
		std::string filename;
		enum sfile_mode mode;
		enum sfile_type filetype;

		static SFile *Create(const std::string&, enum sfile_mode);
		static SFile *Create(const std::string&, enum sfile_mode, enum sfile_type);
		static enum sfile_type GetFileType(const std::string&);
		static enum sfile_type TypeOfFile(const std::string&);

		virtual void *GetMultiArray(const std::string&, const sfilesize_t, sfilesize_t&, sfilesize_t*);

        virtual bool HasVariable(const std::string&) = 0;

		virtual void Close() = 0;
		virtual void CreateStruct(const std::string&) = 0;
		virtual double GetAttributeScalar(const std::string&, const std::string&) = 0;
		virtual std::string GetAttributeString(const std::string&, const std::string&) = 0;
		virtual double *GetMultiArray_linear(const std::string&, const sfilesize_t, sfilesize_t&, sfilesize_t*) = 0;
		virtual std::string GetString(const std::string&) = 0;
		virtual void Open(const std::string&, enum sfile_mode) = 0;
		virtual void WriteAttribute_scalar(const std::string&, const std::string&, double) = 0;
		virtual void WriteAttribute_string(const std::string&, const std::string&, const std::string&) = 0;
		virtual void WriteImage(const std::string&, double**, sfilesize_t) = 0;
		virtual void WriteMultiArray(const std::string&, double*, sfilesize_t, sfilesize_t*) = 0;
		virtual void WriteString(const std::string&, const std::string&) = 0;

		virtual double **GetDoubles(const std::string&, sfilesize_t*) = 0;
        virtual int32_t **GetInt32_2D(const std::string&, sfilesize_t*) = 0;
        virtual int64_t **GetInt64_2D(const std::string&, sfilesize_t*) = 0;
        virtual uint32_t **GetUInt32_2D(const std::string&, sfilesize_t*) = 0;
        virtual uint64_t **GetUInt64_2D(const std::string&, sfilesize_t*) = 0;

        virtual double *GetDoubles1D(const std::string&, sfilesize_t*) = 0;
        virtual int32_t *GetInt32_1D(const std::string&, sfilesize_t*) = 0;
        virtual int64_t *GetInt64_1D(const std::string&, sfilesize_t*) = 0;
        virtual uint32_t *GetUInt32_1D(const std::string&, sfilesize_t*) = 0;
        virtual uint64_t *GetUInt64_1D(const std::string&, sfilesize_t*) = 0;

		virtual void WriteArray(const std::string&, double**, sfilesize_t, sfilesize_t) = 0;
        virtual void WriteInt32Array(const std::string&, int32_t**, sfilesize_t, sfilesize_t) = 0;
        virtual void WriteInt64Array(const std::string&, int64_t**, sfilesize_t, sfilesize_t) = 0;
        virtual void WriteUInt32Array(const std::string&, uint32_t**, sfilesize_t, sfilesize_t) = 0;
        virtual void WriteUInt64Array(const std::string&, uint64_t**, sfilesize_t, sfilesize_t) = 0;

		virtual void WriteList(const std::string&, double*, sfilesize_t) = 0;
        virtual void WriteInt32List(const std::string&, int32_t*, sfilesize_t) = 0;
        virtual void WriteInt64List(const std::string&, int64_t*, sfilesize_t) = 0;
        virtual void WriteUInt32List(const std::string&, uint32_t*, sfilesize_t) = 0;
        virtual void WriteUInt64List(const std::string&, uint64_t*, sfilesize_t) = 0;

		// Functions implemented in the base class
		double *GetList(const std::string&, sfilesize_t*);
        int32_t GetInt32(const std::string&);
        int64_t GetInt64(const std::string&);
        uint32_t GetUInt32(const std::string&);
        uint64_t GetUInt64(const std::string&);
		double GetScalar(const std::string&);
        void WriteScalar(const std::string&, double);
};

#endif/*_SFILE_H*/
