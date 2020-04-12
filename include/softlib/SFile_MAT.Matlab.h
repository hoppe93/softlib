#ifndef _SFILE_MATLAB_H
#define _SFILE_MATLAB_H

#include <mat.h>
#include <softlib/SFile.h>
#include <string>

class SFile_MAT : public SFile {
	private:
		MATFile *mfp;
		std::string GetAttributeName(const std::string&, const std::string&);
	public:
		~SFile_MAT();

		virtual void Close() override;
        virtual bool HasVariable(const std::string&) override;
		virtual void CreateStruct(const std::string&) override;
		virtual double GetAttributeScalar(const std::string&, const std::string&) override;
		virtual std::string GetAttributeString(const std::string&, const std::string&) override;
		virtual double *GetMultiArray_linear(const std::string&, const sfilesize_t, sfilesize_t&, sfilesize_t*) override;
		virtual std::string GetString(const std::string&) override;
		virtual void Open(const std::string&, enum sfile_mode) override;
		virtual void WriteAttribute_scalar(const std::string&, const std::string&, double) override;
		virtual void WriteAttribute_string(const std::string&, const std::string&, const std::string&) override;
		virtual void WriteImage(const std::string&, const double *const*, sfilesize_t) override;
        virtual void WriteMultiArray(const std::string&, const double*, const sfilesize_t, const sfilesize_t*) override;
		virtual void WriteString(const std::string&, const std::string&) override;

		virtual double **GetDoubles(const std::string&, sfilesize_t*) override;
        virtual int32_t **GetInt32_2D(const std::string&, sfilesize_t*) override;
        virtual int64_t **GetInt64_2D(const std::string&, sfilesize_t*) override;
        virtual uint32_t **GetUInt32_2D(const std::string&, sfilesize_t*) override;
        virtual uint64_t **GetUInt64_2D(const std::string&, sfilesize_t*) override;

		virtual double *GetDoubles1D(const std::string&, sfilesize_t*) override;
        virtual int32_t *GetInt32_1D(const std::string&, sfilesize_t*) override;
        virtual int64_t *GetInt64_1D(const std::string&, sfilesize_t*) override;
        virtual uint32_t *GetUInt32_1D(const std::string&, sfilesize_t*) override;
        virtual uint64_t *GetUInt64_1D(const std::string&, sfilesize_t*) override;

		virtual void WriteArray(const std::string&, const double *const*, sfilesize_t, sfilesize_t) override;
        virtual void WriteInt32Array(const std::string&, const int32_t *const*, sfilesize_t, sfilesize_t) override;
        virtual void WriteInt64Array(const std::string&, const int64_t *const*, sfilesize_t, sfilesize_t) override;
        virtual void WriteUInt32Array(const std::string&, const uint32_t *const*, sfilesize_t, sfilesize_t) override;
        virtual void WriteUInt64Array(const std::string&, const uint64_t *const*, sfilesize_t, sfilesize_t) override;

		virtual void WriteList(const std::string&, const double*, sfilesize_t) override;
        virtual void WriteInt32List(const std::string&, const int32_t*, sfilesize_t) override;
        virtual void WriteInt64List(const std::string&, const int64_t*, sfilesize_t) override;
        virtual void WriteUInt32List(const std::string&, const uint32_t*, sfilesize_t) override;
        virtual void WriteUInt64List(const std::string&, const uint64_t*, sfilesize_t) override;

        mxArray *__GetArray(const std::string&, sfilesize_t*);
};

#endif/*_SFILE_MATLAB_H*/
