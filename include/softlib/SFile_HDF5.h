#ifndef _SFILE_HDF5_H
#define _SFILE_HDF5_H

#include "H5Cpp.h"
#include <softlib/SFile.h>

class SFile_HDF5 : public SFile {
    private:
        template<typename T>
		void WriteNumArray(const std::string&, T**, sfilesize_t, sfilesize_t, H5::PredType, H5::PredType);
	protected:
		H5::H5File *file;
	public:
		SFile_HDF5();
		~SFile_HDF5();

		void Close();
        virtual bool HasVariable(const std::string&) override;
		virtual void CreateStruct(const std::string&) override;
		virtual double GetAttributeScalar(const std::string&, const std::string&) override;
		virtual std::string GetAttributeString(const std::string&, const std::string&) override;
		virtual double *GetMultiArray_linear(const std::string&, const sfilesize_t, sfilesize_t&, sfilesize_t*) override;
		virtual std::string GetString(const std::string&) override;
		virtual void Open(const std::string&, enum sfile_mode) override;
		virtual void WriteAttribute_scalar(const std::string&, const std::string&, double) override;
		virtual void WriteAttribute_string(const std::string&, const std::string&, const std::string&) override;
		virtual void WriteImage(const std::string&, double**, sfilesize_t) override;
		virtual void WriteMultiArray(const std::string&, double*, sfilesize_t, sfilesize_t*) override;
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

		virtual void WriteArray(const std::string&, double**, sfilesize_t, sfilesize_t) override;
        virtual void WriteInt32Array(const std::string&, int32_t**, sfilesize_t, sfilesize_t) override;
        virtual void WriteInt64Array(const std::string&, int64_t**, sfilesize_t, sfilesize_t) override;
        virtual void WriteUInt32Array(const std::string&, uint32_t**, sfilesize_t, sfilesize_t) override;
        virtual void WriteUInt64Array(const std::string&, uint64_t**, sfilesize_t, sfilesize_t) override;

		virtual void WriteList(const std::string&, double*, sfilesize_t) override;
        virtual void WriteInt32List(const std::string&, int32_t*, sfilesize_t) override;
        virtual void WriteInt64List(const std::string&, int64_t*, sfilesize_t) override;
        virtual void WriteUInt32List(const std::string&, uint32_t*, sfilesize_t) override;
        virtual void WriteUInt64List(const std::string&, uint64_t*, sfilesize_t) override;

		H5::Group __GetParentGroup(const std::string&, std::string&);
        H5::H5File *__GetFile() { return this->file; }
};

#endif/*_SFILE_HDF5_H*/
