#ifndef _SFILE_MAT_HDF5_H
#define _SFILE_MAT_HDF5_H

#include <H5Cpp.h>
#include <softlib/SFile_HDF5.h>
#include <string>

class SFile_MAT : public SFile_HDF5 {
	private:
		std::string GetAttributeName(const std::string&, const std::string&);
	public:
        /* These are derived from 'SFile_HDF5':
		void Close();
        bool HasVariable(const std::string&);
		double **GetDoubles(const std::string&, sfilesize_t*);
		double *GetDoubles1D(const std::string&, sfilesize_t*);
		std::string *GetString(const std::string&);
        */
		virtual void Open(const std::string&, enum sfile_mode) override;
        H5::H5File *CreateMAT(const std::string&);

		virtual void CreateStruct(const std::string&) override;
		virtual double GetAttributeScalar(const std::string&, const std::string&) override;
		virtual std::string GetAttributeString(const std::string&, const std::string&) override;
        virtual std::string GetString(const std::string&) override;

		virtual void WriteAttribute_scalar(const std::string&, const std::string&, const double) override;
		virtual void WriteAttribute_string(const std::string&, const std::string&, const std::string&) override;
        virtual void WriteMultiArray(const std::string&, const double*, const sfilesize_t, const sfilesize_t*) override;
		virtual void WriteString(const std::string&, const std::string&) override;

		virtual void WriteArray(const std::string&, const double *const*, sfilesize_t, sfilesize_t) override;
        virtual void WriteInt32Array(const std::string&, const int32_t *const*, sfilesize_t, sfilesize_t) override;
        virtual void WriteInt64Array(const std::string&, const int64_t *const*, sfilesize_t, sfilesize_t) override;
        virtual void WriteUInt32Array(const std::string&, const uint32_t *const*, sfilesize_t, sfilesize_t) override;
        virtual void WriteUInt64Array(const std::string&, const uint64_t *const*, sfilesize_t, sfilesize_t) override;

        void WriteMATLAB_class(const std::string&, const std::string&, H5::H5Object *obj=nullptr);
        void WriteMATLAB_int_decode(const std::string&, H5::H5Object *obj=nullptr, const int val=2);
};

#endif/*_SFILE_MAT_HDF5_H*/
