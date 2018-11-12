#ifndef _SFILE_MAT_HDF5_H
#define _SFILE_MAT_HDF5_H

#include <H5Cpp.h>
#include <softlib/SFile_HDF5.h>
#include <string>

class SFile_MAT : public SFile_HDF5 {
	private:
		std::string *GetAttributeName(const std::string&, const std::string&);
	public:
        /* These are derived from 'SFile_HDF5':
		void Close();
		double **GetDoubles(const std::string&, sfilesize_t*);
		double *GetDoubles1D(const std::string&, sfilesize_t*);
		std::string *GetString(const std::string&);
        */
		virtual void Open(const std::string&, enum sfile_mode) override;
        H5::H5File *CreateMAT(const std::string&);

		virtual double GetAttributeScalar(const std::string&, const std::string&) override;
		virtual std::string *GetAttributeString(const std::string&, const std::string&) override;
        virtual std::string *GetString(const std::string&) override;

		virtual void WriteArray(const std::string&, double**, sfilesize_t, sfilesize_t) override;
		virtual void WriteAttribute_scalar(const std::string&, const std::string&, double) override;
		virtual void WriteAttribute_string(const std::string&, const std::string&, const std::string&) override;
		virtual void WriteString(const std::string&, const std::string&) override;

        void WriteMATLAB_class(const std::string&, const std::string&, H5::DataSet *dset=nullptr);
        void WriteMATLAB_int_decode(const std::string&, H5::DataSet *dset=nullptr, const int val=2);
};

#endif/*_SFILE_MAT_HDF5_H*/
