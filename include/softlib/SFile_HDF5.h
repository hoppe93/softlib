#ifndef _SFILE_HDF5_H
#define _SFILE_HDF5_H

#include "H5Cpp.h"
#include <softlib/SFile.h>

class SFile_HDF5 : public SFile {
	protected:
		H5::H5File *file;
	public:
		SFile_HDF5();
		~SFile_HDF5();

		void Close();
		virtual double GetAttributeScalar(const std::string&, const std::string&) override;
		virtual std::string* GetAttributeString(const std::string&, const std::string&) override;
		virtual double **GetDoubles(const std::string&, sfilesize_t*) override;
		virtual double *GetDoubles1D(const std::string&, sfilesize_t*) override;
		virtual std::string *GetString(const std::string&) override;
		virtual void Open(const std::string&, enum sfile_mode) override;
		virtual void WriteArray(const std::string&, double**, sfilesize_t, sfilesize_t) override;
		virtual void WriteAttribute_scalar(const std::string&, const std::string&, double) override;
		virtual void WriteAttribute_string(const std::string&, const std::string&, const std::string&) override;
		virtual void WriteImage(const std::string&, double**, sfilesize_t) override;
		virtual void WriteList(const std::string&, double*, sfilesize_t) override;
		virtual void WriteString(const std::string&, const std::string&) override;
};

#endif/*_SFILE_HDF5_H*/
