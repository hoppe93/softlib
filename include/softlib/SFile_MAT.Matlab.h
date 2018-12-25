#ifndef _SFILE_MATLAB_H
#define _SFILE_MATLAB_H

#include <mat.h>
#include <softlib/SFile.h>
#include <string>

class SFile_MAT : public SFile {
	private:
		MATFile *mfp;
		std::string *GetAttributeName(const std::string&, const std::string&);
	public:
		~SFile_MAT();

		void Close();
        virtual bool HasVariable(const std::string&) override;
		virtual double GetAttributeScalar(const std::string&, const std::string&) override;
		virtual std::string *GetAttributeString(const std::string&, const std::string&) override;
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

#endif/*_SFILE_MATLAB_H*/
