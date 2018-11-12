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
		double GetAttributeScalar(const std::string&, const std::string&);
		std::string *GetAttributeString(const std::string&, const std::string&);
		double **GetDoubles(const std::string&, sfilesize_t*);
		double *GetDoubles1D(const std::string&, sfilesize_t*);
		std::string *GetString(const std::string&);
		void Open(const std::string&, enum sfile_mode);
		void WriteArray(const std::string&, double**, sfilesize_t, sfilesize_t);
		void WriteAttribute_scalar(const std::string&, const std::string&, double);
		void WriteAttribute_string(const std::string&, const std::string&, const std::string&);
		void WriteImage(const std::string&, double**, sfilesize_t);
		void WriteList(const std::string&, double*, sfilesize_t);
		void WriteString(const std::string&, const std::string&);
};

#endif/*_SFILE_MATLAB_H*/
