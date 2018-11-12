#ifndef _SFILE_SDT_H
#define _SFILE_SDT_H

#include <softlib/SFile.h>
#include <fstream>
#include <vector>

class SFile_SDT : public SFile {
	private:
		std::ofstream sdtfile;
		char *GetAttributeName(const char*, const char*);
		std::vector<ssdt_key> keys;
	public:
		void Close();
		double **GetDoubles(const char*, sfilesize_t);
		double *GetDoubles1D(const char*, sfilesize_t);
		char *GetString(const char*);
		void Open(const char*, enum sfile_mode);
		void WriteArray(const char*, double**, int, int);
		void WriteAttribute_scalar(const char*, const char*, double);
		void WriteAttribute_string(const char*, const char*, const char*, int);
		void WriteImage(const char*, double*, int);
		void WriteString(const char*, const char*, int);
};

#endif/*_SFILE_SDT_H*/
