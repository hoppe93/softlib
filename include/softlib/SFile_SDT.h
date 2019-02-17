#ifndef _SFILE_SDT_H
#define _SFILE_SDT_H

#include <softlib/SFile.h>
#include <fstream>
#include <map>
#include <vector>
#include <string>

class SFile_SDT : public SFile {
    public:
        struct str {
            std::string name;
            std::string value;
        };
        struct matrix {
            std::string name;
            sfilesize_t nrows, ncols;
            double **value;
        };

        typedef std::map<std::string, struct matrix> mat_list_t;
        typedef std::map<std::string, struct str> str_list_t;
	private:
		std::fstream sdtfile;
		std::string GetAttributeName(const std::string&, const std::string&);

        mat_list_t matrices;
        str_list_t strings;

        void Load();
        void LoadMatrix(const std::string&, const sfilesize_t, const sfilesize_t);
        void LoadString(const std::string&, const sfilesize_t);
	public:
		void Close();
        virtual bool HasVariable(const std::string&) override;
		virtual double **GetDoubles(const std::string&, sfilesize_t*) override;
		virtual double *GetDoubles1D(const std::string&, sfilesize_t*) override;
        virtual double GetAttributeScalar(const std::string&, const std::string&) override;
        virtual std::string GetAttributeString(const std::string&, const std::string&) override;
		virtual std::string GetString(const std::string&) override;
		virtual void Open(const std::string&, enum sfile_mode) override;
		virtual void WriteArray(const std::string&, double**, sfilesize_t, sfilesize_t) override;
		virtual void WriteAttribute_scalar(const std::string&, const std::string&, double) override;
		virtual void WriteAttribute_string(const std::string&, const std::string&, const std::string&) override;
		virtual void WriteImage(const std::string&, double**, sfilesize_t) override;
		virtual void WriteList(const std::string&, double*, sfilesize_t) override;
		virtual void WriteString(const std::string&, const std::string&) override;
};

#endif/*_SFILE_SDT_H*/
