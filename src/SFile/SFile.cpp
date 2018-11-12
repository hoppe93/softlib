/* SOFTLib file interface */

#include <iostream>
#include <string>
#include <softlib/SFile.h>
#include <softlib/SFileException.h>

#ifdef SOFT_HDF5
#   include <softlib/SFile_HDF5.h>
#endif
#include <softlib/SFile_MAT.h>
//#include <softlib/SFile_SDT.h>

using namespace std;

/**
 * Create and open a new SFile object with access to
 * the file with name 'filename'. The filetype
 * is automatically determined based on the
 * file extension.
 *
 * filename: Name of file to open.
 */
SFile *SFile::Create(const string& filename, enum sfile_mode mode) {
	return Create(filename, mode, TypeOfFile(filename));
}
/**
 * Create and open a new SFile object with access to
 * the file with name 'filename'. The 'type'
 * parameter specifies the way in which to
 * interpret the file, i.e. as a HDF5, MATLAB
 * or SDT file.
 *
 * filename: Name of file to open.
 */
SFile *SFile::Create(const string& filename, enum sfile_mode mode, enum sfile_type type) {
	SFile *sf;

	switch (type) {
		case SFILE_TYPE_HDF5:
#ifdef SOFT_HDF5
			sf = new SFile_HDF5();
			break;
#else
			throw SFileException("This version of softlib was not compiled with HDF5 support.");
#endif/*USE_HDF5*/
		case SFILE_TYPE_MATLAB:
			sf = new SFile_MAT();
			break;
			throw SFileException("This version of softlib was not compiled with Matlab support.");
/*
		case SFILE_TYPE_SDT:
			sf = new SFile_SDT();
			break;
*/
		default: 
			//cerr << "softlib ERROR: Trying to open file of unrecognized format: " << type << "." << endl;
			//exit(-1);
			throw SFileException("Trying to open file of unrecognized format.");
	}

	sf->filetype = type;
	sf->Open(filename, mode);

	return sf;
}

/**
 * Converts a type string to a 'sfile_type'.
 *
 * type: String/name of filetype to convert.
 */
enum sfile_type SFile::GetFileType(const string& type) {
	if (type == "h5") return SFILE_TYPE_HDF5;
	if (type == "hdf5") return SFILE_TYPE_HDF5;
	if (type == "mat") return SFILE_TYPE_MATLAB;
	if (type == "sdt") return SFILE_TYPE_SDT;
	if (type == "out") return SFILE_TYPE_SDT;

	return SFILE_TYPE_UNKNOWN;
}
/**
 * Identifies the filetype of a file based on
 * its filename extension. Assumes filename extension
 * to be separated by dot ('.') from the rest of the
 * filename.
 *
 * filename: Name of file.
 */
enum sfile_type SFile::TypeOfFile(const string& filename) {
	size_t i = filename.find_last_of(".");

	if (i == string::npos) return SFILE_TYPE_UNKNOWN;
	else return GetFileType(filename.substr(i+1));
}

/***************************************
 * HIGH-LEVEL I/O FUNCTIONS
 ***************************************/
/**
 * Reads a list from the file.
 *
 * name: Name of list variable.
 * length: Length of list variable.
 */
double *SFile::GetList(const string& name, sfilesize_t *length) {
	sfilesize_t len[2];
	double **arr = GetDoubles(name, len);
	double *p;

	if (len[0] == 1 || len[1] == 1) {
		if (len[0] == 1) *length = len[1];
		else *length = len[0];

		p = arr[0];
		delete [] arr;
	} else {
		throw new SFileException(filename+": The requested variable is not a list: '"+name+"'.");
	}

	return p;
}

/**
 * Reads a scalar (real) variable from the file.
 *
 * name: Name of scalar variable.
 */
double SFile::GetScalar(const string& name) {
	sfilesize_t length[2];
	double s;
	double **arr = GetDoubles(name, length);

	if (length[0] != 1 || length[1] != 1)
		throw new SFileException(filename+": The requested variable is not a real scalar: '"+name+"'.");
	
	s = arr[0][0];
	delete [] arr[0];
	delete [] arr;

	return s;
}

/**
 * Write a scalar value.
 */
void SFile::WriteScalar(const string& name, double s) {
    WriteList(name, &s, 1);
}
