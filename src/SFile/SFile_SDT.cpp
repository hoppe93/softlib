/**
 * Implementation of the 'Self-Descriptive Text' format,
 * designed specifically for SOFT.
 *
 * The original motivation for this format was to allow
 * SOFT to read magnetic field data exported from systems
 * without any capabilities of generating HDF5 or MAT
 * files (such as those at ASDEX-U).
 *
 * SPECIFICATION
 * An SDT file consists of blocks of variables, separated by
 * empty lines. All data is given as ASCII text.
 *
 * A variable definition consists of two or more lines. The
 * first line defines the size and name of the variable, while
 * the following lines contain the data of the variable. The
 * example below illustrates how a 3-by-2 matrix named 'abc'
 * should be specified:
 *
 * @matrix abc 3 2
 * 1 2
 * 3 4
 * 5 6
 *
 * The available types are
 *   
 *   @matrix   - Matrix (2-dimensional list of double-precision 
 *               floating point values)
 *   @string   - String
 */

#include <cstdlib>
#include <iostream>
#include <fstream>
#include <string>
#include <softlib/SFile.h>
#include <softlib/SFileException.h>
#include <softlib/SFile_SDT.h>

using namespace std;

/**
 * Open SDT file.
 *
 * filename: Name of file to open.
 * openmode: Mode to open the file in (read/write etc.)
 */
void SFile_SDT::Open(const string &filename, enum sfile_mode openmode) {
    mode = openmode;

    switch (mode) {
        case SFILE_MODE_READ:
            sdtfile.open(filename, ios::in);
            Load();
            Close();
            break;
        case SFILE_MODE_UPDATE:
            sdtfile.open(filename, ios::in | ios::out);
            break;
        case SFILE_MODE_WRITE:
            sdtfile.open(filename, ios::out | ios::trunc);
            break;
        default:
            throw SFileException("Unrecognized option for opening the SDT file: %d.", mode);
    }
}

/**
 * Closes the SDT file.
 */
void SFile_SDT::Close() {
    sdtfile.close();
}

/**
 * Checks if a variable of the given name is present
 * in the SDT file.
 *
 * name: Name of variable to search for.
 */
bool SFile_SDT::HasVariable(const string& name) {
    if (matrices.find(name) == matrices.end() &&
        strings.find(name) == strings.end())
        return false;
    else return true;
}

/**
 * HDF5 compatibility function which combines the
 * name of a dataset with the name of a variable to emulate
 * HDF5 dataset attributes.
 *
 * dsetname: Name of dataset.
 * name: Name of variable.
 */
string SFile_SDT::GetAttributeName(const string& dsetname, const string& name) {
	string nname = dsetname + "_" + name;
	
	return nname;
}

/**
 * Loads the matrix with the given name from
 * the SDT file.
 * 
 * name: Name of matrix to load.
 * dims: Contains the matrix dimensions (rows, cols)
 *       upon return.
 */
double **SFile_SDT::GetDoubles(const string& name, sfilesize_t *dims) {
    mat_list_t::iterator it = matrices.find(name);

    if (it == matrices.end())
        throw SFileException("No matrix named '%s' in the SDT file.", name.c_str());

    dims[0] = it->second.nrows;
    dims[1] = it->second.ncols;

    return it->second.value;
}

/**
 * !!! NOT SUPPORTED !!!
 *
 * Loads the multidimensional array with the given
 * name from the SDT file.
 *
 * !!! NOT SUPPORTED !!!
 */
double *SFile_SDT::GetMultiArray_linear(const string& name, const sfilesize_t, sfilesize_t&, sfilesize_t*) {
	throw SFileException("When loading variable '%s': The SDT format has no support for multidimensional arrays.", name.c_str());
}

/**
 * Loads the vector with the given name from
 * the SDT file.
 *
 * name: Name of vector to load.
 * dims: Contains the number of elements in the vector
 *       upon return.
 */
double *SFile_SDT::GetDoubles1D(const string &name, sfilesize_t *dims) {
    mat_list_t::iterator it = matrices.find(name);

    if (it == matrices.end())
        throw SFileException("No vector named '%s' in the SDT file.", name.c_str());

    if (it->second.nrows == 1)
        *dims = it->second.ncols;
    else if (it->second.ncols != 1)
        *dims = it->second.nrows;
    else
        throw SFileException("The requested variable is not a vector: '%s'.", name.c_str());

    return it->second.value[0];
}

/**
 * Loads the scalar with the given name
 * as an attribute.
 *
 * datasetname: Name of dataset owning the attribute.
 * name:        Name of attribute scalar to load.
 */
double SFile_SDT::GetAttributeScalar(const string& datasetname, const string &name) {
    string att_name = GetAttributeName(datasetname, name);

    mat_list_t::iterator it = matrices.find(att_name);

    if (it == matrices.end())
        throw SFileException("No attribute scalar named '%s' for dataset '%s' in the SDT file.", name.c_str(), datasetname.c_str());

    if (it->second.nrows != 1 || it->second.ncols != 1)
        throw SFileException("The requested variable in dataset '%s' is not a scalar: '%s'.", datasetname.c_str(), name.c_str());

    return it->second.value[0][0];
}

/**
 * Loads the string with the given name
 * as an attribute.
 *
 * datasetname: Name of dataset owning the attribute.
 * name:        Name of attribute string to load.
 */
string SFile_SDT::GetAttributeString(const string &datasetname, const string &name) {
    string att_name = GetAttributeName(datasetname, name);

    str_list_t::iterator it = strings.find(att_name);

    if (it == strings.end())
        throw SFileException("No attribute string named '%s' for dataset '%s' in the SDT file.", name.c_str(), datasetname.c_str());

    return it->second.value;
}

/**
 * Loads the string with the given name from
 * the file.
 *
 * name: Name of string to load.
 */
string SFile_SDT::GetString(const string &name) {
    str_list_t::iterator it = strings.find(name);

    if (it == strings.end())
        throw SFileException("No string named '%s' in the SDT file.", name.c_str());

    return it->second.value;
}


/************************
 * ROUTINES FOR WRITING *
 ************************/
void SFile_SDT::WriteArray(const string& name, double **arr, sfilesize_t nrows, sfilesize_t ncols) {
    // Definition string
    sdtfile << "@matrix " << name << " " << nrows << " " << ncols << endl;

    for (sfilesize_t i = 0; i < nrows; i++) {
        for (sfilesize_t j = 0; j < ncols; j++) {
            sdtfile << arr[i][j] << " ";
        }

        sdtfile << endl;
    }

    sdtfile << endl;
}

void SFile_SDT::WriteAttribute_scalar(const string& datasetname, const string& name, double val) {
    string att_name = GetAttributeName(datasetname, name);

    sdtfile << "@matrix " << att_name << " 1 1" << endl << val << endl << endl;
}
void SFile_SDT::WriteAttribute_string(const string& datasetname, const string& name, const string& val) {
    string att_name = GetAttributeName(datasetname, name);
    WriteString(att_name, val);
}

void SFile_SDT::WriteMultiArray(const string& name, double*, sfilesize_t, sfilesize_t*) {
	throw SFileException("Trying to write variable '%s': The SDT format does not support higher than 2-dimensional arrays.", name.c_str());
}

void SFile_SDT::WriteImage(const string& name, double **img, sfilesize_t pixels) {
    WriteArray(name, img, pixels, pixels);
}

void SFile_SDT::WriteList(const string& name, double *list, sfilesize_t length) {
    WriteArray(name, &list, 1, length);
}

void SFile_SDT::WriteString(const string& name, const string& val) {
    sdtfile << "@string " << name << " " << val.size() << endl;

    // Encode and write string value
    for (const char &c : val)
        sdtfile << (int)c << " ";

    sdtfile << endl << endl;
}

