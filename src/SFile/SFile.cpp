/* SOFTLib file interface */

#include <iostream>
#include <string>
#include <functional>
#include <softlib/SFile.h>
#include <softlib/SFileException.h>

#ifdef SOFT_HDF5
#   include <softlib/SFile_HDF5.h>
#endif
#include <softlib/SFile_MAT.h>
#include <softlib/SFile_SDT.h>

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
		case SFILE_TYPE_SDT:
			sf = new SFile_SDT();
			break;
		default: 
			throw SFileException("Trying to open file of unrecognized format: '%s'.", filename.c_str());
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
 * Reads a multidimensional array from the file.
 *
 * name:   Name of variable to read.
 * nndims: Number of elements allowed in the 'dims' array.
 * ndims:  On return, containts the number of dimensions in array.
 * dims:   On return, points to a list of sizes of each dimension.
 *
 * RETURNS an ndims-pointer on success. If ndims > nndims,
 * the function returns a 'nullptr'. Note that all other
 * failures cause an exception.
 *
 * The array is stored contiguously in memory, so dereferencing
 * the pointer 'ndims-1' times gives a (double*) to the
 * contents of the entire array. Note that if you want to access
 * the array linearly, you should use 'GetMultiArray_linear()'
 * instead.
 */
void *SFile::GetMultiArray(const string& name, const sfilesize_t nndims, sfilesize_t &ndims, sfilesize_t *dims) {
	double *ptr = GetMultiArray_linear(name, nndims, ndims, dims);
	
	if (ptr == nullptr)
		return nullptr;

	// Convert to multidimensional pointer
	sfilesize_t nel = 1;
	for (sfilesize_t i = 0; i < ndims; i++)
		nel *= dims[i];

	// Yes, this is horrible abuse of pointers...
	double **arr;
	for (sfilesize_t i = ndims-2; true; i--) {
		sfilesize_t n = dims[i+1];
		nel /= n;

		arr = new double*[nel];
		for (sfilesize_t j = 0; j < nel; j++)
			arr[j] = ptr + j * n;
		
		ptr = (double*)arr;

		if (i == 0) break;
	}

	return (void*)arr;
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
		throw SFileException("%s: The requested variable is not a list: '%s'.", filename.c_str(), name.c_str());
	}

	return p;
}

/**
 * Loads the named array and returns it
 * converted to the desired type.
 *
 * name: Name of variable to load.
 * dims: Dimensions of the output array.
 */
template<typename Tout, typename Tin>
Tout **SFile::__Get2DAs(const string& name, sfilesize_t *dims, Tin** (SFile::*f)(const string&, sfilesize_t*), SFile *sf) {
    Tin **a = (sf->*f)(name, dims);

    Tout **b = new Tout*[dims[0]];
    b[0] = new Tout[dims[0]*dims[1]];
    for (sfilesize_t i = 1; i < dims[0]; i++)
        b[i] = b[i-1] + dims[1];

    for (sfilesize_t i = 0; i < dims[0]; i++) {
        for (sfilesize_t j = 0; j < dims[1]; j++) {
            b[i][j] = (Tout)a[i][j];
        }
    }

    delete [] a[0];
    delete [] a;

    return b;
}
template<typename Tout, typename Tin>
Tout *SFile::__Get1DAs(const string& name, sfilesize_t *dims, Tin* (SFile::*f)(const string&, sfilesize_t*), SFile *sf) {
    Tin *a = (sf->*f)(name, dims);

    Tout *b = new Tout[*dims];
    for (sfilesize_t i = 0; i < *dims; i++)
        b[i] = (Tout)a[i];

    delete [] a;

    return b;
}

template<typename T>
T **SFile::__Get2DAsType(const string& name, sfilesize_t *dims) {
    enum sfile_data_type hint = SFILE_DATA_UNDEFINED;

    if (typeid(T) == typeid(double)) hint = SFILE_DATA_DOUBLE;
    else if (typeid(T) == typeid(int32_t)) hint = SFILE_DATA_INT32;
    else if (typeid(T) == typeid(int64_t)) hint = SFILE_DATA_INT64;
    else if (typeid(T) == typeid(uint32_t)) hint = SFILE_DATA_UINT32;
    else if (typeid(T) == typeid(uint64_t)) hint = SFILE_DATA_UINT64;

    switch (GetDataType(name, hint)) {
        case SFILE_DATA_DOUBLE: return __Get2DAs<T,double>(name, dims, &SFile::GetDoubles, this);
        case SFILE_DATA_INT32:  return __Get2DAs<T,int32_t>(name, dims, &SFile::GetInt32_2D, this);
        case SFILE_DATA_INT64:  return __Get2DAs<T,int64_t>(name, dims, &SFile::GetInt64_2D, this);
        case SFILE_DATA_UINT32: return __Get2DAs<T,uint32_t>(name, dims, &SFile::GetUInt32_2D, this);
        case SFILE_DATA_UINT64: return __Get2DAs<T,uint64_t>(name, dims, &SFile::GetUInt64_2D, this);

        default:
            throw SFileException("%s: Unable to load the variable '%s' as an array of the desired type.", filename.c_str(), name.c_str());
    }
}
template<typename T>
T *SFile::__Get1DAsType(const string& name, sfilesize_t *dims) {
    enum sfile_data_type hint = SFILE_DATA_UNDEFINED;

    if (typeid(T) == typeid(double)) hint = SFILE_DATA_DOUBLE;
    else if (typeid(T) == typeid(int32_t)) hint = SFILE_DATA_INT32;
    else if (typeid(T) == typeid(int64_t)) hint = SFILE_DATA_INT64;
    else if (typeid(T) == typeid(uint32_t)) hint = SFILE_DATA_UINT32;
    else if (typeid(T) == typeid(uint64_t)) hint = SFILE_DATA_UINT64;

    switch (GetDataType(name, hint)) {
        case SFILE_DATA_INT32:  return __Get1DAs<T,int32_t>(name, dims, &SFile::GetInt32_1D, this);
        case SFILE_DATA_INT64:  return __Get1DAs<T,int64_t>(name, dims, &SFile::GetInt64_1D, this);
        case SFILE_DATA_UINT32: return __Get1DAs<T,uint32_t>(name, dims, &SFile::GetUInt32_1D, this);
        case SFILE_DATA_UINT64: return __Get1DAs<T,uint64_t>(name, dims, &SFile::GetUInt64_1D, this);
        case SFILE_DATA_DOUBLE: return __Get1DAs<T,double>(name, dims, &SFile::GetList, this);

        default:
            throw SFileException("%s: Unable to load the variable '%s' as an array of the desired type.", filename.c_str(), name.c_str());
    }
}

template<typename T>
T SFile::__GetScalarAsType(const string& name) {
    enum sfile_data_type hint = SFILE_DATA_UNDEFINED;

    if (typeid(T) == typeid(double)) hint = SFILE_DATA_DOUBLE;
    else if (typeid(T) == typeid(int32_t)) hint = SFILE_DATA_INT32;
    else if (typeid(T) == typeid(int64_t)) hint = SFILE_DATA_INT64;
    else if (typeid(T) == typeid(uint32_t)) hint = SFILE_DATA_UINT32;
    else if (typeid(T) == typeid(uint64_t)) hint = SFILE_DATA_UINT64;

    switch (GetDataType(name, hint)) {
        case SFILE_DATA_INT32:  return (T)GetInt32(name);
        case SFILE_DATA_INT64:  return (T)GetInt64(name);
        case SFILE_DATA_UINT32: return (T)GetUInt32(name);
        case SFILE_DATA_UINT64: return (T)GetUInt64(name);
        case SFILE_DATA_DOUBLE: return (T)GetScalar(name);

        default:
            throw SFileException("%s: Unable to load the variable '%s' as an array of the desired type.", filename.c_str(), name.c_str());
    }
}

/**
 * Returns the given variable as an integer array.
 *
 * name: Name of variable to load.
 * dims: Contains dimensions of array on return.
 */
int64_t **SFile::GetIntArray(const string& name, sfilesize_t *dims) {
    return __Get2DAsType<int64_t>(name, dims);
}

int64_t *SFile::GetIntList(const string& name, sfilesize_t *dims) {
    return __Get1DAsType<int64_t>(name, dims);
}

int64_t SFile::GetInt(const string& name) {
    return __GetScalarAsType<int64_t>(name);
}

/**
 * Reads the variable expecting it to be a
 * scalar value of the given type.
 * 
 * name: Name of variable to read.
 * f:    Function to use to read the variable.
 */
template<typename T>
T SFile::__GetSingle(const string& name, T** (SFile::*f)(const string&, sfilesize_t*), SFile *sf) {
	sfilesize_t length[2] = {0};
	T s;
	T **arr = (sf->*f)(name, length);

	if (length[0] != 1 || length[1] != 1)
		throw SFileException("%s: The requested variable is not a scalar: '%s'.", filename.c_str(), name.c_str());
	
	s = arr[0][0];
	delete [] arr[0];
	delete [] arr;

	return s;
}
double SFile::GetScalar(const string& name) { return __GetSingle<double>(name, &SFile::GetDoubles, this); }
int32_t SFile::GetInt32(const string& name) { return __GetSingle<int32_t>(name, &SFile::GetInt32_2D, this); }
int64_t SFile::GetInt64(const string& name) { return __GetSingle<int64_t>(name, &SFile::GetInt64_2D, this); }
uint32_t SFile::GetUInt32(const string& name) { return __GetSingle<uint32_t>(name, &SFile::GetUInt32_2D, this); }
uint64_t SFile::GetUInt64(const string& name) { return __GetSingle<uint64_t>(name, &SFile::GetUInt64_2D, this); }

/**
 * Write a scalar value.
 */
void SFile::WriteScalar(const string& name, double s) {
    WriteList(name, &s, 1);
}
