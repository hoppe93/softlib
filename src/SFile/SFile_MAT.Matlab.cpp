/**
 * A SOFT MAT interface for simplified I/O of
 * MATLAB .MAT files.
 */

#include <cstring>
#include <string>
#include <mat.h>
#include <matrix.h>
#include <softlib/SFile.h>
#include <softlib/SFile_MAT.Matlab.h>
#include <softlib/SFileException.h>

using namespace std;

SFile_MAT::~SFile_MAT() { Close(); }

/**
 * Closes the MATLAB file.
 */
void SFile_MAT::Close(void) {
	matClose(mfp);
	mfp = NULL;
}

/**
 * Opens a MATLAB '.mat' file for read/write.
 *
 * filename: Name of file to open
 * ot: Purpose for opening file (read or write)
 */
void SFile_MAT::Open(const string& filename, enum sfile_mode openmode) {
	mode = openmode;
	mfp = NULL;

	switch (openmode) {
		case SFILE_MODE_READ:
			mfp = matOpen(filename.c_str(), "r");
			break;
		case SFILE_MODE_UPDATE:
			mfp = matOpen(filename.c_str(), "u");
			break;
		case SFILE_MODE_WRITE:
			mfp = matOpen(filename.c_str(), "w7.3");
			break;
		default:
			throw SFileException("Unrecognized option for opening MATLAB file.");
	}

	if (mfp == NULL) {
		throw SFileException("Unable to create or open MATLAB file: "+filename);
	}

	this->filename = filename;
}

/**
 * Check if this file contains a variable with
 * the given name.
 */
bool SFile_MAT::HasVariable(const string& name) {
	mxArray *arr = matGetVariable(mfp, name.c_str());
    bool exists = (arr != NULL);

    if (exists)
        mxDestroyArray(arr);

    return exists;
}

/******************************
 ************ INPUT ***********
 ******************************/
/**
 * Reads an attribute variable as a real scalar.
 *
 * dsetname: Name of dataset attribute belongs to.
 * name: Name of attribute variable.
 */
double SFile_MAT::GetAttributeScalar(const string& dsetname, const string& name) {
	string aname = GetAttributeName(dsetname, name);
	double s = GetScalar(aname);

	return s;
}
/**
 * Reads an attribute variable as a string.
 *
 * dsetname: Name of dataset attribute belongs to.
 * name: Name of attribute variable.
 */
string SFile_MAT::GetAttributeString(const string& dsetname, const string& name) {
	string aname = GetAttributeName(dsetname, name);
	string s = GetString(aname);

	return s;
}

/**
 * Reads an array of C doubles from the MAT file.
 * The dimensions of the array are returned in
 * "dims".
 *
 * name: Name of variable to read
 * dims: Array with two elements. Contains the
 *  number of rows and columns in the returned
 *  array on return.
 *
 * RETURNS a 1-D array (logically 2-D) which contains
 * the data of the variable. If the named variable
 * does not exist, or is not a 2-D matrix, NULL
 * is returned.
 */
double **SFile_MAT::GetDoubles(const string& name, sfilesize_t *dims) {
	size_t rows, cols;
	mxArray *arr = matGetVariable(mfp, name.c_str());

	if (arr == NULL || mxIsEmpty(arr))
		throw SFileException(filename+": The variable '"+name+"' does not exist in the file.");
	if (!mxIsDouble(arr))
		throw SFileException(filename+": The variable '"+name+"' is not a vector or matrix.");
	
    // NOTE: MATLAB does column major matrices,
    // while we do row major. Therefore our 'rows'
    // correspond to MATLAB's 'cols'.
	rows = mxGetN(arr);
	cols = mxGetM(arr);

	if (dims != NULL) {
		dims[0] = rows;
		dims[1] = cols;
	}

	double *tdata = mxGetPr(arr);

    sfilesize_t len = rows*cols;
	double **data = new double*[rows];
    data[0] = new double[len];

	for (sfilesize_t i = 1; i < rows; i++)
		data[i] = data[i-1] + cols;

    for (sfilesize_t i = 0; i < len; i++)
        data[0][i] = tdata[i];

	mxDestroyArray(arr);
	return data;
}
/*double *SFile_MAT::GetDoubles1D(const string& name, sfilesize_t *dims) {
	size_t rows, cols;
	mxArray *arr = matGetVariable(mfp, name.c_str());

	if (arr == NULL || mxIsEmpty(arr))
		throw SFileException(filename+": The variable '"+name+"' does not exist in the file.");
	if (!mxIsDouble(arr))
		throw SFileException(filename+": The variable '"+name+"' is not a vector or matrix.");
	
    // NOTE: MATLAB does column major matrices,
    // while we do row major. Therefore our 'rows'
    // correspond to MATLAB's 'cols'.
	rows = mxGetN(arr);
	cols = mxGetM(arr);

	if (dims != NULL) {
		dims[0] = rows;
		dims[1] = cols;
	}

	double *data = mxGetPr(arr);
	return data;
}*/
mxArray *SFile_MAT::__GetArray(const string& name, sfilesize_t *dims) {
	size_t rows, cols;
	mxArray *arr = matGetVariable(mfp, name.c_str());

	if (arr == NULL || mxIsEmpty(arr))
		throw SFileException(filename+": The variable '"+name+"' does not exist in the file.");
	if (!mxIsDouble(arr))
		throw SFileException(filename+": The variable '"+name+"' is not a vector or matrix.");
	
    // NOTE: MATLAB does column major matrices,
    // while we do row major. Therefore our 'rows'
    // correspond to MATLAB's 'cols'.
	rows = mxGetN(arr);
	cols = mxGetM(arr);

	if (dims != NULL) {
		dims[0] = rows;
		dims[1] = cols;
	}

    return arr;
}
double *SFile_MAT::GetDoubles1D(const string& name, sfilesize_t *dims) {
    sfilesize_t d[2];
    mxArray *arr = __GetArray(name, dims);
    *dims = d[0]*d[1];
    return mxGetDoubles(arr);
}
int32_t *SFile_MAT::GetInt32_1D(const string& name, sfilesize_t *dims) {
    sfilesize_t d[2];
    mxArray *arr = __GetArray(name, d);
    *dims = d[0]*d[1];
    return mxGetInt32s(arr);
}
int64_t *SFile_MAT::GetInt64_1D(const string& name, sfilesize_t *dims) {
    sfilesize_t d[2];
    mxArray *arr = __GetArray(name, dims);
    *dims = d[0]*d[1];
    return mxGetInt64s(arr);
}
uint32_t *SFile_MAT::GetUInt32_1D(const string& name, sfilesize_t *dims) {
    sfilesize_t d[2];
    mxArray *arr = __GetArray(name, dims);
    *dims = d[0]*d[1];
    return mxGetUInt32s(arr);
}
uint64_t *SFile_MAT::GetUInt64_1D(const string& name, sfilesize_t *dims) {
    sfilesize_t d[2];
    mxArray *arr = __GetArray(name, dims);
    *dims = d[0]*d[1];
    return mxGetUint64s(arr);
}

/**
 * Reads an array of values from the MAT file.
 * The dimensions of the array are returned in
 * "dims".
 *
 * name: Name of variable to read
 * dims: Array with two elements. Contains the
 *  number of rows and columns in the returned
 *  array on return.
 *
 * RETURNS a 1-D array (logically 2-D) which contains
 * the data of the variable. If the named variable
 * does not exist, or is not a 2-D matrix, NULL
 * is returned.
 */
template<typename T>
T **SFile_MAT_GetArray2D(SFile_MAT *sf, const string& name, sfilesize_t *dims) {
    mArray *arr = sf->__GetArray(name, dims);

    sfilesize_t rows = dims[0], cols = dims[1];

    sfilesize_t len = rows*cols;
	T **data = new T*[rows];
    data[0] = new T[len];

	for (sfilesize_t i = 1; i < rows; i++)
		data[i] = data[i-1] + cols;

    for (sfilesize_t i = 0; i < len; i++)
        data[0][i] = tdata[i];

	mxDestroyArray(arr);
	return data;
}
double **SFile_MAT::GetDoubles(const string& name, sfilesize_t *dims) {
    return SFile_MAT_GetArray2D<double>(this, const string& name, dims);
}
int32_t **SFile_MAT::GetInt32_2D(const string& name, sfilesize_t *dims) {
    return SFile_MAT_GetArray2D<int32_t>(this, const string& name, dims);
}
int64_t **SFile_MAT::GetInt64_2D(const string& name, sfilesize_t *dims) {
    return SFile_MAT_GetArray2D<int64_t>(this, const string& name, dims);
}
uint32_t **SFile_MAT::GetUInt32_2D(const string& name, sfilesize_t *dims) {
    return SFile_MAT_GetArray2D<uint32_t>(this, const string& name, dims);
}
uint64_t **SFile_MAT::GetUInt64_2D(const string& name, sfilesize_t *dims) {
    return SFile_MAT_GetArray2D<uint64_t>(this, const string& name, dims);
}

/**
 * Loads a multidimensional array from the given file.
 * Returns the array as a single vector. This vector
 * can be turned into a C++ multidimensional array
 * of the appropriate sizes using 'SFile::GetMultiArray()'.
 *
 * name:   Name of variable to load.
 * nndims: Number of allowed elements in 'dims'.
 * ndims:  On return, contains the number of elements in 'dims'.
 * dims:   On return, contains a list of sizes of each
 *         dimension in the array.
 *
 * RETURNS the data of the array. If 'nndims' < 'ndims', then
 * a 'nullptr' is returned. All other failures lead to an
 * exception being thrown.
 */
double *SFile_MAT::GetMultiArray_linear(const string& name, const sfilesize_t nndims, sfilesize_t& ndims, sfilesize_t *dims) {
	mxArray *ma = matGetVariable(mfp, name.c_str());

	if (ma == NULL || mxIsEmpty(ma))
		throw SFileException(filename+": The variable '"+name+"' does not exist in the file.");
	if (!mxIsDouble(ma))
		throw SFileException(filename+": The variable '"+name+"' is not a vector or matrix.");
	
	ndims = mxGetNumberOfDimensions(ma);
	if (ndims > nndims)
		return nullptr;
	
	const mwSize *mxDims = mxGetDimensions(ma);
	sfilesize_t nel = 1;
	for (sfilesize_t i = 0; i < ndims; i++) {
		dims[i] = mxDims[i];
		nel *= dims[i];
	}
	
	double *tdata = mxGetPr(ma);
	double *data = new double[nel];
	memcpy(data, tdata, sizeof(slibreal_t)*nel);

	mxDestroyArray(ma);

	return data;
}

/**
 * Loads a variable from the MAT file as
 * a string.
 *
 * name: Name of variable to load
 *
 * Returns the variable as a C++ string. If
 * the string does not exist, or if the
 * variable cannot be interpreted as a
 * string, NULL is returned.
 */
string SFile_MAT::GetString(const string& name) {
	mxArray *arr = matGetVariable(mfp, name.c_str());

	if (arr == NULL || mxIsEmpty(arr))
		throw SFileException(filename+": The string '"+name+"' does not exist in the file.");
	if (mxGetClassID(arr) != mxCHAR_CLASS)
		throw SFileException(filename+": The variable '"+name+"' is not a string.");
	
	int length = mxGetN(arr)+1;
	char *str = new char[length];
	mxGetString(arr, str, length);
	
	mxDestroyArray(arr);

	/*string *s = new string(str);
	delete [] str;
	return s;*/
	return string(str);
}

/******************************
 *********** OUTPUT ***********
 ******************************/
/**
 * Create a MATLAB struct with the given name.
 * NOTE: This feature is currently unsupported when linking
 * with the official MATLAB libraries (i.e. when using
 * this implementation of 'SFile_MAT').
 *
 * name: Name of struct to create.
 */
void SFile_MAT::CreateStruct(const string& name) {
	throw SFileException("While attempting to create struct '%s': SFile with official MATLAB libraries does not support creating structs.", name.c_str());
}

void SFile_MAT::WriteString(const string& name, const string& str) {
	int status;
	mxArray *ms;

	ms = mxCreateString(str.c_str());
	if (ms == NULL) throw SFileException("Unable to create MATLAB string '"+name+"'.");

	status = matPutVariable(mfp, name.c_str(), ms);
	mxDestroyArray(ms);

	if (status != 0) throw SFileException("Unable to write string '"+name+"' to MATLAB file.");
}

void SFile_MAT::WriteArray(const string& name, double **arr, sfilesize_t rows, sfilesize_t cols) {
	mxArray *ma = mxCreateDoubleMatrix(cols, rows, mxREAL);
	if (ma == NULL) throw SFileException("Unable to allocate MATLAB array for variable '%s'.", name.c_str());

	double *t = mxGetPr(ma);
	for (sfilesize_t i = 0; i < rows; i++)
		for (sfilesize_t j = 0; j < cols; j++)
			t[i*cols + j] = arr[i][j];

	int status = matPutVariable(mfp, name.c_str(), ma);
	mxDestroyArray(ma);
	if (status != 0) throw SFileException("Unable to write variable '%s' to MATLAB file.", name.c_str());
}
#define DEF_WRITEINTARRAY(NAME, TYPE, CTYPE, FTYPE) \
    void SFile_MAT::NAME(const string& name, TYPE **arr, sfilesize_t rows, sfilesize_t cols) { \
        mxArray *ma = mxCreateNumericMatrix(cols, rows, CTYPE, mxREAL); \
        if (ma == NULL) throw SFileException("Unable to allocate MATLAB array for variable '%s'.", name.c_str()); \
        \
        TYPE *t = FTYPE(ma); \
        for (sfilesize_t i = 0; i < rows; i++) \
            for (sfilesize_t j = 0; j < cols; j++) \
                t[i*cols + j] = arr[i][j]; \
        \
        int status = matPutVariable(mfp, name.c_str(), ma); \
        mxDestroyArray(ma); \
        if (status != 0) throw SFileException("Unable to write variable '%s' to MATLAB file.", name.c_str()); \
    }

DEF_WRITEINTARRAY(WriteInt32Array, int32_t, mxINT32_CLASS, mxGetInt32s)
DEF_WRITEINTARRAY(WriteInt64Array, int64_t, mxINT64_CLASS, mxGetInt64s)
DEF_WRITEINTARRAY(WriteUInt32Array, uint32_t, mxUINT32_CLASS, mxGetUint32s)
DEF_WRITEINTARRAY(WriteUInt64Array, uint64_t, mxUINT64_CLASS, mxGetUint64s)

void SFile_MAT::WriteImage(const string& name, double **image, sfilesize_t n) {
	WriteArray(name, image, n, n);
}
void SFile_MAT::WriteList(const string& name, double *list, sfilesize_t n) {
	WriteArray(name, &list, 1, n);
}
void SFile_MAT::WriteInt32List(const string& name, int32_t *list, sfilesize_t n) {
    WriteInt32Array(name, &list, 1, n);
}
void SFile_MAT::WriteInt32List(const string& name, int64_t *list, sfilesize_t n) {
    WriteInt64Array(name, &list, 1, n);
}
void SFile_MAT::WriteInt32List(const string& name, uint32_t *list, sfilesize_t n) {
    WriteUInt32Array(name, &list, 1, n);
}
void SFile_MAT::WriteInt64List(const string& name, uint64_t *list, sfilesize_t n) {
    WriteUInt64Array(name, &list, 1, n);
}

/**
 * Writes the given multi-dimensional array to
 * the file. Note that the data is given as a single
 * pointer-to-double, meaning that data must be stored
 * contiguously in memory.
 *
 * name:  Name of variable to store.
 * arr:   Array to write (stored contiguously in memory).
 * ndims: Number of dimensions of array.
 * dims:  Array with 'ndims' elements, specifying the
 *        the number of elements in each dimension.
 */
void SFile_MAT::WriteMultiArray(const string& name, double *arr, sfilesize_t ndims, sfilesize_t *dims) {
	mwSize mxDims[ndims];

	for (sfilesize_t i = 0; i < ndims; i++)
		mxDims[i] = dims[i];

	mxArray *ma = mxCreateNumericArray(ndims, mxDims, mxDOUBLE_CLASS, mxREAL);
	if (ma == NULL) throw SFileException("Unable to allocate MATLAB array for variable '%s'.", name.c_str());

	// Calculate total number of elements
	sfilesize_t nel = 1;
	for (sfilesize_t i = 0; i < ndims; i++)
		nel *= dims[i];

	// Copy array
	double *t = (double*)mxGetPr(ma);
	memcpy(t, arr, sizeof(double)*nel);

	// Write to MAT file
	int status = matPutVariable(mfp, name.c_str(), ma);
	mxDestroyArray(ma);
	if (status != 0) throw SFileException("Unable to write variable '%s' to MATLAB file.", name.c_str());
}

/**
 * HDF5 compatibility function which combines the
 * name of a dataset with the name of a variable to emulate
 * HDF5 dataset attributes.
 *
 * dsetname: Name of dataset.
 * name: Name of variable.
 */
string SFile_MAT::GetAttributeName(const string& dsetname, const string& name) {
	string nname = dsetname + "_" + name;
	return nname;
}

/**
 * Since MATLAB files don't have attributes for
 * variables, we instead give attributes names
 * according to 'dsetname_name'.
 *
 * dsetname: Dataset name
 * name: Name of attribute
 * q: Value of scalar to write
 */
void SFile_MAT::WriteAttribute_scalar(const string& dsetname, const string& name, double q) {
	string nname = GetAttributeName(dsetname, name);
	WriteList(nname, &q, 1);
}
/**
 * Since MATLAB files don't have attributes for
 * variables, we instead give attributes names
 * according to 'dsetname_name'.
 *
 * dsetname: Dataset name
 * name: Name of attribute
 * str: Value of attribute
 * len: Length of attribute string
 */
void SFile_MAT::WriteAttribute_string(const string& dsetname, const string& name, const string& str) {
	string nname = GetAttributeName(dsetname, name);
	WriteString(nname, str);
}

