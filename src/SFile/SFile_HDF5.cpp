/**
 * A softlib HDF5 SFile interface
 */

//#include <hdf5.h>
#include <H5Cpp.h>
#include <string>
#include <softlib/SFile.h>
#include <softlib/SFile_HDF5.h>
#include <softlib/SFileException.h>

#include <iostream>

using namespace H5;
using namespace std;

SFile_HDF5::SFile_HDF5() {
    // Disable HDF5 error messages
    H5::Exception::dontPrint();
}
SFile_HDF5::~SFile_HDF5() {
	if (this->file != NULL)
		delete this->file;
}
/**
 * Open HDF5 file.
 *
 * filename: Name of file to open
 * openmode: SFile mode to open file in (read, write or update)
 */
void SFile_HDF5::Open(const string& filename, enum sfile_mode openmode) {
	mode = openmode;

	try {
		switch (openmode) {
			case SFILE_MODE_READ:
				this->file = new H5File(filename, H5F_ACC_RDONLY);
				break;
			case SFILE_MODE_UPDATE:
				this->file = new H5File(filename, H5F_ACC_RDWR);
				break;
			case SFILE_MODE_WRITE:
				this->file = new H5File(filename, H5F_ACC_TRUNC);
				break;
			default:
				throw SFileException("Unrecognized option for opening HDF5 file.");
		}
	} catch (FileIException &error) {
		throw SFileException("%s: Unable to open HDF5 file.\n%s", filename.c_str(), error.getCDetailMsg());
	}

	this->filename = filename;
}

/**
 * Close the currently open HDF5 file.
 */
void SFile_HDF5::Close(void) {
	this->file->close();
	delete this->file;
}

/**
 * Checks if a variable of the given name exists
 * in the file.
 */
bool SFile_HDF5::HasVariable(const string& s) {
    try {
        DataSet dset = file->openDataSet(s);
        dset.close();

        return true;
    } catch (FileIException &ex) {
        return false;
    }
}

/**
 * Returns the type of the given variable.
 *
 * name: Name of variable to check type of.
 * hint: Hint at what data type the user truly desires.
 *       For HDF5, we ignore this parameter, since
 *       HDF5 types are strong (and data cannot be
 *       as easily cast as in, e.g. SDT).
 */
enum SFile::sfile_data_type SFile_HDF5::GetDataType(const string& name, enum sfile_data_type) {
    try {
        DataSet dset = file->openDataSet(name);
        DataType type = dset.getDataType();

        if (type == PredType::STD_I32LE)
            return SFILE_DATA_INT32;
        else if (type == PredType::STD_I64LE)
            return SFILE_DATA_INT64;
        else if (type == PredType::STD_U32LE)
            return SFILE_DATA_UINT32;
        else if (type == PredType::STD_U64LE)
            return SFILE_DATA_UINT64;
        else if (type == PredType::IEEE_F64LE)
            return SFILE_DATA_DOUBLE;
        else if (type == PredType::C_S1)
            return SFILE_DATA_STRING;
        else if (type == PredType::STD_U16LE)
            return SFILE_DATA_STRING;
        else
            return SFILE_DATA_UNDEFINED;
    } catch (FileIException &ex) {
        return SFILE_DATA_UNDEFINED;
    }
}

/*******************************
 ************ INPUT ************
 *******************************/
/**
 * Reads a (real) scalar attribute from the dataset
 * named 'datasetname'.
 *
 * datasetname: Name of dataset to read attribute from.
 * name: Name of attribute to read.
 */
double SFile_HDF5::GetAttributeScalar(const string& datasetname, const string& name) {
    try {
        double s;
        DataSet dataset = file->openDataSet(datasetname);
        Attribute attr = dataset.openAttribute(name);
        DataType type = attr.getDataType();

        attr.read(type, &s);

        return s;
    } catch (FileIException &ex) {
        throw SFileException("Unable to read attribute string '%s'. HDF5 error: %s", datasetname.c_str(), ex.getCDetailMsg());
    }
}
/**
 * Reads a string attribute from the dataset named
 * 'datasetname'.
 *
 * datasetname: Name of dataset to read attribute from.
 * name: Name of attribute to read.
 */
string SFile_HDF5::GetAttributeString(const string& datasetname, const string& name) {
    try {
        string s;

        DataSet dataset = file->openDataSet(datasetname);
        Attribute attr = dataset.openAttribute(name);
        DataType type = attr.getDataType();

        attr.read(type, s);
        return s;
    } catch (FileIException &ex) {
        throw SFileException("Unable to read attribute string '%s'. HDF5 error: %s", datasetname.c_str(), ex.getCDetailMsg());
    }
}
/**
 * Reads a scalar from the dataset with name "name".
 */
/**
 * Reads a string from the dataset with name "name".
 *
 * name: Name of dataset to load string from
 */
string SFile_HDF5::GetString(const string& name) {
    if (!HasVariable(name))
        throw SFileException("A variable with the name '%s' does not exist in the file '%s'.", name.c_str(), filename.c_str());

    string s;

    DataSet dataset = file->openDataSet(name);
    DataType type = dataset.getDataType();
    DataSpace dspace = dataset.getSpace();

    dataset.read(s, type, dspace);

    return s;
}
/**
 * Reads an array of C doubles from the dataset
 * with name "name". The dimensions of the array are
 * returned in "dims".
 *
 * name: Name of dataset to read
 * dims: Array with two elements. Will contain dimensions of returned value on return.
 *
 * Returns a 1-D array (logically 2-D) which contains
 * the data of the dataset. If the named dataset does not exist,
 * NULL is returned.
 */
template<typename T>
T **SFile_HDF5::GetArray2D(const string& name, sfilesize_t *dims, PredType p) {
    if (!HasVariable(name))
        throw SFileException("A variable with the name '%s' does not exist in the file '%s'.", name.c_str(), filename.c_str());

	if (dims == nullptr)
		throw SFileException("Null-pointer given for storing length of vector.");

	T *data, **pointers;
	sfilesize_t ndims;
	unsigned int i;

	DataSet dset = this->file->openDataSet(name);
	DataSpace dspace = dset.getSpace();

	ndims = dspace.getSimpleExtentDims(dims);
	
	if (ndims == 1) {
        dims[1] = dims[0];
        dims[0] = 1;
		data = new T[dims[0]];
		pointers = new T*;
		pointers[0] = data;
	} else {
		data = new T[dims[0]*dims[1]];
		pointers = new T*[dims[0]];
		for (i = 0; i < dims[0]; i++) {
			pointers[i] = data + (i*dims[1]);
		}
	}

	dset.read(data, p, dspace);

	return pointers;
}
double **SFile_HDF5::GetDoubles(const string& name, sfilesize_t *dims) {
    return SFile_HDF5::GetArray2D<double>(name, dims, PredType::IEEE_F64LE);
}
int32_t **SFile_HDF5::GetInt32_2D(const string& name, sfilesize_t *dims) {
    return SFile_HDF5::GetArray2D<int32_t>(name, dims, PredType::STD_I32LE);
}
int64_t **SFile_HDF5::GetInt64_2D(const string& name, sfilesize_t *dims) {
    return SFile_HDF5::GetArray2D<int64_t>(name, dims, PredType::STD_I64LE);
}
uint32_t **SFile_HDF5::GetUInt32_2D(const string& name, sfilesize_t *dims) {
    return SFile_HDF5::GetArray2D<uint32_t>(name, dims, PredType::STD_U32LE);
}
uint64_t **SFile_HDF5::GetUInt64_2D(const string& name, sfilesize_t *dims) {
    return SFile_HDF5::GetArray2D<uint64_t>(name, dims, PredType::STD_U64LE);
}

/**
 * Reads an array of numbers from the dataset
 * with name "name". The dimensions of the array are
 * returned in "dims".
 *
 * name: Name of dataset to read
 * dims: Array with two elements. Will contain dimensions of returned value on return.
 *
 * Returns a 1-D array (logically 1-D) which contains
 * the data of the dataset. If the named dataset does not exist,
 * nullptr is returned.
 */
template<typename T>
T *SFile_HDF5_GetArray1D(SFile_HDF5 *sf, const string& name, sfilesize_t *dims, PredType p) {
    sfilesize_t d[2];
    if (!sf->HasVariable(name))
        throw SFileException("A variable with the name '%s' does not exist in the file '%s'.", name.c_str(), sf->filename.c_str());
	
	if (dims == nullptr)
		throw SFileException("Null-pointer given for storing length of vector.");

	T *data;
	sfilesize_t ndims;

	DataSet dset = sf->__GetFile()->openDataSet(name);
	DataSpace dspace = dset.getSpace();

	ndims = dspace.getSimpleExtentDims(d);
	
	if (ndims == 1) {
        *dims = d[0];
		data = new T[d[0]];
	} else {
        *dims = d[0]*d[1];
		data = new T[*dims];
    }

	dset.read(data, p, dspace);

	return data;
}
double *SFile_HDF5::GetDoubles1D(const string& name, sfilesize_t *dims) {
    return SFile_HDF5_GetArray1D<double>(this, name, dims, PredType::IEEE_F64LE);
}
int32_t *SFile_HDF5::GetInt32_1D(const string& name, sfilesize_t *dims) {
    return SFile_HDF5_GetArray1D<int32_t>(this, name, dims, PredType::STD_I32LE);
}
int64_t *SFile_HDF5::GetInt64_1D(const string& name, sfilesize_t *dims) {
    return SFile_HDF5_GetArray1D<int64_t>(this, name, dims, PredType::STD_I64LE);
}
uint32_t *SFile_HDF5::GetUInt32_1D(const string& name, sfilesize_t *dims) {
    return SFile_HDF5_GetArray1D<uint32_t>(this, name, dims, PredType::STD_U32LE);
}
uint64_t *SFile_HDF5::GetUInt64_1D(const string& name, sfilesize_t *dims) {
    return SFile_HDF5_GetArray1D<uint64_t>(this, name, dims, PredType::STD_U64LE);
}

/**
 * Same as 'GetMultiArray', but returns the array
 * as a "linear" array, i.e. as a single-pointer-to-double.
 */
double *SFile_HDF5::GetMultiArray_linear(const string& name, const sfilesize_t nndims, sfilesize_t &ndims, sfilesize_t *dims) {
    if (!HasVariable(name))
        throw SFileException("A variable with the name '%s' does not exist in the file '%s'.", name.c_str(), filename.c_str());

	DataSet dset = file->openDataSet(name);
	DataSpace dspace = dset.getSpace();

	ndims = dspace.getSimpleExtentNdims();
	if (ndims > nndims)
		return nullptr;

	dspace.getSimpleExtentDims(dims);
	
	sfilesize_t nel = 1;
	for (sfilesize_t i = 0; i < ndims; i++)
		nel *= dims[i];
		
	double *data = new double[nel];
	dset.read(data, PredType::IEEE_F64LE, dspace);

	return data;
}

/*******************************
 *********** OUTPUT ************
 *******************************/
 /**
  * Creates a struct in the output file with the given name.
  *
  * name: Name of struct to create.
  */
void SFile_HDF5::CreateStruct(const string& name) {
	Group *grp = new Group(file->createGroup(name));
	delete grp;
}

/**
 * Write a dataset consisting of a string
 *
 * fileid: hid_t id of file to write to
 * name: Name of dataset to create
 * str: String to write
 * length: Length of string to write, or <= 0
 *   to automatically determine the length.
 */
void SFile_HDF5::WriteString(const string& name, const string& str) {
	hsize_t dim[] = {str.size()};
	DataSpace dspace(1, dim);
	//StrType tp(0, H5T_VARIABLE);
	
	DataSet dset = file->createDataSet(name, PredType::C_S1, dspace);
	dset.write(str, PredType::C_S1);
	dset.close();
}

/**
 * Write a 2-D array to the target file.
 *
 * fileid: hid_t id of the HDF5 file
 * name: Name of the dataset to create
 * arr: 2-D array data to write
 * rows: Number of rows of array
 * cols: Number of columns of data
 */
template<typename T>
void SFile_HDF5::WriteNumArray(
    const string& name, T **arr, sfilesize_t rows,
    sfilesize_t cols, PredType op, PredType ip
) {
	hsize_t dims[] = {rows, cols};
	DataSpace dspace(2, dims);
	DataSet dset = file->createDataSet(name, op, dspace);
	dset.write(arr[0], ip);
	dset.close();
}

void SFile_HDF5::WriteArray(const string& name, double **arr, sfilesize_t rows, sfilesize_t cols) {
    return SFile_HDF5::WriteNumArray<double>(name, arr, rows, cols, PredType::IEEE_F64LE, PredType::NATIVE_DOUBLE);
}
void SFile_HDF5::WriteInt32Array(const string& name, int32_t **arr, sfilesize_t rows, sfilesize_t cols) {
    return SFile_HDF5::WriteNumArray<int32_t>(name, arr, rows, cols, PredType::STD_I32LE, PredType::NATIVE_INT32);
}
void SFile_HDF5::WriteInt64Array(const string& name, int64_t **arr, sfilesize_t rows, sfilesize_t cols) {
    return SFile_HDF5::WriteNumArray<int64_t>(name, arr, rows, cols, PredType::STD_I64LE, PredType::NATIVE_INT64);
}
void SFile_HDF5::WriteUInt32Array(const string& name, uint32_t **arr, sfilesize_t rows, sfilesize_t cols) {
    return SFile_HDF5::WriteNumArray<uint32_t>(name, arr, rows, cols, PredType::STD_U32LE, PredType::NATIVE_UINT32);
}
void SFile_HDF5::WriteUInt64Array(const string& name, uint64_t **arr, sfilesize_t rows, sfilesize_t cols) {
    return SFile_HDF5::WriteNumArray<uint64_t>(name, arr, rows, cols, PredType::STD_U64LE, PredType::NATIVE_UINT64);
}

/**
 * Write an image to the HDF5 file. This
 * function is just a wrapper for
 * 'shdf5_write_array()'.
 *
 * fileid: hid_t id of file to write to
 * name: Name of dataset to create
 * image: Image data to write
 * n: Number of pixels to write (image is assumed square)
 */
void SFile_HDF5::WriteImage(const string& name, double **image, sfilesize_t n) {
	WriteArray(name, image, n, n);
}
/**
 * Write a simple list (1-D array) of data
 * to the HDF5 file. This function is just
 * a wrapper for 'shdf5_write_array()'.
 *
 * fileid: hid_t id of file to write to
 * name: Name of dataset to create
 * list: Data to write
 * n: Number of elements in list
 */
void SFile_HDF5::WriteList(const string& name, double *list, sfilesize_t n) {
	WriteArray(name, &list, 1, n);
}
void SFile_HDF5::WriteInt32List(const string& name, int32_t *list, sfilesize_t n) {
    WriteInt32Array(name, &list, 1, n);
}
void SFile_HDF5::WriteInt64List(const string& name, int64_t *list, sfilesize_t n) {
    WriteInt64Array(name, &list, 1, n);
}
void SFile_HDF5::WriteUInt32List(const string& name, uint32_t *list, sfilesize_t n) {
    WriteUInt32Array(name, &list, 1, n);
}
void SFile_HDF5::WriteUInt64List(const string& name, uint64_t *list, sfilesize_t n) {
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
void SFile_HDF5::WriteMultiArray(const string& name, double *arr, sfilesize_t ndims, sfilesize_t *dims) {
	hsize_t *h_dims = new hsize_t[ndims];
	for (sfilesize_t i = 0; i < ndims; i++)
		h_dims[i] = dims[i];

	DataSpace dspace(ndims, dims);
	DataSet dset = file->createDataSet(name, PredType::IEEE_F64LE, dspace);
	dset.write(arr, PredType::NATIVE_DOUBLE);
	dset.close();
}

/**
 * Add a scalar attribute to a HDF5 dataset.
 *
 * fileid: hid_t id of file to write to
 * dsetname: Name of dataset to apply this attribute to
 * name: Name of the attribute to create
 * val: Value to give the attribute
 */
void SFile_HDF5::WriteAttribute_scalar(const string& dsetname, const string& name, double val) {
	DataSet dset = file->openDataSet(dsetname);
	DataSpace dspace(H5S_SCALAR);
	Attribute att = dset.createAttribute(name, PredType::IEEE_F64LE, dspace);
	att.write(PredType::NATIVE_DOUBLE, &val);
}

/**
 * Add a string attribute to a HDF5 dataset.
 *
 * fileid: hid_t id of file to write to
 * dsetname: Name of dataset to apply this attribute to
 * name: Name of the attribute to create
 * str: String to write to the attribute
 * length: Length of string to write, or <= 0
 *   to automatically determine length
 */
void SFile_HDF5::WriteAttribute_string(const string& dsetname, const string& name, const string& str) {
	sfilesize_t dims[] = {str.size()};

	DataSet dset = file->openDataSet(dsetname);
	DataSpace dspace(1, dims);

	Attribute att = dset.createAttribute(name, PredType::C_S1, dspace);
	att.write(PredType::C_S1, str);
}

