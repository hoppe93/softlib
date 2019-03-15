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
    //H5::Exception::dontPrint();
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
				//identifier = H5Fopen(filename, H5F_ACC_RDONLY, H5P_DEFAULT);
				this->file = new H5File(filename, H5F_ACC_RDONLY);
				break;
			case SFILE_MODE_UPDATE:
				this->file = new H5File(filename, H5F_ACC_RDWR);
				break;
			case SFILE_MODE_WRITE:
				this->file = new H5File(filename, H5F_ACC_TRUNC);
				break;
			default:
				//fprintf(stderr, "Unrecognized option for opening HDF5 file: %d.\n", mode);
				throw SFileException("Unrecognized option for opening HDF5 file.");
		}
	} catch (FileIException &error) {
		throw SFileException("Unable to open HDF5 file.\n"+error.getDetailMsg());
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
    //} catch (GroupIException &ex) {
        return false;
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
double **SFile_HDF5::GetDoubles(const string& name, sfilesize_t *dims) {
    if (!HasVariable(name))
        throw SFileException("A variable with the name '%s' does not exist in the file '%s'.", name.c_str(), filename.c_str());

	if (dims == nullptr)
		throw SFileException("Null-pointer given for storing length of vector.");

	double *data, **pointers;
	sfilesize_t ndims;
	unsigned int i;
	DataSet dset = file->openDataSet(name);
	DataSpace dspace = dset.getSpace();

	ndims = dspace.getSimpleExtentDims(dims);
	
	if (ndims == 1) {
		data = new double[dims[0]];
		pointers = new double*;
		pointers[0] = data;
	} else {
		data = new double[dims[0]*dims[1]];
		pointers = new double*[dims[0]];
		for (i = 0; i < dims[0]; i++) {
			pointers[i] = data + (i*dims[1]);
		}
	}

	dset.read(data, PredType::IEEE_F64LE, dspace);

	return pointers;
}
double *SFile_HDF5::GetDoubles1D(const string& name, sfilesize_t *dims) {
    if (!HasVariable(name))
        throw SFileException("A variable with the name '%s' does not exist in the file '%s'.", name.c_str(), filename.c_str());
	
	if (dims == nullptr)
		throw SFileException("Null-pointer given for storing length of vector.");

	double *data;
	sfilesize_t ndims;
	DataSet dset = file->openDataSet(name);
	DataSpace dspace = dset.getSpace();

	ndims = dspace.getSimpleExtentDims(dims);
	
	if (ndims == 1)
		data = new double[dims[0]];
	else
		data = new double[dims[0]*dims[1]];

	dset.read(data, PredType::IEEE_F64LE, dspace);

	return data;
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
void SFile_HDF5::WriteArray(const string& name, double **arr, sfilesize_t rows, sfilesize_t cols) {
	hsize_t dims[] = {rows, cols};
	DataSpace dspace(2, dims);
	DataSet dset = file->createDataSet(name, PredType::IEEE_F64LE, dspace);
	dset.write(arr[0], PredType::NATIVE_DOUBLE);
	dset.close();
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

