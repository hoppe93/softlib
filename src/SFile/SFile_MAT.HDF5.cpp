/**
 * Implementation of support for Matlab MAT files
 * using only the HDF5 library.
 *
 * Since no specification of the HDF5-based MAT file format is
 * published, this library relies on empirical studies. It is
 * however able to make assumptions on the data that greatly
 * improves the API performance, and allows simultaneous support
 * for MAT and HDF5.
 */

#include <chrono>
#include <ctime>
#include <H5Cpp.h>
#include <string>
#include <softlib/SFile_MAT.HDF5.h>
#include <softlib/SFileException.h>

using namespace std;
using namespace H5;
using namespace std::chrono;

/**
 * Open MAT file.
 *
 * filename: Name of file to open
 * openmode: SFile mode to open file in (read, write or update)
 */
void SFile_MAT::Open(const string& filename, enum sfile_mode openmode) {
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
				//this->file = new H5File(filename, H5F_ACC_TRUNC);
                this->file = CreateMAT(filename);
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
 * Create a new MAT file.
 * This implementation is based on the one found in
 * the Python library 'hdf5storage'. Unofficial
 * documentation for the MATLAB header is available at:
 *
 *   https://pythonhosted.org/hdf5storage/storage_format.html#matlab-file-header
 *
 * filename: Name of MAT file to create.
 */
H5File *SFile_MAT::CreateMAT(const string& filename) {
    #define MATLAB_STRING_LENGTH 116
    #define MATLAB_USER_BLOCK_SIZE 512
    #define NMAGIC 3
    char *s = new char[MATLAB_STRING_LENGTH+1];
    int32_t magic[NMAGIC] = {0x00000000, 0x00000000, 0x4D490200};
    time_point<system_clock> mom = system_clock::now();
    time_t tt = system_clock::to_time_t(mom);
    tm local_tm = *localtime(&tt);

    int
        dayw = local_tm.tm_wday,
        imon = local_tm.tm_mon-1,
        iday = local_tm.tm_mday,
        hour = local_tm.tm_hour,
        minute = local_tm.tm_min,
        second = local_tm.tm_sec,
        year = 1900 + local_tm.tm_year;

    const char days[7][4] = {
        "Sun", "Mon", "Tue", "Wed",
        "Thu", "Fri", "Sat"
    };
    const char months[12][4] = {
        "Jan", "Feb", "Mar", "Apr", "May", "Jun",
        "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"
    };

    int n = snprintf(
        s, MATLAB_STRING_LENGTH+1,
        "MATLAB 7.3 MAT-file, "
        "Platform: libsoft, "
        "Created on: %s %s %02d %02d:%02d:%02d %d "
        "HDF5 schema 1.00 .",
        days[dayw], months[imon], iday,
        hour, minute, second, year
    );

    if (n < 0)
        throw SFileException("Encoding error occured when generating MAT file header.");
    else if (n > MATLAB_STRING_LENGTH)
        throw SFileException("Error when writing MAT file header: string too long (%d bytes too many).", n-MATLAB_STRING_LENGTH);

    // Pad with spaces
    for (int i = n; i < MATLAB_STRING_LENGTH; i++)
        s[i] = ' ';
    s[MATLAB_STRING_LENGTH] = 0;

    // Create initial HDF5 file
    FileCreatPropList fcpl(FileCreatPropList::DEFAULT);
    fcpl.setUserblock(MATLAB_USER_BLOCK_SIZE);
    H5File *hf = new H5File(filename, H5F_ACC_TRUNC, fcpl);
    hf->close();

    // Write file
    FILE *f = fopen(filename.c_str(), "r+b");
    if (!f) {
        perror("ERROR");
        throw SFileException("Unable to create MAT file '%s'.", filename.c_str());
    }

    fwrite(s, 1, MATLAB_STRING_LENGTH, f);
    fwrite(magic, sizeof(int32_t), NMAGIC, f);
    fclose(f);

    delete [] s;

    // Re-open using HDF5 API
    hf = new H5File(filename, H5F_ACC_RDWR);

    return hf;
}

/************************
 ******** INPUT *********
 ************************/
/**
 * Most of the routines in the input section are
 * derived from the implementation of HDF5, and
 * so are not reimplemented here.
 */

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
 * Reads an attribute variable as a real scalar.
 *
 * dsetname: Name of dataset attribute belongs to.
 * name: Name of attribute variable.
 */
double SFile_MAT::GetAttributeScalar(const string& dsetname, const string& name) {
	string aname = GetAttributeName(dsetname, name);
	return GetScalar(aname);
}

/**
 * Reads an attribute variable as a string.
 *
 * dsetname: Name of dataset attribute belongs to.
 * name: Name of attribute variable.
 */
string SFile_MAT::GetAttributeString(const string& dsetname, const string& name) {
	string aname = GetAttributeName(dsetname, name);
	return GetString(aname);
}

/**
 * Reads a string.
 *
 * name: Name of variable to read.
 */
string SFile_MAT::GetString(const string& name) {
    if (!HasVariable(name))
        throw SFileException("A variable with the name '%s' does not exist in the file '%s'.", name.c_str(), filename.c_str());

    string s;
    DataSet dataset = file->openDataSet(name);
    DataType type = dataset.getDataType();
    DataSpace dspace = dataset.getSpace();

    sfilesize_t dims[2]={0,0};
    if (dspace.getSimpleExtentDims(dims) != 2)
        throw SFileException("Invalid number of dimensions for string.");
    
    if (type == PredType::STD_U16LE) {
        sfilesize_t len = (dims[1]==1?dims[0]:dims[1]);
        unsigned short arr[len];
        dataset.read(arr, type, dspace);

        // Convert to proper string
        char *str = new char[len+1];
        for (sfilesize_t i = 0; i < len; i++)
            str[i] = (char)arr[i];
        str[len] = 0;

        s = str;
    } else if (type == PredType::C_S1) {
        dataset.read(s, type, dspace);
    } else
        throw SFileException("Unrecognized data type for string.");

    return s;
}

/************************
 ******** OUTPUT ********
 ************************/
/**
 * Creates a MATLAB 'struct' with the given name
 * in the file.
 *
 * name: Name of the struct to create.
 */
void SFile_MAT::CreateStruct(const string& name) {
	Group *grp = new Group(file->createGroup(name));
	WriteMATLAB_class(name, "struct", grp);

	delete grp;
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
void SFile_MAT::WriteArray(const string& name, double **arr, sfilesize_t rows, sfilesize_t cols) {
    SFile_HDF5::WriteArray(name, arr, rows, cols);
    WriteMATLAB_class(name, "double");
}
void SFile_MAT::WriteInt32Array(const string& name, int32_t **arr, sfilesize_t rows, sfilesize_t cols) {
    SFile_HDF5::WriteInt32Array(name, arr, rows, cols);
    WriteMATLAB_class(name, "int32");
}
void SFile_MAT::WriteInt64Array(const string& name, int64_t **arr, sfilesize_t rows, sfilesize_t cols) {
    SFile_HDF5::WriteInt64Array(name, arr, rows, cols);
    WriteMATLAB_class(name, "int64");
}
void SFile_MAT::WriteUInt32Array(const string& name, uint32_t **arr, sfilesize_t rows, sfilesize_t cols) {
    SFile_HDF5::WriteUInt32Array(name, arr, rows, cols);
    WriteMATLAB_class(name, "uint32");
}
void SFile_MAT::WriteUInt64Array(const string& name, uint64_t **arr, sfilesize_t rows, sfilesize_t cols) {
    SFile_HDF5::WriteUInt64Array(name, arr, rows, cols);
    WriteMATLAB_class(name, "uint64");
}

void SFile_MAT::WriteAttribute_scalar(const string& dsetname, const string& name, double q) {
    string nname = GetAttributeName(dsetname, name);
    WriteList(nname, &q, 1);
}
void SFile_MAT::WriteAttribute_string(const string& dsetname, const string& name, const string &str) {
    string nname = GetAttributeName(dsetname, name);
    WriteString(nname, str);
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
	SFile_HDF5::WriteMultiArray(name, arr, ndims, dims);
	WriteMATLAB_class(name, "double");
}

void SFile_MAT::WriteString(const string& name, const string& str) {
    sfilesize_t dims[] = {str.size(),1};
    DataSpace dspace(2, dims);

    DataSet dset = file->createDataSet(name, PredType::STD_U16LE, dspace);
    // Convert string to array of shorts
    unsigned short arr[dims[0]];
    for (sfilesize_t i = 0; i < dims[0]; i++)
        arr[i] = (unsigned short)str.at(i);

    dset.write(arr, PredType::STD_U16LE);

    // Write attributes required by MATLAB
    WriteMATLAB_class(name, "char", &dset);
    WriteMATLAB_int_decode(name, &dset);

    dset.close();
}

/**
 * Set the "MATLAB_class" attribute of a given dataset.
 *
 * name: Name of dataset.
 * cls:  MATLAB class.
 * obj:  H5Object to apply the attribute to. If nullptr,
 *       opens the dataset.
 */
void SFile_MAT::WriteMATLAB_class(const string& name, const string& cls, H5Object *obj) {
    if (obj == nullptr) {
        DataSet ds = file->openDataSet(name);

        DataSpace dspace_class(H5S_SCALAR);
        StrType dtype(PredType::C_S1, cls.length());
        Attribute att_class = ds.createAttribute(
            "MATLAB_class", dtype, dspace_class
        );
        att_class.write(dtype, cls);
        ds.close();
    } else {
        DataSpace dspace_class(H5S_SCALAR);
        StrType dtype(PredType::C_S1, cls.length());
        Attribute att_class = obj->createAttribute(
            "MATLAB_class", dtype, dspace_class
        );
        att_class.write(dtype, cls);
    }
}

/**
 * Set the "MATLAB_int_decode" attribute of a given dataset.
 *
 * name: Name of dataset.
 * obj:  H5Object to apply the attribute to. If nullptr,
 *       opens the dataset.
 */
void SFile_MAT::WriteMATLAB_int_decode(const string& name, H5Object *obj, const int val) {
    if (obj == nullptr) {
        DataSet ds = file->openDataSet(name);

        DataSpace dspace_int_decode(H5S_SCALAR);
        Attribute att_int_decode = ds.createAttribute(
            "MATLAB_int_decode", PredType::STD_I32LE, dspace_int_decode
        );
        att_int_decode.write(PredType::NATIVE_INT, &val);

        ds.close();
    } else {
        DataSpace dspace_int_decode(H5S_SCALAR);
        Attribute att_int_decode = obj->createAttribute(
            "MATLAB_int_decode", PredType::STD_I32LE, dspace_int_decode
        );
        int val = 2;
        att_int_decode.write(PredType::NATIVE_INT, &val);
    }
}

