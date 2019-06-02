
#include <cmath>
#include <cstdio>
#include <iostream>
#include <limits>
#include <string>
#include <softlib/SFile.h>
#include <softlib/SOFTLibException.h>

using namespace std;

/**
 * Check if two arrays are equal.
 */
template<typename T>
bool sfile_compareArray(T **arr1, T **arr2, sfilesize_t *size1, sfilesize_t *size2) {
	if (size1[0] != size2[0] || size1[1] != size2[1]) return false;

	sfilesize_t i, j;
	for (i = 0; i < size1[0]; i++)
		for (j = 0; j < size1[1]; j++) {
            if (typeid(T) == typeid(double)) {
                if (fabs(arr1[i][j] - arr2[i][j]) > 5.0*(numeric_limits<double>::epsilon()))
                    return false;
            } else {
                if (arr1[i][j] != arr2[i][j])
                    return false;
            }
        }
	
	return true;
}

/**
 * Check if two 3D arrays are equal.
 */
template<typename T>
bool sfile_compare3Array(T ***arr1, T ***arr2, sfilesize_t *size1, sfilesize_t *size2) {
	if (size1[0] != size2[0] ||
		size1[1] != size2[1] ||
		size1[2] != size2[2])
		return false;
	
	for (sfilesize_t i = 0; i < size1[0]; i++)
		for (sfilesize_t j = 0; j < size1[1]; j++)
			for (sfilesize_t k = 0; k < size1[2]; k++) {
                if (typeid(T) == typeid(double)) {
                    if (fabs(arr1[i][j][k] - arr2[i][j][k]) > 5.0*(numeric_limits<double>::epsilon()))
                        return false;
                } else {
                    if (arr1[i][j][k] != arr2[i][j][k])
                        return false;
                }
            }
	
	return true;
}

/**
 * Check if two lists are equal.
 */
template<typename T>
bool sfile_compareLists(T *list1, T *list2, sfilesize_t l1, sfilesize_t l2) {
	if (l1 != l2) return false;

	sfilesize_t i;
	for (i = 0; i < l1; i++) {
        if (typeid(T) == typeid(double)) {
            if (fabs(list1[i] - list2[i]) > 5.0*(numeric_limits<double>::epsilon()))
                return false;
        } else {
            if (list1[i] != list2[i])
                return false;
        }
    }
	
	return true;
}

template<typename T>
T *sfile_construct_list(sfilesize_t n) {
    T *a = new T[n];
    for (sfilesize_t i = 0; i < n; i++)
        a[i] = i+1;

    return a;
}
template<typename T>
T **sfile_construct_array(sfilesize_t nr, sfilesize_t nc) {
    T **a = new T*[nr];
    a[0] = new T[nr*nc];

    for (sfilesize_t i = 1; i < nr; i++)
        a[i] = a[i-1] + nc;

    for (sfilesize_t i = 0; i < nr; i++)
        for (sfilesize_t j = 0; j < nc; j++)
            a[i][j] = (T)(i*nc + j);

    return a;
}

/**
 * Do a file test on the given SFile object.
 * Write output to the file named 'testname'.
 *
 * sf:         SFile object to use for testing.
 * testname:   Name of test (used for generating file name).
 * multisupp:  Run the multi-dimensional array tests
 * structsupp: Run the 'struct'-writing/reading tests.
 */
bool sfile_test(SFile *sf, const string& testname, const bool multisupp, const bool structsupp) {
	sfilesize_t i, j, k;

	// Construct the array
	sfilesize_t array_rows=3, array_cols=4;
    double **array = sfile_construct_array<double>(array_rows, array_cols);
	
	// Construct the multi-dimensional array
	sfilesize_t mularr_ni=3, mularr_nj=4, mularr_nk=5;
	double ***mularr = new double**[mularr_ni];
	sfilesize_t mularr_dims[3] = {mularr_ni, mularr_nj, mularr_nk};
	
	for (i = 0; i < mularr_ni; i++) {
		if (i == 0) {
			mularr[0] = new double*[mularr_ni*mularr_nj];
			mularr[0][0] = new double[mularr_ni*mularr_nj*mularr_nk];
		} else {
			mularr[i] = mularr[i-1] + mularr_nj;
			mularr[i][0] = mularr[i-1][0] + mularr_nj*mularr_nk;
		}

		for (j = 1; j < mularr_nj; j++)
			mularr[i][j] = mularr[i][j-1] + mularr_nk;
	}

	for (i = 0; i < mularr_ni; i++)
		for (j = 0; j < mularr_nj; j++)
			for (k = 0; k < mularr_nk; k++)
				mularr[i][j][k] = (double)((mularr_nj*i + j)*mularr_nk + k);

    // Construct the 1D integer array
    sfilesize_t int64list_n = 5;
    int64_t *int64list = sfile_construct_list<int64_t>(int64list_n);

    // Construct the 2D integer array
    sfilesize_t int64arr_rows=4, int64arr_cols=3;
    int64_t **int64arr = sfile_construct_array<int64_t>(int64arr_rows, int64arr_cols);

	// Construct the scalar
	double att_scalar = 3.14159;
	// Construct the attribute string
	string att_string = "sfile_test_string1";
	// Construct the image
	sfilesize_t image_size=3;
	double **image = new double*[image_size];
	image[0] = new double[image_size*image_size];
	for (i = 1; i < image_size; i++)
		image[i] = image[i-1] + image_size;
	for (i = 0; i < image_size; i++)
		for (j = 0; j < image_size; j++)
			image[i][j] = (double)(i*image_size + j);

	// Construct the list
	double list[] = {1.0,2.0,3.0,4.0,5.0,6.0,7.0,8.0,9.0}; sfilesize_t list_size=9;
	// Construct the regular string
	string sfile_string = "sfile_test_string2";

	double **arrbuf, ***mularrbuf, *listbuf, scalbuf, scalar = 10;
	string strbuf;
	sfilesize_t lenbuf1, lenbuf2[2], mularrlenbuf[3];
	bool success = true;

	/* Write file */
	sf->Open(testname, SFILE_MODE_WRITE);

	sf->WriteArray("array", array, array_rows, array_cols);
	sf->WriteAttribute_scalar("array", "att_scalar", att_scalar);
	sf->WriteAttribute_string("array", "att_string", att_string);
	sf->WriteImage("image", image, image_size);
	sf->WriteList("list", list, list_size);
    sf->WriteScalar("scalar", scalar);
	sf->WriteString("string", sfile_string);

    sf->WriteInt64List("int64list", int64list, int64list_n);
    sf->WriteInt64Array("int64arr", int64arr, int64arr_rows, int64arr_cols);

	if (multisupp)
		sf->WriteMultiArray("multiArray", **mularr, 3, mularr_dims);
	
	if (structsupp) {
		sf->CreateStruct("struct");
		sf->WriteString("struct/struct_string", sfile_string);
	}

	sf->Close();

	/* Read back file */
	sf->Open(testname, SFILE_MODE_READ);

	// Read array
	arrbuf = sf->GetDoubles("array", lenbuf2);
	sfilesize_t tlen[] = {array_rows, array_cols};
	if (!sfile_compareArray(arrbuf, array, lenbuf2, tlen))
		throw SOFTLibException("Reading/writing arrays did not work.");

	delete [] arrbuf[0];
	delete [] arrbuf;

    // Read int64 list
    int64_t *int64bufl = sf->GetIntList("int64list", &lenbuf1);
    if (!sfile_compareLists(int64bufl, int64list, lenbuf1, int64list_n))
        throw SOFTLibException("Reading/writing int64 lists did not work.");

    delete [] int64bufl;

    // Read int64 array
    int64_t **int64buf = sf->GetIntArray("int64arr", lenbuf2);
    tlen[0] = int64arr_rows; tlen[1] = int64arr_cols;
    if (!sfile_compareArray(int64buf, int64arr, lenbuf2, tlen))
        throw SOFTLibException("Reading/writing int64 arrays did not work.");

    delete [] int64buf[0];
    delete [] int64buf;

	// Read scalar attribute
	scalbuf = sf->GetAttributeScalar("array", "att_scalar");
	if (fabs(scalbuf - att_scalar) > 5.0*(numeric_limits<double>::epsilon()))
		throw SOFTLibException("Reading/writing scalar attributes did not work.");

	// Read string attribute
	strbuf = sf->GetAttributeString("array", "att_string");
	if (strbuf != att_string)
		throw SOFTLibException("Reading/writing string attributes did not work.");

	// Read image
	arrbuf = sf->GetDoubles("image", lenbuf2);
	tlen[0] = tlen[1] = image_size;
	if (!sfile_compareArray(arrbuf, image, lenbuf2, tlen))
		throw SOFTLibException("Reading/writing images did not work.");

	delete [] arrbuf[0];
	delete [] arrbuf;

	// Read list
	listbuf = sf->GetList("list", &lenbuf1);
	if (!sfile_compareLists(listbuf, list, lenbuf1, list_size))
		throw SOFTLibException("Reading/writing lists did not work.");

	// Read multi-dimensional array
	if (multisupp) {
		sfilesize_t mularrndims;
		mularrbuf = (double***)sf->GetMultiArray("multiArray", 3, mularrndims, mularrlenbuf);
		if (mularrbuf == nullptr)
			throw SOFTLibException("Reading/writing multi-dimensional arrays did not work. Error: 1.");
		
		if (mularrndims != 3)
			throw SOFTLibException("Reading/writing multi-dimensional arrays did not work. The loaded array had %llu dimensions, not 3.", mularrndims);
		
		if (!sfile_compare3Array(mularrbuf, mularr, mularrlenbuf, mularr_dims))
			throw SOFTLibException("Reading/writing multi-dimensional arrays did not work. Error: 3.");

		delete [] mularrbuf[0][0];
		delete [] mularrbuf[0];
		delete [] mularrbuf;

		delete [] mularr[0][0];
		delete [] mularr[0];
		delete [] mularr;
	}

	// Read scalar
    scalbuf = sf->GetScalar("scalar");
    if (scalbuf != scalar)
        throw SOFTLibException("Reading/writing scalar did not work.");

	// Read string
	strbuf = sf->GetString("string");
	if (strbuf != sfile_string)
		throw SOFTLibException("Reading/writing strings did not work.");

	if (structsupp) {
		strbuf = sf->GetString("struct/struct_string");
		if (strbuf != sfile_string)
			throw SOFTLibException("Reading/writing struct did not work.");
	}

    sf->Close();

    // Remove test file (if successfull)
    if (success)
        remove(testname.c_str());
	
	return success;
}

