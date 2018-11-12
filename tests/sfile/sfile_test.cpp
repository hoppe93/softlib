
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
bool sfile_compareArray(double **arr1, double **arr2, sfilesize_t *size1, sfilesize_t *size2) {
	if (size1[0] != size2[0] || size1[1] != size2[1]) return false;

	sfilesize_t i, j;
	for (i = 0; i < size1[0]; i++)
		for (j = 0; j < size1[1]; j++)
			if (fabs(arr1[i][j] - arr2[i][j]) > 5.0*(numeric_limits<double>::epsilon()))
				return false;
	
	return true;
}
/**
 * Check if two lists are equal.
 */
bool sfile_compareLists(double *list1, double *list2, sfilesize_t l1, sfilesize_t l2) {
	if (l1 != l2) return false;

	sfilesize_t i;
	for (i = 0; i < l1; i++)
		if (fabs(list1[i] - list2[i]) > 5.0*(numeric_limits<double>::epsilon()))
			return false;
	
	return true;
}
/**
 * Do a file test on the given SFile object.
 * Write output to the file named 'testname'.
 */
bool sfile_test(SFile *sf, const string& testname) {
	sfilesize_t array_rows=3, array_cols=4, i, j;
	double **array = new double*[array_rows];
	array[0] = new double[array_cols*array_rows];
	for (i = 1; i < array_rows; i++)
		array[i] = array[i-1] + array_cols;
	for (i = 0; i < array_rows; i++)
		for (j = 0; j < array_cols; j++)
			array[i][j] = (double)(i*array_cols + j);

	double att_scalar = 3.14159;
	string att_string = "sfile_test_string1";
	sfilesize_t image_size=3;
	double **image = new double*[image_size];
	image[0] = new double[image_size*image_size];
	for (i = 1; i < image_size; i++)
		image[i] = image[i-1] + image_size;
	for (i = 0; i < image_size; i++)
		for (j = 0; j < image_size; j++)
			image[i][j] = (double)(i*image_size + j);

	double list[] = {1.0,2.0,3.0,4.0,5.0,6.0,7.0,8.0,9.0}; sfilesize_t list_size=9;
	string sfile_string = "sfile_test_string2";

	double **arrbuf, *listbuf, scalbuf, scalar = 10;
	string *strbuf;
	sfilesize_t lenbuf1, lenbuf2[2];
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

	// Read scalar attribute
	scalbuf = sf->GetAttributeScalar("array", "att_scalar");
	if (fabs(scalbuf - att_scalar) > 5.0*(numeric_limits<double>::epsilon()))
		throw SOFTLibException("Reading/writing scalar attributes did not work.");

	// Read string attribute
	strbuf = sf->GetAttributeString("array", "att_string");
	if (*strbuf != att_string)
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

    scalbuf = sf->GetScalar("scalar");
    if (scalbuf != scalar)
        throw SOFTLibException("Reading/writing scalar did not work.");

	// Read string
	strbuf = sf->GetString("string");
	if (*strbuf != sfile_string)
		throw SOFTLibException("Reading/writing strings did not work.");

    sf->Close();

    // Remove test file
    //remove(testname.c_str());
	
	return success;
}

