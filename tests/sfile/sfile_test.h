#ifndef _SFILE_TEST_H
#define _SFILE_TEST_H

#include <softlib/SFile.h>
#include <string>

using namespace std;

bool sfile_compareArray(double**,double**,sfilesize_t*,sfilesize_t*);
bool sfile_compareLists(double*, double*, sfilesize_t, sfilesize_t);
bool sfile_test(SFile*, const string&);

#endif/*_SFILE_TEST_H*/
