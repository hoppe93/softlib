#ifndef _SFILE_TEST_H
#define _SFILE_TEST_H

#include <softlib/SFile.h>
#include <string>

using namespace std;

template<typename T>
bool sfile_compareLists(T*, T*, sfilesize_t, sfilesize_t);
template<typename T>
bool sfile_compareArray(T**,T**,sfilesize_t*,sfilesize_t*);

bool sfile_test(SFile*, const string&, const bool multisupp=true, const bool structsupp=true);

#endif/*_SFILE_TEST_H*/
