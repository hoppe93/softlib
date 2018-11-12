#ifndef _SFILE_MATLAB_TEST_H
#define _SFILE_MATLAB_TEST_H

#include <runtest.h>
#include <softlib/SFile.h>

class Test_SFile_MAT : public UnitTest {
	public:
		Test_SFile_MAT(const string&);
		bool Run();
};

#endif/*_SFILE_MATLAB_TEST_H*/
