#ifndef _SFILE_HDF5_TEST_H
#define _SFILE_HDF5_TEST_H

#include <runtest.h>
#include <softlib/SFile.h>

class Test_SFile_HDF5 : public UnitTest {
	public:
		Test_SFile_HDF5(const string&);
		bool Run();
};

#endif/*_SFILE_HDF5_TEST_H*/
