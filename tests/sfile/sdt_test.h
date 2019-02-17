#ifndef _SFILE_SDT_TEST_H
#define _SFILE_SDT_TEST_H

#include <runtest.h>
#include <softlib/SFile.h>

class Test_SFile_SDT : public UnitTest {
	public:
		Test_SFile_SDT(const string&);
		bool Run();
};

#endif/*_SFILE_SDT_TEST_H*/
