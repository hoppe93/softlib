#ifndef _CONFIGURATION_TEST_H
#define _CONFIGURATION_TEST_H

#include <runtest.h>

class Test_Configuration : public UnitTest {
	public:
		Test_Configuration(const string&);
		bool Run();
};

#endif/*_CONFIGURATION_TEST_H*/
