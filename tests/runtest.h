#ifndef _RUNTEST_H
#define _RUNTEST_H

#include <string>
#include <softlib/config.h>

using namespace std;

class UnitTest {
	protected:
		string name;
	public:
		UnitTest(const string&);
		string& GetName();
		bool HasName(const string&);
		void PrintError(const string&, ...);
		void PrintOK(const string&, ...);
		void PrintStatus(const string&, ...);
		void PrintWarning(const string&, ...);
		virtual bool Run() = 0;

        void InitRand();
        slibreal_t Rand();
        slibreal_t Rand(slibreal_t, slibreal_t);
};

#endif/*_RUNTEST_H*/
