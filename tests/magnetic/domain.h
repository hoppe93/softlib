#ifndef _TEST_DOMAIN_H
#define _TEST_DOMAIN_H

#include <softlib/MagneticField/MagneticField2D.h>
#include <softlib/MagneticField/MagneticFieldNumeric2D.h>
#include "runtest.h"

class Test_Domain : public UnitTest {
	public:
		Test_Domain(const string&);

		MagneticField2D *GenerateMF();
		bool Run();
};

#endif/*_TEST_DOMAIN_H*/
