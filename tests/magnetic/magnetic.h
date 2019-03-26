#ifndef _TEST_MAGNETIC_H
#define _TEST_MAGNETIC_H

#include <softlib/MagneticField/MagneticField2D.h>
#include "magfield_points.h"
#include "runtest.h"

class Test_MagneticField : public UnitTest {
	private:
		slibreal_t threshold;
	public:
		Test_MagneticField(const string&);
		Test_MagneticField(const string&, const slibreal_t);
		bool ComparePoints(const slibreal_t[MAGNETIC_FIELD_TEST_NPOINTS][12], MagneticField2D*, slibreal_t, const string&, bool);
		bool CompareDerivatives(const slibreal_t[MAGNETIC_FIELD_TEST_NPOINTS][12], MagneticField2D*, slibreal_t, const string&);
		void SetThreshold(const slibreal_t);
		virtual bool Run() = 0;
};

#endif/*_TEST_MAGNETIC_H*/
