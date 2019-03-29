#ifndef _MAGNETIC_FIELD_NUMERIC2D_TEST_H
#define _MAGNETIC_FIELD_NUMERIC2D_TEST_H

#include <softlib/MagneticField/MagneticFieldAnalytical2D.h>
#include <softlib/MagneticField/MagneticFieldNumeric2D.h>
#include <runtest.h>
#include "magfield_points.h"
#include "magnetic.h"

class Test_MagneticFieldNumeric2D : public Test_MagneticField {
    private:
        const slibreal_t
            MAXERROR = 3e-3;
        const unsigned int
            NPOINTS = 1000;
	public:
		Test_MagneticFieldNumeric2D(const string&);
		bool CheckMagneticFieldSave(MagneticFieldNumeric2D*, const slibreal_t);
		MagneticFieldNumeric2D *GenerateMF(MagneticFieldAnalytical2D*, slibreal_t, slibreal_t, unsigned int, unsigned int, unsigned int);
        bool CheckNumericDerivatives(MagneticFieldAnalytical2D*, MagneticFieldNumeric2D*);
		bool Run();
};

#endif/*_MAGNETIC_FIELD_NUMERIC2D_TEST_H*/
