#ifndef _MAGNETIC_FIELD_ANALYTICAL2D_TEST_H
#define _MAGNETIC_FIELD_ANALYTICAL2D_TEST_H

#include <softlib/MagneticField/MagneticFieldAnalytical2D.h>
#include <runtest.h>
#include "magfield_points.h"
#include "magnetic.h"

class Test_MagneticFieldAnalytical2D : public Test_MagneticField {
    private:
        const unsigned NR=200,          // No. points in R for numerical evaluation of magnetic field.
                       NZ=200,          // No. points in Z for numerical evaluation of magnetic field.
                       NRANDOM=200;     // No. random points to test jacobian in.
        const slibreal_t 
            JAC1TOL = 100.0*REAL_EPSILON,
            JAC2TOL = 5e-2,             // This tolerance is quite high, but necessary since we use random sample points
            CONVTOL = 5e-4;
	public:
		Test_MagneticFieldAnalytical2D(const string&);
		bool Run();
        bool TestConversion(MagneticFieldAnalytical2D*, const string&);
        bool TestJacobian1(MagneticFieldAnalytical2D*);
        bool TestJacobian2(MagneticFieldAnalytical2D*);
};

#endif/*_MAGNETIC_FIELD_ANALYTICAL2D_TEST_H*/
