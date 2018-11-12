#ifndef _TEST_ANALYTICAL_AVALANCHE_H
#define _TEST_ANALYTICAL_AVALANCHE_H

#include "runtest.h"
#include "Test_DistributionFunction.h"

class Test_AnalyticalAvalanche : public Test_DistributionFunction {
    private:
        /* Test configuration */
        const unsigned int
            NEHAT       =10,
            NLOGLAMBDA  =3,
            NZEFF       =4,
            NPOINTS     =20;

        const slibreal_t
            EHAT0       =1.3,
            EHAT1       =15.0,
            LOGLAMBDA0  =14.0,
            LOGLAMBDA1  =18.0,
            ZEFF0       =1.0,
            ZEFF1       =5.0,

            PMAX        =100.0,

            TOLERANCE  =100.0*REAL_EPSILON;
        
        slibreal_t *r=nullptr, *p=nullptr, *xi=nullptr, *f=nullptr;

    public:
        Test_AnalyticalAvalanche(const string& name) : Test_DistributionFunction(name) {}
        ~Test_AnalyticalAvalanche();
        slibreal_t AA(slibreal_t, slibreal_t, slibreal_t, slibreal_t, slibreal_t);
        bool Run();
};

#endif/*_TEST_ANALYTICAL_AVALANCHE_H*/
