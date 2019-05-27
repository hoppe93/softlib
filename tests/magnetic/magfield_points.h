/* This file was auto-generated */
#ifndef _MAGFIELD_POINTS_H
#define _MAGFIELD_POINTS_H

#include <softlib/config.h>

#define MAGNETIC_FIELD_TEST_GUARANTEED_PRECISION (max(100*1e-15,20*REAL_EPSILON))

const unsigned int MAGNETIC_FIELD_TEST_NPOINTS=25;
const slibreal_t magnetic_field_test_data_B0      =5,
                 magnetic_field_test_data_Rm      =0.68,
                 magnetic_field_test_data_zm      =0,
                 magnetic_field_test_data_rminor  =0.22,
                 magnetic_field_test_data_qa_const=0.5,
                 magnetic_field_test_data_qa_lin  =2,
                 magnetic_field_test_data_qa_quad =2,
                 magnetic_field_test_data_qa_exp  =0.693147180559945;
extern const slibreal_t magnetic_field_test_data_const[25][12];
extern const slibreal_t magnetic_field_test_data_linear[25][12];
extern const slibreal_t magnetic_field_test_data_quadratic[25][12];
extern const slibreal_t magnetic_field_test_data_exponential[25][12];

#endif/*_MAGFIELD_POINTS_H*/