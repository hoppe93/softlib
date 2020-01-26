
set(core_tests
	"${PROJECT_SOURCE_DIR}/tests/runtest.cpp"
	"${PROJECT_SOURCE_DIR}/tests/configuration/configuration_test.cpp"
	"${PROJECT_SOURCE_DIR}/tests/sfile/sfile_test.cpp"
	"${PROJECT_SOURCE_DIR}/tests/sfile/sdt_test.cpp"
)
set(distribution_tests
    "${PROJECT_SOURCE_DIR}/tests/distribution/Test_AnalyticalAvalanche.cpp"
    "${PROJECT_SOURCE_DIR}/tests/distribution/Test_CODEDistributionFunction.cpp"
    "${PROJECT_SOURCE_DIR}/tests/distribution/Test_DistributionFunction.cpp"
    "${PROJECT_SOURCE_DIR}/tests/distribution/Test_LUKEDistributionFunction.cpp"
    "${PROJECT_SOURCE_DIR}/tests/distribution/Test_NumericMomentumSpaceDistribution.cpp"
    "${PROJECT_SOURCE_DIR}/tests/distribution/Test_NumericDistributionFunction.cpp"
    "${PROJECT_SOURCE_DIR}/tests/distribution/Test_SOFTDistributionFunction.cpp"
)
set(hdf5_test "${PROJECT_SOURCE_DIR}/tests/sfile/hdf5_test.cpp" "${PROJECT_SOURCE_DIR}/tests/sfile/matlab_test.cpp")
set(integrator_tests
	"${PROJECT_SOURCE_DIR}/tests/integrator/integrator.cpp"
	"${PROJECT_SOURCE_DIR}/tests/integrator/rkdp45.cpp"
)
set(magnetic_test
	"${PROJECT_SOURCE_DIR}/tests/magnetic/analytical2d.cpp"
	"${PROJECT_SOURCE_DIR}/tests/magnetic/domain.cpp"
	"${PROJECT_SOURCE_DIR}/tests/magnetic/luke.cpp"
	"${PROJECT_SOURCE_DIR}/tests/magnetic/magfield_points.cpp"
	"${PROJECT_SOURCE_DIR}/tests/magnetic/magnetic.cpp"
	"${PROJECT_SOURCE_DIR}/tests/magnetic/numeric2d.cpp"
)
set(matlab_test "${PROJECT_SOURCE_DIR}/tests/sfile/matlab_test.cpp")

if (OFFICIAL_MATLAB)
	set(core_tests ${core_tests} ${matlab_test})
else (OFFICIAL_MATLAB)
	set(core_tests ${core_tests} ${hdf5_test})
endif (OFFICIAL_MATLAB)

set(target_tests ${core_tests} ${distribution_tests} ${integrator_tests} ${magnetic_test})
add_executable(softlib_tests ${target_tests})

target_include_directories(softlib_tests PUBLIC "${PROJECT_BINARY_DIR}/include" "${PROJECT_SOURCE_DIR}/tests" "${PROJECT_SOURCE_DIR}/include")
target_link_libraries(softlib_tests PUBLIC softlib)

# Require C++17
set_target_properties(softlib_tests PROPERTIES CXX_STANDARD 17 CXX_STANDARD_REQUIRED YES)

