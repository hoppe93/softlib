include_directories("${PROJECT_SOURCE_DIR}/include")
include_directories("${PROJECT_SOURCE_DIR}/tests")

set(CMAKE_CXX_STANDARD 17)

option(BUILD_TESTS "Build unit tests" OFF)
option(DEBUG "Compile with debug symbols and no optimizations" OFF)
option(INTERP_SPLINTER "Use SPLINTER for interpolation" OFF)
option(PRECISION_DOUBLE "Use double precision floating-point arithmetic" ON)
option(OFFICIAL_MATLAB "Link against official Matlab libraries (disables HDF5 support)" OFF)

# Dependent options
#CMAKE_DEPENDENT_OPTION(INTERP_GSL "Use GSL for interpolation" ON
#					   "INTERP_SPLINTER" OFF)
#CMAKE_DEPENDENT_OPTION(PRECISION_SINGLE "Use single precision floating-point arithmetic" OFF
#					   "PRECISION_DOUBLE" ON)

# Set floating point precision
if (PRECISION_SINGLE AND PRECISION_DOUBLE)
	message(FATAL_ERROR "Can only compile with either single or double precision, not both simultaneously.")
elseif (PRECISION_SINGLE)
	set(REALTYPE "float")
	set(REALEPSILON "FLT_EPSILON")
elseif (PRECISION_DOUBLE)
	set(REALTYPE "double")
	set(REALEPSILON "DBL_EPSILON")
else (PRECISION_SINGLE AND PRECISION_DOUBLE)
	message(FATAL_ERROR "Both single and double precision cannot be turned off.")
endif (PRECISION_SINGLE AND PRECISION_DOUBLE)

# Set RK solver default tolerance
set(RK_TOLERANCE "1e-6")

set(core
	"${PROJECT_SOURCE_DIR}/src/Configuration/ConfigBlock.cpp"
	"${PROJECT_SOURCE_DIR}/src/Configuration/ConfigToken.cpp"
	"${PROJECT_SOURCE_DIR}/src/Configuration/Configuration.cpp"
	"${PROJECT_SOURCE_DIR}/src/Configuration/Configuration.interpreter.cpp"
	"${PROJECT_SOURCE_DIR}/src/Configuration/Configuration.lexer.cpp"
	"${PROJECT_SOURCE_DIR}/src/Configuration/Setting.cpp"
	"${PROJECT_SOURCE_DIR}/src/General/SOFTLibException.cpp"
	"${PROJECT_SOURCE_DIR}/src/General/ProgressTracker.cpp"
	"${PROJECT_SOURCE_DIR}/src/General/Timer.cpp"
	"${PROJECT_SOURCE_DIR}/src/General/Vector.cpp"
	"${PROJECT_SOURCE_DIR}/src/SFile/SFile.cpp"
)
set(core_tests
	"${PROJECT_SOURCE_DIR}/tests/runtest.cpp"
	"${PROJECT_SOURCE_DIR}/tests/configuration/configuration_test.cpp"
	"${PROJECT_SOURCE_DIR}/tests/sfile/sfile_test.cpp"
)
set(distribution
    "${PROJECT_SOURCE_DIR}/src/DistributionFunction/AnalyticalAvalanche.cpp"
    "${PROJECT_SOURCE_DIR}/src/DistributionFunction/CODEDistributionFunction.cpp"
    "${PROJECT_SOURCE_DIR}/src/DistributionFunction/LinearRadialProfile.cpp"
    "${PROJECT_SOURCE_DIR}/src/DistributionFunction/MomentumSpaceDistributionFunction.cpp"
    "${PROJECT_SOURCE_DIR}/src/DistributionFunction/NumericDistributionFunction.cpp"
    "${PROJECT_SOURCE_DIR}/src/DistributionFunction/NumericMomentumSpaceDistributionFunction.cpp"
    "${PROJECT_SOURCE_DIR}/src/DistributionFunction/PowerRadialProfile.cpp"
    "${PROJECT_SOURCE_DIR}/src/DistributionFunction/RadialDistributionFunction.cpp"
    "${PROJECT_SOURCE_DIR}/src/DistributionFunction/SOFTDistributionFunction.cpp"
    "${PROJECT_SOURCE_DIR}/src/DistributionFunction/UnitDistributionFunction.cpp"
)
set(distribution_tests
    "${PROJECT_SOURCE_DIR}/tests/distribution/Test_AnalyticalAvalanche.cpp"
    "${PROJECT_SOURCE_DIR}/tests/distribution/Test_CODEDistributionFunction.cpp"
    "${PROJECT_SOURCE_DIR}/tests/distribution/Test_DistributionFunction.cpp"
    "${PROJECT_SOURCE_DIR}/tests/distribution/Test_NumericMomentumSpaceDistribution.cpp"
    "${PROJECT_SOURCE_DIR}/tests/distribution/Test_SOFTDistributionFunction.cpp"
)
set(hdf5 "${PROJECT_SOURCE_DIR}/src/SFile/SFile_HDF5.cpp" "${PROJECT_SOURCE_DIR}/src/SFile/SFile_MAT.HDF5.cpp")
set(hdf5_test "${PROJECT_SOURCE_DIR}/tests/sfile/hdf5_test.cpp" "${PROJECT_SOURCE_DIR}/tests/sfile/matlab_test.cpp")
set(integration
	"${PROJECT_SOURCE_DIR}/src/Integration/PeakedIntegration.cpp"
)
#set(integrator
#	"${PROJECT_SOURCE_DIR}/src/Integrator/GeneralEquation.cpp"
#	"${PROJECT_SOURCE_DIR}/src/Integrator/Integrator.cpp"
#	"${PROJECT_SOURCE_DIR}/src/Integrator/RKDP45.cpp"
#	"${PROJECT_SOURCE_DIR}/src/Integrator/RungeKutta.cpp"
#)
set(integrator_tests
	"${PROJECT_SOURCE_DIR}/tests/integrator/integrator.cpp"
	"${PROJECT_SOURCE_DIR}/tests/integrator/rkdp45.cpp"
)
set(magnetic
	"${PROJECT_SOURCE_DIR}/src/MagneticField/MagneticField2D.cpp"
	"${PROJECT_SOURCE_DIR}/src/MagneticField/MagneticFieldNumeric2D.cpp"
	"${PROJECT_SOURCE_DIR}/src/MagneticField/MagneticFieldAnalytical2D.cpp"
	"${PROJECT_SOURCE_DIR}/src/MagneticField/MagneticFieldAnalytical2D.q.cpp"
)
set(magnetic_test
	"${PROJECT_SOURCE_DIR}/tests/magnetic/analytical2d.cpp"
	"${PROJECT_SOURCE_DIR}/tests/magnetic/domain.cpp"
	"${PROJECT_SOURCE_DIR}/tests/magnetic/magfield_points.cpp"
	"${PROJECT_SOURCE_DIR}/tests/magnetic/magnetic.cpp"
	"${PROJECT_SOURCE_DIR}/tests/magnetic/numeric2d.cpp"
)
set(matlab "${PROJECT_SOURCE_DIR}/src/SFile/SFile_MAT.Matlab.cpp")
set(matlab_test "${PROJECT_SOURCE_DIR}/tests/sfile/matlab_test.cpp")

if (DEBUG)
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -Wno-misleading-indentation -fopenmp -O0 -g -pg -D_FILE_OFFSET_BITS=64")
else (DEBUG)
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -Wno-misleading-indentation -fopenmp -O3 -march=native -D_FILE_OFFSET_BITS=64")
endif (DEBUG)

if (OFFICIAL_MATLAB)
	set(core ${core} ${matlab})
	set(core_tests ${core_tests} ${matlab_test})
else (OFFICIAL_MATLAB)
	set(core ${core} ${hdf5})
	set(core_tests ${core_tests} ${hdf5_test})
endif (OFFICIAL_MATLAB)

set(target ${core} ${distribution} ${integration} ${magnetic})
add_library(softlib ${target})

if (BUILD_TESTS)
	set(target_tests ${core_tests} ${distribution_tests} ${integrator_tests} ${magnetic_test})
	add_executable(softlib_tests ${target} ${target_tests})
endif (BUILD_TESTS)

# Find GSL!
find_package(GSL REQUIRED)
if (GSL_FOUND)
	if ((GSL_VERSION VERSION_GREATER "2.0") OR (GSL_VERSION VERSION_EQUAL "2.0"))
		include_directories(${GSL_INCLUDE_DIRS})
		target_link_libraries(softlib ${GSL_LIBRARIES})
		target_link_libraries(softlib ${GSL_CBLAS_LIBRARY})

		if (BUILD_TESTS)
			target_link_libraries(softlib_tests ${GSL_LIBRARIES})
			target_link_libraries(softlib_tests ${GSL_CBLAS_LIBRARY})
		endif (BUILD_TESTS)
	else ((GSL_VERSION VERSION_GREATER "2.0") OR (GSL_VERSION VERSION_EQUAL "2.0"))
		message(FATAL_ERROR "softlib needs GSL >= 2.0")
	endif ((GSL_VERSION VERSION_GREATER "2.0") OR (GSL_VERSION VERSION_EQUAL "2.0"))
endif (GSL_FOUND)

# Compile with official Matlab OR HDF5 support?
if (OFFICIAL_MATLAB)
	find_package(Matlab COMPONENTS MAT_LIBRARY MX_LIBRARY)
	if (Matlab_FOUND)
		include_directories(${Matlab_INCLUDE_DIRS})
		target_link_libraries(softlib ${Matlab_MAT_LIBRARY} ${Matlab_MX_LIBRARY})

		if (BUILD_TESTS)
			target_link_libraries(softlib_tests ${Matlab_MAT_LIBRARY} ${Matlab_MX_LIBRARY})
		endif (BUILD_TESTS)
	else (Matlab_FOUND)
		message(FATAL_ERROR "Official MATLAB support was requested, but no MATLAB installation was found")
	endif (Matlab_FOUND)
else (OFFICIAL_MATLAB)
	find_package(HDF5 COMPONENTS CXX)
	if (HDF5_FOUND)
		include_directories(${HDF5_INCLUDE_DIRS})
		target_link_libraries(softlib ${HDF5_LIBRARIES})

		if (BUILD_TESTS)
			target_link_libraries(softlib_tests ${HDF5_LIBRARIES})
		endif (BUILD_TESTS)
	else (HDF5_FOUND)
		message(FATAL_ERROR "HDF5 support was requested, but no HDF5 library was found on this system")
	endif (HDF5_FOUND)
endif (OFFICIAL_MATLAB)

configure_file(
	"${PROJECT_SOURCE_DIR}/include/softlib/config.h.in"
	"${PROJECT_SOURCE_DIR}/include/softlib/config.h"
)

