
set(core
	"${PROJECT_SOURCE_DIR}/src/Configuration/ConfigBlock.cpp"
	"${PROJECT_SOURCE_DIR}/src/Configuration/Script/ConfigToken.cpp"
	"${PROJECT_SOURCE_DIR}/src/Configuration/Configuration.cpp"
	"${PROJECT_SOURCE_DIR}/src/Configuration/Script/ConfigurationScript.cpp"
	"${PROJECT_SOURCE_DIR}/src/Configuration/Script/Interpreter.cpp"
	"${PROJECT_SOURCE_DIR}/src/Configuration/Script/Lexer.cpp"
	"${PROJECT_SOURCE_DIR}/src/Configuration/Script/SettingScript.cpp"
	"${PROJECT_SOURCE_DIR}/src/Configuration/SFile/SettingSFile.cpp"
	"${PROJECT_SOURCE_DIR}/src/Configuration/Setting.cpp"
	"${PROJECT_SOURCE_DIR}/src/General/SOFTLibException.cpp"
	"${PROJECT_SOURCE_DIR}/src/General/ProgressTracker.cpp"
	"${PROJECT_SOURCE_DIR}/src/General/Timer.cpp"
	#"${PROJECT_SOURCE_DIR}/src/General/Vector.cpp"
	"${PROJECT_SOURCE_DIR}/src/SFile/SFile.cpp"
	"${PROJECT_SOURCE_DIR}/src/SFile/SFile_SDT.cpp"
	"${PROJECT_SOURCE_DIR}/src/SFile/SFile_SDT.load.cpp"
)
set(distribution
    "${PROJECT_SOURCE_DIR}/src/DistributionFunction/AnalyticalAvalanche.cpp"
    "${PROJECT_SOURCE_DIR}/src/DistributionFunction/BesselRadialProfile.cpp"
    "${PROJECT_SOURCE_DIR}/src/DistributionFunction/CODEDistributionFunction.cpp"
    "${PROJECT_SOURCE_DIR}/src/DistributionFunction/ConnorHastieDistribution.cpp"
    "${PROJECT_SOURCE_DIR}/src/DistributionFunction/ExponentialPitch.cpp"
    "${PROJECT_SOURCE_DIR}/src/DistributionFunction/GOCODEDistributionFunction.cpp"
    "${PROJECT_SOURCE_DIR}/src/DistributionFunction/LinearRadialProfile.cpp"
    "${PROJECT_SOURCE_DIR}/src/DistributionFunction/LUKEDistributionFunction.cpp"
    "${PROJECT_SOURCE_DIR}/src/DistributionFunction/MomentumSpaceDistributionFunction.cpp"
    "${PROJECT_SOURCE_DIR}/src/DistributionFunction/NumericDistributionFunction.cpp"
    "${PROJECT_SOURCE_DIR}/src/DistributionFunction/NumericMomentumSpaceDistributionFunction.cpp"
    "${PROJECT_SOURCE_DIR}/src/DistributionFunction/PowerRadialProfile.cpp"
    "${PROJECT_SOURCE_DIR}/src/DistributionFunction/RadialDistributionFunction.cpp"
    "${PROJECT_SOURCE_DIR}/src/DistributionFunction/deprSOFTDistributionFunction.cpp"
    "${PROJECT_SOURCE_DIR}/src/DistributionFunction/SOFTDistributionFunction.cpp"
    "${PROJECT_SOURCE_DIR}/src/DistributionFunction/UnitDistributionFunction.cpp"
)
set(hdf5 "${PROJECT_SOURCE_DIR}/src/SFile/SFile_HDF5.cpp" "${PROJECT_SOURCE_DIR}/src/SFile/SFile_MAT.HDF5.cpp")
set(imagegenerator
    "${PROJECT_SOURCE_DIR}/src/ImageGenerator/ImageGenerator.cpp"
    "${PROJECT_SOURCE_DIR}/src/ImageGenerator/ImageGeneratorPPM.cpp"
)
set(imagegeneratorpng
    "${PROJECT_SOURCE_DIR}/src/ImageGenerator/ImageGeneratorPNG.cpp"
)
set(integration
	"${PROJECT_SOURCE_DIR}/src/Integration/PeakedIntegration.cpp"
)
set(magnetic
	"${PROJECT_SOURCE_DIR}/src/MagneticField/MagneticField2D.cpp"
	"${PROJECT_SOURCE_DIR}/src/MagneticField/MagneticFieldNumeric2D.cpp"
	"${PROJECT_SOURCE_DIR}/src/MagneticField/MagneticFieldAnalytical2D.cpp"
	"${PROJECT_SOURCE_DIR}/src/MagneticField/MagneticFieldAnalytical2D.q.cpp"
	"${PROJECT_SOURCE_DIR}/src/MagneticField/MagneticFieldLUKE.cpp"
)
set(matlab "${PROJECT_SOURCE_DIR}/src/SFile/SFile_MAT.Matlab.cpp")

#if (DEBUG)
#	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -Wno-misleading-indentation -fopenmp -O0 -g -pg -D_FILE_OFFSET_BITS=64")
#else (DEBUG)
#	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -Wno-misleading-indentation -fopenmp -O3 -march=native -D_FILE_OFFSET_BITS=64")
#endif (DEBUG)

if (OFFICIAL_MATLAB)
	set(core ${core} ${matlab})
else (OFFICIAL_MATLAB)
	set(core ${core} ${hdf5})
endif (OFFICIAL_MATLAB)

set(target ${core} ${distribution} ${imagegenerator} ${integration} ${magnetic})

# Check for libpng (if it exists => add PNG image generation support)
find_package(PNG)
if (PNG_FOUND)
    set(target ${target} ${imagegeneratorpng})
endif (PNG_FOUND)

# Define library target
add_library(softlib STATIC ${target})
target_include_directories(softlib PUBLIC "${PROJECT_BINARY_DIR}/include" "${PROJECT_SOURCE_DIR}/include")

# Ensure we are able to open files that are larger than 4 GiB
target_compile_definitions(softlib PUBLIC _FILE_OFFSET_BITS=64)

# Require C++17
set_target_properties(softlib PROPERTIES CXX_STANDARD 17 CXX_STANDARD_REQUIRED YES)

# Set floating point precision
if (PRECISION_SINGLE AND PRECISION_DOUBLE)
	message(FATAL_ERROR "Can only compile with either single or double precision, not both simultaneously.")
endif (PRECISION_SINGLE AND PRECISION_DOUBLE)

########################
# MODULES
########################
# Math library
find_library(MATH_LIBRARY m)
if (MATH_LIBRARY)
    target_link_libraries(softlib PUBLIC ${MATH_LIBRARY})
endif()

# Interprocedural optimizations
#if (NOT CMAKE_BUILD_TYPE STREQUAL "Debug")
#    include(CheckIPOSupport)
#    check_ipo_supported(RESULT result)
#    if (result)
#        set_target_properties(softlib PROPERTIES INTERPROCEDURAL_OPTIMIZATION TRUE)
#    endif()
#endif()

# Find OpenMP
find_package(OpenMP)
if (OpenMP_CXX_FOUND)
    target_link_libraries(softlib PUBLIC OpenMP::OpenMP_CXX)
endif()

# Find GSL!
find_package(GSL REQUIRED)
if (GSL_FOUND)
	if ((GSL_VERSION VERSION_GREATER "2.0") OR (GSL_VERSION VERSION_EQUAL "2.0"))
		target_include_directories(softlib PUBLIC "${GSL_INCLUDE_DIRS}")
		target_link_libraries(softlib PUBLIC "${GSL_LIBRARIES}")
		target_link_libraries(softlib PUBLIC "${GSL_CBLAS_LIBRARY}")
	else ((GSL_VERSION VERSION_GREATER "2.0") OR (GSL_VERSION VERSION_EQUAL "2.0"))
		message(FATAL_ERROR "softlib needs GSL >= 2.0")
	endif ((GSL_VERSION VERSION_GREATER "2.0") OR (GSL_VERSION VERSION_EQUAL "2.0"))
endif (GSL_FOUND)

# Find libpng (optional) -- checked earlier above
#find_package(PNG QUIET)
if (PNG_FOUND)
    target_include_directories(softlib PUBLIC "${PNG_INCLUDE_DIRS}")
    target_link_libraries(softlib PUBLIC "${PNG_LIBRARIES}")
else (PNG_FOUND)
    message(INFO "No libpng was found. Configuring without PNG image support...")
endif (PNG_FOUND)

# Compile with official Matlab OR HDF5 support?
if (OFFICIAL_MATLAB)
	find_package(Matlab 9.4 COMPONENTS MAT_LIBRARY MX_LIBRARY)
	if (Matlab_FOUND)
		target_include_directories(softlib PUBLIC "${Matlab_INCLUDE_DIRS}")
		target_link_libraries(softlib PUBLIC "${Matlab_MAT_LIBRARY}" "${Matlab_MX_LIBRARY}")
	else (Matlab_FOUND)
		message(FATAL_ERROR "Official MATLAB support was requested, but no MATLAB installation was found")
	endif (Matlab_FOUND)
else (OFFICIAL_MATLAB)
	find_package(HDF5 COMPONENTS CXX)
	if (HDF5_FOUND)
		target_include_directories(softlib PUBLIC "${HDF5_INCLUDE_DIRS}")
		target_link_libraries(softlib PUBLIC "${HDF5_LIBRARIES}")
	else (HDF5_FOUND)
		message(FATAL_ERROR "HDF5 support was requested, but no HDF5 library was found on this system")
	endif (HDF5_FOUND)
endif (OFFICIAL_MATLAB)

