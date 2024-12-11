cmake_minimum_required(VERSION 3.10..3.14)

if (${CMAKE_VERSION} VERSION_LESS 3.13)
    cmake_policy(VERSION ${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION})
else()
    cmake_policy(VERSION 3.13)
endif()

# Declare project
project(softlib VERSION 1.0
                DESCRIPTION "Library in support of the SOFT2 code"
                LANGUAGES CXX)

option(BUILD_TESTS "Build unit tests" OFF)
option(DEBUG "Compile with debug symbols and no optimizations" OFF)
option(INTERP_SPLINTER "Use SPLINTER for interpolation" OFF)
option(PRECISION_DOUBLE "Use double precision floating-point arithmetic" ON)
option(OFFICIAL_MATLAB "Link against official Matlab libraries (disables HDF5 support)" OFF)
option(POS_INDP_CODE "Enable position-independent code generation (decreases performance but allows dynamic linking)" OFF)

configure_file(
	"${PROJECT_SOURCE_DIR}/include/softlib/config.h.in"
	"${PROJECT_BINARY_DIR}/include/softlib/config.h"
)

# Require out-of-source builds
file(TO_CMAKE_PATH "${PROJECT_BINARY_DIR}/CMakeLists.txt" LOC_PATH)
if (EXISTS "${LOC_PATH}")
    message(FATAL_ERROR "You cannot build in a source directory")
endif()

add_subdirectory(src)
if (BUILD_TESTS)
    add_subdirectory(tests)
endif()

#subdirs(src)

