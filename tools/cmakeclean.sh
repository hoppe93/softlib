#!/usr/bin/env bash

pwd=${PWD##*/}

if [[ $pwd == "build" ]]; then
    rm CMakeCache.txt
    rm -rf CMakeFiles
    rm cmake_install.cmake
    rm Makefile
    rm -rf src
else
    echo "ERROR: This script must be run from the 'build' directory."
fi
