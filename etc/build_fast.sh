#!/bin/bash
# > ./etc/build.sh

# exit when any command fails
set -e

# Clean build directory
cd build

# Build LEAP
cmake ..
cmake --build . -j8
