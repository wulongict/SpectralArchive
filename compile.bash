#!/bin/bash
set -e
echo "Usage: "$0 "[with_GPU] [build_type]"
echo "A script for compile spectroscape with CMake"
echo ""
echo "if with_GPU is TRUE, compile GPU version, otherwise, CPU version; Default to CPU."
echo "if build_type is release, compile release version, otherwise, debug version, Default to release. "
echo ""
echo "............................. "
echo $0 $@
echo "............................. "
echo ""
build_type=release
if (( $# > 1)); then 
    build_type=$2
    if [ ${build_type} = "release" ]; then
        echo "compile release version";
    elif [ ${build_type} = "debug" ]; then
        echo "compile debug version";
    else 
        echo "invalid option of build_type. Default to release. "
    fi
fi

releasePath=cmake-build-${build_type}-cpu
with_GPU=FALSE
if (( $# > 0 )); then
    with_GPU=$1
    if [ ${with_GPU} = "TRUE" ]; then
        echo "compile GPU version";
        releasePath=cmake-build-${build_type}-gpu
    elif [ ${with_GPU} = "FALSE" ]; then 
        echo "compile CPU version";
        releasePath=cmake-build-${build_type}-cpu
    else 
        echo "invalid option of WITH_GPU. Default to CPU. "
    fi
fi

currentPath=$(pwd)
echo $currentPath

mkdir -p ${releasePath}
cd ${currentPath}/${releasePath}/ 
# rm -f CMakeCache.txt

# echo which python `which python`
# add -pg option to CXX_FLAGS, LINKER_FLAGS AND SHARED_LINKER_FLAGS.
# THANKS: https://stackoverflow.com/questions/26491948/how-to-use-gprof-with-cmake
# -DCMAKE_CXX_FLAGS=-pg -DCMAKE_EXE_LINKER_FLAGS=-pg -DCMAKE_SHARED_LINKER_FLAGS=-pg
#--debug-output -DCMAKE_VERBOSE_MAKEFILE:BOOL=ON



# create cmake build options. 
cmake_build_options="-DWITH_GPU=${with_GPU}"

# if build type is release, add CXX flag -g
if [ ${build_type} = "debug" ]; then
    # to use gdb, add -DCMAKE_CXX_FLAGS=-g
    cmake_build_options="${cmake_build_options} -DCMAKE_CXX_FLAGS=-g"
fi
cmake ${cmake_build_options} ..

# -DCMAKE_INSTALL_PREFIX=${currentPath}/build
# to use gdb, add -DCMAKE_CXX_FLAGS=-g
# cmake -DWITH_GPU=${with_GPU} -DCMAKE_CXX_FLAGS=-g ..
#/usr/local/spectralarchive ..
# cmake  --graphviz=foo.dot ..
# cmake ..

# compile libfcgi will multiple cores will crash on macos+docker(ubuntu2204) environment. 
# it turns out that this is a resource issue on macos.
cores=$((`nproc` / 16 + 1)) 
cmake  --build ../${releasePath}  --target spectroscape msmstest using_all_cpu_cores -- -j $cores

# run tests added into CMakeLists.txt file with add_test(Name XXX Command YYY)
ctest --verbose

cmake --install . --prefix ${currentPath}/build 


cpack --config CPackConfig.cmake
cpack -D "CPACK_PACKAGING_INSTALL_PREFIX=" --config CPackSourceConfig.cmake


