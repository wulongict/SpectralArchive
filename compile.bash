#!/bin/bash
set -e


if (( $# < 1 )); then
    echo "compile CPU version"
    with_GPU=FALSE
else 
    with_GPU=$1
fi

echo "compile with GPU: ${with_GPU}"

# to continue.
releasePath=cmake-build-release-cpu
if [ ${with_GPU} = "TRUE" ]; then
    echo "compile GPU version";
    releasePath=cmake-build-release-gpu
elif [ ${with_GPU} = "FALSE" ]; then 
    echo "compile CPU version";
    releasePath=cmake-build-release-cpu
else 
    echo "invalid option of WITH_GPU. It must be TRUE or FALSE. "
fi

currentPath=$(pwd)
echo $currentPath

mkdir -p ${releasePath}
cd ${currentPath}/${releasePath}/ 
rm -f CMakeCache.txt

# echo which python `which python`
# add -pg option to CXX_FLAGS, LINKER_FLAGS AND SHARED_LINKER_FLAGS.
# THANKS: https://stackoverflow.com/questions/26491948/how-to-use-gprof-with-cmake
# -DCMAKE_CXX_FLAGS=-pg -DCMAKE_EXE_LINKER_FLAGS=-pg -DCMAKE_SHARED_LINKER_FLAGS=-pg
#--debug-output -DCMAKE_VERBOSE_MAKEFILE:BOOL=ON

# -DCMAKE_INSTALL_PREFIX=${currentPath}/build
cmake -DWITH_GPU=${with_GPU}  ..
#/usr/local/spectralarchive ..
# cmake  --graphviz=foo.dot ..
# cmake ..


cmake  --build ../${releasePath}   --target  boost  -- -j `nproc`
cmake  --build ../${releasePath}  --target spectroscape msmstest -- -j `nproc`

# run tests added into CMakeLists.txt file with add_test(Name XXX Command YYY)
ctest --verbose

cmake --install . --prefix ${currentPath}/build 


cpack --config CPackConfig.cmake
cpack -D "CPACK_PACKAGING_INSTALL_PREFIX=" --config CPackSourceConfig.cmake


