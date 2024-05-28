#!/bin/bash
set -e

show_help() {
    echo "Usage: $0 [option...] " 
    echo "A script for compile spectroscape with CMake"
    echo ""
    echo "Options:"
    echo "  -g, --gpu           compile GPU version, otherwise, CPU version; Default to CPU."
    echo "  -b, --build-type    compile release version, otherwise, debug version, Default to release. "
    echo "  -c, --cores         number of cores to use for compilation. Default to nproc / 16 + 1"
    echo "  -h, --help          display this help and exit"
    echo ""
    echo "Examples:"
    echo "  $0 --gpu TRUE --build-type release"
    echo "  $0 --gpu FALSE --build-type debug"
    echo "  $0 --gpu TRUE --build-type release --cores 8"
    exit 1
}

is_docker(){
    if [ -f /.dockerenv ] || [ -f /run/.containerenv ]; then
        echo "Running inside Docker"
        # Adjust script behavior for Docker
    else
        echo "Running outside Docker"
        # Normal script behavior
    fi
}

is_docker

POSITIONAL_ARGS=()
with_GPU=FALSE
build_type=release
cores=`nproc`
path_suffix="cpu"
while [[ $# -gt 0 ]]; do
  case $1 in
    -g|--gpu)
      with_GPU="$2"
      if [ ${with_GPU} = "TRUE" ]; then
          echo "compile GPU version";
          path_suffix="gpu"
      elif [ ${with_GPU} = "FALSE" ]; then 
          echo "compile CPU version";
          path_suffix="cpu"
      else 
          echo "invalid option of with_GPU. Default to CPU. "
          with_GPU=FALSE
          path_suffix="cpu"
      fi
      shift # past argument
      shift # past value
      ;;
    -b|--build-type)
      build_type="$2"
        if [ ${build_type} = "release" ]; then
            echo "compile release version";
        elif [ ${build_type} = "debug" ]; then
            echo "compile debug version";
        else 
            echo "invalid option of build_type. Default to release. "
            build_type=release
        fi
      shift # past argument
      shift # past value
      ;;
    -c|--cores)
      cores="$2"
      # limit core number to 1 to 128
        if ((cores < 1)); then
            cores=1
        elif ((cores > 128)); then
            cores=128
        fi
      shift # past argument
      shift # past value
      ;;
    -*|--*)
      echo "Unknown option $1"
      # show help
      show_help
      exit 1
      ;;
    *)
      POSITIONAL_ARGS+=("$1") # save positional arg
      shift # past argument
      ;;
  esac
done

set -- "${POSITIONAL_ARGS[@]}" # restore positional parameters

echo "............................. "
echo "with_GPU  = ${with_GPU}"
echo "build_type= ${build_type}"
echo "cores     = ${cores}"
echo "POSITIONAL_ARGS = ${POSITIONAL_ARGS}"
echo "............................."


releasePath=cmake-build-${build_type}-${path_suffix}
echo "releasePath = ${releasePath}"


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

cmake  --build ../${releasePath}  --target spectroscape msmstest using_all_cpu_cores -- -j $cores

# run tests added into CMakeLists.txt file with add_test(Name XXX Command YYY)
ctest --verbose

cmake --install . --prefix ${currentPath}/build 


cpack --config CPackConfig.cmake
cpack -D "CPACK_PACKAGING_INSTALL_PREFIX=" --config CPackSourceConfig.cmake


