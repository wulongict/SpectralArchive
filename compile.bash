#!/usr/bin/bash

releasePath=$(pwd)
echo $releasePath
mkdir cmake-build-release
cd ${releasePath}/cmake-build-release/ 
 rm CMakeCache.txt

#source ~/anaconda3/etc/profile.d/conda.sh
#conda deactivate

# echo which python `which python`
# add -pg option to CXX_FLAGS, LINKER_FLAGS AND SHARED_LINKER_FLAGS.
# THANKS: https://stackoverflow.com/questions/26491948/how-to-use-gprof-with-cmake
# -DCMAKE_CXX_FLAGS=-pg -DCMAKE_EXE_LINKER_FLAGS=-pg -DCMAKE_SHARED_LINKER_FLAGS=-pg
cmake  -DCMAKE_INSTALL_PREFIX=${releasePath}/build ..
#/usr/local/spectralarchive ..
# cmake -DCMAKE_INSTALL_PREFIX=/usr/local/spectralarchive --graphviz=foo.dot ..
# cmake ..
cmake  --build ../cmake-build-release   --target  boost  -- -j 30
cmake  --build ../cmake-build-release  --target spectroscape  -- -j 30

cmake --install .

make package -j `nproc`
make package_source  -j `nproc`
cd ../build/bin && ln -s spectroscape fastcgi_similarity.fcgi
