#!/usr/bin/bash

releasePath=$(pwd)
echo $releasePath
mkdir cmake-build-release
cd ${releasePath}/cmake-build-release/ 
# rm CMakeCache.txt

source ~/anaconda3/etc/profile.d/conda.sh
conda deactivate 
#echo which python `which python`
cmake -DCMAKE_INSTALL_PREFIX=${releasePath}/build ..
#/usr/local/spectralarchive ..
# cmake -DCMAKE_INSTALL_PREFIX=/usr/local/spectralarchive --graphviz=foo.dot ..
# cmake ..
cmake  --build ../cmake-build-release   --target  boost  -- -j 30 
cmake  --build ../cmake-build-release   --target  fastcgi_similarity.fcgi  -- -j 30 
# ctest 
cmake --install .  
#--prefix /tools/archive 
#cpack -G ZIP -V 
make package -j `nproc`
make package_source  -j `nproc`