# clean cpu version
currentPath=$(pwd)
echo $currentPath
releasePath=cmake-build-release-cpu
cd ${currentPath}/$releasePath/
cmake --build ../$releasePath  --target clean
# rm CMakeCache.txt
rm -rf ../build 
rm -rf ../$releasePath/*

# clean gpu version
releasePath=cmake-build-release-gpu
cd ${currentPath}/$releasePath/
cmake --build ../$releasePath  --target clean
# rm CMakeCache.txt
rm -rf ../build 
rm -rf ../$releasePath/*