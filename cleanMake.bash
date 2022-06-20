releasePath=$(pwd)
echo $releasePath
cd ${releasePath}/cmake-build-release/
cmake --build ../cmake-build-release  --target clean
# rm CMakeCache.txt
rm -rf ../build 
rm -rf ../cmake-build-release/*
