# clean cpu version
currentPath=$(pwd)
echo $currentPath
clean_build (){
    currentPath=$1
    releasePath=$2
    if [ -d ${currentPath}/$releasePath/ ]; then
        cd ${currentPath}/$releasePath/
        # if empty folder, do nothing
        if [ "$(ls -A $currentPath/$releasePath)" ]; then
            echo "cleaning $currentPath/$releasePath"
        else
            echo "$currentPath/$releasePath is empty. Skip cleaning."
            return
        fi
        cmake --build ../$releasePath  --target clean
        # rm CMakeCache.txt
        rm -rf ../build 
        rm -rf ../$releasePath/*
    fi
}

# clean gpu version
for device in cpu gpu; do 
    for build_type in release debug; do 
        releasePath=cmake-build-${build_type}-${device}
        echo $releasePath
        clean_build ${currentPath} ${releasePath}
    done

done
