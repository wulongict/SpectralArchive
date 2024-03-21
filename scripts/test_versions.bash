#!/bin/bash

# Test different versions of Spectroscape. 
# The test folder should be prepared as ~/data/spectroscape
# three subfolders: mass_spectra, new_data, and spectral_archives

# Example:
# ./test_versions.bash 1.1.8 1.1.9
# ./test_versions.bash latest

# Prerequisite:
# 1. download the test data from GoogleDrive link:
#    go to
#    https://drive.google.com/drive/folders/1ZaHdkzxclzYboYf_67pPxgkmFxpnVFAG?usp=
#    download the data.zip file
# 2. unzip it to ~/data/spectroscape/
#    Then, there will be three subfolders: mass_spectra, new_data, and spectral_archives
# 3. try run the following example to test the version 1.1.8
#    example: ./test_versions.bash 1.1.8 1.1.9

# 4. The test will be run in the spectral_archives_${version} folder.
#    The test will be run with a clean start, i.e., the spectral_archives_${version} folder will be removed before the test.

function verify_folder_settings(){
    if [ ! -d ~/data/spectroscape ]; then
        echo "The test folder should be prepared as ~/data/spectroscape"
        echo "three subfolders: mass_spectra, new_data, and spectral_archives"
        exit 1
    fi
    if [ ! -d ~/data/spectroscape/mass_spectra ]; then
        echo "The test folder should be prepared as ~/data/spectroscape"
        echo "three subfolders: mass_spectra, new_data, and spectral_archives"
        exit 1
    fi
    if [ ! -d ~/data/spectroscape/new_data ]; then
        echo "The test folder should be prepared as ~/data/spectroscape"
        echo "three subfolders: mass_spectra, new_data, and spectral_archives"
        exit 1
    fi
    if [ ! -d ~/data/spectroscape/spectral_archives ]; then
        echo "The test folder should be prepared as ~/data/spectroscape"
        echo "three subfolders: mass_spectra, new_data, and spectral_archives"
        exit 1
    fi

}


version=1.1.8

function test_version() {
    version=$1
    echo "testing version: ${version}"
    pushd ~/data/spectroscape
    wget https://github.com/wulongict/SpectralArchive/releases/download/v${version}/Spectroscape_CPU-${version}.deb

    # install
    # remove
    sudo apt remove -y spectroscape_cpu
    sudo apt install -y ./Spectroscape_CPU-${version}.deb
    rm ./Spectroscape_CPU-${version}.deb
    # download the data
    
    test_path=spectral_archives_${version}
    mkdir -p ${test_path}

    cd ${test_path}
    # test with a clean start
    rm -rf *
    # run test
    spectroscape --init --datasearchpath ../mass_spectra
    spectroscape --run 
    spectroscape --add --datasearchpath ../new_data
    popd
}


function test_latest() {
    version=$1
    echo "testing version: ${version}"
    pushd ~/data/spectroscape
    (cd ~/code/SpectralArchive && ./compile.bash)
    test_path=spectral_archives_${version}
    mkdir -p ${test_path}

    cd ${test_path}
    # test with a clean start
    rm -rf *
    # run test
    ~/code/SpectralArchive/build/bin/spectroscape --init --datasearchpath ../mass_spectra
    ~/code/SpectralArchive/build/bin/spectroscape --run 
    ~/code/SpectralArchive/build/bin/spectroscape --add --datasearchpath ../new_data
    popd
}



# print help information, if no parameter is provided.
if (( $# < 1 )); then
    echo "Usage: test_versions.bash version1 version2 ..."
    echo "version1, version2, ... are the versions to be tested."
    echo -e "\nPrerequisite:" 
    echo -e "1. download the test data from GoogleDrive link:"
    echo -e "\tgo to"
    echo -e "\thttps://drive.google.com/drive/folders/1ZaHdkzxclzYboYf_67pPxgkmFxpnVFAG?usp="
    echo -e "\tdownload the data.zip file"
    echo -e "\n2. unzip it to ~/data/spectroscape/"
    echo -e "\tThen, there will be three subfolders: mass_spectra, new_data, and spectral_archives"
    echo -e "\n3. try run the following example to test the version 1.1.8"
    echo -e "\texample: $0 1.1.8 1.1.9"
    exit 1
else 
    verify_folder_settings
    for version in $@; do
        if [ $version = "latest" ]; then
            test_latest $version
        else
            test_version $version
        fi
        # test_version $version
    done
fi