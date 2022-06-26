# the following bash commands is from:
# https://www.linuxcapable.com/install-cmake-on-ubuntu-20-04-lts/
threads=`nproc`

sudo apt install build-essential checkinstall zlib1g-dev libssl-dev -y
wget -nc https://github.com/Kitware/CMake/releases/download/v3.22.2/cmake-3.22.2.tar.gz
tar -zxvf cmake-3.22.2.tar.gz
cd cmake-3.22.2
./bootstrap --parallel=$threads
make -j $threads
sudo make install
cd ..

rm -rf cmake-3.22.2
rm cmake-3.22.2.tar.gz

