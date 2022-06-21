# Building and searching spectral archive

## Compile from source code

Using the following command to compile the code. 

### Prerequisites

```bash
./install_prerequisites.bash
```

Install cmake version >= 3.20. The default version of Ubuntu 20.04 is cmake 3.16. We can update it using the following script.
```bash
# the following bash commands is from:
# https://www.linuxcapable.com/install-cmake-on-ubuntu-20-04-lts/

sudo apt install build-essential checkinstall zlib1g-dev libssl-dev -y
wget https://github.com/Kitware/CMake/releases/download/v3.22.2/cmake-3.22.2.tar.gz
tar -zxvf cmake-3.22.2.tar.gz
cd cmake-3.22.2
sudo ./bootstrap
sudo make
sudo make install
cd ..
```

Install faiss-1.7.1. This faiss library requires cmake version >=3.20.
```bash
cd ./faiss-1.7.1
./compile.bash
cd ..
```

### Compile spectral archive tool
From the source code folder, run the following scripts. 
```bash
# to have a clean start
./cleanMake.bash
# compile spectral archive tool.
./compile.bash
```

## Usage
### Build archive
#### Add MS data files

#### Add search results

### Search archive
#### Search a given mzXML file




