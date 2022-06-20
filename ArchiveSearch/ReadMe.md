# Create Spectral Archive from A list of Data files

## Introduction
This tool is implemented to build spectral archive from a list of data files, in mzXML or mzML format. The spectra are indexed with library [FAISS](https://github.com/facebookresearch/faiss). The index enables us to retrieve the approximate nearest neighbors (ANNs) of a query spectrum in real-time. For a spectral archive of 1 million spectra, it takes only 0.5ms on average to retrieve the top 1024 ANNs of a query spectrum for a batch of 1000 query spectra. 

The spectral archive contains three core components, a indexer, a SQL-database for annotation, a mz file for peak list. We can refine the 1024 ANNs by calculating the true dot product similarity score with peaks retrieved from peak list file. Then the refined list of nearest neighbors can be visualized as a network using force directed graph from [D3.js](https://d3js.org). The SQL-database records the peptide identifications of each spectrum and the corresponding scores reported by search engines and post-processing tools. 

## User Manual
### Installation
#### requirements
- FAISS
- BOOST
- COMET
- SPECTRAST
- NGINX

```
sudo apt-get install -y spawn-fcgi
sudo apt-get install -y libarmadillo-dev
sudo apt-get install libgomp1
sudo apt-get install make
sudo apt-get install build-essential

sudo dpkg --configure -a


 8810  sudo dpkg -i cuda-repo-ubuntu1404-8-0-local-ga2_8.0.61-1_amd64.deb
 8811  sudo dpkg -i ./cuda-repo-ubuntu1404-8-0-local-ga2_8.0.61-1_amd64.deb
 8815  sudo dpkg -i "cuda-repo-ubuntu1404-8-0-local-ga2_8.0.61-1_amd64-deb"
 8869  dpkg -i code_1.49.3-1601661857_amd64.deb
 8870  sudo dpkg -i code_1.49.3-1601661857_amd64.deb
 8899  sudo dpkg -i libboost1.58-dev_1.58.0+dfsg-5ubuntu3.1_amd64.deb
 8900  sudo dpkg -i libstdc++-4.8-dev_4.8.5-4ubuntu2_amd64.deb
 8909  sudo dpkg -i libboost1.58-dev_1.58.0+dfsg-5ubuntu3.1_i386.deb
 8910  sudo dpkg -i libstdc++-4.8-dev_4.8.5-4ubuntu2_i386.deb
 8929  sudo dpkg -i code_1.49.3-1601661857_amd64.deb
 8987  sudo dpkg -i libboost-system1.58-dev_1.58.0+dfsg-5ubuntu3_amd64.deb
 8989  sudo dpkg -i libboost-system1.58-dev_1.58.0+dfsg-5ubuntu3_i386.deb
 8992  sudo dpkg -r libboost1.54-dev
 8993  sudo dpkg -r libboost-system1.54-dev
 8994  sudo dpkg -r libboost-system-dev:amd64
 8995  sudo dpkg -r libboost-system1.54-dev:amd64
 8996  sudo dpkg -i libboost1.58-dev_1.58.0+dfsg-5ubuntu3.1_amd64.deb
 8997  sudo dpkg -r libboost1.58-dev:i386
 8998  sudo dpkg --configure -a
 9000  sudo dpkg -i libboost1.58-dev_1.58.0+dfsg-5ubuntu3.1_i386.deb
 9298  sudo dpkg -i gcc-4.8-base_4.8.4-2ubuntu1~14.04.4_amd64.deb
 9299  sudo dpkg -i gcc-4.8-base_4.8.4-2ubuntu1~14.04.4_i386.deb
 9301  sudo dpkg -i libstdc++6_4.8.4-2ubuntu1~14.04.4_i386.deb
 9302  sudo dpkg -i libstdc++6_4.8.4-2ubuntu1~14.04.4_amd64.deb
 9328  sudo dpkg -i libstdc++6_4.8.2-19ubuntu1_i386.deb
 9331  sudo dpkg -i libstdc++6_4.8.2-19ubuntu1_amd64.deb
 9332  sudo dpkg -i libstdc++6_4.8.2-19ubuntu1_i386.deb
 9370  sudo dpkg -i gcc-4.8-base_4.8.5-4ubuntu8~14.04.2_amd64.deb
 9372  sudo dpkg -i gcc-4.8-base_4.8.5-4ubuntu2_i386.deb
 9373  sudo dpkg -i gcc-4.8-base_4.8.5-4ubuntu2_amd64.deb
 9377  sudo dpkg -i libstdc++-4.8-dev_4.8.5-4ubuntu2_amd64.deb
 9378  sudo dpkg -i libstdc++-4.8-dev_4.8.5-4ubuntu2_i386.deb
 9381  sudo dpkg -i libgcc-4.8-dev_4.8.5-4ubuntu2_amd64.deb
 9383  sudo dpkg -i libgcc-4.8-dev_4.8.5-4ubuntu2_i386.deb
 9395  sudo dpkg -i libgcc-4.8-dev_4.8.2-19ubuntu1_i386.deb
 9407  sudo dpkg -i apt_1.0.1ubuntu2_arm64.deb
 9411  sudo dpkg -i libapt-pkg4.12_1.0.1ubuntu2_amd64.deb
 9412  sudo dpkg -i libapt-pkg4.12_1.0.1ubuntu2_i386.deb
 9413  dpkg --config
 9416  dpkg --configure -a
 9417  sudo dpkg --configure -a
 9418  sudo dpkg -r libapt-pkg5.0
 9432  sudo dpkg –configure -a
 9433  sudo dpkg –-configure -a
 9441  dpkg -l | grep ^iU | awk '{print $2}' | xargs sudo dpkg --purge
 9445  sudo dpkg -i apt_1.0.1ubuntu2_i386.deb
 9446  sudo dpkg -i apt_1.0.1ubuntu2_amd64.deb
 9447  sudo dpkg -r lintian
 9448  sudo dpkg -r libapt-pkg-perl
 9449  sudo dpkg -r libapt-pkg5.0:amd64
 9450  sudo dpkg -r apt-utils
 9569  sudo dpkg -i libstdc++6_4.8.2-19ubuntu1_i386.deb
 9570  sudo dpkg -i libstdc++-4.8-dev_4.8.5-4ubuntu2_amd64.deb
 9572  sudo dpkg -i libstdc++6_4.8.2-19ubuntu1_amd64.deb
 9573  sudo dpkg --force-overwrite  -i libstdc++6_4.8.2-19ubuntu1_amd64.deb
 9574  sudo dpkg --force-overwrite  -i libstdc++6_4.8.2-19ubuntu1_i386.deb
 9575  sudo dpkg -r gcc-4.8-base
 9576  sudo dpkg -i --force-overwrite  gcc-4.8-base_4.8.2-19ubuntu1_amd64.deb
 9577  sudo dpkg -r gcc-4.8-base:i386
 9578  sudo dpkg -r libstdc++6:i386
 9581  sudo dpkg -i gcc-4.8-base_4.8.2-19ubuntu1_amd64.deb
 9794  sudo dpkg -i gcc-4.8-base_4.8.5-4ubuntu2_i386.deb
 9795  sudo dpkg -i gcc-4.8-base_4.8.5-4ubuntu2_amd64.deb
 9796  sudo dpkg -r gcc-4.8-base_4.8.5-4ubuntu2_amd64.deb
 9797  sudo dpkg -r gcc-4.8-base
 9798  sudo dpkg -r cpp-4.8
 9799  sudo dpkg -r gcc-4.8
 9802  sudo dpkg -r gcc
 9803  sudo dpkg -r dkms
 9831  sudo dpkg -i gcc-4.8-base_4.8.2-19ubuntu1_amd64.deb
 9832  sudo dpkg -i gcc-4.8-base_4.8.2-19ubuntu1_i386.deb
 9870  sudo dpkg -i gedit_3.10.4-0ubuntu4_amd64.deb
 9935  sudo dpkg -i libc6-amd64_2.19-0ubuntu6_i386.deb
 9936  sudo dpkg -r libc6:i386
 9975  sudo dpkg -i binutils-2.26_2.26.1-1ubuntu1~14.04_i386.deb
 9977  sudo dpkg -i binutils-2.26_2.26.1-1ubuntu1~14.04_amd64.deb
 9992  history | grep dpkg



sudo apt-get install -y gcc-4.8 g++-4.8 gcc-4.9 g++-4.9 gcc-5 g++-5 gcc-6 g++-6 gcc-7 g++-7 gcc-8 g++-8 gcc-9 g++-9 

sudo apt-get install cmake
sudo apt-get install libcgi-session-perl 
sudo apt-get install cuda

sudo apt-get install libopenblas-dev liblapack-dev
sudo apt-get install -y libarmadillo-dev
```

#### compile from source code


```
cd ArchiveSearch/
cmake .
make -j 3
```
### Quick One
Prepare a list of mzxml files. 
```
/path/to/a.mzxml
/path/to/b.mzxml
/path/to/c.mzxml
```
Create the spectral archive with the following command:

