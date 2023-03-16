# Spectroscape: searching for similar PSMs in spectral archives

![example workflow](https://github.com/wulongict/SpectralArchive/actions/workflows/cmake.yml/badge.svg)

## Binary installation 
This following command has been teted on Ubuntu 22.04. The .deb file of latest version of spectroscape can be found in the following [link](https://github.com/wulongict/SpectralArchive/releases/latest). Spectroscape comes with CPU and GPU versions. If CUDA environment is not available, please use CPU version.  

Spectroscape (CPU version) can be installed on Ubuntu 22.04 using following command lines. 
```bash
wget https://github.com/wulongict/SpectralArchive/releases/download/v1.0.8/Spectroscape_CPU-1.0.8-Linux.deb
sudo apt update
sudo apt install ./Spectroscape_CPU-1.0.8-Linux.deb
```

The [GPU version](https://github.com/wulongict/SpectralArchive/releases/download/v1.0.8/Spectroscape_GPU-1.0.8-Linux.deb) can be installed similarly. 


```bash
wget https://github.com/wulongict/SpectralArchive/releases/download/v1.0.8/Spectroscape_GPU-1.0.8-Linux.deb
sudo apt update
sudo apt install ./Spectroscape_GPU-1.0.8-Linux.deb
```

However, users should first make sure CUDA enviroment avaiable.  . Otheriwise, users may see error information as follows. 

```bash
spectroscape: error while loading shared libraries: libcudart.so.11.0: cannot open shared object file: No such file or directory
```

## Uninstallation
use the following command line to remove spectroscape (both GPU and CPU versions) from Ubuntu system. 
```bash
sudo apt remove spectroscape_cpu spectroscape_gpu
```

## Source code installation

### Prerequisites

CMake and gcc are reqired to compilation of C++ code.  
```bash
sudo apt update
sudo apt install cmake build-essential 
```

The source code requires two extra libraries, libfcgi and liblapack. 

```bash
sudo apt install libfcgi-dev liblapack-dev 
```

To make the web interface work, two more tools should be installed, spawn-fcgi and nginx. 
```bash
sudo apt install spawn-fcgi nginx
```

Finally, to compile GPU version, CUDA environment is required. 

### Compile

First, get the latest source code of spectroscape from GitHub. 

```bash
 git clone --recurse-submodules  https://github.com/wulongict/SpectralArchive.git
```

Start from here, all the command should be excuted under the source code folder, namely, SpectralArchive. 

Run the following scripts to remove any intermediate files generated from previous failed compilation and have a clean start. 
```bash
./cleanMake.bash
```

Users can compile a CPU or GPU version using option FALSE or TRUE. 

```bash
# CPU version
./compile.bash FALSE
```

```bash
# GPU version 
./compile.bash TRUE
```

After the compilation, the excutable files are under the build/bin folder inside the source code directory. 
```bash
build/
├── bin
├── include
├── lib
└── share
```

## Usage
### Build archive
#### Add MS data files
First create a text file, mzxmllist, which contains a list of mzXML files to initialize a spectral archive. The raw files corresponding to the mzXML files below can be downloaded from pride archive [PXD000561](http://ftp.ebi.ac.uk/pride-archive/2014/04/PXD000561/).
```bash
$ cat mzxmllist
Adult_Adrenalgland_Gel_Elite_49_f01.mzXML
Adult_Adrenalgland_Gel_Elite_49_f02.mzXML
Adult_Adrenalgland_Gel_Elite_49_f03.mzXML
Adult_Adrenalgland_Gel_Elite_49_f04.mzXML
Adult_Adrenalgland_Gel_Elite_49_f05.mzXML
Adult_Adrenalgland_Gel_Elite_49_f06.mzXML
Adult_Adrenalgland_Gel_Elite_49_f07.mzXML
Adult_Adrenalgland_Gel_Elite_49_f08.mzXML
Adult_Adrenalgland_Gel_Elite_49_f09.mzXML
Adult_Adrenalgland_Gel_Elite_49_f10.mzXML
Adult_Adrenalgland_Gel_Elite_49_f11.mzXML
Adult_Adrenalgland_Gel_Elite_49_f12.mzXML
Adult_Adrenalgland_Gel_Elite_49_f13.mzXML
Adult_Adrenalgland_Gel_Elite_49_f14.mzXML
Adult_Adrenalgland_Gel_Elite_49_f15.mzXML
Adult_Adrenalgland_Gel_Elite_49_f16.mzXML
Adult_Adrenalgland_Gel_Elite_49_f17.mzXML
Adult_Adrenalgland_Gel_Elite_49_f18.mzXML
Adult_Adrenalgland_Gel_Elite_49_f19.mzXML
Adult_Adrenalgland_Gel_Elite_49_f20.mzXML
Adult_Adrenalgland_Gel_Elite_49_f21.mzXML
Adult_Adrenalgland_Gel_Elite_49_f22.mzXML
Adult_Adrenalgland_Gel_Elite_49_f23.mzXML
Adult_Adrenalgland_Gel_Elite_49_f24.mzXML

```
Then run the following command to build a spectral archive
```bash
SpectralArchive/build/bin/spectroscape  -m mzxmllist
```
#### Add search results
The spectral archive should be properly annotated. Currently, it supports the following input format.
- .pep.xml file generated by xinteract or search engine (e.g. Comet)
- spectral library, text foramt, sptxt

Run the following command to update the annotation of spectra in the 24 mzXML used above. One can get the interact-Adult_Adrenalgland_Gel_Elite_49.ipro.pep.xml file from a Comet+xinteract database searching pipeline in TPP.
```
SpectralArchive/build/bin/spectroscape  -m mzxmllist --update --updategt interact-Adult_Adrenalgland_Gel_Elite_49.ipro.pep.xml
```

### Add new MS data file
The spectral archive can be expanded to include more MS data file. Run the following command to add a new mzXML file.
```bash
SpectralArchive/build/bin/spectroscape  -m mzxmllist --update --updaterawdata <input>.mzXML
```
Currently, it supports the following input formats of MS data file.
- mzXML
- mzML
- sptxt

### Search archive
#### Search a given mzXML file
Run the following command to search a data file. Before searching against an archvie, make sure the spectra are annotated by search results under FDR control, e.g. annotated by pepXML files of iProphet/PeptidePropeht. 

```bash
SpectralArchive/build/bin/spectroscape  -m mzxmllist --inputsource cmd --datafile <input>.mzXML
```



## Issues
- When running archive tool, I got an error said "libdpgpu.so: cannot open shared object file: No such file or directory"?  
    If the binary of archive tool is called /path/to/archive/bin, then try add the library path, /path/to/archive/lib to LD_LIBRARY_PATH variable.
    One can check the RUNPATH of a binary or library file with following command: 
    ```bash
    objdump -x /path/to/the/binary-or-library-file | grep RUNPATH
    ```

