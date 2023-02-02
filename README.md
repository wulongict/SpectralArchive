# Query for nearest neighbors against a spectral archive with Spectroscape

## Install with a deb file
Users could install spectroscape from compiled binary files. This following command has been teted on Ubuntu 22.04. The .deb file of version v1.0.5 can be found in the following [link](https://github.com/wulongict/SpectralArchive/releases/download/v1.0.5/Spectroscape-1.0.5-Linux.deb). More newer version will be released on [this page](https://github.com/wulongict/SpectralArchive/releases). 

To use the following deb package, users have to install libopenblas-dev first with the following command. 
```bash
sudo apt-get install libopenblas-dev
```

The CUDA enviroment is also required. If not, users could try use CPU-only version, which can be downloaded from release v1.0.6

```bash 
# the following wget command downloads the v1.0.5 version. Change the url to get a newer version.
wget https://github.com/wulongict/SpectralArchive/releases/download/v1.0.5/Spectroscape-1.0.5-Linux.deb
sudo dpkg -i Spectroscape-1.0.5-Linux.deb
```

User may encounter the following error: 
```bash
spectroscape: error while loading shared libraries: libopenblas.so.0: cannot open shared object file: No such file or directory
```

Try the following command to get the missing package.

```bash
sudo apt-get install libopenblas-dev
```

User may also encounter the following error:
```bash
spectroscape: error while loading shared libraries: libcudart.so.11.0: cannot open shared object file: No such file or directory
```

This means that there is no CUDA available on this computer, we can only use CPU. 

# Jan 26, 2022, try to install  in binary, fail, try to fix it with two different binary package. 

## Compile from source code

Using the following command to compile the code. 

### Prerequisites

To compile the source code, you need to have cmake, gcc installed. 
```bash
sudo apt install cmake build-essential 
```

The source code requires two extra libraries, fcgi and lapack. 

```bash
sudo apt install libfcgi-dev liblapack-dev
```

To make the web service work, two more tools should be installed, spawn-fcgi and nginx. 
```bash
sudo apt install spawn-fcgi nginx
```


### Compile

First you can run the following command to get the source code from GitHub. 

```bash
 git clone --recurse-submodules  https://github.com/wulongict/SpectralArchive.git
```

From the source code folder, run the following scripts. 
```bash
# to have a clean start
./cleanMake.bash

```

Users can compile a GPU version or CPU version using option TRUE or FALSE. 

```bash
# compile spectral archive tool, compile GPU version with 
./compile.bash TRUE
```

```bash
# compile spectral archive tool, compile CPU version with 
./compile.bash FALSE
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
SpectralArchive/build/bin/fastcgi_similarity.fcgi  -m mzxmllist
```
#### Add search results
The spectral archive should be properly annotated. Currently, it supports the following input format.
- .pep.xml file generated by xinteract or search engine (e.g. Comet)
- spectral library, text foramt, sptxt

Run the following command to update the annotation of spectra in the 24 mzXML used above. One can get the interact-Adult_Adrenalgland_Gel_Elite_49.ipro.pep.xml file from a Comet+xinteract database searching pipeline in TPP.
```
SpectralArchive/build/bin/fastcgi_similarity.fcgi  -m mzxmllist --update --updategt interact-Adult_Adrenalgland_Gel_Elite_49.ipro.pep.xml
```

### Add new MS data file
The spectral archive can be expanded to include more and more MS data file. Run the following command to add a new mzXML file.
```bash
SpectralArchive/build/bin/fastcgi_similarity.fcgi  -m mzxmllist --update --updaterawdata <input>.mzXML
```
Currently, it supports the following input formats of MS data file.
- mzXML
- mzML
- sptxt

### Search archive
#### Search a given mzXML file
Run the following command to search a data file. Before searching against an archvie, make sure the spectra are annotated by search results under FDR control, e.g. annotated by pepXML files of iProphet/PeptidePropeht. 

```bash
SpectralArchive/build/bin/fastcgi_similarity.fcgi  -m mzxmllist --inputsource cmd --datafile <input>.mzXML
```



## Issues
- When running archive tool, I got an error said "libdpgpu.so: cannot open shared object file: No such file or directory"?  
    If the binary of archive tool is called /path/to/archive/bin, then try add the library path, /path/to/archive/lib to LD_LIBRARY_PATH variable.
    One can check the RUNPATH of a binary or library file with following command: 
    ```bash
    objdump -x /path/to/the/binary-or-library-file | grep RUNPATH
    ```

