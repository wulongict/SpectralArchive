# Query for nearest neighbors against a spectral archive with Spectroscape

## Compile from source code

Using the following command to compile the code. 

### Prerequisites

```bash
./install_prerequisites.bash
```

Install cmake version >= 3.20. The default version of Ubuntu 20.04 is cmake 3.16. We can update it using the following script.
```bash
./install_cmake_3.22.bash
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

