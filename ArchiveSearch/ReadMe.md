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

