#ifndef SPECTRASTCONSTANTS_HPP_
#define SPECTRASTCONSTANTS_HPP_

#define SPECTRAST_VERSION 5
#define SPECTRAST_SUB_VERSION 0

#define MAX_LINE 8192

// longest query name - will really just affect output in a fixed-width format
#define MAX_NAME_LEN 64

// m/z limits of any peaks in any spectra
#define MAX_MZ 2010
#define MIN_MZ 10

#define DEFAULT_SEARCH_PARAMS_FILE "spectrast.params"
#define DEFAULT_CREATE_PARAMS_FILE "spectrast_create.params"

// #define ANNOTATE_MZ_TOLERANCE 1.0
// #define ISOTOPIC_PEAK_DELTA_TOLERANCE 0.3

// #define ALIGN_MZ_TOLERANCE 0.5

// this is the maximum number of *mzXML* files opened by SpectraST at any given time
// if your system allows the opening of N files simultaneously, set this number not to exceed approximately N/2 - 10
// this is to allow for opening the same number of output files during searching, plus the opening of library files and log files
// but opening too many files can bog down the file system, and may not be efficient. 
#define MAX_NUM_OPEN_FILES 50

//#define DECOY_BATCH_SIZE 100
//#define DECOY_PIECE_SIZE 200

#ifdef STANDALONE_LINUX
#define szTPPVersionInfo "STANDALONE"
#define MAX_CHARGE 7
#else
#include "common/TPPVersion.h"
#endif

#endif /* SPECTRASTCONSTANTS_HPP_ */
