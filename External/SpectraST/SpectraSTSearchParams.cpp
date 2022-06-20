#include "SpectraSTSearchParams.hpp"
#include "SpectraSTLog.hpp"
#include "SpectraSTConstants.hpp"
#include "FileUtils.hpp"
#include <fstream>
#include <iostream>
#include <stdlib.h>


/*

Program       : Spectrast
Author        : Henry Lam <hlam@systemsbiology.org>                                                       
Date          : 03.06.06 


Copyright (C) 2006 Henry Lam

This library is free software; you can redistribute it and/or
modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation; either
version 2.1 of the License, or (at your option) any later version.

This library is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public
License along with this library; if not, write to the Free Software
Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307, USA

Henry Lam
Institute for Systems Biology
401 Terry Avenue North 
Seattle, WA  98109  USA
hlam@systemsbiology.org

*/

/* Class: SpectraSTSearchParams
 * 
 * Manages the search parameters.
 *  
 */



using namespace std;

extern bool g_quiet;
extern SpectraSTLog *g_log;

// constructor
SpectraSTSearchParams::SpectraSTSearchParams() :
        m_options() {
    setDefault();

}

// copy constructor
SpectraSTSearchParams::SpectraSTSearchParams(SpectraSTSearchParams &s) {
    (*this) = s;

}

// assignment operator
SpectraSTSearchParams &SpectraSTSearchParams::operator=(SpectraSTSearchParams &s) {

    this->paramsFileName = s.paramsFileName;
    this->libraryFile = s.libraryFile;
    this->databaseFile = s.databaseFile;
    this->databaseType = s.databaseType;
    this->indexCacheAll = s.indexCacheAll;
    this->filterSelectedListFileName = s.filterSelectedListFileName;

    this->precursorMzTolerance = s.precursorMzTolerance;
    // this->indexRetrievalMzTolerance = s.indexRetrievalMzTolerance;
    this->precursorMzUseAverage = s.precursorMzUseAverage;
    // this->indexRetrievalUseAverage = s.indexRetrievalUseAverage;
    // this->expectedCysteineMod = s.expectedCysteineMod;
    this->detectHomologs = s.detectHomologs;
    // this->ignoreChargeOneLibSpectra = s.ignoreChargeOneLibSpectra;
    // this->ignoreAbnormalSpectra = s.ignoreAbnormalSpectra;
    // this->ignoreSpectraWithUnmodCysteine = s.ignoreSpectraWithUnmodCysteine;
    this->searchAllCharges = s.searchAllCharges;
    this->useSp4Scoring = s.useSp4Scoring;
    this->usePValue = s.usePValue;
    this->useTierwiseOpenModSearch = s.useTierwiseOpenModSearch;
    this->useRankTransformWithQuota = s.useRankTransformWithQuota;
    this->useRankTransformWithQuotaNumberOfPeaks = s.useRankTransformWithQuotaNumberOfPeaks;
    this->useRankTransformWithQuotaWindowSize = s.useRankTransformWithQuotaWindowSize;

    this->outputExtension = s.outputExtension;
    this->outputDirectory = s.outputDirectory;
    this->hitListTopHitFvalThreshold = s.hitListTopHitFvalThreshold;
    this->hitListLowerHitsFvalThreshold = s.hitListLowerHitsFvalThreshold;
    this->hitListShowMaxRank = s.hitListShowMaxRank;
    this->hitListExcludeNoMatch = s.hitListExcludeNoMatch;
    this->hitListShowHomologs = s.hitListShowHomologs;
    this->enzymeForPepXMLOutput = s.enzymeForPepXMLOutput;
    this->printFingerprintingSummary = s.printFingerprintingSummary;

    this->peakScalingMzPower = s.peakScalingMzPower;
    this->peakScalingIntensityPower = s.peakScalingIntensityPower;
    this->peakScalingUnassignedPeaks = s.peakScalingUnassignedPeaks;
    this->peakBinningNumBinsPerMzUnit = s.peakBinningNumBinsPerMzUnit;
    this->peakBinningFractionToNeighbor = s.peakBinningFractionToNeighbor;
    this->peakNoBinning = s.peakNoBinning;

    this->filterAllPeaksBelowMz = s.filterAllPeaksBelowMz;
    this->filterMinPeakCount = s.filterMinPeakCount;
    this->filterCountPeakIntensityThreshold = s.filterCountPeakIntensityThreshold;
    this->filterRemovePeakIntensityThreshold = s.filterRemovePeakIntensityThreshold;
    // this->filterRemoveHuge515Threshold = s.filterRemoveHuge515Threshold;
    this->filterMinMzRange = s.filterMinMzRange;
    this->filterMaxPeaksUsed = s.filterMaxPeaksUsed;
    this->filterMaxDynamicRange = s.filterMaxDynamicRange;
    this->filterMaxIntensityBelow = s.filterMaxIntensityBelow;
    this->filterLibMaxPeaksUsed = s.filterLibMaxPeaksUsed;
    this->filterITRAQReporterPeaks = s.filterITRAQReporterPeaks;
    this->filterTMTReporterPeaks = s.filterTMTReporterPeaks;
    this->filterLightIonsMzThreshold = s.filterLightIonsMzThreshold;

    this->fvalFractionDelta = s.fvalFractionDelta;
    this->fvalUseDotBias = s.fvalUseDotBias;

    this->m_options = s.m_options;

    return (*this);

}

// destructor
SpectraSTSearchParams::~SpectraSTSearchParams() {

}

// finalizeOptions - go through the options stored in m_options and actually set them.
void SpectraSTSearchParams::finalizeOptions() {


    for (vector<string>::iterator i = m_options.begin(); i != m_options.end(); i++) {

        double f = -1.0;
        int k = 0;
        bool valid = false;
        string extension;

        char optionType = (*i)[0];
        string optionValue = (*i).length() > 1 ? (*i).substr(1) : "";

        // all the command line options. note that if the option F (read from params file) is specified,
        // it should be the first option to be processed (see addOption() method). this way, the options
        // set in the params file can be overridden by other command line options.

        switch (optionType) {

            case 'F' :
                if (optionValue.empty()) {
                    // no params file name specified. use default
                    paramsFileName = DEFAULT_SEARCH_PARAMS_FILE;
                    readFromFile();
                    valid = true;
                } else {
                    //fixpath(optionValue);
                    paramsFileName = optionValue;
                    readFromFile();
                    valid = true;

                }
                if (!g_quiet) {
                    cout << "Search parameter file loaded: \"" << paramsFileName << "\"." << endl;
                }
                break;

            case 'L' :
                if (!optionValue.empty()) {
                    //fixpath(optionValue);
                    getExtension(optionValue, extension);
                    if (extension == ".splib") {
                        libraryFile = optionValue;
                        valid = true;
                    }
                }
                break;

            case 'D' :
                if (!optionValue.empty()) {
                    //fixpath(optionValue);
                    databaseFile = optionValue;
                    valid = true;
                }
                break;

            case 'T' :
                if (databaseType == "AA" || databaseType == "DNA") {
                    databaseType = optionValue;
                    valid = true;
                }
                break;


            case 'R' :
                if (optionValue.empty()) {
                    indexCacheAll = true;
                    valid = true;
                } else if (optionValue == "!") {
                    indexCacheAll = false;
                    valid = true;
                }
                break;

            case 'S' :
                if (!optionValue.empty()) {
                    //fixpath(optionValue);
                    filterSelectedListFileName = optionValue;
                    valid = true;
                }
                break;



                // CANDIDATE SELECTION AND SCORING

            case 'M' :

                if (!optionValue.empty()) {
                    f = atof(optionValue.c_str());
                    if (f >= 0.0) {
                        precursorMzTolerance = f;
                        // indexRetrievalMzTolerance = f;
                        valid = true;
                    }
                }
                break;

            case 'A' :
                if (optionValue.empty()) {
                    precursorMzUseAverage = true;
                    // indexRetrievalUseAverage = true;
                    valid = true;
                } else if (optionValue == "!") {
                    precursorMzUseAverage = false;
                    // indexRetrievalUseAverage = false;
                    valid = true;
                }
                break;


/*      
    case 'C' :
      if (optionValue == "ICAT_cl" || optionValue == "ICAT_uc" || optionValue == "CAM") {
        expectedCysteineMod = optionValue;
        valid = true;
      } else if (optionValue.empty() || optionValue == "!") {
        expectedCysteineMod = "";
        valid = true;					
      }
      break;
		
      
    case 'c' :
      if (optionValue.empty()) {
        ignoreSpectraWithUnmodCysteine = true;
        valid = true;
      } else if (optionValue == "!") {
        ignoreSpectraWithUnmodCysteine = false;
        valid = true;
      }       
      break;
  */

            case 'z' :
                if (optionValue.empty()) {
                    searchAllCharges = true;
                    valid = true;
                } else if (optionValue == "!") {
                    searchAllCharges = false;
                    valid = true;
                }
                break;

                // OUTPUT DISPLAY
            case 'E' :
                if (optionValue == "xls" || optionValue == "nxls" || optionValue == "txt" || optionValue == "ntxt" ||
                    optionValue == "pep.xml" || optionValue == "xml" ||
                    optionValue == "pepXML" || optionValue == "html") {
                    outputExtension = optionValue;
                    valid = true;
                }
                break;

            case 'O' :
                if (!optionValue.empty()) {
                    outputDirectory = optionValue;
                    if (outputDirectory[outputDirectory.length() - 1] != '/') {
                        outputDirectory += "/";
                    }
                    valid = true;
                }
                break;

                // ADVANCED OPTIONS
            case '_' :
                finalizeAdvancedOption(optionValue);
                valid = true;
                break;

            case 'H' :
            case '2' :
            case 'G' :
            case 'V' :
            case 'v' :
            case '1' :
            case '0' :
            case 'h' :
            case 'k' :
            case 'm' :
            case 'i' :
            case 'u' :
            case 'b' :
            case 'f' :
            case 'w' :
            case 'n' :
            case 'P' :
            case 'p' :
            case '5' :
            case 'Q' :
            case 'y' :
            case 'x' :
            case 'C' :
            case 'c' :
                if (!g_quiet) {
                    cout << "Option \"-s" << optionType << "\" is deprecated.  Ignored. ";
                    cout << "Type \"spectrast -s_\" to see if an equivalent advanced option is available." << endl;
                }
                valid = true;
                break;

            default :
                if (!g_quiet) cout << "Option \"-s" << optionType << "\" is undefined. Ignored. ";
                valid = true;
                break;

        }


        if (!valid && !g_quiet) {
            cout << "Invalid value for the \"-s" << optionType << "\" in the command line. Ignored." << endl;
        }

    }

    // resolve conflicting options
    if (outputExtension == ".xml" || outputExtension == ".pepXML") {
        hitListExcludeNoMatch = true;
    }

    if (!useSp4Scoring) {
        fvalUseDotBias = false;
    }

    if (useTierwiseOpenModSearch) {
        useSp4Scoring = false;
        peakNoBinning = true;
        fvalFractionDelta = 0.0;
        fvalUseDotBias = false;
    }

    // Fingerprinting
    if (!(printFingerprintingSummary.empty())) {
        indexCacheAll = true;
    }
    // END Fingerprinting

}

void SpectraSTSearchParams::finalizeAdvancedOption(string option) {

    if (option.length() < 3) {
        cout << "Advanced option \"-s_" << option << " is undefined. Ignored." << endl;
        return;
    }

    double f = -1.0;
    int k = 0;
    bool valid = false;

    string optionType = option.substr(0, 3);
    string optionValue = option.substr(3);

    if (optionType == "HOM") {

        if (!optionValue.empty()) {
            k = atoi(optionValue.c_str());
            if (k >= 0) {
                detectHomologs = k;
                valid = true;
            }
        }

    } else if (optionType == "RTQ") {

        if (optionValue.empty()) {
            useRankTransformWithQuota = true;
            valid = true;
        } else if (optionValue == "!") {
            useRankTransformWithQuota = false;
            valid = true;
        }

    } else if (optionType == "RTN") {

        if (!optionValue.empty()) {
            k = atoi(optionValue.c_str());
            if (k >= 0) {
                useRankTransformWithQuotaNumberOfPeaks = k;
                useRankTransformWithQuota = true;
                valid = true;
            }
        }

    } else if (optionType == "RTW") {

        if (!optionValue.empty()) {
            k = atoi(optionValue.c_str());
            if (k >= 0) {
                useRankTransformWithQuotaWindowSize = k;
                useRankTransformWithQuota = true;
                valid = true;
            }
        }
/*
  } else if (optionType == "NO1") {
    if (optionValue.empty()) {
      ignoreChargeOneLibSpectra = true;
      valid = true;
    } else if (optionValue == "!") {
      ignoreChargeOneLibSpectra = false;
      valid = true;
    }

  } else if (optionType == "NOS") {
    if (optionValue.empty()) {
      ignoreAbnormalSpectra = true;
      valid = true;
    } else if (optionValue == "!") {
      ignoreAbnormalSpectra = false;
      valid = true;
    }
*/
    } else if (optionType == "FV1") {

        if (!optionValue.empty()) {
            f = atof(optionValue.c_str());
            hitListTopHitFvalThreshold = f;
            valid = true;
        }

    } else if (optionType == "FV2") {

        if (!optionValue.empty()) {
            f = atof(optionValue.c_str());
            hitListLowerHitsFvalThreshold = f;
            valid = true;
        }

    } else if (optionType == "SHM") {
        if (optionValue.empty()) {
            hitListExcludeNoMatch = true;
            valid = true;
        } else if (optionValue == "!") {
            hitListExcludeNoMatch = false;
            valid = true;
        }

    } else if (optionType == "SHH") {
        if (optionValue.empty()) {
            hitListShowHomologs = true;
            valid = true;
        } else if (optionValue == "!") {
            hitListShowHomologs = false;
            valid = true;
        }

    } else if (optionType == "SHR") {

        if (!optionValue.empty()) {
            k = atoi(optionValue.c_str());
            if (k >= 1) {
                hitListShowMaxRank = (unsigned int) k;
                valid = true;
            }
        }

    } else if (optionType == "FIN") {

        if (!optionValue.empty()) {
            fixpath(optionValue);
            printFingerprintingSummary = optionValue;
        }
        valid = true;

    } else if (optionType == "ENZ") {
        if (!optionValue.empty()) {
            enzymeForPepXMLOutput = optionValue;
            valid = true;
        }

        // SPECTRUM FILTERING AND PROCESSING

    } else if (optionType == "SP4") {
        if (optionValue.empty()) {
            useSp4Scoring = true;
            valid = true;
        } else if (optionValue == "!") {
            useSp4Scoring = false;
            valid = true;
        }

    } else if (optionType == "PVL") {
        if (optionValue.empty()) {
            usePValue = true;
            valid = true;
        } else if (optionValue == "!") {
            usePValue = false;
            valid = true;
        }

    } else if (optionType == "OMT") {
        if (optionValue.empty()) {
            useTierwiseOpenModSearch = true;
            valid = true;
        } else if (optionValue == "!") {
            useTierwiseOpenModSearch = false;
            valid = true;
        }

    } else if (optionType == "MZS") {

        if (!optionValue.empty()) {
            f = atof(optionValue.c_str());
            if (f >= 0.0) {
                peakScalingMzPower = f;
                valid = true;
            }
        }

    } else if (optionType == "INS") {

        if (!optionValue.empty()) {
            f = atof(optionValue.c_str());
            if (f >= 0.0) {
                peakScalingIntensityPower = f;
                valid = true;
            }
        }

    } else if (optionType == "UAS") {

        if (!optionValue.empty()) {
            f = atof(optionValue.c_str());
            if (f >= 0.0) {
                peakScalingUnassignedPeaks = f;
                valid = true;
            }
        }

    } else if (optionType == "BIN") {

        if (!optionValue.empty()) {
            k = atoi(optionValue.c_str());
            if (k > 0) {
                peakBinningNumBinsPerMzUnit = k;
                valid = true;
            }
        }

    } else if (optionType == "NEI") {

        if (!optionValue.empty()) {
            f = atof(optionValue.c_str());
            if (f >= 0.0) {
                peakBinningFractionToNeighbor = f;
                valid = true;
            }
        }

    } else if (optionType == "NOB") {
        if (optionValue.empty()) {
            peakNoBinning = true;
            valid = true;
        } else if (optionValue == "!") {
            peakNoBinning = false;
            valid = true;
        }

    } else if (optionType == "XMZ") {

        if (!optionValue.empty()) {
            f = atof(optionValue.c_str());
            if (f >= 0.0) {
                filterAllPeaksBelowMz = f;
                valid = true;
            }
        }

    } else if (optionType == "XNP") {

        if (!optionValue.empty()) {
            k = atoi(optionValue.c_str());
            if (k >= 0) {
                filterMinPeakCount = k;
                valid = true;
            }
        }

    } else if (optionType == "CNT") {

        if (!optionValue.empty()) {
            f = atof(optionValue.c_str());
            if (f >= 0.0) {
                filterCountPeakIntensityThreshold = f;
                valid = true;
            }

        }
    } else if (optionType == "RNT") {

        if (!optionValue.empty()) {
            f = atof(optionValue.c_str());
            if (f >= 0.0) {
                filterRemovePeakIntensityThreshold = f;
                valid = true;
            }
        }

/*
  } else if (optionType == "R51") {
 
    if (!optionValue.empty()) {
      f = atof(optionValue.c_str());
      if (f >= 0.0) {
	filterRemoveHuge515Threshold = f;
	valid = true;
      }
    }
  */

    } else if (optionType == "XMR") {

        if (!optionValue.empty()) {
            f = atof(optionValue.c_str());
            if (f >= 0.0) {
                filterMinMzRange = f;
                valid = true;
            }
        }

    } else if (optionType == "XMR") {

        if (!optionValue.empty()) {
            f = atof(optionValue.c_str());
            if (f >= 0.0) {
                filterMinMzRange = f;
                valid = true;
            }
        }

    } else if (optionType == "RNP") {

        if (!optionValue.empty()) {
            k = atoi(optionValue.c_str());
            if (k >= 1) {
                filterMaxPeaksUsed = k;
                valid = true;
            }
        }

    } else if (optionType == "RDR") {

        if (!optionValue.empty()) {
            f = atof(optionValue.c_str());
            if (f >= 1.0) {
                filterMaxDynamicRange = f;
                valid = true;
            }
        }

    } else if (optionType == "XIN") {

        if (!optionValue.empty()) {
            f = atof(optionValue.c_str());
            if (f >= 0.0) {
                filterMaxIntensityBelow = f;
                valid = true;
            }
        }

    } else if (optionType == "ITQ") {

        if (optionValue.empty()) {
            filterITRAQReporterPeaks = true;
            valid = true;
        } else if (optionValue == "!") {
            filterITRAQReporterPeaks = false;
            valid = true;
        }

    } else if (optionType == "TMT") {

        if (optionValue.empty()) {
            filterTMTReporterPeaks = true;
            valid = true;
        } else if (optionValue == "!") {
            filterTMTReporterPeaks = false;
            valid = true;
        }

    } else if (optionType == "RLI") {

        if (!optionValue.empty()) {
            f = atof(optionValue.c_str());
            if (f >= 0.0) {
                filterLightIonsMzThreshold = f;
                valid = true;
            }
        }

    } else if (optionType == "FDL") {

        if (!optionValue.empty()) {
            f = atof(optionValue.c_str());
            if (f >= 0.0 && f <= 1.0) {
                fvalFractionDelta = f;
                valid = true;
            }
        }

    } else if (optionType == "FBI") {
        if (optionValue.empty()) {
            fvalUseDotBias = true;
            valid = true;
        } else if (optionValue == "!") {
            fvalUseDotBias = false;
            valid = true;
        }


    } else if (optionType == "LNP") {

        if (!optionValue.empty()) {
            k = atoi(optionValue.c_str());
            if (k >= 1) {
                filterLibMaxPeaksUsed = (unsigned int) k;
                valid = true;
            }
        }

    } else {
        if (!g_quiet)
            cout << "Advanced option \"-s_" << option << " is undefined. Ignored." << endl;
        valid = true;
    }

    if (!valid && !g_quiet) {
        cout << "Invalid value for advanced option \"-s_" << option << ". Ignored." << endl;
    }

}

// addOption - add an option to the list of options
bool SpectraSTSearchParams::addOption(string option) {

    // to allow command-line override of options set by reading a params file (-sF option),
    // options are first read from the file, and then the command-line options are allowed
    // to overwrite them. Hence, we want to make sure the -sF option (if any) is the first
    // to be executed.

    // check if it's the -sF option, in which case, put as the first option
    if (option.length() == 0) {
        return (true);
    } else if (option[0] == 'F') {
        m_options.insert(m_options.begin(), option);
        return (true);
    } else {
        if (isExpectingArg(option)) {
            return (false);
        } else {
            // else, put at the back
            m_options.push_back(option);
            return (true);
        }
    }
}

bool SpectraSTSearchParams::isExpectingArg(string option) {

    return (option == "L" || option == "D" || option == "T" || option == "S" || option == "E" || option == "O");
}

// setDefault - sets the defaults for all options
void SpectraSTSearchParams::setDefault() {

    // GENERAL

    // the name of the params file
    paramsFileName = "";

    // the library file - actually, this is not strictly an "option" because if it's
    // not specified (on command-line or in the params file) the program will not run.
    libraryFile = "";

    // the database file and type - won't affect the search at all, but will show up
    // in the output (if in .pepXML format) for downstream processing
    databaseFile = "";
    databaseType = "AA";

    // name of file containing all query to be searched (will skip all others)
    filterSelectedListFileName = "";

    // whether or not to load the entire library in memory (faster but needs lots of RAM)
    indexCacheAll = false;

    // CANDIDATE SELECTION AND SCORING

    // precursor m/z tolerance -- m/z tolerance for candidates that are actually compared
    precursorMzTolerance = 3.0;

    // m/z tolerance for index retrieval, i.e.
    // guarantee to retrieve all library spectra of m/z within this tolerance of the query's
    // indexRetrievalMzTolerance = 3.0; // HARD-CODE FOR NOW

    // use isotopically average m/z. what actually happens is, the m/z window will be expanded to the left (lighter) by 1 Th.
    // e.g. if one searches AVERAGE +/- 3.0 Th, actually MONOISOTOPIC -4.0 Th to +3.0 Th will be searched
    precursorMzUseAverage = false;
    // indexRetrievalUseAverage = false;

    // expected cysteine modification: ICAT_cl, ICAT_uc or CAM. Search will ignore those library spectra
    // that have a different modification (but will still consider those without ANY modification)
    // i.e. if ICAT_cl is specified, all library spectra with ICAT_uc or CAM will be ignored, but those
    // with unmodified C will still be considered (unless if ignoreSpectraWithUnmodCystein is also turned on).
    // expectedCysteineMod = "";

    // whether to detect if the lower hits are homologous or identical to the top hit. It will then calculate delta dot by
    // taking the difference between the top hit and the first non-homologous lower hit, so that peptide ions with more than
    // one spectra in the library will not get penalized.
    detectHomologs = 4;

    // whether or not to ignore all +1 spectra in the library (they don't seem to be too reliable at this point)
    // ignoreChargeOneLibSpectra = false;

    // whether or not to ignore spectra that do not have a status of "Normal". These abnormal spectra could be those that
    // fail some quality tests (e.g. Inquorate and impure spectra)
    // ignoreAbnormalSpectra = false;

    // whether or not to ignore spectra of peptides with an unmodified cysteines. (If an alkylation agent like an ICAT agent or
    // iodoacetamide is used, one typically expects most, if not all, cysteines to be modified. Sequence searching typically uses a
    // static modification for alkylation which will not allow an unmodified cysteine, although it is possible that some cysteines will // remain unmodified.)
    // ignoreSpectraWithUnmodCysteine = false;

    // whether or not to ignore the charge state specified in the query spectrum (if any) and search all library spectra regardless of
    // charge state
    searchAllCharges = false;

    // Use SpectraST 4.0 scoring (sqrt intensity dot product, with dot bias)
    useSp4Scoring = false;

    // compute P-value and use it for scoring
    usePValue = false;

    // use tierwise open modifications search
    useTierwiseOpenModSearch = false;

    // use peak quota in a sliding window for rank transform
    useRankTransformWithQuota = false;
    useRankTransformWithQuotaNumberOfPeaks = 8;
    useRankTransformWithQuotaWindowSize = 30;

    // OUTPUT DISPLAY

    // the type of output
    outputExtension = "pep.xml";

    // the directory to which the result files are written
    // if empty, the result file is written to the same directory as its corresponding search file.
    outputDirectory = "";

    // threshold of the F value of the top hit, below which the found spectrum will be
    // removed from the hit list (the value of 0.0 means all top hits will be displayed
    // no matter how good they are)
    hitListTopHitFvalThreshold = 0.03;

    // threshold of the F value of the lower hit, above which the lower hit will also be
    // displayed (provided hitListOnlyTopHit is set to false).
    hitListLowerHitsFvalThreshold = 0.45;

    // shows the homologous/identical lower hits (if detected - see detectHomologs above)
    // in the hit list regardless of their f-values
    hitListShowHomologs = true;

    // only prints the top hit
    // hitListOnlyTopHit = true; // deprecated, replaced by hitListShowMaxRank = 1

    // shows the top N hits
    hitListShowMaxRank = 1;

    // exclude all NO_MATCH's - that query will not show up at all in the output file
    hitListExcludeNoMatch = true;

    // print fingerprinting summary to file specified by string. Empty string = not performing fingerprinting
    printFingerprintingSummary = "";

    // the enzyme to use for outputting to PepXML, and for calculating NTT and NMC (WON'T AFFECT WHAT CANDIDATES ARE SEARCHED!)
    enzymeForPepXMLOutput = "";


    // SPECTRUM FILTERING AND PROCESSING

    // Peaks are scaled by: scaled intensity = (mass)^(peakScalingMzPower) * (intensity)^(peakScalingIntensityPower)
    peakScalingMzPower = 0.0;
    peakScalingIntensityPower = 0.5;

    // Scaling for unassigned peaks (those having '?' annotation in the library spectrum)
    // This is applied after the intensity is scaled by peakScalingIntensityPower
    peakScalingUnassignedPeaks = 1.0; // 0.1;

    // This specifies the bin width.
    peakBinningNumBinsPerMzUnit = 1;

    // This specifies the "spread" of the peak when binning.
    // The bin corresponding to the mass of the peak gets its full scaled intensity,
    // whereas the two neighboring bins each get a fraction of the scaled intensities.
    peakBinningFractionToNeighbor = 0.5;

    // This disables binning and instead perform direct peak to peak comparison
    // Only valid for SP5
    peakNoBinning = false;

    // minimum number of significant peaks required in the query spectra
    filterMinPeakCount = 10;

    // minimum intensity for a peak to be considered significant
    filterCountPeakIntensityThreshold = 2.01;

    // if intensity is smaller than this value, it will be removed from the spectra completely
    filterRemovePeakIntensityThreshold = 2.01;

    // if there is a peak near 515.3 Th, remove it if its intensity is greater than this fraction of the
    // total intensity of all peaks (0.0 means don't remove)
    // filterRemoveHuge515Threshold = 0.0;

    // if the m/z range is smaller than this value, filter out this query spectrum
    filterMinMzRange = 350.0;

    // if almost all peaks (95% of total intensity) are below this m/z value, filter out this query spectrum
    filterAllPeaksBelowMz = 520;

    // removes all but the top X peaks
    filterMaxPeaksUsed = 150;

    // removes all peaks smaller than 1/X of the highest-intensity peak
    filterMaxDynamicRange = 1000.0;

    // if the highest intensity is below this value, filter out this query spectrum (0.0 means don't filter)
    filterMaxIntensityBelow = 0.0;

    // remove all but the top X peaks IN THE LIBRARY SPECTRA
    filterLibMaxPeaksUsed = 50;

    // remove the iTRAQ reporter peaks at 112-122 Th
    filterITRAQReporterPeaks = false;

    // remove the TMT reporter peaks at 126-132 Th
    filterTMTReporterPeaks = false;

    // remove light ions below the m/z threshold set below (y1 of R-ending peptide is at 175.0 Th)
    filterLightIonsMzThreshold = 180;

    // fraction of fval that is delta/dot (the rest is dot)
    fvalFractionDelta = 0.4;

    // use dot bias to penalize high-scoring matches with huge noise and dominant peak
    fvalUseDotBias = true;

}

// readFromFile - loads the options from the params file
// the params file has the format
// <param name> = <param value> (the '=' can be replaced by ':' or just blank space)
void SpectraSTSearchParams::readFromFile() {

    if (paramsFileName.empty()) {
        return;
    }

    ifstream fin;
    if (!myFileOpen(fin, paramsFileName)) {
        g_log->error("SEARCH", "Cannot open PARAMS file \"" + paramsFileName +
                               "\" for reading Search params. Using all defaults.");
        setDefault();
        return;
    }

    string line("");
    while (nextLine(fin, line, "", "")) {
        if (line == "_EOF_") {
            return;
        }
        string::size_type pos = 0;
        string param = nextToken(line, 0, pos, " #\t\r\n=:");
        string value = nextToken(line, pos, pos, " #\t\r\n", " \t\r\n=:");
        if (param.empty()) continue;

        double f;
        int k;
        bool valid = false;

        if (param == "libraryFile") {
            if (!value.empty()) {
                //fixpath(value);
                string extension;
                getExtension(value, extension);
                if (extension == ".splib") {
                    libraryFile = value;
                    valid = true;
                }
            }

        } else if (param == "databaseFile") {
            if (!value.empty()) {
                //fixpath(value);
                databaseFile = value;
                makeFullPath(databaseFile); // absolute path, save trouble later when parsed by RefreshParser
                valid = true;
            }

        } else if (param == "databaseType") {
            if (value == "AA" || value == "DNA") {
                databaseType = value;
                valid = true;
            }

        } else if (param == "indexCacheAll") {
            indexCacheAll = (value == "true");
            valid = true;

        } else if (param == "filterSelectedListFileName") {
            if (!value.empty()) {
                //fixpath(value);
                filterSelectedListFileName = value;
                valid = true;
            }


            // CANDIDATE SELECTION AND SCORING
        } else if (param == "precursorMzTolerance") {
            //    } else if (param == "indexRetrievalMzTolerance") {
            if (!value.empty()) {
                f = atof(value.c_str());
                if (f >= 0.0) {
                    precursorMzTolerance = f;
                    // indexRetrievalMzTolerance = f;
                    valid = true;
                }
            }

        } else if (param == "precursorMzUseAverage") {
            //    } else if (param == "indexRetrievalUseAverage") {
            //  indexRetrievalUseAverage = (value == "true");
            precursorMzUseAverage = (value == "true");
            valid = true;

            /*
            } else if (param == "expectedCysteineMod") {
              if (value == "ICAT_cl" || value == "ICAT_uc" || value == "CAM" || value.empty()) {
                expectedCysteineMod = value;
                valid = true;
              }
              */

        } else if (param == "detectHomologs") {
            if (!value.empty()) {
                k = atoi(value.c_str());
                if (k >= 0) {
                    detectHomologs = (unsigned int) k;
                    valid = true;
                }
            }

            /*
        } else if (param == "ignoreChargeOneLibSpectra") {
          ignoreChargeOneLibSpectra = (value == "true");
          valid = true;

        } else if (param == "ignoreAbnormalSpectra") {
          ignoreAbnormalSpectra = (value == "true");
          valid = true;

        } else if (param == "ignoreSpectraWithUnmodCysteine") {
          ignoreSpectraWithUnmodCysteine = (value == "true");
          valid = true;
          */

        } else if (param == "useSp4Scoring") {
            useSp4Scoring = (value == "true");
            valid = true;

        } else if (param == "useTierwiseOpenModSearch") {
            useTierwiseOpenModSearch = (value == "true");
            valid = true;

        } else if (param == "usePValue") {
            usePValue = (value == "true");
            valid = true;

            // OUTPUT DISPLAY

        } else if (param == "outputExtension") {
            if (value == "xls" || value == "txt" || value == "pep.xml" || value == "xml" || value == "pepXML" ||
                value == "html") {
                outputExtension = value;
                valid = true;
            }
        } else if (param == "outputDirectory") {
            if (!value.empty()) {
                outputDirectory = value;
                if (outputDirectory[outputDirectory.length() - 1] != '/') {
                    outputDirectory += "/";
                }
                valid = true;
            }


        } else if (param == "hitListTopHitFvalThreshold") {
            if (!value.empty()) {
                f = atof(value.c_str());
                hitListTopHitFvalThreshold = f;
                valid = true;
            }

        } else if (param == "hitListLowerHitsFvalThreshold") {
            if (!value.empty()) {
                f = atof(value.c_str());
                hitListLowerHitsFvalThreshold = f;
                valid = true;
            }

        } else if (param == "hitListShowMaxRank") {
            if (!value.empty()) {
                k = atoi(value.c_str());
                if (k >= 1) {
                    hitListShowMaxRank = (unsigned int) k;
                    valid = true;
                }
            }

        } else if (param == "hitListExcludeNoMatch") {
            hitListExcludeNoMatch = (value == "true");
            valid = true;

        } else if (param == "hitListShowHomologs") {
            hitListShowHomologs = (value == "true");
            valid = true;

        } else if (param == "enzymeForPepXMLOutput") {
            if (!value.empty()) {
                enzymeForPepXMLOutput = value;
                valid = true;
            }

        } else if (param == "printFingerprintingSummary") {
            if (!value.empty()) {
                fixpath(value);
                printFingerprintingSummary = value;
            }
            valid = true;

            // SPECTRUM FILTERING AND PROCESSING

        } else if (param == "peakScalingMzPower") {
            if (!value.empty()) {
                f = atof(value.c_str());
                if (f >= 0.0) {
                    peakScalingMzPower = f;
                    valid = true;
                }
            }

        } else if (param == "peakScalingIntensityPower") {
            if (!value.empty()) {
                f = atof(value.c_str());
                if (f >= 0.0) {
                    peakScalingIntensityPower = f;
                    valid = true;
                }
            }

        } else if (param == "peakScalingUnassignedPeaks") {
            if (!value.empty()) {
                f = atof(value.c_str());
                if (f >= 0.0) {
                    peakScalingUnassignedPeaks = f;
                    valid = true;
                }
            }

        } else if (param == "peakBinningNumBinsPerMzUnit") {
            if (!value.empty()) {
                k = atoi(value.c_str());
                if (k > 0) {
                    peakBinningNumBinsPerMzUnit = (unsigned int) k;
                    valid = true;
                }
            }

        } else if (param == "peakBinningFractionToNeighbor") {
            if (!value.empty()) {
                f = atof(value.c_str());
                if (f >= 0.0) {
                    peakBinningFractionToNeighbor = f;
                    valid = true;
                }
            }

        } else if (param == "peakNoBinning") {
            peakNoBinning = (value == "true");
            valid = true;

        } else if (param == "filterAllPeaksBelowMz") {
            if (!value.empty()) {
                f = atof(value.c_str());
                if (f >= 0.0) {
                    filterAllPeaksBelowMz = f;
                    valid = true;
                }
            }

        } else if (param == "filterMinPeakCount") {
            if (!value.empty()) {
                k = atoi(value.c_str());
                if (k >= 0) {
                    filterMinPeakCount = (unsigned int) k;
                    valid = true;
                }
            }

        } else if (param == "filterCountPeakIntensityThreshold") {
            if (!value.empty()) {
                f = atof(value.c_str());
                if (f >= 0.0) {
                    filterCountPeakIntensityThreshold = f;
                    valid = true;
                }
            }

        } else if (param == "filterRemovePeakIntensityThreshold") {
            if (!value.empty()) {
                f = atof(value.c_str());
                if (f >= 0.0) {
                    filterRemovePeakIntensityThreshold = f;
                    valid = true;
                }
            }

/*
    } else if (param == "filterRemoveHuge515Threshold") {
      if (!value.empty()) {
	f = atof(value.c_str());
	if (f >= 0.0 && f <= 1.0) {
	  filterRemoveHuge515Threshold = f;
	  valid = true;
	}
      }
*/

        } else if (param == "filterMinMzRange") {
            if (!value.empty()) {
                f = atof(value.c_str());
                if (f >= 0.0) {
                    filterMinMzRange = f;
                    valid = true;
                }
            }

        } else if (param == "filterMinMzRange") {
            if (!value.empty()) {
                f = atof(value.c_str());
                if (f >= 0.0) {
                    filterMinMzRange = f;
                    valid = true;
                }
            }

        } else if (param == "filterMaxPeaksUsed") {
            if (!value.empty()) {
                k = atoi(value.c_str());
                if (k >= 1) {
                    filterMaxPeaksUsed = (unsigned int) k;
                    valid = true;
                }
            }

        } else if (param == "filterMaxDynamicRange") {
            if (!value.empty()) {
                f = atof(value.c_str());
                if (f >= 1.0) {
                    filterMaxDynamicRange = f;
                    valid = true;
                }
            }

        } else if (param == "filterMaxIntensityBelow") {
            if (!value.empty()) {
                f = atof(value.c_str());
                if (f >= 0.0) {
                    filterMaxIntensityBelow = f;
                    valid = true;
                }
            }

        } else if (param == "filterLibMaxPeaksUsed") {
            if (!value.empty()) {
                k = atoi(value.c_str());
                if (k >= 1) {
                    filterLibMaxPeaksUsed = (unsigned int) k;
                    valid = true;
                }
            }

        } else if (param == "filterITRAQReporterPeaks") {
            filterITRAQReporterPeaks = (value == "true");
            valid = true;

        } else if (param == "filterTMTReporterPeaks") {
            filterTMTReporterPeaks = (value == "true");
            valid = true;

        } else if (param == "filterLightIonsMzThreshold") {
            if (!value.empty()) {
                f = atof(value.c_str());
                if (f >= 0.0) {
                    filterLightIonsMzThreshold = f;
                    valid = true;
                }
            }

        } else if (param == "fvalFractionDelta") {
            if (!value.empty()) {
                f = atof(value.c_str());
                if (f >= 0.0 && f <= 1.0) {
                    fvalFractionDelta = f;
                    valid = true;
                }
            }

        } else if (param == "fvalUseDotBias") {
            fvalUseDotBias = (value == "true");
            valid = true;


        } else {
            if (!g_quiet) {
                cout << "Unknown option in " << paramsFileName << " : \"" << param << "\". Ignored. " << endl;
            }
            valid = true;
        }

        if (!valid && !g_quiet) {
            cout << "Invalid value of \"" << param << "\" in " << paramsFileName << ". Ignored." << endl;
        }

    }
}


void SpectraSTSearchParams::printPepXMLSearchParams(ofstream &fout) {

    string fullLibraryFile(libraryFile);
    makeFullPath(fullLibraryFile);

    fout << "<parameter name=\"spectral_library\" value=\"" << fullLibraryFile << "\"/>" << endl;

    fout << "<parameter name=\"precursor_mz_tolerance\" value=\"" << precursorMzTolerance << "\"/>" << endl;
//  fout << "<parameter name=\"precursor_mz_tolerance\" value=\"" << indexRetrievalMzTolerance << "\"/>" << endl;

    fout << "<parameter name=\"precursor_mz_use_average\" value=\"" << (precursorMzUseAverage ? "true" : "false")
         << "\"/>" << endl;
//  fout << "<parameter name=\"precursor_mz_use_average\" value=\"" << (indexRetrievalUseAverage ? "true" : "false") << "\"/>" << endl;


    fout << "<parameter name=\"params_file\" value=\"" << paramsFileName << "\"/>" << endl;

    fout << "<parameter name=\"detect_homologs_max_rank\" value=\"" << detectHomologs << "\"/>" << endl;

    fout << "<parameter name=\"use_sp4_scoring\" value=\"" << (useSp4Scoring ? "true" : "false") << "\"/>" << endl;

    fout << "<parameter name=\"use_P_value\" value=\"" << (usePValue ? "true" : "false") << "\"/>" << endl;

    fout << "<parameter name=\"use_tierwise_open_modification_search\" value=\""
         << (useTierwiseOpenModSearch ? "true" : "false") << "\"/>" << endl;

    fout << "<parameter name=\"peak_scaling_mz_power\" value=\"" << peakScalingMzPower << "\"/>" << endl;
    fout << "<parameter name=\"peak_scaling_intensity_power\" value=\"" << peakScalingIntensityPower << "\"/>" << endl;

    fout << "<parameter name=\"peak_scaling_unassigned_peaks\" value=\"" << peakScalingUnassignedPeaks << "\"/>"
         << endl;

    fout << "<parameter name=\"peak_binning_num_bins_per_mz\" value=\"" << peakBinningNumBinsPerMzUnit << "\"/>"
         << endl;

    fout << "<parameter name=\"peak_binning_fraction_to_neighbor\" value=\"" << peakBinningFractionToNeighbor << "\"/>"
         << endl;

    fout << "<parameter name=\"peak_no_binning\" value=\"" << (peakNoBinning ? "true" : "false") << "\"/>" << endl;

    fout << "<parameter name=\"filter_min_peak_count\" value=\"" << filterMinPeakCount << "\"/>" << endl;

    fout << "<parameter name=\"filter_min_mz_range\" value=\"" << filterMinMzRange << "\"/>" << endl;

    fout << "<parameter name=\"filter_all_peaks_below_mz\" value=\"" << filterAllPeaksBelowMz << "\"/>" << endl;

    fout << "<parameter name=\"filter_max_peaks_used\" value=\"" << filterMaxPeaksUsed << "\"/>" << endl;

    fout << "<parameter name=\"filter_max_dynamic_range\" value=\"" << filterMaxDynamicRange << "\"/>" << endl;

    fout << "<parameter name=\"filter_max_intensity_below\" value=\"" << filterMaxIntensityBelow << "\"/>" << endl;

    fout << "<parameter name=\"filter_lib_max_peaks_used\" value=\"" << filterLibMaxPeaksUsed << "\"/>" << endl;

    fout << "<parameter name=\"filter_iTRAQ_reporter_peaks\" value=\"" << (filterITRAQReporterPeaks ? "true" : "false")
         << "\"/>" << endl;

    fout << "<parameter name=\"filter_TMT_reporter_peaks\" value=\"" << (filterTMTReporterPeaks ? "true" : "false")
         << "\"/>" << endl;

    fout << "<parameter name=\"filter_light_ions_mz_threshold\" value=\"" << filterLightIonsMzThreshold << "\"/>"
         << endl;

    fout << "<parameter name=\"fval_fraction_delta\" value=\"" << fvalFractionDelta << "\"/>" << endl;

    fout << "<parameter name=\"fval_use_dot_bias\" value=\"" << (fvalUseDotBias ? "true" : "false") << "\"/>" << endl;

}

void SpectraSTSearchParams::printAdvancedOptions(ostream &out) {

    out << "Spectrast (version " << SPECTRAST_VERSION << "." << SPECTRAST_SUB_VERSION << ", " << szTPPVersionInfo
        << ") by Henry Lam." << endl;
    out << endl;

    out << "SEARCH MODE ADVANCED OPTIONS" << endl;
    out << endl;

    out << "         CANDIDATE SELECTION AND SCORING OPTIONS" << endl;
    out << "         -s_HOM<rank>    Detect homologous lower hits up to <rank>." << endl;
    out << "                           Looks for lower hits homologous to the first one and adjust delta accordingly."
        << endl;
    out << "                           No detection will be done if <rank> < 2." << endl;
//  out << "         -s_NO1          Ignore all +1 spectra in the library. (Turn off with -sNO1!)" << endl; 
//  out << "         -s_NOS          Ignore all spectra which have non-Normal status. (Turn off with -s_NOS!)" << endl; 
    out << "         -s_FDL<frac>    Specify fraction of f-value that is delta/dot. (The rest is dot.)" << endl;
    out
            << "         -s_FBI          Use dot bias to penalize high-scoring matches with massive noise and/or dominant peak. (Turn off with -s_FBI!)"
            << endl;
    out
            << "         -s_SP4          Use old (SpectraST 4.0) scoring -- square-rooted intensity dot product with dot bias."
            << endl;
    out
            << "         -s_PVL          Compute P-value by fitting score distribution of lower hits, and use it for scoring. (Turn off with -s_PVL!)"
            << endl;
    out
            << "                           NOTE: Only applicable to new (SpectraST 5.0) scoring. Tested for low-resolution CID spectra only."
            << endl;
    out
            << "         -s_OMT          Perform tier-wise open modification search for modifications within precursor m/z window. (Turn off with -s_OMT!)"
            << endl;
    out << endl;

    out << "         OUTPUT AND DISPLAY OPTIONS" << endl;
    out << "         -s_FV1<thres>   Minimum F value threshold for the top hit. " << endl;
    out
            << "                           Only search results with top hits having F value greater than <thres> will be displayed."
            << endl;
    out << "         -s_FV2<thres>   Minimum F value threshold for the lower hits. " << endl;
    out
            << "         -s_SHH          Always displays homologous lower hits (only valid if -s_HOM is on. Turn off with -s_SHH!). "
            << endl;
    out << "                           Lower hits having F value greater than <thres> will also be displayed." << endl;
    out << "         -s_SHR<rank>    Maximum rank of hits shown (e.g. -s_SHR3 will show the top 3 hits). " << endl;
    out << "         -s_SHM          Exclude NO_MATCH's. (Turn off with -s_SHM!)" << endl;
    out
            << "                           Do not display queries for which there is no match above the dot product threshold."
            << endl;
    out << "         -s_ENZ<enzyme>  Specify proteolytic enzyme used for the purpose of pepXML report. " << endl;
    out << "                           <enzyme> can be trypsin, lysc, lysn, chymotrypsin, etc. " << endl;
    out
            << "                           NOTE: This does not affect searching; it only affects how the results are written to pepXML."
            << endl;
    out << "         -s_FIN<file>    Print a text file of name <file> summarizing fingerprinting results." << endl;
    out << endl;

    out << "         SPECTRUM FILTERING OPTIONS" << endl;
    out
            << "         -s_XNP<thres>   Discard query spectra with fewer than <thres> peaks above threshold set with -s_CNT."
            << endl;
    out << "         -s_XMZ<m/z>     Discard query spectra with (almost) no peaks above a certain m/z value." << endl;
    out << "                           All query spectra with 95%+ intensity below <m/z> will be removed." << endl;
    out << "         -s_XIN<inten>   Discard query spectra with no peaks with intensity above <inten>." << endl;
    out << "         -s_XMR<range>   Discard query spectra with m/z range narrower than <range>." << endl;
    out << "         -s_CNT<thres>   Minimum peak intensity for peaks to be counted in query spectra." << endl;
    out << "                           Only peaks in query spectra with intensity above <thres> will be counted "
        << endl;
    out << "                           to meet the requirement for minimum number of peaks." << endl;

    out << "         SPECTRUM PROCESSING OPTIONS" << endl;
    out << "         -s_RNT<thres>   Noise peak threshold for query spectra." << endl;
    out << "                           All peaks in query spectra with intensities below <thres> will be zeroed."
        << endl;
    // out << "         -s_R51<thres>   Remove dominant peak at 515.3 Th" << endl;
    // out << "                           All dominant peaks in query spectra near 515.3 Th (with intensity greater than <thres> of the total " << endl;
    // out << "                           intensity of the spectrum) will be zeroed." << endl;
    // out << "                           NOTE: Dominant 515.3 Th peaks are a common impurity artifact in cleavable ICAT experiments." << endl;
    out << "         -s_RNP<num>     Remove all but the top <num> peaks in query spectra." << endl;
    out << "                           NOTE: In SP5, the query spectra will be simplified to the same number" << endl;
    out << "                                 of peaks as the library spectra, specified by -s_LNP." << endl;
    out << "         -s_RDR<num>     Remove all peaks smaller than 1/<num> of the base (highest) peak in query spectra."
        << endl;
    out << "         -s_MZS<mzpow>, " << endl;
    out << "         -s_INS<intpow>  Intensity scaling powers." << endl;
    out << "                           All peaks are scaled to (m/z)^<mzpow> * (intensity)^<intpow>." << endl;
    out << "         -s_UAS<factor>  Scaling factor for unassigned peaks in library spectra." << endl;
    out << "                           Unassigned peaks in the library spectra will be scaled by <factor>." << endl;
    out << "         -s_BIN<num>     Number of bins per Th." << endl;
    out << "         -s_NOB          Disable binning and instead perform peak-to-peak matching. (Turn off with -s_NOB!)"
        << endl;
    out
            << "                           Only applicable for SpectraST 5.0. Specify fragment m/z tolerance by -s_BIN<num> (Tolerance = 1/<num>)."
            << endl;
    out << "         -s_NEI<frac>    Fraction of the scaled intensity assigned to neighboring bins." << endl;
    out
            << "         -s_LNP<num>     Remove all but the top <num> peaks in the LIBRARY spectra. (For best results, <num> should be at least 50.)"
            << endl;
    out << "         -s_RLI<thres>   Remove all light ions with m/z lower than <thres> Th for both library and query."
        << endl;
    out << "         -s_ITQ          Remove iTRAQ reporter peaks in the range 112-122 Th. (Turn off with -s_ITQ!)"
        << endl;
    out << "         -s_TMT          Remove TMT reporter peaks in the range 126-132 Th. (Turn off with -s_TMT!)"
        << endl;

    out << endl;

}


// printOptions - prints to console the usage of the options
void SpectraSTSearchParams::printUsage(ostream &out) {

    out << "(II) SEARCH MODE " << endl;
    out << "Usage: spectrast [ options ] <SearchFileName1> [ <SearchFileName2> ... <SearchFileNameN> ]" << endl;
    out << "where: SearchFileNameX = Name(s) of file containing unknown spectra to be searched." << endl;
    out << "                         (Extension specifies format of file. Supports .mzXML, .mzData, .dta and .msp)"
        << endl;

    out << endl;
    out << "Options: GENERAL OPTIONS" << endl;
    out << "         -sF<file>    Read search options from file. " << endl;
    out << "                           If <file> is not given, \"spectrast.params\" is assumed." << endl;
    out
            << "                           NOTE: All options set in the file will be overridden by command-line options, if specified."
            << endl;
    out << "         -sL<file>    Specify library file." << endl;
    out
            << "                           <file> must have .splib extension. The existence of the corresponding .spidx file of the same name"
            << endl;
    out << "                           in the same directory is assumed." << endl;
    out << "         -sD<file>    Specify a sequence database file." << endl;
    out << "                           <file> must be in .fasta format. This will not affect the search in any way,"
        << endl;
    out
            << "                           but this information will be included in the output for any downstream data processing."
            << endl;
    out << "         -sT<type>    Specify the type of the sequence database file." << endl;
    out << "                           <type> must be either \"AA\" or \"DNA\"." << endl;
    out << "         -sR          Cache all entries in RAM. (Turn off with -sR!)" << endl;
    out
            << "                           Requires a lot of memory (the library will usually be loaded almost in its entirety), but speeds up search for unsorted queries."
            << endl;
    out << "         -sS<file>    Only search a subset of the query spectra in the search file." << endl;
    out << "                           Only query spectra with names matching a line of <file> will be searched."
        << endl;

    out << endl;

    out << "         CANDIDATE SELECTION AND SCORING OPTIONS" << endl;
    out << "         -sM<tol>     Specify precursor m/z tolerance. " << endl;
    out
            << "                           <tol> is the m/z tolerance within which candidate entries are compared to the query."
            << endl;
    out << "         -sA          Use isotopically averaged mass instead of monoisotopic mass. (Turn of with -sA!)"
        << endl;
//  out << "         -sC<type>    Specify the expected kind of cysteine modification for the query spectra." << endl;
//  out << "                           <type> must be \"ICAT_cl\" for cleavable ICAT, \"ICAT_uc\" for uncleavable ICAT, or \"CAM\" for CarbAmidoMethyl." << endl;
//  out << "                           Those library spectra with a different kind of cysteine modification will be ignored." << endl;
//  out << "                           The ICAT type, if any, will also be included in the pepXML output for validation by PeptideProphet." << endl;
//  out << "         -sc          Ignore any library spectra of peptide ions with an unmodified cysteine. (Turn off with -sc!)" << endl;
//  out << "                           A suitable co-option when -sC is non-empty. Typically one expects almost all cysteines to be modified" << endl;
//  out << "                           if an ICAT agent or an alkylation agent such as iodoacetamide is used." << endl; 
    out
            << "         -sz          Search library spectra of all charge states, i.e., ignore specified charge state (if any) of the query spectrum."
            << endl;
    out << endl;

    out << "         OUTPUT AND DISPLAY OPTIONS" << endl;
    out << "         -sE<ext>     Output format." << endl;
    out
            << "                           The search result will be written to a file with the same base name as the search file, with extension <ext>"
            << endl;
    out << "                           <ext> = txt : Fixed-width text format." << endl;
    out << "                           <ext> = xls : Tab-delimited text format." << endl;
    out << "                           <ext> = pep.xml or xml or pepXML : PepXML format." << endl;
    out << "                           <ext> = html : HTML format." << endl;
    out << "         -sO<dir>     Write output file to <dir> instead of the same directory as the query data file."
        << endl;
    out << endl;

    out << "         OTHER ADVANCED OPTIONS" << endl;
    out << "         Type \"spectrast -s_\" for a full list of advanced (and obscure and not-so-useful) options."
        << endl;
    out << endl;

    out << endl;

}

