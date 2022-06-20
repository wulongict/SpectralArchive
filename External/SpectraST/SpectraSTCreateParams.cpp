#include "SpectraSTCreateParams.hpp"
#include "SpectraSTLog.hpp"
#include "SpectraSTConstants.hpp"
#include "FileUtils.hpp"
#include <fstream>
#include <iostream>
#include <sstream>
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

Institute for Systems Biology, hereby disclaims all copyright interest 
in Spectrast written by Henry Lam.

*/

/* Class: SpectraSTCreateParams
 * 
 * Manages the create parameters.
 *  
 */

using namespace std;

extern bool g_quiet;
extern SpectraSTLog *g_log;

// constructor
SpectraSTCreateParams::SpectraSTCreateParams() :
        m_options() {
    setDefault();

}

// copy constructor
SpectraSTCreateParams::SpectraSTCreateParams(SpectraSTCreateParams &s) {
    (*this) = s;

}

// assignment operator
SpectraSTCreateParams &SpectraSTCreateParams::operator=(SpectraSTCreateParams &s) {


    this->paramsFileName = s.paramsFileName;
    this->outputFileName = s.outputFileName;
    this->binaryFormat = s.binaryFormat;
    this->reannotatePeaks = s.reannotatePeaks;
    this->remark = s.remark;
    this->writeDtaFiles = s.writeDtaFiles;
    this->writeMgfFile = s.writeMgfFile;
    this->writePAIdent = s.writePAIdent;
    this->filterCriteria = s.filterCriteria;
    this->useProbTable = s.useProbTable;
    this->useProteinList = s.useProteinList;
    this->printMRMTable = s.printMRMTable;
    this->minimumMRMQ3MZ = s.minimumMRMQ3MZ;
    this->maximumMRMQ3MZ = s.maximumMRMQ3MZ;
    this->removeDecoyProteins = s.removeDecoyProteins;

    this->minimumProbabilityToInclude = s.minimumProbabilityToInclude;
    this->maximumFDRToInclude = s.maximumFDRToInclude;
    this->rawSpectraNoiseThreshold = s.rawSpectraNoiseThreshold;
    this->rawSpectraMaxDynamicRange = s.rawSpectraMaxDynamicRange;
    this->setDeamidatedNXST = s.setDeamidatedNXST;
    this->datasetName = s.datasetName;
    this->addMzXMLFileToDatasetName = s.addMzXMLFileToDatasetName;
    this->minimumDeltaCnToInclude = s.minimumDeltaCnToInclude;
    this->minimumNumAAToInclude = s.minimumNumAAToInclude;
    this->minimumNumPeaksToInclude = s.minimumNumPeaksToInclude;
    this->maximumMassDiffToInclude = s.maximumMassDiffToInclude;
    this->centroidPeaks = s.centroidPeaks;
    this->setFragmentation = s.setFragmentation;
    this->skipRawAnnotation = s.skipRawAnnotation;
    this->evaluatePhosphoSiteAssignment = s.evaluatePhosphoSiteAssignment;
    this->keepRawIntensities = s.keepRawIntensities;
    this->bracketSpectra = s.bracketSpectra;
    this->mergeBracket = s.mergeBracket;

    this->combineAction = s.combineAction;
    this->buildAction = s.buildAction;
    this->plotSpectra = s.plotSpectra;
    this->reduceSpectrum = s.reduceSpectrum;


    this->minimumNumReplicates = s.minimumNumReplicates;
    this->peakQuorum = s.peakQuorum;
    this->maximumNumReplicates = s.maximumNumReplicates;
    this->maximumNumPeaksUsed = s.maximumNumPeaksUsed;
    this->removeDissimilarReplicates = s.removeDissimilarReplicates;
    this->replicateWeight = s.replicateWeight;
    this->maximumNumPeaksKept = s.maximumNumPeaksKept;
    this->recordRawSpectra = s.recordRawSpectra;

    this->useBayesianDenoiser = s.useBayesianDenoiser;
    this->trainBayesianDenoiser = s.trainBayesianDenoiser;
    this->denoiserMinimumSignalProb = s.denoiserMinimumSignalProb;
    this->denoiserParamFile = s.denoiserParamFile;


    this->qualityLevelRemove = s.qualityLevelRemove;
    this->qualityLevelMark = s.qualityLevelMark;
    this->qualityPenalizeSingletons = s.qualityPenalizeSingletons;
    this->qualityImmuneProbThreshold = s.qualityImmuneProbThreshold;
    this->qualityImmuneMultipleEngines = s.qualityImmuneMultipleEngines;

    this->decoyConcatenate = s.decoyConcatenate;
    this->decoySizeRatio = s.decoySizeRatio;
    this->decoyPrecursorSwap = s.decoyPrecursorSwap;

    this->normalizeRTWithLandmarks = s.normalizeRTWithLandmarks;
    this->normalizeRTLinearRegression = s.normalizeRTLinearRegression;

    this->allowableModTokens = s.allowableModTokens;

    this->refreshDatabase = s.refreshDatabase;
    this->refreshDeleteUnmapped = s.refreshDeleteUnmapped;
    this->refreshDeleteMultimapped = s.refreshDeleteMultimapped;
    this->refreshTrypticOnly = s.refreshTrypticOnly;

    this->unidentifiedClusterIndividualRun = s.unidentifiedClusterIndividualRun;
    this->unidentifiedClusterMinimumDot = s.unidentifiedClusterMinimumDot;
    this->unidentifiedRemoveSinglyCharged = s.unidentifiedRemoveSinglyCharged;
    this->unidentifiedMinimumNumPeaksToInclude = s.unidentifiedMinimumNumPeaksToInclude;
    this->unidentifiedSingletonXreaThreshold = s.unidentifiedSingletonXreaThreshold;

    this->m_options.clear();
    for (vector<string>::iterator i = s.m_options.begin(); i != s.m_options.end(); i++) {
        this->m_options.push_back(*i);
    }
    return (*this);

}

// destructor
SpectraSTCreateParams::~SpectraSTCreateParams() {

}

// finalizeOptions - go through the options stored in m_options and actually set them.
void SpectraSTCreateParams::finalizeOptions() {


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

            // GENERAL
            case 'F' :
                if (optionValue.empty()) {
                    // no params file name specified. use default
                    paramsFileName = DEFAULT_CREATE_PARAMS_FILE;
                    readFromFile();
                    valid = true;
                } else {
                    // fixpath(optionValue);
                    paramsFileName = optionValue;
                    readFromFile();
                    valid = true;
                }

                if (!g_quiet) {
                    cout << "Create parameter file loaded: \"" << paramsFileName << "\"." << endl;
                }
                break;

            case 'N' :
                if (!optionValue.empty()) {
                    // fixpath(optionValue);
                    outputFileName = optionValue;
                    valid = true;
                }
                break;

            case 'm' :
                if (!optionValue.empty()) {
                    remark = optionValue;
                    valid = true;
                }
                break;


            case 'f' :
                if (!optionValue.empty()) {
                    filterCriteria = optionValue;
                    valid = true;
                }
                break;

            case 'T' :
                if (!optionValue.empty()) {
                    // fixpath(optionValue);
                    useProbTable = optionValue;
                    makeFullPath(useProbTable);
                    valid = true;
                } else {
                    valid = false;
                }
                break;

            case 'O' :
                if (!optionValue.empty()) {
                    // fixpath(optionValue);
                    useProteinList = optionValue;
                    makeFullPath(useProteinList);
                    valid = true;
                } else {
                    valid = false;
                }
                break;

            case 'M' :
                if (optionValue.empty()) {
                    printMRMTable = "DEFAULT";
                    valid = true;
                } else if (optionValue == "!") {
                    printMRMTable = "";
                    valid = true;
                } else {
                    printMRMTable = optionValue;
                    valid = true;
                }
                break;

                // PEPXML
            case 'n' :
                if (!optionValue.empty()) {
                    datasetName = optionValue;
                    valid = true;
                }
                break;

            case 'o' :
                if (optionValue.empty()) {
                    addMzXMLFileToDatasetName = true;
                    valid = true;
                } else if (optionValue == "!") {
                    addMzXMLFileToDatasetName = false;
                    valid = true;
                }
                break;

            case 'P' :
                f = atof(optionValue.c_str());
                if (f >= 0.0 && f <= 1.000001) {
                    minimumProbabilityToInclude = f;
                    valid = true;
                }
                break;

            case 'q' :
                f = atof(optionValue.c_str());
                if (f >= 0.0 && f <= 1.0) {
                    maximumFDRToInclude = f;
                    valid = true;
                }
                break;

            case 'g' :
                if (optionValue.empty()) {
                    setDeamidatedNXST = true;
                    valid = true;
                } else if (optionValue == "!") {
                    setDeamidatedNXST = false;
                    valid = true;
                }
                break;

            case 'I' :
                if (!optionValue.empty()) {
                    setFragmentation = optionValue;
                    valid = true;
                }
                break;

                // LIBRARY MANIPULATION

            case 'J' :
                if (optionValue[0] == 'U') {
                    combineAction = "UNION";
                    valid = true;
                } else if (optionValue[0] == 'I') {
                    combineAction = "INTERSECT";
                    valid = true;
                } else if (optionValue[0] == 'S') {
                    combineAction = "SUBTRACT";
                    valid = true;
                } else if (optionValue[0] == 'H') {
                    combineAction = "SUBTRACT_HOMOLOGS";
                    valid = true;
                } else if (optionValue[0] == 'A') {
                    combineAction = "APPEND";
                    valid = true;
                }
                break;

            case 'A' :
                if (optionValue[0] == 'B') {
                    buildAction = "BEST_REPLICATE";
                    valid = true;
                } else if (optionValue[0] == 'C') {
                    buildAction = "CONSENSUS";
                    valid = true;
                } else if (optionValue[0] == 'Q') {
                    buildAction = "QUALITY_FILTER";
                    valid = true;
                } else if (optionValue[0] == 'D') {
                    buildAction = "DECOY";
                    valid = true;
                } else if (optionValue[0] == 'N') {
                    buildAction = "SORT_BY_NREPS";
                    valid = true;
                } else if (optionValue[0] == 'M') {
                    buildAction = "USER_SPECIFIED_MODS";
                    valid = true;
                } else if (optionValue[0] == 'S') {
                    buildAction = "SIMILARITY_CLUSTERING";
                    valid = true;
                } else if (optionValue[0] == 'E') {
                    buildAction = "SEMI_EMPIRICAL_SPLIB";
                    valid = true;
                }
                break;

            case 'Q' :
                k = atoi(optionValue.c_str());
                if (k >= 1) {
                    reduceSpectrum = k;
                    valid = true;
                }
                break;



                // CONSENSUS
            case 'r' :
                k = atoi(optionValue.c_str());
                if (k > 0) {
                    minimumNumReplicates = k;
                    valid = true;
                }
                break;


                // QUALITY FILTER
            case 'L' :
                k = atoi(optionValue.c_str());
                if (k >= 0 && k <= 5) {
                    qualityLevelRemove = k;
                    valid = true;
                }
                break;

            case 'l' :
                k = atoi(optionValue.c_str());
                if (k >= 0 && k <= 5) {
                    qualityLevelMark = k;
                    valid = true;
                }
                break;

                // DECOY
            case 'c' :
                if (optionValue.empty()) {
                    decoyConcatenate = true;
                    valid = true;
                } else if (optionValue == "!") {
                    decoyConcatenate = false;
                    valid = true;
                }
                break;

            case 'y' :
                k = atoi(optionValue.c_str());
                if (k >= 0) {
                    decoySizeRatio = k;
                    valid = true;
                }
                break;

                // REFRESH PROTEIN MAPPINGS
            case 'u' :
                if (optionValue.empty()) {
                    refreshDeleteUnmapped = true;
                    valid = true;
                } else if (optionValue == "!") {
                    refreshDeleteUnmapped = false;
                    valid = true;
                }
                break;

            case 'd' :
                if (optionValue.empty()) {
                    refreshDeleteMultimapped = true;
                    valid = true;
                } else if (optionValue == "!") {
                    refreshDeleteMultimapped = false;
                    valid = true;
                }
                break;

            case 'D' :
                if (!optionValue.empty()) {
                    // fixpath(optionValue);
                    refreshDatabase = optionValue;
                    makeFullPath(refreshDatabase);
                    valid = true;
                } else {
                    valid = false;
                }
                break;

            case 'x' :
                allowableModTokens = optionValue;
                valid = true;
                break;

            case '_' :
                finalizeAdvancedOption(optionValue);
                valid = true;
                break;

            case 'a' :
            case 'b' :
                //    case 'D' :
            case 't' :
            case 'v' :
                // case 'u' :
            case 's' :
            case 'j' :
            case 'z' :
            case 'w' :
            case 'k' :
            case 'p' :
                //    case 'd' :
            case '1' :
            case 'i' :
            case 'e' :
                if (!g_quiet) {
                    cout << "Option \"-c" << optionType << "\" is deprecated.  Ignored. ";
                    cout << "Type \"spectrast -c_\" to see if an equivalent advanced option is available." << endl;
                }
                valid = true;
                break;

            default :
                if (!g_quiet) cout << "Option \"-c" << optionType << "\" is undefined. Ignored. ";
                valid = true;
                break;

        }

        if (!valid && !g_quiet) {
            cout << "Invalid value for the \"-c" << optionType << "\" in the command line. Ignored." << endl;
        }

    }

}

void SpectraSTCreateParams::finalizeAdvancedOption(string option) {

    if (option.length() < 3) {
        cout << "Advanced option \"-c_" << option << " is undefined. Ignored." << endl;
        return;
    }

    double f = -1.0;
    int k = 0;
    bool valid = false;

    string optionType = option.substr(0, 3);
    string optionValue = option.substr(3);

    if (optionType == "ANN") {

        if (optionValue.empty()) {
            reannotatePeaks = true;
            valid = true;
        } else if (optionValue == "!") {
            reannotatePeaks = false;
            valid = true;
        }

    } else if (optionType == "BIN") {

        if (optionValue.empty()) {
            binaryFormat = true;
            valid = true;
        } else if (optionValue == "!") {
            binaryFormat = false;
            valid = true;
        }

    } else if (optionType == "DTA") {

        if (optionValue.empty()) {
            writeDtaFiles = true;
            valid = true;
        } else if (optionValue == "!") {
            writeDtaFiles = false;
            valid = true;
        }

    } else if (optionType == "MGF") {

        if (optionValue.empty()) {
            writeMgfFile = true;
            valid = true;
        } else if (optionValue == "!") {
            writeMgfFile = false;
            valid = false;
        }

    } else if (optionType == "PAI") {

        if (optionValue.empty()) {
            writePAIdent = true;
            valid = true;
        } else if (optionValue == "!") {
            writePAIdent = false;
            valid = false;
        }

    } else if (optionType == "Q3L") {

        if (!optionValue.empty()) {
            f = atof(optionValue.c_str());
            if (f >= 0.0) {
                minimumMRMQ3MZ = f;
                valid = true;
            }
        }

    } else if (optionType == "Q3H") {

        if (!optionValue.empty()) {
            f = atof(optionValue.c_str());
            if (f >= 0.0) {
                maximumMRMQ3MZ = f;
                valid = true;
            }
        }

    } else if (optionType == "RDY") {

        removeDecoyProteins = optionValue;
        valid = true;

    } else if (optionType == "RNT") {

        if (!optionValue.empty()) {
            f = atof(optionValue.c_str());
            if (f >= 0.0) {
                rawSpectraNoiseThreshold = f;
                valid = true;
            }
        }

    } else if (optionType == "RDR") {

        if (!optionValue.empty()) {
            f = atof(optionValue.c_str());
            if (f >= 1.0) {
                rawSpectraMaxDynamicRange = f;
                valid = true;
            }
        }

    } else if (optionType == "DCN") {

        if (!optionValue.empty()) {
            f = atof(optionValue.c_str());
            if (f >= 0.0) {
                minimumDeltaCnToInclude = f;
                valid = true;
            }
        }

    } else if (optionType == "NAA") {

        if (!optionValue.empty()) {
            k = atoi(optionValue.c_str());
            if (k >= 0) {
                minimumNumAAToInclude = k;
                valid = true;
            }
        }

    } else if (optionType == "NPK") {

        if (!optionValue.empty()) {
            k = atoi(optionValue.c_str());
            if (k >= 0) {
                minimumNumPeaksToInclude = k;
                valid = true;
            }
        }

    } else if (optionType == "MDF") {

        if (!optionValue.empty()) {
            f = atof(optionValue.c_str());
            if (f >= 0.0) {
                maximumMassDiffToInclude = f;
                valid = true;
            }
        }

    } else if (optionType == "RWI") {

        if (optionValue.empty()) {
            keepRawIntensities = true;
            valid = true;
        } else if (optionValue == "!") {
            keepRawIntensities = false;
            valid = true;
        }


    } else if (optionType == "BRK") {

        if (optionValue.empty()) {
            bracketSpectra = true;
            valid = true;
        } else if (optionValue == "!") {
            bracketSpectra = false;
            valid = true;
        }

    } else if (optionType == "BRM") {

        if (optionValue.empty()) {
            mergeBracket = true;
            valid = true;
        } else if (optionValue == "!") {
            mergeBracket = false;
            valid = true;
        }


    } else if (optionType == "PLT") {

        if (optionValue.empty()) {
            plotSpectra = "ALL";
            valid = true;
        } else if (optionValue == "!") {
            plotSpectra = "";
            valid = true;
        } else {
            plotSpectra = optionValue;
            valid = true;
        }

    } else if (optionType == "CEN") {

        if (optionValue.empty()) {
            centroidPeaks = true;
            valid = true;
        } else if (optionValue == "!") {
            centroidPeaks = false;
            valid = true;
        }

    } else if (optionType == "XAN") {

        if (optionValue.empty()) {
            skipRawAnnotation = true;
            valid = true;
        } else if (optionValue == "!") {
            skipRawAnnotation = false;
            valid = true;
        }

    } else if (optionType == "WGT") {

        if (!optionValue.empty()) {
            if (optionValue[0] == 'X') {
                replicateWeight = "XCORR";
                valid = true;
            } else if (optionValue[0] == 'P') {
                replicateWeight = "PROB";
                valid = true;
            } else if (optionValue[0] == 'S') {
                replicateWeight = "SN";
                valid = true;
            } else if (optionValue[0] == 'N') {
                replicateWeight = "NONE";
                valid = true;
            } else if (optionValue[0] == 'I') {
                replicateWeight = "INTP";
                valid = true;
            }
        }

    } else if (optionType == "DIS") {

        if (optionValue.empty()) {
            removeDissimilarReplicates = true;
            valid = true;
        } else if (optionValue == "!") {
            removeDissimilarReplicates = false;
            valid = true;
        }

    } else if (optionType == "QUO") {

        if (!optionValue.empty()) {
            f = atof(optionValue.c_str());
            if (k >= 0.0 && k <= 1.0) {
                peakQuorum = f;
                valid = true;
            }
        }

    } else if (optionType == "XPU") {

        if (!optionValue.empty()) {
            k = atoi(optionValue.c_str());
            if (k > 0) {
                maximumNumPeaksUsed = k;
                valid = true;
            }
        }

    } else if (optionType == "XNR") {

        if (!optionValue.empty()) {
            k = atoi(optionValue.c_str());
            if (k > 1) {
                maximumNumReplicates = k;
                valid = true;
            }
        }

    } else if (optionType == "XPK") {

        if (!optionValue.empty()) {
            k = atoi(optionValue.c_str());
            if (k >= 0) {
                maximumNumPeaksKept = k;
                valid = true;
            }
        }

    } else if (optionType == "RRS") {

        if (optionValue.empty()) {
            recordRawSpectra = true;
            valid = true;
        } else if (optionValue == "!") {
            recordRawSpectra = false;
            valid = true;
        }

    } else if (optionType == "BDU") {
        if (optionValue.empty()) {
            useBayesianDenoiser = true;
            valid = true;
        } else if (optionValue == "!") {
            useBayesianDenoiser = false;
            valid = true;
        }

    } else if (optionType == "BDT") {
        if (optionValue.empty()) {
            trainBayesianDenoiser = true;
            valid = true;
        } else if (optionValue == "!") {
            trainBayesianDenoiser = false;
            valid = true;
        }

    } else if (optionType == "BDP") {

        if (!optionValue.empty()) {
            f = atof(optionValue.c_str());
            if (f >= 0.0 && f <= 1.0) {
                denoiserMinimumSignalProb = f;
                valid = true;
            }
        }

    } else if (optionType == "BDF") {

        if (!optionValue.empty()) {
            // fixpath(optionValue);
            denoiserParamFile = optionValue;
            makeFullPath(denoiserParamFile);
        }
        valid = true;

    } else if (optionType == "QP1") {

        if (optionValue.empty()) {
            qualityPenalizeSingletons = true;
            valid = true;
        } else if (optionValue == "!") {
            qualityPenalizeSingletons = false;
            valid = true;
        }

    } else if (optionType == "QIP") {

        if (!optionValue.empty()) {
            f = atof(optionValue.c_str());
            if (f >= 0.0) {
                qualityImmuneProbThreshold = f;
                valid = true;
            }
        }

    } else if (optionType == "QIE") {

        if (optionValue.empty()) {
            qualityImmuneMultipleEngines = true;
            valid = true;
        } else if (optionValue == "!") {
            qualityImmuneMultipleEngines = false;
            valid = true;
        }

    } else if (optionType == "RTO") {

        if (optionValue.empty()) {
            refreshTrypticOnly = true;
            valid = true;
        } else if (optionValue == "!") {
            refreshTrypticOnly = false;
            valid = true;
        }

    } else if (optionType == "IRT") {

        if (!optionValue.empty()) {
            normalizeRTWithLandmarks = optionValue;
            makeFullPath(normalizeRTWithLandmarks);
            valid = true;
        } else {
            valid = false;
        }

    } else if (optionType == "IRR") {

        if (optionValue.empty()) {
            normalizeRTLinearRegression = true;
            valid = true;
        } else if (optionValue == "!") {
            normalizeRTLinearRegression = false;
            valid = true;
        }

    } else if (optionType == "PHO") {

        if (optionValue.empty()) {
            evaluatePhosphoSiteAssignment = true;
            valid = true;
        } else if (optionValue == "!") {
            evaluatePhosphoSiteAssignment = false;
            valid = true;
        }

    } else if (optionType == "UCR") {

        if (optionValue.empty()) {
            unidentifiedClusterIndividualRun = true;
            valid = true;
        } else if (optionValue == "!") {
            unidentifiedClusterIndividualRun = false;
            valid = true;
        }

    } else if (optionType == "UCD") {

        if (!optionValue.empty()) {
            f = atof(optionValue.c_str());
            if (f >= 0.0 && f <= 1.0) {
                unidentifiedClusterMinimumDot = f;
                valid = true;
            }
        }

    } else if (optionType == "UX1") {

        if (optionValue.empty()) {
            unidentifiedRemoveSinglyCharged = true;
            valid = true;
        } else {
            unidentifiedRemoveSinglyCharged = false;
            valid = true;
        }

    } else if (optionType == "UNP") {

        if (!optionValue.empty()) {
            k = atoi(optionValue.c_str());
            if (k >= 1) {
                unidentifiedMinimumNumPeaksToInclude = k;
                valid = true;
            }
        }

    } else if (optionType == "USX") {

        if (!optionValue.empty()) {
            f = atof(optionValue.c_str());
            if (f >= 0.0 && f <= 1.0) {
                unidentifiedSingletonXreaThreshold = f;
                valid = true;
            }
        }

    } else if (optionType == "DPS") {

        if (optionValue.empty()) {
            decoyPrecursorSwap = true;
            valid = true;
        } else if (optionValue == "!") {
            decoyPrecursorSwap = false;
            valid = true;
        }


    } else if (optionType == "OTL") {
        predictionOrigTargetFileName = optionValue;
        valid = true;
    } else if (optionType == "PEP") {
        predictionTargetPeptidesFileName = optionValue;
        valid = true;
    } else if (optionType == "SNP") {
        allowableSNP = optionValue;
        valid = true;
    } else if (optionType == "PTM") {
        allowableMutatedModifications = optionValue;
        valid = true;
    } else if (optionType == "DST") {
        allowablePredictionDistance = atoi(optionValue.c_str());
        valid = true;


    } else {
        if (!g_quiet) cout << "Advanced option \"-c_" << option << " is undefined. Ignored." << endl;
        valid = true;
    }

    if (!valid && !g_quiet) {
        cout << "Invalid value for advanced option \"-c_" << option << ". Ignored." << endl;
    }


}

// addOption - add an option to the list of options
bool SpectraSTCreateParams::addOption(string option) {

    // to allow command-line override of options set by reading a params file (-cF option),
    // options are first read from the file, and then the command-line options are allowed
    // to overwrite them. Hence, we want to make sure the -cF option (if any) is the first
    // to be executed.

    // check if it's the -cF option, in which case, put as the first option

    if (option.length() == 0) {
        return (true);
    } else if (option[0] == 'F') {
        m_options.insert(m_options.begin(), option);
        return (true);
    } else {
        // else, put at the back
        if (isExpectingArg(option)) {
            return (false);
        } else {
            m_options.push_back(option);
            return (true);
        }
    }
}

bool SpectraSTCreateParams::isExpectingArg(string option) {

    return (option == "N" || option == "n" || option == "m" || option == "f" || option == "T" || option == "O" ||
            option == "D" || option == "I");
}

// setDefault - sets the defaults for all options
void SpectraSTCreateParams::setDefault() {


    // Defaults for all options - see SpectraSTCreateParams::printUsage() for more information

    // GENERAL
    paramsFileName = DEFAULT_CREATE_PARAMS_FILE;
    outputFileName = "";
    reannotatePeaks = false;
    binaryFormat = true;
    writeDtaFiles = false;
    writeMgfFile = false;
    writePAIdent = false;
    remark = "";
    useProbTable = "";
    useProteinList = "";
    printMRMTable = "";
    minimumMRMQ3MZ = 200.0;
    maximumMRMQ3MZ = 1400.0;
    removeDecoyProteins = "";
    setFragmentation = "";

    // PEPXML
    minimumProbabilityToInclude = 0.9;
    maximumFDRToInclude = 9999.0; // not using by default
    datasetName = "";
    addMzXMLFileToDatasetName = false;
    rawSpectraNoiseThreshold = 0.0;
    rawSpectraMaxDynamicRange = 100000.0;
    setDeamidatedNXST = false;
    minimumDeltaCnToInclude = 0.0;
    minimumNumAAToInclude = 6;
    minimumNumPeaksToInclude = 10;
    maximumMassDiffToInclude = 9999.0;
    centroidPeaks = false;
    skipRawAnnotation = false;
    evaluatePhosphoSiteAssignment = false;
    keepRawIntensities = false;
    bracketSpectra = false;
    mergeBracket = false;

    // LIBRARY MANIPULATION
    filterCriteria = "";
    combineAction = "UNION";
    buildAction = "";
    plotSpectra = "";
    reduceSpectrum = 0;

    // CONSENSUS
    minimumNumReplicates = 1;
    peakQuorum = 0.6;
    maximumNumPeaksUsed = 300;
    maximumNumReplicates = 100;
    removeDissimilarReplicates = true;
    replicateWeight = "SN";
    maximumNumPeaksKept = 150;
    recordRawSpectra = false;

    // DENOISER
    useBayesianDenoiser = false;
    trainBayesianDenoiser = false;
    denoiserMinimumSignalProb = 0.0;
    denoiserParamFile = "";

    // QUALITY FILTER
    qualityLevelRemove = 2;
    qualityLevelMark = 5;
    qualityPenalizeSingletons = true;
    qualityImmuneProbThreshold = 0.999;
    qualityImmuneMultipleEngines = true;

    // DECOY
    decoyConcatenate = false;
    decoySizeRatio = 1;
    decoyPrecursorSwap = false;

    // REFRESH PROTEIN MAPPINGS
    refreshDatabase = "";
    refreshDeleteUnmapped = false;
    refreshDeleteMultimapped = false;
    refreshTrypticOnly = false;

    // USER SPECIFIED MODIFICATIONS
    allowableModTokens = "";
    predictionOrigTargetFileName = "";
    predictionTargetPeptidesFileName = "";
    allowableSNP = "ALL";
    allowableMutatedModifications = "ALL";
    allowablePredictionDistance = 1;

    // UNIDENTIFIED LIBRARIES
    unidentifiedClusterIndividualRun = false;
    unidentifiedClusterMinimumDot = 0.7;
    unidentifiedRemoveSinglyCharged = true;
    unidentifiedMinimumNumPeaksToInclude = 35;
    unidentifiedSingletonXreaThreshold = 0.6;

    // new features



}

// readFromFile - loads the options from the params file
// the params file has the format
// <param name> = <param value> (the '=' can be replaced by ':' or just blank space)
void SpectraSTCreateParams::readFromFile() {

    if (paramsFileName.empty()) {
        return;
    }

    ifstream fin;
    if (!myFileOpen(fin, paramsFileName)) {
        g_log->error("CREATE", "Cannot open PARAMS file \"" + paramsFileName +
                               "\" for reading Create params. Using all defaults.");
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

        double f;
        int k;
        bool valid = false;

        if (param.empty()) continue;

        // GENERAL
        if (param == "outputFileName") {
            if (!value.empty()) {
                // fixpath(value);
                outputFileName = value;
                valid = true;
            }
        } else if (param == "reannotatePeaks") {
            reannotatePeaks = (value == "true");
            valid = true;
        } else if (param == "binaryFormat") {
            binaryFormat = (value == "true");
            valid = true;
        } else if (param == "remark") {
            if (!value.empty()) {
                remark = value;
                valid = true;
            }
        } else if (param == "writeDtaFiles") {
            writeDtaFiles = (value == "true");
            valid = true;
        } else if (param == "writeMgfFile") {
            writeMgfFile = (value == "true");
            valid = true;
        } else if (param == "writePAIdent") {
            writePAIdent = (value == "true");
            valid = true;
        } else if (param == "filterCriteria") {
            if (value[0] == '\"' && value[value.length() - 1] == '\"') {
                filterCriteria = value.substr(1, value.length() - 2);
                valid = true;
            } else if (value[0] == '\'' && value[value.length() - 1] == '\'') {
                filterCriteria = value.substr(1, value.length() - 2);
                valid = true;
            } else {
                valid = false;
            }

        } else if (param == "useProbTable") {
            if (!value.empty()) {
                // fixpath(value);
                useProbTable = value;
                valid = true;
            }
        } else if (param == "useProteinList") {
            if (!value.empty()) {
                // fixpath(value);
                useProteinList = value;
                valid = true;
            }
        } else if (param == "printMRMTable") {
            printMRMTable = value;
            valid = true;
        } else if (param == "minimumMRMQ3MZ") {
            if (!value.empty()) {
                f = atof(value.c_str());
                if (f >= 0.0) {
                    minimumMRMQ3MZ = f;
                    valid = true;
                }
            }
        } else if (param == "maximumMRMQ3MZ") {
            if (!value.empty()) {
                f = atof(value.c_str());
                if (f >= 0.0) {
                    maximumMRMQ3MZ = f;
                    valid = true;
                }
            }
        } else if (param == "removeDecoyProteins") {
            if (!value.empty()) {
                removeDecoyProteins = value;
                valid = true;
            }
        } else if (param == "setFragmentation") {
            if (!value.empty()) {
                setFragmentation = value;
                valid = true;
            }

            // PEPXML
        } else if (param == "minimumProbabilityToInclude") {
            if (!value.empty()) {
                f = atof(value.c_str());
                if (f >= 0.0 && f <= 1.0) {
                    minimumProbabilityToInclude = f;
                    valid = true;
                }
            }
        } else if (param == "maximumFDRToInclude") {
            if (!value.empty()) {
                f = atof(value.c_str());
                if (f >= 0.0 && f <= 1.0) {
                    maximumFDRToInclude = f;
                    valid = true;
                }
            }
        } else if (param == "rawSpectraNoiseThreshold") {
            if (!value.empty()) {
                f = atof(value.c_str());
                if (f >= 0.0) {
                    rawSpectraNoiseThreshold = f;
                    valid = true;
                }
            }
        } else if (param == "rawSpectraMaxDynamicRange") {
            if (!value.empty()) {
                f = atof(value.c_str());
                if (f >= 1.0) {
                    rawSpectraMaxDynamicRange = f;
                    valid = true;
                }
            }
        } else if (param == "datasetName") {
            if (!value.empty()) {
                datasetName = value;
                valid = true;
            }
        } else if (param == "setDeamidatedNXST") {
            setDeamidatedNXST = (value == "true");
            valid = true;
        } else if (param == "minimumDeltaCnToInclude") {
            if (!value.empty()) {
                f = atof(value.c_str());
                if (f >= 0.0) {
                    minimumDeltaCnToInclude = f;
                    valid = true;
                }
            }
        } else if (param == "minimumNumAAToInclude") {
            if (!value.empty()) {
                k = atoi(value.c_str());
                if (k >= 0) {
                    minimumNumAAToInclude = k;
                    valid = true;
                }
            }
        } else if (param == "minimumNumPeaksToInclude") {
            if (!value.empty()) {
                k = atoi(value.c_str());
                if (k >= 0) {
                    minimumNumPeaksToInclude = k;
                    valid = true;
                }
            }
        } else if (param == "maximumMassDiffToInclude") {
            if (!value.empty()) {
                f = atof(value.c_str());
                if (f >= 0.0) {
                    maximumMassDiffToInclude = f;
                    valid = true;
                }
            }
        } else if (param == "centroidPeaks") {
            centroidPeaks = (value == "true");
            valid = true;
        } else if (param == "skipRawAnnotation") {
            skipRawAnnotation = (value == "true");
            valid = true;
        } else if (param == "evaluatePhosphoSiteAssignment") {
            evaluatePhosphoSiteAssignment = (value == "true");
            valid = true;
        } else if (param == "keepRawIntensities") {
            keepRawIntensities = (value == "true");
            valid = true;
        } else if (param == "bracketSpectra") {
            bracketSpectra = (value == "true");
            valid = true;
        } else if (param == "mergeBracket") {
            mergeBracket = (value == "true");
            valid = true;



            // LIBRARY MANIPULATION

        } else if (param == "combineAction") {
            if (value == "UNION" || value == "INTERSECT" || value == "SUBTRACT" || value == "SUBTRACT_HOMOLOG" ||
                value == "APPEND") {
                combineAction = value;
                valid = true;
            }
        } else if (param == "buildAction") {
            if (value == "BEST_REPLICATE" || value == "CONSENSUS" || value == "QUALITY_FILTER" ||
                value == "DECOY" || value == "SORT_BY_NREPS" || value == "USER_SPECIFIED_MODS" ||
                value == "SIMILARITY_CLUSTERING") {
                buildAction = value;
                valid = true;
            }
        } else if (param == "plotSpectra") {
            plotSpectra = value;
            valid = true;
        } else if (param == "reduceSpectrum") {
            if (!value.empty()) {
                k = atoi(value.c_str());
                if (k >= 0) {
                    reduceSpectrum = k;
                    valid = true;
                }
            }


            // CONSENSUS
        } else if (param == "minimumNumReplicates") {
            if (!value.empty()) {
                k = atoi(value.c_str());
                if (k > 0) {
                    minimumNumReplicates = k;
                    valid = true;
                }
            }
        } else if (param == "peakQuorum") {
            if (!value.empty()) {
                f = atof(value.c_str());
                if (f >= 0.0 && f <= 1.0) {
                    peakQuorum = f;
                    valid = true;
                }
            }
        } else if (param == "maximumNumPeaksKept") {
            if (!value.empty()) {
                k = atoi(value.c_str());
                if (k >= 0) {
                    maximumNumPeaksKept = k;
                    valid = true;
                }
            }
        } else if (param == "replicateWeight") {
            if (value == "XCORR" || value == "PROB" || value == "SN" || value == "INTP" || value == "NONE") {
                replicateWeight = value;
                valid = true;
            }
        } else if (param == "removeDissimilarReplicates") {
            removeDissimilarReplicates = (value == "true");
            valid = true;
        } else if (param == "maximumNumReplicates") {
            if (!value.empty()) {
                k = atoi(value.c_str());
                if (k > 1) {
                    maximumNumReplicates = k;
                    valid = true;
                }
            }
        } else if (param == "maximumNumPeaksUsed") {
            if (!value.empty()) {
                k = atoi(value.c_str());
                if (k > 0) {
                    maximumNumPeaksUsed = k;
                    valid = true;
                }
            }
        } else if (param == "recordRawSpectra") {
            recordRawSpectra = (value == "true");
            valid = true;

            // DENOISER
        } else if (param == "useBayesianDenoiser") {
            useBayesianDenoiser = (value == "true");
            valid = true;
        } else if (param == "trainBayesianDenoiser") {
            trainBayesianDenoiser = (value == "true");
            valid = true;
        } else if (param == "denoiserMinimumSignalProb") {
            if (!value.empty()) {
                f = atof(value.c_str());
                if (f >= 0.0 && f <= 1.0) {
                    denoiserMinimumSignalProb = f;
                    valid = true;
                }
            }
        } else if (param == "denoiserParamFile") {
            if (!value.empty()) {
                // fixpath(value);
                denoiserParamFile = value;
                makeFullPath(denoiserParamFile);
            }
            valid = true;

            // QUALITY FILTER
        } else if (param == "qualityLevelRemove") {
            if (!value.empty()) {
                k = atoi(value.c_str());
                if (k >= 0 && k <= 5) {
                    qualityLevelRemove = k;
                    valid = true;
                }
            }
        } else if (param == "qualityLevelMark") {
            if (!value.empty()) {
                k = atoi(value.c_str());
                if (k >= 0 && k <= 5) {
                    qualityLevelMark = k;
                    valid = true;
                }
            }
        } else if (param == "qualityPenalizeSingletons") {
            qualityPenalizeSingletons = (value == "true");
            valid = true;
        } else if (param == "qualityImmuneProbThreshold") {
            if (!value.empty()) {
                f = atof(value.c_str());
                if (f >= 0.0) {
                    qualityImmuneProbThreshold = f;
                    valid = true;
                }
            }
        } else if (param == "qualityImmuneMultipleEngines") {
            qualityImmuneMultipleEngines = (value == "true");
            valid = true;

            // DECOY
        } else if (param == "decoyConcatenate") {
            decoyConcatenate = (value == "true");
            valid = true;
        } else if (param == "decoySizeRatio") {
            if (!value.empty()) {
                k = atoi(value.c_str());
                if (k >= 0) {
                    decoySizeRatio = k;
                    valid = true;
                }
            }
        } else if (param == "decoyPrecursorSwap") {
            decoyPrecursorSwap = (value == "true");
            valid = true;


            // REFRESH PROTEIN MAPPINGS
        } else if (param == "refreshDatabase") {
            if (!value.empty()) {
                // fixpath(value);
                refreshDatabase = value;
                valid = true;
            }
        } else if (param == "refreshDeleteUnmapped") {
            refreshDeleteUnmapped = (value == "true");
            valid = true;
        } else if (param == "refreshDeleteMultimapped") {
            refreshDeleteMultimapped = (value == "true");
            valid = true;
        } else if (param == "refreshTrypticOnly") {
            refreshTrypticOnly = (value == "true");
            valid = true;

            // USER SPECIFIED MODIFICATIONS
        } else if (param == "allowableModTokens") {
            allowableModTokens = value;
            valid = true;

            // SIMILARITY CLUSTERING
        } else if (param == "unidentifiedClusterIndividualRun") {
            unidentifiedClusterIndividualRun = (value == "true");
            valid = true;
        } else if (param == "unidentifiedClusterMinimumDot") {
            if (!value.empty()) {
                f = atof(value.c_str());
                if (f >= 0.0 && f <= 1.0) {
                    unidentifiedClusterMinimumDot = f;
                    valid = true;
                }
            }
        } else if (param == "unidentifiedRemoveSinglyCharged") {
            unidentifiedRemoveSinglyCharged = (value == "true");
            valid = true;
        } else if (param == "unidentifiedMinimumNumPeaksToInclude") {
            if (!value.empty()) {
                k = atoi(value.c_str());
                if (k >= 0) {
                    unidentifiedMinimumNumPeaksToInclude = k;
                    valid = true;
                }
            }
        } else if (param == "unidentifiedSingletonXreaThreshold") {
            if (!value.empty()) {
                f = atof(value.c_str());
                if (f >= 0.0 && f <= 1.0) {
                    unidentifiedSingletonXreaThreshold = f;
                    valid = true;
                }
            }


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

void SpectraSTCreateParams::printUsage(ostream &out) {

    out << "(I) CREATE MODE " << endl;
    out << "Usage: spectrast [ options ] <FileName1> [ <FileName2> ... <FileNameN> ]" << endl;
    out << "where: FileNameX = Name of file containing spectra from which library is to be created." << endl;
    out
            << "                          Extension specifies format of file. Supports .msp, .hlf, .pepXML (or .pep.xml or .xml), .ms2, and .splib."
            << endl;
    out << endl;

    out << "Options: GENERAL OPTIONS" << endl;
    out << "         -cF<file>    Read create options from file <file>. " << endl;
    out << "                           If <file> is not given, \"spectrast_create.params\" is assumed." << endl;
    out
            << "                           NOTE: All options set in the file will be overridden by command-line options, if specified."
            << endl;
    out << "         -cN<name>    Specify output file name for .splib, .spidx and .pepidx files. " << endl;
    out << "         -cm<remark>  Remark. Add a Remark=<remark> comment to all library entries created. " << endl;
    out
            << "         -cM<format>  Write all library spectra as MRM transition tables. Leave <format> blank for default. (Turn off with -cM!) "
            << endl;
    out
            << "         -cT<file>    Use probability table in <file>. Only those peptide ions included in the table will be imported. "
            << endl;
    out
            << "                           A probability table is a text file with one peptide ion in the format AC[160]DEFGHIK/2 per line. "
            << endl;
    out
            << "                           If a probability is supplied following the peptide ion separated by a tab, it will be used to replace the original probability of that library entry."
            << endl;
    out
            << "         -cO<file>    Use protein list in <file>. Only those peptide ions associated with proteins in the list will be imported. "
            << endl;
    out << "                           A protein list is a text file with one protein identifier per line. " << endl;
    out
            << "                           If a number X is supplied following the protein separated by a tab, then at most X peptide ions associated with that protein will be imported."
            << endl;

    out << endl;
    out << "         PEPXML IMPORT OPTIONS (Applicable with .pepXML files)" << endl;
    out << "         -cP<prob>    Include all spectra identified with probability no less than <prob> in the library."
        << endl;
    out
            << "         -cq<fdr>     (Only PepXML import) Only include spectra with global FDR no greater than <fdr> in the library."
            << endl;
    out << "         -cn<name>    Specify a dataset identifier for the file to be imported." << endl;
    out
            << "         -co          Add the originating mzXML file name to the dataset identifier. Good for keeping track of in which"
            << endl;
    out << "                           MS run the peptide is observed. (Turn off with -co!)" << endl;
    out
            << "         -cg          Set all asparagines (N) in the motif NX(S/T) as deamidated (N[115]). Use for glycocaptured peptides. (Turn off with -cg!)."
            << endl;
    out
            << "         -cI          Set the instrument and acquisition settings of the spectra (in case not specified in data files)."
            << endl;
    out
            << "                           Examples: -cICID, -cIETD, -cICID-QTOF, -cIHCD. The latter two are treated as high-mass accuracy spectra."
            << endl;

    out << endl;
    out << "         LIBRARY MANIPULATION OPTIONS (Applicable with .splib files)" << endl;
    out << "         -cf<pred>    Filter library. Keep only those entries satisfying the predicate <pred>. " << endl;
    out << "                           <pred> should be a C-style predicate in quotes. " << endl;
    out << "         -cJU         Union. Include all the peptide ions in all the files. " << endl;
    out << "         -cJI         Intersection. Only include peptide ions that are present in all the files. " << endl;
    out
            << "         -cJS         Subtraction. Only include peptide ions in the first file that are not present in any of the other files."
            << endl;
    out << "         -cJH         Subtraction of homologs. Only include peptide ions in the first file " << endl;
    out
            << "                           that do not have any homologs with same charge and similar m/z in any of the other files."
            << endl;
    out
            << "         -cJA         Appending. Each peptide ion is added from only one library: the first file in the argument list that contains that peptide ion."
            << endl;
    out
            << "                           Useful for keeping existing consensus spectra unchanged while adding only previously unseen peptide ions."
            << endl;
    out << "         -cAB         Best replicate. Pick the best replicate of each peptide ion. " << endl;
    out
            << "         -cAC         Consensus. Create the consensus spectrum of all replicate spectra of each peptide ion. "
            << endl;
    out << "         -cAQ         Quality filter. Apply quality filters to library." << endl;
    out
            << "                           IMPORTANT: Quality filter can only be applied on a SINGLE .splib file with no peptide ion represented by more than one spectrum."
            << endl;
    out << "         -cAD         Create artificial decoy spectra. " << endl;
    out
            << "         -cAN         Sort library entries by descending number of replicates used (tie-breaking by probability). "
            << endl;
    out
            << "         -cAM         Create semi-empirical spectra based on allowable modifications specified by -cx option. "
            << endl;
    // HIDDEN FOR NOW: out << "         -cAS         Cluster spectra by similarity and merge clusters into consensus spectra. " << endl;
    out << "         -cQ<num>     Produce reduced spectra of at most <num> peaks. Inactive with -cAQ and -cAD." << endl;
    out
            << "         -cD<file>    Refresh protein mappings of each library entry against the protein database <file> (Must be in .fasta format)."
            << endl;
    out
            << "         -cu          Delete entries whose peptide sequences do not map to any protein during refreshing with -cD option."
            << endl;
    out
            << "                           When off, unmapped entries will be marked with Protein=0/UNMAPPED but retained in library. (Turn off with -cu!)."
            << endl;
    out
            << "         -cd          Delete entries whose peptide sequences map to multiple proteins during refreshing with -cD option. (Turn off with -cd!)."
            << endl;

    out << endl;

    out << "         CONSENSUS/BEST-REPLICATE OPTIONS (Applicable with -cAC and -cAB options)" << endl;
    out << "         -cr<num>     Minimum number of replicates required for each library entry." << endl;
    out << "                           Peptide ions failing to have originated from enough replicates" << endl;
    out << "                           will be excluded from library when creating consensus/best-replicate library."
        << endl;
    out << endl;

    out << "         QUALITY FILTER OPTIONS (Applicable with -cAQ option)" << endl;
    out << "         -cr<num>     Replicate quorum. Its value affects behavior of quality filter (see below)." << endl;
    out << "         -cL<level>   Specify the stringency of the quality filter." << endl;
    out << "         -cl<level>        -cL specifies the level for removal, -cl specifies the level for marking."
        << endl;
    out << "                           <level> = 0: No filter." << endl;
    out << "                           <level> = 1: Remove/mark impure spectra." << endl;
    out
            << "                           <level> = 2: Also remove/mark spectra with a spectrally similar counterpart in the library that is better."
            << endl;
    out
            << "                           <level> = 3: Also remove/mark inquorate entries (defined with -cr) that share no peptide sub-sequences with any other entries in the library."
            << endl;
    out << "                           <level> = 4: Also remove/mark all singleton entries." << endl;
    out << "                           <level> = 5: Also remove/mark all inquorate entries (defined with -cr)." << endl;

    out << endl;
    out << "         DECOY CREATION OPTIONS (Applicable with -cAD option)" << endl;
    out << "         -cc          Concatenate real and decoy libraries. (Turn off with -cc!)" << endl;
    out << "         -cy<num>     Specify the (decoy / real) size ratio. Must be an integer." << endl;

    out << endl;
    out << "         SEMI-EMPIRICAL SPECTRUM CREATION OPTIONS (Applicable with -cAM option)" << endl;
    out
            << "         -cx<str>     Specify allowable modification tokens used to generate new peptide ions for which semi-empirical spectra are to be created."
            << endl;
    out
            << "                           (e.g. -cx\"K[136]\" for static mod of +6 heavy lysines, -cx\"MM[147]\" for variable mod of methionine oxidation.)"
            << endl;

    out << endl;
    out << "         OTHER ADVANCED OPTIONS" << endl;
    out << "         Type \"spectrast -c_\" for a full list of advanced (and obscure and not-so-useful) options."
        << endl;
    out << endl;

    out << endl;

}

void SpectraSTCreateParams::printAdvancedOptions(ostream &out) {

    out << "Spectrast (version " << SPECTRAST_VERSION << "." << SPECTRAST_SUB_VERSION << ", " << szTPPVersionInfo
        << ") by Henry Lam." << endl;
    out << endl;

    out << "CREATE MODE ADVANCED OPTIONS" << endl;
    out << endl;
    out << "GENERAL OPTIONS:" << endl;
    out << "         -c_BIN          Write library in binary format (Enables quicker search). (Turn off with -c_BIN!) "
        << endl;
    out << "                           A human-readable text-format library file will also be created." << endl;
    out << "         -c_DTA          Write all library spectra as .dta files. (Turn off with -c_DTA!) " << endl;
    out << "         -c_MGF          Write all library spectra as .mgf files. (Turn off with -c_MGF!) " << endl;
    out
            << "         -c_RDY<prefix>  Remove spectra of decoys, for which all proteins have names starting with <prefix>."
            << endl;
    out
            << "                           Also remove decoy proteins from Protein field for peptides mapped to both target and decoy proteins."
            << endl;

    out << "LIBRARY IMPORT OPTIONS (Applicable with .pep.xml, .tsv, .msp, .hlf, .ms2, .mz(X)ML)" << endl;
    out << "         -c_CEN          Centroid peaks." << endl;
    out << "         -c_RNT<thres>   Absolute noise filter. Filter out noise peaks with intensity below <thres>."
        << endl;
    out
            << "         -c_RDR<range>   Relative noise filter. Filter out noise peaks with intensity below 1/<range> of that of the highest peak."
            << endl;
    out << "         -c_NAA<num>     Exclude spectra with IDs of fewer than <num> amino acids." << endl;
    out << "         -c_NPK<num>     Exclude spectra with fewer than <num> peaks." << endl;
    out << "         -c_XAN          Skip (re-)annotation of the imported spectra." << endl;
    out << "         -c_DCN<thres>   (Only PepXML import) Exclude spectra with deltaCn smaller than <thres>." << endl;
    out << "                           Useful for excluding spectra with indiscriminate modification sites." << endl;
    out
            << "         -c_MDF<thres>   (Only PepXML import) Exclude spectra with precursor mass difference (absolute value) of over <thres> Daltons."
            << endl;
    out
            << "         -c_BRK          (Only PepXML import) Bracket import: for each confident ID, also search neighboring scan numbers for repeated scans to import. (Turn off with -c_BRK!)"
            << endl;
    out
            << "         -c_BRM          (Only PepXML import) Merge bracketed spectra: merge repeated scans of a bracket into one consensus spectrum for import. (Turn off with -c_BRM!)"
            << endl;

    out << "LIBRARY MANIPULATION OPTIONS (Applicable with .splib files)" << endl;
    out << "         -c_ANN          Re-annotate peaks (Turn off with -c_ANN!). " << endl;
    out << "         -c_Q3L          Specify the lower m/z limit for Q3 in MRM table generation." << endl;
    out << "         -c_Q3H          Specify the upper m/z limit for Q3 in MRM table generation. " << endl;
    out << "         -c_NAA<num>     Exclude spectra with IDs of fewer than <num> amino acids." << endl;
    out << "         -c_NPK<num>     Exclude spectra with fewer than <num> peaks." << endl;
    out
            << "         -c_RTO          With -cD option, only map peptide to protein when the peptide is tryptic in that particular protein. (Turn off with -c_RTO!)"
            << endl;

    out << "CONSENSUS/BEST-REPLICATE OPTIONS (Applicable with -cAC and -cAB options)" << endl;
    out
            << "         -c_DIS          Remove dissimilar replicates before creating consensus spectrum. (Turn off with -c_DIS!)"
            << endl;
    out << "         -c_QUO<frac>    Peak quorum: the fraction of all replicates required to contain a certain peak."
        << endl;
    out << "                           Peaks not present in enough replicates will be deleted." << endl;
    out << "         -c_XPU<num>     Maximum number of peaks in each replicate to be considered in creating consensus."
        << endl;
    out << "                           Only the top <num> peaks by intensity will be considered." << endl;
    out << "         -c_XNR<num>     Maximum number of replicates used to build consensus spectrum. " << endl;
    out
            << "                           Note: This is not an absolute hard cap for consensus spectrum created from consensus spectra."
            << endl;
    out << "         -c_XPK<num>     Maximum number of peaks kept in the final library spectrum. " << endl;
    out << "                           The most intense <num> peaks will be kept. <num> = 0 means keeping all peaks."
        << endl;
    out << "         -c_WGT<score>   Select the type of score to weigh and rank the replicates." << endl;
    out << "                           <score> = 'X': will use a function of the SEQUEST xcorr score as the weight."
        << endl;
    out << "                           <score> = 'S': will use a measure of signal-to-noise ratio as the weight."
        << endl;
    out
            << "                           <score> = 'P': will use a function of the PeptideProphet probability as the weight."
            << endl;
    out << "                           <score> = 'I': will use a function of the precursor intensity as the weight."
        << endl;
    out
            << "                           <score> = everything else: all replicates will be weighted equally and ranked randomly."
            << endl;
    out
            << "         -c_RRS          Record all raw spectra (in the format file.scan.scan) used to build the consensus in the Comment."
            << endl;

    out << "BAYESIAN DENOISER OPTIONS" << endl;
    out
            << "         -c_BDU          Use Bayesian denoiser. Default parameters are used unless trained on the fly with -c_BDT. (Turn off with -c_BDU!)"
            << endl;
    out << "         -c_BDT          Train Bayesian denoiser. Only active in consensus mode (-cAC)." << endl;
    out << "         -c_BDP<thres>   Minimum signal probability to retain a peak when denoiser is used." << endl;
    out << "         -c_BDF<file>    Specify parameter file for Bayesian denoiser for writing and reading." << endl;

    out << "QUALITY FILTER OPTIONS" << endl;
    out
            << "         -c_QP1          Apply stricter thresholds to singleton spectra during quality filters. (Turn off with -c_QP1!)"
            << endl;
    out << "         -c_QIP<thres>   Specify a probability above which library spectra are immune to quality filters."
        << endl;
    out
            << "                           Set <thres> to above 1.0 to subject all spectra to quality filters regardless of probability."
            << endl;
    out
            << "         -c_QIE          Make spectra identified by multiple sequence search engines immune to quality filters. (Turn off with -c_QIE!)"
            << endl;
    out << endl;

    out << "DECOY CREATION OPTIONS" << endl;
    out << "         -c_DPS          Use the precursor swap method for generating decoys. (Turn off with -c_DPS!)"
        << endl;

    out << "RETENTION TIME NORMALIZATION OPTIONS (Applicable with .pep.xml)" << endl;
    out << "         -c_IRT          Use landmark peptides in <file> to normalize retention times to iRT's." << endl;
    out
            << "         -c_IRR          Regress the real RTs of landmark peptides (i.e. assume they form a straight line). (Turn off with -c_IRR!)"
            << endl;

    out << "UNIDENTIFIED LIBRARY/CLUSTERING OPTIONS" << endl;
    out
            << "         -c_UCR          Cluster spectra in each run as they are imported from data files. (Turn off with -c_UCR!)"
            << endl;
    out << "         -c_UCD<thres>   Specify minimum dot products for two spectra to be clustered." << endl;
    out << "         -c_UX1          Remove spectra that appear to be singly charged. (Turn off with -c_UX1!)" << endl;
    out << "         -c_UNP<num>     Remove spectra that have fewer than <num> peaks." << endl;
    out << "         -c_USX<thres>   Apply an Xrea (quality measure) filter to singleton spectra after clustering."
        << endl;


}

string SpectraSTCreateParams::constructDescrStr(string fileList, string fileType) {


    stringstream ss;
    if (fileType == ".splib") {

        ss << "COMPILE FROM " << fileList;


        if (buildAction == "CONSENSUS") {
            ss << " [" << buildAction;
            ss << ";r=" << minimumNumReplicates;
            ss << ";I=" << setFragmentation;
            ss << ";_QUO=" << peakQuorum;
            ss << ";_XNR=" << maximumNumReplicates;
            ss << ";_WGT=" << replicateWeight;
            ss << ";_DIS=" << (removeDissimilarReplicates ? "TRUE" : "FALSE");
            ss << ";_XPK=" << maximumNumPeaksKept;
            ss << ";_XPU=" << maximumNumPeaksUsed;
            ss << ";_BDN=" << (useBayesianDenoiser ? "TRUE" : "FALSE");
            ss << ";_TBD=" << (trainBayesianDenoiser ? "TRUE" : "FALSE");
            if (useBayesianDenoiser) ss << ";_MSP=" << denoiserMinimumSignalProb;
            ss << "]";

        } else if (buildAction == "BEST_REPLICATE") {
            ss << " [" << buildAction;
            ss << ";r=" << minimumNumReplicates;
            ss << ";I=" << setFragmentation;
            ss << ";_DIS=" << (removeDissimilarReplicates ? "TRUE" : "FALSE");
            ss << ";_BDN=" << (useBayesianDenoiser ? "TRUE" : "FALSE");
            if (useBayesianDenoiser) ss << ";_MSP=" << denoiserMinimumSignalProb;
            ss << "]";

        } else if (buildAction == "QUALITY_FILTER") {
            ss << " [" << buildAction;
            ss << ";L=" << qualityLevelRemove;
            ss << ";l=" << qualityLevelMark;
            ss << ";r=" << minimumNumReplicates;
            ss << ";I=" << setFragmentation;
            ss << ";_QP1=" << (qualityPenalizeSingletons ? "TRUE" : "FALSE");
            ss << ";_QIP=" << qualityImmuneProbThreshold;
            ss << ";_QIE=" << (qualityImmuneMultipleEngines ? "TRUE" : "FALSE");
            ss << "]";

        } else if (buildAction == "DECOY") {
            ss << " [" << buildAction;
            ss << ";c=" << (decoyConcatenate ? "TRUE" : "FALSE");
            ss << ";y=" << decoySizeRatio;
            ss << ";DPS=" << (decoyPrecursorSwap ? "TRUE" : "FALSE");
            ss << ";I=" << setFragmentation;
            ss << "]";

        } else if (buildAction == "SORT_BY_NREPS") {
            ss << " [" << buildAction;
            ss << ";I=" << setFragmentation;
            ss << "]";

        } else if (buildAction == "USER_SPECIFIED_MODS") {
            ss << " [" << buildAction;
            ss << ";x=" << allowableModTokens;
            ss << ";I=" << setFragmentation;
            ss << "]";

        } else if (buildAction == "SIMILARITY_CLUSTERING") {
            ss << " [" << buildAction;
            ss << ";r=" << minimumNumReplicates;
            ss << ";I=" << setFragmentation;
            ss << ";_QUO=" << peakQuorum;
            ss << ";_XNR=" << maximumNumReplicates;
            ss << ";_WGT=" << replicateWeight;
            ss << ";_DIS=" << (removeDissimilarReplicates ? "TRUE" : "FALSE");
            ss << ";_XPK=" << maximumNumPeaksKept;
            ss << ";_XPU=" << maximumNumPeaksUsed;
            ss << ";_UCD=" << unidentifiedClusterMinimumDot;
            ss << ";_USX=" << unidentifiedSingletonXreaThreshold;
            ss << ";_BDN=" << (useBayesianDenoiser ? "TRUE" : "FALSE");
            if (useBayesianDenoiser) ss << ";_MSP=" << denoiserMinimumSignalProb;
            ss << "]";
        }

    } else if (fileType == ".msp") {
        ss << "IMPORT FROM MSP " << fileList;
    } else if (fileType == ".hlf") {
        ss << "IMPORT FROM HLF " << fileList;
    } else if (fileType == ".ms2") {
        ss << "IMPORT FROM MS2 " << fileList;

    } else if (fileType == ".pepXML") {
        ss << "IMPORT FROM PepXML " << fileList;
        ss << "[P=" << minimumProbabilityToInclude;
        ss << ";q=" << maximumFDRToInclude;
        ss << ";n=" << datasetName;
        ss << ";g=" << (setDeamidatedNXST ? "TRUE" : "FALSE");
        ss << ";o=" << (addMzXMLFileToDatasetName ? "TRUE" : "FALSE");
        ss << ";I=" << setFragmentation;
        ss << ";_RNT=" << rawSpectraNoiseThreshold;
        ss << ";_RDR=" << rawSpectraMaxDynamicRange;
        ss << ";_DCN=" << minimumDeltaCnToInclude;
        ss << ";_NAA=" << minimumNumAAToInclude;
        ss << ";_NPK=" << minimumNumPeaksToInclude;
        ss << ";_MDF=" << maximumMassDiffToInclude;
        ss << ";_CEN=" << (centroidPeaks ? "TRUE" : "FALSE");
        ss << ";_XAN=" << (skipRawAnnotation ? "TRUE" : "FALSE");
        ss << ";_BRK=" << (bracketSpectra ? "TRUE" : "FALSE");
        ss << ";_BRM=" << (mergeBracket ? "TRUE" : "FALSE");
        ss << ";_IRT=" << normalizeRTWithLandmarks;
        ss << ";_IRR=" << (normalizeRTLinearRegression ? "TRUE" : "FALSE");
//    ss << ";_PHO=" << (evaluatePhosphoSiteAssignment ? "TRUE" : "FALSE");
        ss << "]";

    } else if (fileType == ".tsv") {
        ss << "IMPORT FROM TSV " << fileList;
        ss << "[P=" << minimumProbabilityToInclude;
        ss << ";n=" << datasetName;
        ss << ";g=" << (setDeamidatedNXST ? "TRUE" : "FALSE");
        ss << ";o=" << (addMzXMLFileToDatasetName ? "TRUE" : "FALSE");
        ss << ";I=" << setFragmentation;
        ss << ";_RNT=" << rawSpectraNoiseThreshold;
        ss << ";_RDR=" << rawSpectraMaxDynamicRange;
        ss << ";_DCN=" << minimumDeltaCnToInclude;
        ss << ";_NAA=" << minimumNumAAToInclude;
        ss << ";_NPK=" << minimumNumPeaksToInclude;
        ss << ";_CEN=" << (centroidPeaks ? "TRUE" : "FALSE");
        ss << ";_XAN=" << (skipRawAnnotation ? "TRUE" : "FALSE");
        ss << "]";

    } else if (fileType == ".mzXML") {
        ss << "IMPORT FROM MZXML " << fileList;
        ss << ";n=" << datasetName;
        ss << ";o=" << (addMzXMLFileToDatasetName ? "TRUE" : "FALSE");
        ss << ";I=" << setFragmentation;
        ss << ";_RNT=" << rawSpectraNoiseThreshold;
        ss << ";_RDR=" << rawSpectraMaxDynamicRange;
        ss << ";_CEN=" << (centroidPeaks ? "TRUE" : "FALSE");
        ss << ";_UCR=" << (unidentifiedClusterIndividualRun ? "TRUE" : "FALSE");
        ss << ";_UCD=" << unidentifiedClusterMinimumDot;
        ss << ";_UX1=" << (unidentifiedRemoveSinglyCharged ? "TRUE" : "FALSE");
        ss << ";_UNP=" << unidentifiedMinimumNumPeaksToInclude;
        ss << "]";


    }


    // add filter string
    if (!filterCriteria.empty()) {
        ss << "; FILTER for " << filterCriteria;
    }
    if (!useProbTable.empty()) {
        ss << "; FILTER for entries in probability table \"" << useProbTable << "\"";
    }
    if (!useProteinList.empty()) {
        ss << "; FILTER for entries in protein list \"" << useProteinList << "\"";
    }

    // add refresh string (if .splib)
    if (!refreshDatabase.empty() && fileType == ".splib") {

        ss << "; REFRESH against " << refreshDatabase;
        if (refreshTrypticOnly) {
            ss << " (Tryptic Only)";
        }
        if (refreshDeleteUnmapped || refreshDeleteMultimapped) {
            ss << " (Delete";
            if (refreshDeleteUnmapped) {
                ss << " Unmapped";
            }
            if (refreshDeleteMultimapped) {
                ss << " Multimapped";
            }
            ss << ")";

        }
    }

    return (ss.str());


}
