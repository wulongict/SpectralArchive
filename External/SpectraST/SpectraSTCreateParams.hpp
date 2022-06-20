#ifndef SPECTRASTCREATEPARAMS_HPP_
#define SPECTRASTCREATEPARAMS_HPP_

#include <vector>
#include <string>
#include <map>


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

class SpectraSTCreateParams {

public:
    SpectraSTCreateParams();

    SpectraSTCreateParams(SpectraSTCreateParams &s);

    SpectraSTCreateParams &operator=(SpectraSTCreateParams &s);

    virtual ~SpectraSTCreateParams();

    // Options - see SpectraSTCreateParams::printUsage() for information

    // GENERAL
    string paramsFileName; // -cF
    string outputFileName; // -cN
    bool binaryFormat; // -cb  -c_BIN
    bool writeDtaFiles; // -cD  -c_DTA
    bool writeMgfFile; // -c_MGF
    bool writePAIdent; // -c_PAI
    string remark; // -cm
    string filterCriteria; // -cf
    string useProbTable; // -cT
    string useProteinList; // -cO (Oh)
    string printMRMTable; // -cM
    double minimumMRMQ3MZ; // -c_Q3L
    double maximumMRMQ3MZ; // -c_Q3H
    string removeDecoyProteins; // -c_RDY

    // LIBRARY IMPORT
    double minimumProbabilityToInclude; // -cP
    double maximumFDRToInclude; // -cq
    string datasetName; // -cn
    bool addMzXMLFileToDatasetName; // -co
    double rawSpectraNoiseThreshold; // -ct  -c_RNT
    double rawSpectraMaxDynamicRange; // -cv  -c_RDR
    bool setDeamidatedNXST; // -cg  -c_GLY
    double minimumDeltaCnToInclude; // -cu  -c_DCN
    unsigned int minimumNumAAToInclude; // -cs  -c_NAA
    unsigned int minimumNumPeaksToInclude; // -cj  -c_NPK
    double maximumMassDiffToInclude; // -c_MDF
    bool centroidPeaks; // -c_CEN
    string setFragmentation; // -cI
    bool skipRawAnnotation; // -c_XAN
    bool evaluatePhosphoSiteAssignment; // -c_PHO
    bool keepRawIntensities; // -c_RWI
    bool bracketSpectra; // -c_BRK
    bool mergeBracket; // -c_BRM

    // LIBRARY MANIPULATION
    string combineAction; // -cJ
    string buildAction; // -cA
    string plotSpectra; // -cz
    int reduceSpectrum; // -cQ
    bool reannotatePeaks; // -ca  -c_ANN

    // CONSENSUS
    unsigned int minimumNumReplicates; // -cr
    bool removeDissimilarReplicates; // -ck  -c_DIS
    double peakQuorum; // -cq  -c_QUO
    string replicateWeight; // -cw  -c_WGT
    unsigned int maximumNumReplicates; // -cx  -c_XNR
    unsigned int maximumNumPeaksUsed; // -cp  -c_XPU
    unsigned int maximumNumPeaksKept; // -cd  -c_XPK
    bool recordRawSpectra; // -c_RRS

    // DENOISER
    bool useBayesianDenoiser; // -c_BDU (this will be done for consensus, best-replicate, and similarity clustering)
    bool trainBayesianDenoiser; // -c_BDT (does nothing unless it's building consensus)
    double denoiserMinimumSignalProb; // -c_BDP
    string denoiserParamFile; // -c_BDF

    // QUALITY FILTER
    int qualityLevelRemove; // -cL
    int qualityLevelMark; // -cl
    bool qualityPenalizeSingletons; // -c1  -c_QP1
    double qualityImmuneProbThreshold; // -ci  -c_QIP
    bool qualityImmuneMultipleEngines; // -ce  -c_QIE

    // DECOY CREATION
    bool decoyConcatenate; // -cc
    int decoySizeRatio; // -cy
    bool decoyPrecursorSwap; // -c_DPS

    // REFRESH PROTEIN MAPPINGS
    string refreshDatabase; // -cD
    bool refreshDeleteUnmapped; // -cu
    bool refreshDeleteMultimapped; // -cd
    bool refreshTrypticOnly; // -c_RTO

    // RETENTION TIME NORMALIZATION
    string normalizeRTWithLandmarks; // -c_IRT
    bool normalizeRTLinearRegression; // -c_IRR

    // SEMI-EMPIRICAL SPECTRUM PREDICTION
    string allowableModTokens; // -cx
    string predictionOrigTargetFileName;
    string predictionTargetPeptidesFileName;
    string allowableSNP;
    string allowableMutatedModifications;
    int allowablePredictionDistance;


    // UNIDENTIFIED LIBRARIES
    bool unidentifiedClusterIndividualRun; // -c_UCR
    double unidentifiedClusterMinimumDot; // -c_UCD
    bool unidentifiedRemoveSinglyCharged; // -c_UX1
    int unidentifiedMinimumNumPeaksToInclude; // -c_UNP
    double unidentifiedSingletonXreaThreshold; // -c_USX

    // methods
    bool addOption(string option);

    void finalizeOptions();

    void finalizeAdvancedOption(string option);

    void readFromFile();

    string constructDescrStr(string fileList, string fileType);

    static void printUsage(ostream &out);

    static void printAdvancedOptions(ostream &out);


private:
    // m_options - stores the options as strings, wait until finalizeOptions() is called to actually set them
    vector<string> m_options;

    void setDefault();

    bool isExpectingArg(string option);


};

#endif /*SPECTRASTCREATEPARAMS_HPP_*/
