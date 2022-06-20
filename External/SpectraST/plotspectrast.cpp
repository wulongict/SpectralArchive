#include <string>
#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>
#include <algorithm>
#include <math.h>
#include <map>
#include "Peptide.hpp"
#include "FileUtils.hpp"

#ifdef STANDALONE_LINUX

#include "SpectraST_cramp.hpp"

#else
#include "cramp.hpp"
#include "common/hooks_tpp.h"
#include "common/TPPVersion.h"
#include "common/constants.h"
#endif

#include "SpectraSTConstants.hpp"
#include "SpectraSTLibEntry.hpp"
#include "SpectraSTPeakList.hpp"
#include "SpectraSTQuery.hpp"
#include "SpectraSTLog.hpp"

#define BGCOLOR     "#42D4FD"
#define DATACELLCOLOR   "#FFDDDD"
#define HEADERCELLCOLOR   "#42D4FD"
#define SUBHEADERCELLCOLOR   "#FF9999"

#define LIB_LABEL_BKGRND "#FFDDDD"
#define QUERY_LABEL_BKGRND "#FF9999"

#define UNLABEL_BKGRND "#FFFFFF"

#define IMAGEWIDTH  1.0
#define IMAGEHEIGHT 1.0

#define LIB_BASE_COLOR 3 /* blue */
#define LIB_LABEL_COLOR 1 /* red */
#define QUERY_BASE_COLOR -1 /* black */
#define QUERY_LABEL_COLOR 1 /* red */
#define AXIS_COLOR -1

#define ATLAS_WEBSERVER_ROOT "/net/dblocal/wwwspecial/phosphopep/"
// #define ATLAS_WEBSERVER_ROOT "/net/dblocal/wwwspecial/peptideatlas/"

// #define PLOT_EPS

#ifdef PLOT_EPS
// grey and black only
#define LIB_BASE_COLOR "rgb \"#888888\""
#define LIB_LABEL_COLOR -1 
#define QUERY_BASE_COLOR "rgb \"#888888\"" 
#define QUERY_LABEL_COLOR -1 
#endif

#define INTENSITY_PAD 0.02
#define DEFAULT_MATCH_TOLERANCE 1.0
#define LABEL_MZ_TOLERANCE 1.0

//#define MAX_CHARGE 5

using namespace std;

// environment variables
typedef struct _plotspectrastEnv {
    string splibFileName;
    fstream::off_type splibFileOffset;
    string mzXMLFileName;
    fstream::off_type mzXMLScanNum;

    double minMz;
    double maxMz;
    double matchTolerance;
    double intensityZoom;

    int labelType;

    int showA[MAX_CHARGE];
    int showB[MAX_CHARGE];
    int showC[MAX_CHARGE];
    int showY[MAX_CHARGE];
    int showZ[MAX_CHARGE];

/*  
  int showA1;
  int showA2;
  int showA3;
  int showB1;
  int showB2;
  int showB3;
  int showY1;
  int showY2;
  int showY3;
  */

    int showPLoss;
    int showCommonLoss;
    int showAll;
    double showMinIntensity;
    int showNumPeaks;
    int blankNearParentRegion;

    string outputPath;
    string pngFileName;

} PlotSpectraSTEnv;

void readLib(FileName &splibFile, fstream::off_type &splibFileOffset, SpectraSTLibEntry **lib, PlotSpectraSTEnv &env);

void readQuery(FileName &mzXMLFile, int mzXMLScanNum, SpectraSTQuery **query);

void plotLibQuerySpecs(SpectraSTLibEntry *lib, SpectraSTQuery *query, map<string, Peak> &labels, PlotSpectraSTEnv &env);

void
plotLibLibSpecs(SpectraSTLibEntry *upLib, SpectraSTLibEntry *downLib, map<string, Peak> &labels, PlotSpectraSTEnv &env);

void plotSpecs(SpectraSTPeakList *up, SpectraSTPeakList *down, string &outputBaseName, map<string, Peak> &labels,
               PlotSpectraSTEnv &env);

bool labelLib(Peak &peak, map<string, Peak> &labels, PlotSpectraSTEnv &env);

bool labelQuery(Peak &peak, map<string, Peak> &labels, PlotSpectraSTEnv &env);

void plotAllLabels(ofstream &plotFout, map<string, Peak> &labels, double upMaxIntensity, double downMaxIntensity,
                   PlotSpectraSTEnv &env);

void plotLabels(ofstream &plotFout, vector<pair<string, Peak> > &labelsVtr, double maxIntensity, PlotSpectraSTEnv &env,
                bool isDown);

bool sortLabelsByMzAsc(pair<string, Peak> a, pair<string, Peak> b);

void printHTML(SpectraSTLibEntry *upLib, SpectraSTQuery *query, SpectraSTLibEntry *downLib, map<string, Peak> &labels,
               PlotSpectraSTEnv &env);

void printIonTable(Peptide *pep, map<string, Peak> &labels, PlotSpectraSTEnv &env);

void printComments(SpectraSTLibEntry *lib);

string translateScoreName(string scoreName, string engine);

#ifdef RUN_AS_CGI
string processCGIStr(string s);
void extractCGI(PlotSpectraSTEnv& env);
#endif

SpectraSTLog *g_log;

int main(int argc, char **argv) {

    // command structure: plotspectrast <.splib> offset <.mzXML> scanNum
    //                    plotspectrast <.splib> offset <.dta or .msp or .none>
    //                    plotspectrast <.msp> <.dta or .msp or .none>

#ifndef STANDALONE_LINUX
    hooks_tpp handler(argc,argv); // set up install paths etc
#endif

    Peptide::defaultTables();
    PlotSpectraSTEnv env;

#ifndef RUN_AS_CGI


    if (argc < 3 || argc > 5) {
        cerr << "need 2 - 4 arguments!" << endl;
        return (1);
    }

    if (argc == 3) {
        env.splibFileName = argv[1];
        env.splibFileOffset = 0;
        env.mzXMLFileName = argv[2];
        env.mzXMLScanNum = 0;
    } else if (argc == 4) {
        env.splibFileName = argv[1];
        env.splibFileOffset = strtoull(argv[2], NULL, 10);
        env.mzXMLFileName = argv[3];
        env.mzXMLScanNum = 0;
    } else if (argc == 5) {
        env.splibFileName = argv[1];
        env.splibFileOffset = strtoull(argv[2], NULL, 10);
        env.mzXMLFileName = argv[3];
        env.mzXMLScanNum = strtoull(argv[4], NULL, 10);
    }

    env.labelType = 0;
    env.minMz = 0.0;
    env.maxMz = 0.0;
    env.matchTolerance = 1.0;
    env.intensityZoom = 1.0;

    for (int c = 0; c < MAX_CHARGE; c++) {
        env.showA[c] = 1;
        env.showB[c] = 1;
        env.showC[c] = 0;
        env.showY[c] = 1;
        env.showZ[c] = 0;
    }
    /*
    env.showA1 = 1;
    env.showA2 = 1;
    env.showA3 = 1;
    env.showB1 = 1;
    env.showB2 = 1;
    env.showB3 = 1;
    env.showY1 = 1;
    env.showY2 = 1;
    env.showY3 = 1;
    */
    env.showCommonLoss = 1;
    env.showPLoss = 1;
    env.showAll = 1;
    env.pngFileName = "";
    env.outputPath = "";
    env.showMinIntensity = 0.0;
    env.showNumPeaks = 99999;
    env.blankNearParentRegion = 0;

#endif

#ifdef RUN_AS_CGI
    // print out the http/html headers before anything else, so we'll at
    // least be able to see any error messages:
    cout << "Content-type: text/html" << endl << endl; // be sure to flush this
    // right away (endl)
    cout << "<HTML>" << endl;
    cout << "  <HEAD>" << endl;
    cout << "    <TITLE>PlotSpectraST by Henry Lam (c) ISB 2006 - TPP Version " << szTPPVersionInfo << "</TITLE>" << endl;
    cout << "  </HEAD>" << endl;
    cout << endl;
    cout << "<BODY BGCOLOR=\"#EEEEEE\" OnLoad=\"self.focus();\">" << endl;

    extractCGI(env);

#endif

    makeFullPath(env.splibFileName);
    makeFullPath(env.mzXMLFileName);

    FileName splibFile;
    parseFileName(env.splibFileName, splibFile);
    if (splibFile.ext != ".splib" && splibFile.ext != ".msp") {
        g_log->error("PLOTSPECTRAST", "Library file specified is not of supported format. Exiting.");
        return (1);
    }

    FileName mzXMLFile;
    parseFileName(env.mzXMLFileName, mzXMLFile);
    if (mzXMLFile.ext != ".mzML" && mzXMLFile.ext != ".mzXML" && mzXMLFile.ext != ".dta" && mzXMLFile.ext != ".msp"
        && mzXMLFile.ext != ".none" && mzXMLFile.ext != ".splib") {
        g_log->error("PLOTSPECTRAST", "Query spectrum file specified is not of supported format. Exiting.");
        return (1);
    }

    g_log = new SpectraSTLog(mzXMLFile.path + "plotspectrast.log");

    if (env.outputPath.empty()) {
        env.outputPath = mzXMLFile.path;
    }

    SpectraSTLibEntry *upLib = NULL;
    SpectraSTLibEntry *downLib = NULL;
    SpectraSTQuery *query = NULL;

    readLib(splibFile, env.splibFileOffset, &upLib, env);

    if (mzXMLFile.ext != ".splib") {
        readQuery(mzXMLFile, env.mzXMLScanNum, &query);
    } else {
        readLib(mzXMLFile, env.mzXMLScanNum, &downLib, env);
    }

    // minMz and maxMz unset, use default (minMz = 1.0, maxMz enough to cover both spectra)
    if (env.minMz < 0.001 || env.maxMz < 0.001) {

        double upMaxMz = upLib->getPeakList()->getMaxMz();
        double downMaxMz = query ? query->getPeakList()->getMaxMz() : downLib->getPeakList()->getMaxMz();

        env.minMz = 1.0;
        env.maxMz = upMaxMz > downMaxMz ? upMaxMz : downMaxMz;
        env.maxMz = 100.0 * ((int) (env.maxMz / 100) + 1); // round up to the next 100
    }

    map<string, Peak> labels;
    // put all plot files, including the png in the same directory as the mzXML file

    if (query) {
        plotLibQuerySpecs(upLib, query, labels, env);
    } else {
        plotLibLibSpecs(upLib, downLib, labels, env);
    }

#ifdef RUN_AS_CGI
    // print HTML file
    printHTML(upLib, query, downLib, labels, env);

#endif

    delete (upLib);
    if (query) delete (query);
    if (downLib) delete (downLib);

    Peptide::deleteTables();

    delete (g_log);

    return (0);

}

// readLib - reads the library entry
void readLib(FileName &splibFile, fstream::off_type &splibFileOffset, SpectraSTLibEntry **lib, PlotSpectraSTEnv &env) {

    ifstream fin;
    if (!myFileOpen(fin, splibFile.path + splibFile.name + splibFile.ext, true)) {
        g_log->error("PLOTSPECTRAST",
                     "Cannot open library file \"" + splibFile.path + splibFile.name + splibFile.ext + "\". Exiting.");
        g_log->crash();
    }

    bool binaryLib = true;
    char firstChar = fin.peek();
    if (firstChar == '#' || firstChar == 'N') {
        binaryLib = false;
    }

    if (splibFile.ext == ".msp") {
        fin.seekg(0);
    } else {
        if (splibFileOffset < 0) {
            g_log->error("PLOTSPECTRAST",
                         "No file offset given for \"" + splibFile.path + splibFile.name + splibFile.ext +
                         "\". Exiting.");
            g_log->crash();
        }
        fin.seekg(splibFileOffset);
    }

    *lib = new SpectraSTLibEntry(fin, binaryLib);
    (*lib)->annotatePeaks();

    if ((*lib)->getFragType() == "ETD") {
        int charge = (*lib)->getCharge();
        for (int c = 0; c < charge - 1; c++) {
            env.showC[c] = 1;
            env.showZ[c] = 1;
            env.showA[c] = 0;
            env.showB[c] = 0;
        }
        env.showY[charge - 1] = 0;
        env.showB[charge - 1] = 0;
    }

}

// readQuery - reads the query spectrum
void readQuery(FileName &mzXMLFile, int mzXMLScanNum, SpectraSTQuery **query) {

    if (mzXMLFile.ext == ".none") {
        // no query file supplied - plot library spectrum onlyInputStruct
        *query = new SpectraSTQuery(mzXMLFile.name, 0.0, 1, "", new SpectraSTPeakList(0.0, 0));
        return;
    }


    if (mzXMLScanNum < 0) {
        // this should be a dta file or a msp file
        if (mzXMLFile.ext != ".dta" && mzXMLFile.ext != ".msp") {
            g_log->error("PLOTSPECTRAST",
                         "No scan number given for \"" + mzXMLFile.path + mzXMLFile.name + mzXMLFile.ext +
                         "\". Only library spectrum is plotted.");
            *query = new SpectraSTQuery(mzXMLFile.name, 0.0, 1, "", new SpectraSTPeakList(0.0, 0));
            return;
        }

        ifstream fin;
        if (!myFileOpen(fin, mzXMLFile.path + mzXMLFile.name + mzXMLFile.ext)) {
            g_log->error("PLOTSPECTRAST",
                         "Cannot open query spectrum file \"" + mzXMLFile.path + mzXMLFile.name + mzXMLFile.ext +
                         "\". Only library spectrum is plotted.");
            *query = new SpectraSTQuery(mzXMLFile.name, 0.0, 1, "", new SpectraSTPeakList(0.0, 0));
            return;
        }

        string line("");

        string name;
        double precursorMH = 0.0;
        double precursorMz = 0.0;
        int charge = 0;
        string comments;
        int numPeaks = 0;

        string::size_type pos = 0;

        if (!nextLine(fin, line, "", "")) {
            // nothing in the file!
            return;
        }

        if (line.compare(0, 5, "Name:") == 0) {
            // msp style
            do {
                if (line.compare(0, 5, "Name:") == 0) {
                    name = nextToken(line, 5, pos, "\r\n", " \t\r\n");
                } else if (line.compare(0, 3, "MW:") == 0) {
                    precursorMH = atof(nextToken(line, 3, pos, "\r\n", " \t\r\n").c_str());
                } else if (line.compare(0, 12, "PrecursorMZ:") == 0) {
                    precursorMz = atof(nextToken(line, 12, pos, "\r\n", " \t\r\n").c_str());
                } else if (line.compare(0, 8, "Comment:") == 0) {
                    comments = nextToken(line, 8, pos, "\r\n", " \t\r\n");
                } else if (line.compare(0, 10, "Num peaks:") == 0) {
                    numPeaks = atoi(nextToken(line, 10, pos, "\r\n", " \t\r\n").c_str());
                } else {
                    break;
                }
            } while (nextLine(fin, line, "", ""));

            if (precursorMH > 0.001 && precursorMz > 0.001) {
                charge = (int) (precursorMH / precursorMz + 0.5);
            }

            if (line == "_EOF_") {
                // no peaks
                return;
            }

        } else {
            // dta style

            pos = 0;
            name = mzXMLFile.name;

            precursorMH = atof(nextToken(line, pos, pos, " \t\r\n", " \t\r\n").c_str());
            charge = atoi(nextToken(line, pos, pos, " \t\r\n", " \t\r\n").c_str());
            precursorMz = 0.0;
            if (charge == 0) {
                precursorMz = precursorMH;
            } else {
                precursorMz = (precursorMH + charge * (*Peptide::AAAverageMassTable)['+']) / (double) charge;
            }

            if (!nextLine(fin, line, "", "")) {
                // no peaks
                return;
            }
        }

        // line should point to first line of peaks now

        SpectraSTPeakList *peakList = new SpectraSTPeakList(precursorMz, charge);
        *query = new SpectraSTQuery(name, precursorMz, charge, "", peakList);

        do {
            pos = 0;
            double mz = atof((nextToken(line, 0, pos)).c_str());
            float intensity = atof((nextToken(line, pos, pos)).c_str());
            peakList->insert(mz, intensity, "", "");

        } while (nextLine(fin, line, "", ""));


    } else {

        // normal mzXML file operation

        string mzXMLFileName(mzXMLFile.path + mzXMLFile.name + mzXMLFile.ext);
        cRamp cr(mzXMLFileName.c_str());

        rampScanInfo *scanInfo = cr.getScanHeaderInfo(mzXMLScanNum);

        stringstream logss;
        logss << "Cannot read scan #" << mzXMLScanNum << " from \"" << mzXMLFile.path << mzXMLFile.name << mzXMLFile.ext
              << "\". Only library spectrum is plotted.";

        if (!scanInfo) {
            g_log->error("PLOTSPECTRAST", logss.str());
            *query = new SpectraSTQuery(mzXMLFile.name, 0.0, 1, "", new SpectraSTPeakList(0.0, 0));
            return;
        }

        // Ignore MS1 scans
        if (scanInfo->m_data.msLevel == 1) {
            g_log->error("PLOTSPECTRAST", logss.str());
            *query = new SpectraSTQuery(mzXMLFile.name, 0.0, 1, "", new SpectraSTPeakList(0.0, 0));
            delete scanInfo;
            return;
        }
        if (scanInfo->m_data.acquisitionNum != mzXMLScanNum) {
            g_log->error("PLOTSPECTRAST", logss.str());
            *query = new SpectraSTQuery(mzXMLFile.name, 0.0, 1, "", new SpectraSTPeakList(0.0, 0));
            delete scanInfo;
            return;
        }


        int precursorCharge = scanInfo->m_data.precursorCharge;
        if (precursorCharge < 1) precursorCharge = 0;

        rampPeakList *peakList = cr.getPeakList(mzXMLScanNum);

        if (!peakList) {
            g_log->error("PLOTSPECTRAST", logss.str());
            *query = new SpectraSTQuery(mzXMLFile.name, 0.0, 1, "", new SpectraSTPeakList(0.0, 0));
            delete scanInfo;
            return;
        }

        int peakCount = peakList->getPeakCount();

        SpectraSTPeakList *pl = new SpectraSTPeakList(scanInfo->m_data.precursorMZ, peakCount);
        *query = new SpectraSTQuery(mzXMLFile.name, mzXMLScanNum, mzXMLScanNum, scanInfo->m_data.precursorMZ,
                                    precursorCharge, "", pl);

        for (int j = 0; j < peakCount; j++) {
            double mz = peakList->getPeak(j)->mz;
            float intensity = (float) (peakList->getPeak(j)->intensity);
            pl->insert(mz, intensity, "", "");
        }

        delete peakList;

        delete scanInfo;

    }
}

// plotSpecs - given the library and query spectra, the labels, writes the .gp and .spec file and calls gnuplot to plot
void
plotLibQuerySpecs(SpectraSTLibEntry *lib, SpectraSTQuery *query, map<string, Peak> &labels, PlotSpectraSTEnv &env) {

    string outputBaseName(query->getName());
    // added according to henry's version: L WU Sep 20
    outputBaseName += "_" + lib->getSafeName();

    plotSpecs(lib->getPeakList(), query->getPeakList(), outputBaseName, labels, env);

}

void plotLibLibSpecs(SpectraSTLibEntry *upLib, SpectraSTLibEntry *downLib, map<string, Peak> &labels,
                     PlotSpectraSTEnv &env) {

    string outputBaseName(downLib->getSafeName());
    plotSpecs(upLib->getPeakList(), downLib->getPeakList(), outputBaseName, labels, env);

}

void plotSpecs(SpectraSTPeakList *up, SpectraSTPeakList *down, string &outputBaseName, map<string, Peak> &labels,
               PlotSpectraSTEnv &env) {

    string specFileName(env.outputPath + outputBaseName + ".spec");
    ofstream specFout;
    if (!myFileOpen(specFout, specFileName)) {
        g_log->error("PLOTSPECTRAST", "Cannot open \"" + specFileName + "\" for outputting spectrum files. Exiting.");
        g_log->crash();
    }

    string plotFileName(env.outputPath + outputBaseName + ".gp");
    ofstream plotFout;
    if (!myFileOpen(plotFout, plotFileName)) {
        g_log->error("PLOTSPECTRAST", "Cannot open \"" + plotFileName + "\" for writing gnuplot plot file. Exiting.");
        g_log->crash();
    }

#ifdef PLOT_EPS
    env.pngFileName = env.outputPath + outputBaseName + ".eps";
#else
    env.pngFileName = env.outputPath + outputBaseName + ".png";
#endif

    //up->rankTransformWithQuota(100, true, 0, 30, 5);
    //down->rankTransformWithQuota(100, true, 0, 30, 5);

    double upMaxIntensity = up->getOrigMaxIntensity();
    double downMaxIntensity = down->getOrigMaxIntensity();
    double upPrecursorMz = up->getParentMz();
    double downPrecursorMz = down->getParentMz();

    if (env.blankNearParentRegion) {
        double newUpMaxIntensity = 0.0;
        for (unsigned int i = 0; i < up->getNumPeaks(); i++) {
            Peak p;
            up->getPeak(i, p);
            if (up->isNearPrecursor(p.mz)) {
                continue;
            }
            if (newUpMaxIntensity < p.intensity) newUpMaxIntensity = p.intensity;
        }
        upMaxIntensity = newUpMaxIntensity;

        double newDownMaxIntensity = 0.0;
        for (unsigned int i = 0; i < down->getNumPeaks(); i++) {
            Peak p;
            down->getPeak(i, p);
            if (down->isNearPrecursor(p.mz)) {
                continue;
            }
            if (newDownMaxIntensity < p.intensity) newDownMaxIntensity = p.intensity;
        }
        downMaxIntensity = newDownMaxIntensity;


    }

    bool label = true;
    for (unsigned int i = 0; i < up->getNumPeaks(); i++) {
        if ((int) i >= env.showNumPeaks) {
            label = false;
        }
        Peak p;
        up->getNthLargestPeak(i + 1, p);
        if (p.intensity < env.showMinIntensity) {
            label = false;
        }

        if (env.blankNearParentRegion && up->isNearPrecursor(p.mz)) {
            continue;
        }

        if (label && labelLib(p, labels, env)) {
            // label
            specFout << p.mz << ' ' << p.intensity / upMaxIntensity * env.intensityZoom << ' ' << 0.0 << ' ' << 0.0
                     << ' ' << 0.0 << ' ' << 0.0 << endl;
        } else {
            // don't label
            specFout << p.mz << ' ' << 0.0 << ' ' << p.intensity / upMaxIntensity * env.intensityZoom << ' ' << 0.0
                     << ' ' << 0.0 << ' ' << 0.0 << endl;
        }

    }

    if (down->isAnnotated()) {

        label = true;
        for (unsigned int i = 0; i < down->getNumPeaks(); i++) {

            if ((int) i >= env.showNumPeaks) {
                label = false;
            }
            Peak p;
            down->getNthLargestPeak(i + 1, p);
            if (p.intensity < env.showMinIntensity) {
                label = false;
            }

            if (env.blankNearParentRegion && down->isNearPrecursor(p.mz)) {
                continue;
            }

            p.intensity *= -1.0;

            if (label && labelLib(p, labels, env)) {
                // label
                specFout << p.mz << ' ' << 0.0 << ' ' << 0.0 << ' '
                         << p.intensity / downMaxIntensity * env.intensityZoom << ' ' << 0.0 << ' ' << 0.0 << endl;
            } else {
                // don't label
                specFout << p.mz << ' ' << 0.0 << ' ' << 0.0 << ' ' << 0.0 << ' '
                         << p.intensity / downMaxIntensity * env.intensityZoom << ' ' << 0.0 << endl;

            }
        }

    } else {

        for (unsigned int i = 0; i < down->getNumPeaks(); i++) {

            Peak p;
            down->getPeak(i, p);

            if (env.blankNearParentRegion && down->isNearPrecursor(p.mz)) {
                continue;
            }

            if (labelQuery(p, labels, env)) {
                // label
                specFout << p.mz << ' ' << 0.0 << ' ' << 0.0 << ' '
                         << p.intensity * (-1.0) / downMaxIntensity * env.intensityZoom << ' ' << 0.0 << ' ' << 0.0
                         << endl;
            } else {
                // don't label
                specFout << p.mz << ' ' << 0.0 << ' ' << 0.0 << ' ' << 0.0 << ' '
                         << p.intensity * (-1.0) / downMaxIntensity * env.intensityZoom << ' ' << 0.0 << endl;

            }
        }

    }

#ifdef PLOT_EPS
    plotFout << "set terminal postscript eps" << endl;
#else
    plotFout << "set terminal png" << endl;
#endif

    plotFout << "set output \"" << env.pngFileName << "\"" << endl;
    plotFout << "set size " << IMAGEWIDTH << "," << IMAGEHEIGHT << endl;
    plotFout << "set nokey" << endl;
    plotFout << "set border 1" << endl;
    plotFout << "set xtics border nomirror" << endl;
    plotFout << "set xzeroaxis linetype -1 linewidth 1.0" << endl;
    plotFout << "set yzeroaxis linetype -1 linewidth 1.0" << endl;
    plotFout << "set mxtics" << endl;
    plotFout << "set noytics" << endl;

    plotFout << "set origin 0.0,0.0" << endl;

    if (down->getNumPeaks() > 0) {
        plotFout << "set yrange [-1.2:1.2]" << endl;
    } else {
        plotFout << "set yrange [0.0:1.2]" << endl;
    }

    plotFout.precision(3);
    plotFout << "set label \"" << upMaxIntensity / env.intensityZoom << "\" at " << env.minMz << ",1.04 front" << endl;
    plotFout << "set label \"Library\" at " << env.maxMz << ",1.04 right front" << endl;

    if (down->getNumPeaks() > 0) {
        plotFout.precision(3);
        plotFout << "set label \"" << downMaxIntensity / env.intensityZoom << "\" at " << env.minMz << ",-1.04 front"
                 << endl;
        plotFout << "set label \"Query\" at " << env.maxMz << ",-1.04 right front" << endl;
    }

    plotAllLabels(plotFout, labels, upMaxIntensity / env.intensityZoom, downMaxIntensity / env.intensityZoom, env);

    plotFout.precision(3);
    plotFout << "plot [" << env.minMz << ":" << env.maxMz << "] \"" << specFileName
             << "\" using 1:3 with impulse lt 1 lc " << LIB_BASE_COLOR;
    plotFout << ", \"" << specFileName << "\" using 1:2 with impulse lt 1 lc " << LIB_LABEL_COLOR;
    plotFout << ", \"" << specFileName << "\" using 1:5 with impulse lt 1 lc " << QUERY_BASE_COLOR;
    plotFout << ", \"" << specFileName << "\" using 1:4 with impulse lt 1 lc " << QUERY_LABEL_COLOR;
    plotFout << ", \"" << specFileName << "\" using 1:6 with impulse lt 1 lc " << AXIS_COLOR << endl;


    // PLOT!
#ifdef ATLAS
    myGnuplot(plotFileName, "/usr/local/bin/gnuplot");
#else
    myGnuplot(plotFileName);
#endif

    // delete the plot file and spec file
    plotFout.close();
    specFout.close();
    removeFile(plotFileName);
    removeFile(specFileName);


}

// labelLib - create labels for the library spectrum
bool labelLib(Peak &peak, map<string, Peak> &labels, PlotSpectraSTEnv &env) {

    bool COLOR_ALL_ASSIGNED = true;

    string::size_type pos = 0;
    string::size_type dummyPos = 0;
    bool assigned = false;

    stringstream ss;
    ss.precision(1);
    ss << fixed << peak.mz;
    string mzStr(ss.str());

    bool isDown = (peak.intensity < 0.0);

    bool matched = false;

    if (isDown) {
        for (map<string, Peak>::iterator i = labels.begin(); i != labels.end(); i++) {
            if (i->first[0] != '~' && fabs(peak.mz - i->second.mz) < env.matchTolerance) {
                i->second.info = "M";
                matched = true;
                break;
            }
        }
    }

    while (true) {
        string ionSlashSthg = nextToken(peak.annotation, pos, pos, ",\t\r\n");
        if (ionSlashSthg.empty()) {
            break;
        }

        if (COLOR_ALL_ASSIGNED && ionSlashSthg[0] != '?') {
            assigned = true;
        }

        string labelText("");

        if (ionSlashSthg[0] != 'b' &&
            ionSlashSthg[0] != 'y' &&
            ionSlashSthg[0] != 'a' &&
            ionSlashSthg[0] != 'c' &&
            ionSlashSthg[0] != 'z' &&
            ionSlashSthg[0] != 'p' &&
            ionSlashSthg[0] != 'I') {
            // no ion
            pos++;
            continue;
        }

        string::size_type slashPos;

        string ion = nextToken(ionSlashSthg, 0, slashPos, "/\t\r\n");
        labelText = ion;

        bool show = false;

        if (ion[0] == 'I' && ion.find('i', 0) == string::npos) {
            labelText = ion;
            show = (env.showAll == 1);

        } else {

            string::size_type minusPos = ion.find('-', 0);
            string::size_type plusPos = ion.find('+', 0);
            string::size_type iPos = ion.find('i', 0);
            string::size_type xPos = ion.find('x', 0);
            string::size_type hatPos = ion.find('^', 0);

            if (plusPos != string::npos || iPos != string::npos || xPos != string::npos) {
                pos++;
                continue;
            }

            string ionType = ion.substr(0, 1); // first char of ion must be ion type
            string ionAA = nextToken(ion, 1, dummyPos,
                                     "ix+-/^ \t\r\n"); // AA number must start at second position of ion

            string neutralLoss = "";
            if (minusPos != string::npos) {
                neutralLoss = nextToken(ion, minusPos + 1, dummyPos, "ix/^ \t\r\n");
            }

            int charge = 1;
            string chargeStr("1");
            if (hatPos != string::npos) {
                chargeStr = nextToken(ion, hatPos + 1, dummyPos, "/ \t\r\n");
                charge = atoi(chargeStr.c_str());
            }

            labelText = ionType + ionAA;
            if (!neutralLoss.empty()) {
                labelText += "-";
                if (neutralLoss == "17") {    // to decongest, -17 and -18 will only show up as one (-18)
                    labelText += "18";
                } else if (neutralLoss == "44") { // to decongest, -44, -45 and -46 will only show up as one (-45)
                    labelText += "45";
                } else if (neutralLoss == "46") {
                    labelText += "45";
                } else if (neutralLoss == "34") { // to decongest, -34, -35 and -36 will only show up as one (-35)
                    labelText += "35";
                } else if (neutralLoss == "36") {
                    labelText += "35";
                } else {
                    labelText += neutralLoss;
                }
            }
            if (charge > 1) labelText += "^" + chargeStr;

            if (env.showAll == 1) {
                show = true;
            } else if (charge <= MAX_CHARGE &&
                       ((env.showA[charge - 1] == 1 && ionType == "a") ||
                        (env.showB[charge - 1] == 1 && ionType == "b") ||
                        (env.showY[charge - 1] == 1 && ionType == "y") ||
                        (env.showC[charge - 1] == 1 && ionType == "c") ||
                        (env.showZ[charge - 1] == 1 && ionType == "z"))) {
                if (neutralLoss.empty()) {
                    show = true;
                } else if (env.showCommonLoss == 1 &&
                           (neutralLoss == "17" ||
                            neutralLoss == "18" ||
                            neutralLoss == "98" ||
                            neutralLoss == "80")) {
                    show = true;
                } else {
                    show = false;
                }

            } else if (ionType == "p") {
                if (neutralLoss.empty()) {
                    show = true;
                } else if (env.showPLoss == 1) {
                    show = true;
                } else {
                    show = false;
                }
            }
        }

        if (isDown && !labelText.empty()) {
            labelText = "~" + labelText;
        }

        if (show && !labelText.empty()) {
            map<string, Peak>::iterator found;

            if ((found = labels.find(labelText)) != labels.end()) {
                // duplicate label, only replace if this one has higher intensity
                if (fabs((*found).second.intensity) < fabs(peak.intensity)) {
                    (*found).second.mz = peak.mz;
                    (*found).second.intensity = peak.intensity;
                }
                assigned = true;
            } else {
                assigned = true;
                labels[labelText].mz = peak.mz;
                labels[labelText].intensity = peak.intensity;

                if (matched) {
                    labels[labelText].info = "M";
                } else {
                    labels[labelText].info = "";
                }

            }
        }

        // onto next annotation of the same peak
        pos++;
    }
    return (assigned && env.labelType != 2);
}

// labelQuery - given the labels in the library spectrum, mark the labeled peaks in the query spectra
bool labelQuery(Peak &peak, map<string, Peak> &labels, PlotSpectraSTEnv &env) {

    bool matched = false;
    for (map<string, Peak>::iterator i = labels.begin(); i != labels.end(); i++) {
        if (fabs(peak.mz - (*i).second.mz) < env.matchTolerance) {
            // mark them so that they can be displayed in a different color
            (*i).second.info = "M";
            matched = true;
        }
    }
    return (matched && env.labelType != 2);
}

// plotLabels - create "set label" commands for gnuplots given the labels

void plotAllLabels(ofstream &plotFout, map<string, Peak> &labels, double upMaxIntensity, double downMaxIntensity,
                   PlotSpectraSTEnv &env) {

    vector<pair<string, Peak> > upLabelsVtr;
    vector<pair<string, Peak> > downLabelsVtr;
    for (map<string, Peak>::iterator i = labels.begin(); i != labels.end(); i++) {
        if (i->first[0] != '~') {
            upLabelsVtr.push_back(*i);
        } else {
            pair<string, Peak> p;
            p.first = i->first.substr(1);
            p.second = i->second;
            downLabelsVtr.push_back(p);
        }
    }

    std::sort(upLabelsVtr.begin(), upLabelsVtr.end(), sortLabelsByMzAsc);
    plotLabels(plotFout, upLabelsVtr, upMaxIntensity, env, false);

    if (!(downLabelsVtr.empty())) {
        std::sort(downLabelsVtr.begin(), downLabelsVtr.end(), sortLabelsByMzAsc);
        plotLabels(plotFout, downLabelsVtr, downMaxIntensity, env, true);
    }
}

void plotLabels(ofstream &plotFout, vector<pair<string, Peak> > &labelsVtr, double maxIntensity, PlotSpectraSTEnv &env,
                bool isDown) {

    double labelMz = 0.0;
    double labelIntensity = 0.0;
    bool labelMatched = false;
    string labelText("");

    vector<pair<string, Peak> >::iterator j;
    for (j = labelsVtr.begin(); j != labelsVtr.end(); j++) {
        if (!labelText.empty() && ((*j).second.mz - labelMz < 1.0)) {
            // too close to the last one, will overlap.
            // instead create one combined label with all the overlapping ones
            if (env.labelType == 0) {
                labelText += ",";
                labelText += (*j).first;
            } else if (env.labelType == 1) {
                stringstream ss;
                ss.precision(1);
                ss << fixed << (*j).second.mz;
                if (labelText != ss.str()) {
                    labelText += ",";
                    labelText += ss.str();
                }
            }
            labelMz = 0.5 * (labelMz + (*j).second.mz);
            labelIntensity = labelIntensity > (*j).second.intensity ? labelIntensity : (*j).second.intensity;
            labelMatched = labelMatched || ((*j).second.info == "M");
        } else {
            // print what's already created
            if (!labelText.empty()) {
                plotFout << "set label \"" << labelText << "\" at ";
                plotFout.precision(3);
                plotFout << labelMz << ",";
                plotFout.precision(3);
                if (!isDown) {
                    plotFout << labelIntensity / maxIntensity + INTENSITY_PAD << " left rotate front";
                } else {
                    plotFout << labelIntensity / maxIntensity - 0.035 * labelText.length() - INTENSITY_PAD
                             << " left rotate front";
                }
                if (labelMatched) {
                    plotFout << " tc lt 1";
                } else {
                    plotFout << " tc lt -1";
                }
                plotFout << endl;

            }
            // create new label
            if (env.labelType == 0) {
                labelText = (*j).first;
            } else if (env.labelType == 1) {
                stringstream ss;
                ss.precision(1);
                ss << fixed << (*j).second.mz;
                labelText = ss.str();
            }
            labelMz = (*j).second.mz;
            labelIntensity = (*j).second.intensity;
            labelMatched = ((*j).second.info == "M");
        }


    }

    // print the last label
    if (!labelText.empty()) {
        plotFout << "set label \"" << labelText << "\" at ";
        plotFout.precision(3);
        plotFout << labelMz << ",";
        plotFout.precision(3);
        if (!isDown) {
            plotFout << labelIntensity / maxIntensity + INTENSITY_PAD << " left rotate front";
        } else {
            plotFout << labelIntensity / maxIntensity - INTENSITY_PAD << " right rotate front";
        }
        if (labelMatched) {
            plotFout << " tc lt 1";
        } else {
            plotFout << " tc lt -1";
        }
        plotFout << endl;
    }

}

// printHTML - prints the HTML file (in CGI operation) with the plotted spectra and the various tables
void printHTML(SpectraSTLibEntry *upLib, SpectraSTQuery *query, SpectraSTLibEntry *downLib, map<string, Peak> &labels,
               PlotSpectraSTEnv &env) {

    int charge = upLib->getCharge();

    if (downLib && downLib->getCharge() > charge) {
        charge = downLib->getCharge();
    }

    /*
     /// debug
     cout << lib->getFullName() << "<BR>" << endl;
     cout << query->getName() << "<BR>" << endl;
     cout << "</BODY></HTML>" << endl;

     return;
     /// debug
     */

    char *scriptPtr = getenv("SCRIPT_NAME");
    string script;
    if (!scriptPtr) {
        cerr << "Not run as a CGI script! HTML form won't work!" << endl;
        script = "SCRIPT";
    } else {
        script = scriptPtr;
    }

    cout << "<TABLE CELLPADDING=\"3\"><TR ALIGN=MIDDLE BGCOLOR=\"" << BGCOLOR << "\">" << endl;
    cout << "<TD ALIGN=\"LEFT\" NOWRAP><FORM ACTION=\"" << script << "\" METHOD=\"GET\">" << endl;

    cout << "<TT><B>X-range:</B><BR>" << endl;
    cout << " &nbsp;<INPUT TYPE=\"TEXT\" NAME=\"MinMz\" VALUE=\"" << env.minMz << "\" SIZE=\"5\">&nbsp;-&nbsp;";
    cout << "<INPUT TYPE=\"TEXT\" NAME=\"MaxMz\" VALUE=\"" << env.maxMz << "\" SIZE=\"5\"><BR>" << endl;

    cout << "<B>MassTol: &nbsp;Y-zoom:</B><BR>" << endl;
    cout << " &nbsp;<INPUT TYPE=\"TEXT\" NAME=\"MatchTol\" VALUE=\"" << env.matchTolerance << "\" SIZE=\"5\">" << endl;
    cout << " &nbsp;&nbsp;<INPUT TYPE=\"TEXT\" NAME=\"IntensityZoom\" VALUE=\"" << env.intensityZoom
         << "\" SIZE=\"5\"><BR>" << endl;

    cout << "<B>BlankPrecRegion</B><INPUT TYPE=\"CHECKBOX\" NAME=\"BlankNearParentRegion\" VALUE=\"1\""
         << (env.blankNearParentRegion == 1 ? " CHECKED" : "") << "><BR>" << endl;


    cout << "<BR><BR><B>ANNOTATION OPTIONS</B><BR>" << endl;

    cout << "<B>LabelType:</B><BR>" << endl;
    cout << "   <INPUT TYPE=\"RADIO\" NAME=\"LabelType\" VALUE=\"0\"" << (env.labelType == 0 ? " CHECKED" : "")
         << ">Ion";
    cout << "   <INPUT TYPE=\"RADIO\" NAME=\"LabelType\" VALUE=\"1\"" << (env.labelType == 1 ? " CHECKED" : "")
         << ">M/Z";
    cout << "   <INPUT TYPE=\"RADIO\" NAME=\"LabelType\" VALUE=\"2\"" << (env.labelType == 2 ? " CHECKED" : "")
         << ">Off<BR>" << endl;

    cout << "<B>NumPeaks: &nbsp;&nbsp;MinInten:</B><BR>" << endl;
    cout << " &nbsp;<INPUT TYPE=\"TEXT\" NAME=\"ShowNumPeaks\" VALUE=\"" << env.showNumPeaks << "\" SIZE=\"5\">"
         << endl;
    cout << " &nbsp;&nbsp;<INPUT TYPE=\"TEXT\" NAME=\"ShowMinIntensity\" VALUE=\"" << env.showMinIntensity
         << "\" SIZE=\"5\"><BR>" << endl;

    cout << "<B>Ions:</B><BR>" << endl;
    for (int c = 0; c < MAX_CHARGE && c < charge; c++) {
        cout << "<B>&nbsp;&nbsp;a<SUP>" << c + 1 << "+</SUP></B><INPUT TYPE=\"CHECKBOX\" NAME=\"ShowA" << c + 1
             << "\" VALUE=\"1\"" << (env.showA[c] == 1 ? " CHECKED" : "") << ">";
    }
    cout << "<BR>" << endl;

    for (int c = 0; c < MAX_CHARGE && c < charge; c++) {
        cout << "<B>&nbsp;&nbsp;b<SUP>" << c + 1 << "+</SUP></B><INPUT TYPE=\"CHECKBOX\" NAME=\"ShowB" << c + 1
             << "\" VALUE=\"1\"" << (env.showB[c] == 1 ? " CHECKED" : "") << ">";
    }
    cout << "<BR>" << endl;

    for (int c = 0; c < MAX_CHARGE && c < charge - 1; c++) {
        cout << "<B>&nbsp;&nbsp;c<SUP>" << c + 1 << "+</SUP></B><INPUT TYPE=\"CHECKBOX\" NAME=\"ShowC" << c + 1
             << "\" VALUE=\"1\"" << (env.showC[c] == 1 ? " CHECKED" : "") << ">";
    }
    cout << "<BR>" << endl;

    for (int c = 0; c < MAX_CHARGE && c < charge; c++) {
        cout << "<B>&nbsp;&nbsp;y<SUP>" << c + 1 << "+</SUP></B><INPUT TYPE=\"CHECKBOX\" NAME=\"ShowY" << c + 1
             << "\" VALUE=\"1\"" << (env.showY[c] == 1 ? " CHECKED" : "") << ">";
    }
    cout << "<BR>" << endl;

    for (int c = 0; c < MAX_CHARGE && c < charge - 1; c++) {
        cout << "<B>&nbsp;&nbsp;z<SUP>" << c + 1 << "+</SUP></B><INPUT TYPE=\"CHECKBOX\" NAME=\"ShowZ" << c + 1
             << "\" VALUE=\"1\"" << (env.showZ[c] == 1 ? " CHECKED" : "") << ">";
    }
    cout << "<BR>" << endl;

    cout << "<B>&nbsp;&nbsp;Prec losses</B> <INPUT TYPE=\"CHECKBOX\" NAME=\"ShowPLoss\" VALUE=\"1\""
         << (env.showPLoss == 1 ? " CHECKED" : "") << "><BR>" << endl;

    cout
            << "<B>&nbsp;&nbsp;-H<sub>2</sub>O/-NH<sub>3</sub>/-P</B> <INPUT TYPE=\"CHECKBOX\" NAME=\"ShowCommonLoss\" VALUE=\"1\""
            << (env.showCommonLoss == 1 ? " CHECKED" : "") << "><BR>" << endl;

    cout << "<B>&nbsp;&nbsp;All</B> <INPUT TYPE=\"CHECKBOX\" NAME=\"ShowAll\" VALUE=\"1\""
         << (env.showAll == 1 ? " CHECKED" : "") << ">" << endl;



    // the hidden fields
    cout << "<INPUT TYPE=\"HIDDEN\" NAME=\"LibFile\" VALUE=\"" << env.splibFileName << "\">" << endl;
    cout << "<INPUT TYPE=\"HIDDEN\" NAME=\"LibFileOffset\" VALUE=\"" << env.splibFileOffset << "\">" << endl;
    cout << "<INPUT TYPE=\"HIDDEN\" NAME=\"QueryFile\" VALUE=\"" << env.mzXMLFileName << "\">" << endl;
    cout << "<INPUT TYPE=\"HIDDEN\" NAME=\"QueryScanNum\" VALUE=\"" << env.mzXMLScanNum << "\">" << endl;
    cout << "<INPUT TYPE=\"HIDDEN\" NAME=\"OutputPath\" VALUE=\"" << env.outputPath << "\">" << endl;


    cout << "<CENTER>" << endl;
    cout << "    <INPUT TYPE=\"SUBMIT\" VALUE=\"GO\"><BR>" << endl;
    cout << "</CENTER>" << endl;
    cout << "</FORM></TD>" << endl;

    cout << "<TD>";

    Peptide *upPep = upLib->getPeptidePtr();

    if (upPep) {
        cout << "<TT><B>" << upPep->htmlStyle() << '/' << upPep->charge << "</B>  (M/Z = ";
    } else {
        cout << "<TT><B>" << upLib->getFullName() << "</B>  (M/Z = ";
    }

    cout.precision(3);
    cout << fixed << upLib->getPrecursorMz();

    double upProb = upLib->getProb(-1.0);
    if (upProb >= 0.0) {
        cout.precision(4);
        cout << ", P = " << fixed << upProb;
    }

    unsigned int upNreps = upLib->getNrepsUsed(0);
    if (upNreps > 0) {
        cout << ", N = " << upNreps;
    }

    string specStr("");
    if (upLib->getOneComment("Spec", specStr)) {
        cout << ", " << specStr;
    }

    if (upLib->getStatus() != "Normal") {
        cout << " **" << upLib->getStatus() << "**";
    }

    cout << ")";

    int proteinCount = 0;
    string protein = upLib->getFirstProtein(proteinCount);
    if (!(protein.empty())) {
        cout << " from " << protein;
        if (proteinCount > 1) {
            cout << " (+" << proteinCount - 1 << ")";
        }
    }

    cout << "</TT><BR>" << endl;


    string webroot("");

    const char *try_webroot = getWebserverRoot();
    if (try_webroot) {
        webroot = try_webroot;
    }

#ifdef ATLAS
    webroot = ATLAS_WEBSERVER_ROOT;
    string relpath(env.pngFileName);
    if (!webroot.empty() && env.pngFileName.substr(0, webroot.length()) == webroot) {
      if (webroot[webroot.length()-1] == '/') {
        relpath = env.pngFileName.substr(webroot.length()-1);
      } else {
        relpath = env.pngFileName.substr(webroot.length());
      }
    }


#else

    std::string relpath(env.pngFileName);
    translate_absolute_filesystem_path_to_relative_webserver_root_path(relpath);

#endif

    cout << "<IMG SRC=\"" << relpath << "\" BORDER=0><BR>" << endl;

    if (query) {
        cout << "<TT><B>" << query->getName() << "</B>  (M/Z = " << query->getPrecursorMz() << ")</TT><BR>" << endl;
    } else if (downLib) {

        Peptide *downPep = downLib->getPeptidePtr();

        if (downPep) {
            cout << "<TT><B>" << downPep->htmlStyle() << '/' << downPep->charge << "</B>  (M/Z = ";
        } else {
            cout << "<TT><B>" << downLib->getFullName() << "</B>  (M/Z = ";
        }

        cout.precision(3);
        cout << fixed << downLib->getPrecursorMz();

        double downProb = downLib->getProb(-1.0);
        if (downProb >= 0.0) {
            cout.precision(4);
            cout << ", P = " << fixed << downProb;
        }

        unsigned int downNreps = downLib->getNrepsUsed(0);
        if (downNreps > 0) {
            cout << ", N = " << downNreps;
        }

        string specStr("");
        if (downLib->getOneComment("Spec", specStr)) {
            cout << ", " << specStr;
        }

        if (downLib->getStatus() != "Normal") {
            cout << " **" << downLib->getStatus() << "**";
        }

        cout << ")";

        int proteinCount = 0;
        string protein = downLib->getFirstProtein(proteinCount);
        if (!(protein.empty())) {
            cout << " from " << protein;
            if (proteinCount > 1) {
                cout << " (+" << proteinCount - 1 << ")";
            }
        }

        cout << "</TT><BR>" << endl;

    }

    cout << "</TD></TR>" << endl;

    if (upPep) {
        cout << "<TR><TD ALIGN=MIDDLE COLSPAN=2>" << endl;
        printIonTable(upPep, labels, env);
        cout << "</TD></TR>" << endl;
    }

    cout << "<TR><TD ALIGN=MIDDLE COLSPAN=2>" << endl;
    printComments(upLib);
    cout << "</TD></TR>" << endl;

    if (downLib) {
        if (downLib->getPeptidePtr()) {
            map<string, Peak> downLabels;
            for (map<string, Peak>::iterator l = labels.begin(); l != labels.end(); l++) {
                if (l->first[0] == '~') {
                    pair<string, Peak> p;
                    p.first = l->first.substr(1);
                    p.second = l->second;
                    downLabels.insert(p);
                }
            }

            cout << "<TR><TD ALIGN=MIDDLE COLSPAN=2>" << endl;
            printIonTable(downLib->getPeptidePtr(), downLabels, env);
            cout << "</TD></TR>" << endl;
        }

        cout << "<TR><TD ALIGN=MIDDLE COLSPAN=2>" << endl;
        printComments(downLib);
        cout << "</TD></TR>" << endl;
    }

    cout << "</TABLE>" << endl;
    cout << "</BODY></HTML>" << endl;

}

// printIonTable - prints the ion table
void printIonTable(Peptide *pep, map<string, Peak> &labels, PlotSpectraSTEnv &env) {

    int charge = pep->charge;

    cout << endl;
    cout << "   <TABLE BORDER=0 CELLPADDING=\"2\">" << endl;
    cout << "      <TR BGCOLOR=\"" << BGCOLOR << "\" ALIGN=\"CENTER\">";

    for (int c = 0; c < MAX_CHARGE && c < charge; c++) {
        if (env.showA[c] || env.showAll) cout << "<TH><B>a<SUP>" << c + 1 << "+</SUP></B></TH>" << endl;
    }
    for (int c = 0; c < MAX_CHARGE && c < charge; c++) {
        if (env.showB[c] || env.showAll) cout << "<TH><B>b<SUP>" << c + 1 << "+</SUP></B></TH>" << endl;
    }
    for (int c = 0; c < MAX_CHARGE && c < charge - 1; c++) {
        if (env.showC[c] || env.showAll) cout << "<TH><B>c<SUP>" << c + 1 << "+</SUP></B></TH>" << endl;
    }

    cout << "<TH><B>#</B></TH>" << endl;
    cout << "<TH><B>AA</B></TH>" << endl;
    cout << "<TH><B>#</B></TH>" << endl;

    for (int c = 0; c < MAX_CHARGE && c < charge; c++) {
        if (env.showY[c] || env.showAll) cout << "<TH><B>y<SUP>" << c + 1 << "+</SUP></B></TH>" << endl;
    }
    for (int c = 0; c < MAX_CHARGE && c < charge - 1; c++) {
        if (env.showZ[c] || env.showAll) cout << "<TH><B>z<SUP>" << c + 1 << "+</SUP></B></TH>" << endl;
    }
    cout << "</TR>" << endl;

    string::size_type len = pep->stripped.length();

    for (string::size_type i = 0; i < len; i++) {
        cout << "      <TR ALIGN=\"RIGHT\">" << endl;

        if (i == len - 1) {
            // last line, just print empty space
            for (int c = 0; c < MAX_CHARGE && c < charge; c++) {
                if (env.showA[c] || env.showAll) {
                    cout << "<TD BGCOLOR=\"" << UNLABEL_BKGRND << "\"><TT>&nbsp;</TT></TD>" << endl;
                }
            }
            for (int c = 0; c < MAX_CHARGE && c < charge; c++) {
                if (env.showB[c] || env.showAll) {
                    cout << "<TD BGCOLOR=\"" << UNLABEL_BKGRND << "\"><TT>&nbsp;</TT></TD>" << endl;
                }
            }
            for (int c = 0; c < MAX_CHARGE && c < charge - 1; c++) {
                if (env.showC[c] || env.showAll) {
                    cout << "<TD BGCOLOR=\"" << UNLABEL_BKGRND << "\"><TT>&nbsp;</TT></TD>" << endl;
                }
            }

        } else {
            stringstream aion;
            aion << 'a' << (i + 1);

            stringstream bion;
            bion << 'b' << (i + 1);

            stringstream cion;
            cion << 'c' << (i + 1);

            map<string, Peak>::iterator found;

            for (int c = 0; c < MAX_CHARGE && c < charge; c++) {
                if (env.showA[c] || env.showAll) {
                    double mz = pep->monoisotopicMZFragment('a', i + 1, c + 1);
                    stringstream labelss;
                    labelss << aion.str();
                    if (c > 0) labelss << '^' << c + 1;
                    if ((found = labels.find(labelss.str())) != labels.end()) {
                        if (found->second.info == "M") {
                            cout << "<TD BGCOLOR=\"" << QUERY_LABEL_BKGRND << "\"><TT><B>";
                        } else {
                            cout << "<TD BGCOLOR=\"" << LIB_LABEL_BKGRND << "\"><TT><B>";
                        }
                        cout.precision(4);
                        cout << fixed << mz << "</B></TT></TD>" << endl;

                    } else {
                        cout << "<TD BGCOLOR=\"" << UNLABEL_BKGRND << "\"><TT>";
                        cout.precision(4);
                        cout << fixed << mz << "</TT></TD>" << endl;
                    }
                }
            }

            for (int c = 0; c < MAX_CHARGE && c < charge; c++) {

                if (env.showB[c] || env.showAll) {
                    double mz = pep->monoisotopicMZFragment('b', i + 1, c + 1);
                    stringstream labelss;
                    labelss << bion.str();
                    if (c > 0) labelss << '^' << c + 1;
                    if ((found = labels.find(labelss.str())) != labels.end()) {
                        if (found->second.info == "M") {
                            cout << "<TD BGCOLOR=\"" << QUERY_LABEL_BKGRND << "\"><TT><B>";
                        } else {
                            cout << "<TD BGCOLOR=\"" << LIB_LABEL_BKGRND << "\"><TT><B>";
                        }
                        cout.precision(4);
                        cout << fixed << mz << "</B></TT></TD>" << endl;

                    } else {
                        cout << "<TD BGCOLOR=\"" << UNLABEL_BKGRND << "\"><TT>";
                        cout.precision(4);
                        cout << fixed << mz << "</TT></TD>" << endl;
                    }
                }
            }

            for (int c = 0; c < MAX_CHARGE && c < charge - 1; c++) {

                if (env.showC[c] || env.showAll) {
                    double mz = pep->monoisotopicMZFragment('c', i + 1, c + 1);
                    stringstream labelss;
                    labelss << cion.str();
                    if (c > 0) labelss << '^' << c + 1;
                    if ((found = labels.find(labelss.str())) != labels.end()) {
                        if (found->second.info == "M") {
                            cout << "<TD BGCOLOR=\"" << QUERY_LABEL_BKGRND << "\"><TT><B>";
                        } else {
                            cout << "<TD BGCOLOR=\"" << LIB_LABEL_BKGRND << "\"><TT><B>";
                        }
                        cout.precision(4);
                        cout << fixed << mz << "</B></TT></TD>" << endl;

                    } else {
                        cout << "<TD BGCOLOR=\"" << UNLABEL_BKGRND << "\"><TT>";
                        cout.precision(4);
                        cout << fixed << mz << "</TT></TD>" << endl;
                    }
                }
            }
        }
        // middle part
        cout << "<TD ALIGH=\"CENTER\" BGCOLOR=\"" << BGCOLOR << "\"><TT><B>" << (i + 1) << "</B></TT></TD>";
        cout << "<TD ALIGN=\"CENTER\" BGCOLOR=\"" << BGCOLOR << "\"><TT><B>" << pep->interactStyleAA(i)
             << "</B></TT></TD>";
        cout << "<TD ALIGN=\"CENTER\" BGCOLOR=\"" << BGCOLOR << "\"><TT><B>" << len - i << "</B></TT></TD>";

        // right part

        if (i == 0) {
            // first line, just print empty space
            for (int c = 0; c < MAX_CHARGE && c < charge; c++) {
                if (env.showY[c] || env.showAll) {
                    cout << "<TD BGCOLOR=\"" << UNLABEL_BKGRND << "\"><TT>&nbsp;</TT></TD>" << endl;
                }
            }
            for (int c = 0; c < MAX_CHARGE && c < charge - 1; c++) {
                if (env.showZ[c] || env.showAll) {
                    cout << "<TD BGCOLOR=\"" << UNLABEL_BKGRND << "\"><TT>&nbsp;</TT></TD>" << endl;
                }
            }

        } else {
            stringstream yion;
            yion << 'y' << len - i;
            stringstream zion;
            zion << 'z' << len - i;


            map<string, Peak>::iterator found;

            for (int c = 0; c < MAX_CHARGE && c < charge; c++) {
                if (env.showY[c] || env.showAll) {
                    double mz = pep->monoisotopicMZFragment('y', len - i, c + 1);
                    stringstream labelss;
                    labelss << yion.str();
                    if (c > 0) labelss << '^' << c + 1;
                    if ((found = labels.find(labelss.str())) != labels.end()) {
                        if (found->second.info == "M") {
                            cout << "<TD BGCOLOR=\"" << QUERY_LABEL_BKGRND << "\"><TT><B>";
                        } else {
                            cout << "<TD BGCOLOR=\"" << LIB_LABEL_BKGRND << "\"><TT><B>";
                        }
                        cout.precision(4);
                        cout << fixed << mz << "</B></TT></TD>" << endl;

                    } else {
                        cout << "<TD BGCOLOR=\"" << UNLABEL_BKGRND << "\"><TT>";
                        cout.precision(4);
                        cout << fixed << mz << "</TT></TD>" << endl;
                    }
                }
            }

            for (int c = 0; c < MAX_CHARGE && c < charge - 1; c++) {
                if (env.showZ[c] || env.showAll) {
                    double mz = pep->monoisotopicMZFragment('z', len - i, c + 1);
                    stringstream labelss;
                    labelss << zion.str();
                    if (c > 0) labelss << '^' << c + 1;
                    if ((found = labels.find(labelss.str())) != labels.end()) {
                        if (found->second.info == "M") {
                            cout << "<TD BGCOLOR=\"" << QUERY_LABEL_BKGRND << "\"><TT><B>";
                        } else {
                            cout << "<TD BGCOLOR=\"" << LIB_LABEL_BKGRND << "\"><TT><B>";
                        }
                        cout.precision(4);
                        cout << fixed << mz << "</B></TT></TD>" << endl;

                    } else {
                        cout << "<TD BGCOLOR=\"" << UNLABEL_BKGRND << "\"><TT>";
                        cout.precision(4);
                        cout << fixed << mz << "</TT></TD>" << endl;
                    }
                }
            }
        }

        cout << "</TR>" << endl;
    }

    cout << "    </TABLE>" << endl << endl;
}

// printComments - prints the Sequence Search Info table and Sample Source Info table
void printComments(SpectraSTLibEntry *lib) {


    map<char, map<string, pair<double, double> > *> seq;
    lib->getSeqInfo(seq);

    if (!seq.empty()) {

        cout << "   <TABLE BORDER=0 CELLPADDING=\"2\">" << endl;
        cout << "      <THEAD><TR BGCOLOR=\"" << BGCOLOR << "\" ALIGN=\"CENTER\">";
        cout << "<TH BGCOLOR=\"" << HEADERCELLCOLOR << "\" COLSPAN=\"3\"><TT><B>"
             << "Library Spectrum Sequence Search Information (Found by " << seq.size()
             << " search engine(s))</B></TT></TH>";
        cout << "</TR></THEAD>" << endl;

        cout << "       <TBODY>";
        cout << "<TR BGCOLOR=\"" << SUBHEADERCELLCOLOR << "\" ALIGN=\"CENTER\">";
        cout << "<TH BGCOLOR=\"" << SUBHEADERCELLCOLOR << "\"><TT>" << "Engine" << "</TT></TH>";
        cout << "<TH BGCOLOR=\"" << SUBHEADERCELLCOLOR << "\"><TT>" << "Score" << "</TT></TH>";
        cout << "<TH BGCOLOR=\"" << SUBHEADERCELLCOLOR << "\"><TT>" << "Ave/StDev" << "</TT></TH>";
        cout << "</TR>" << endl;

        for (map<char, map<string, pair<double, double> > *>::iterator en = seq.begin(); en != seq.end(); en++) {
            string engine("Unknown");
            if (en->first == 'S') engine = "SEQUEST";
            if (en->first == 'M') engine = "Mascot";
            if (en->first == 'X') engine = "X!Tandem (Native)";
            if (en->first == 'K') engine = "X!Tandem (K-score)";
            if (en->first == 'O') engine = "OMSSA";
            if (en->first == 'P') engine = "Phenyx";


            map<string, pair<double, double> > *scs = en->second;

            if (!scs) continue;

            int numSpectra = 0;
            map<string, pair<double, double> >::iterator sc = scs->find("nr");
            if (sc != scs->end()) {
                numSpectra = (int) (sc->second.first + 0.5);
            }

            cout << "<TR BGCOLOR=\"" << DATACELLCOLOR << "\" ALIGN=\"CENTER\">";
            cout << "<TD BGCOLOR=\"" << DATACELLCOLOR << "\"><TT>" << engine << "(" << numSpectra
                 << " spectra)</TT></TD>";

            vector<string> scoreNames;
            vector<pair<double, double> > scoreValues;
            for (sc = scs->begin(); sc != scs->end(); sc++) {
                if (sc->first == "nr") continue;
                if (sc->first[0] == 'b') continue;
                scoreNames.push_back(sc->first);
                scoreValues.push_back(sc->second);
            }

            cout << "<TD BGCOLOR=\"" << DATACELLCOLOR << "\"><TT>";
            for (vector<string>::iterator s = scoreNames.begin(); s != scoreNames.end(); s++) {
                cout << translateScoreName((*s), engine) << "<BR>";
            }
            cout << "</TT></TD>";

            cout << "<TD BGCOLOR=\"" << DATACELLCOLOR << "\"><TT>";
            for (vector<pair<double, double> >::iterator v = scoreValues.begin(); v != scoreValues.end(); v++) {
                double mean = 0.0;
                if (v->first > 0.00001 || v->first < 0.00001) mean = v->first / (double) numSpectra;
                double variance = 0.0;
                if (v->second > 0.00001) variance = v->second / (double) numSpectra - mean * mean;
                double stdev = 0.0;
                if (variance > 0.00001) stdev = sqrt(variance);
                cout.precision(2);
                cout << mean << '/' << stdev << "<BR>";
            }
            cout << "</TT></TD>";


            cout << "</TR>" << endl;

        }
        cout << "</TBODY></TABLE><BR>" << endl;

    }

    map<string, pair<unsigned int, unsigned int> > sample;
    lib->getSampleInfo(sample);

    if (!sample.empty()) {

        cout << "   <TABLE BORDER=0 CELLPADDING=\"2\">" << endl;

        cout << "      <THEAD><TR BGCOLOR=\"" << BGCOLOR << "\" ALIGN=\"CENTER\">";
        cout << "<TH BGCOLOR=\"" << HEADERCELLCOLOR << "\" COLSPAN=\"3\"><TT><B>" << "Library Spectrum Source Samples ("
             << sample.size() << ")</B></TT></TH>";
        cout << "</TR></THEAD>" << endl;

        cout << "      <TBODY>" << endl;
        cout << "<TR BGCOLOR=\"" << SUBHEADERCELLCOLOR << "\" ALIGN=\"CENTER\">";
        cout << "<TH BGCOLOR=\"" << SUBHEADERCELLCOLOR << "\"><TT>" << "Sample" << "</TT></TH>";
        cout << "<TH BGCOLOR=\"" << SUBHEADERCELLCOLOR << "\"><TT>" << "# replicates used for consensus"
             << "</TT></TH>";
        cout << "<TH BGCOLOR=\"" << SUBHEADERCELLCOLOR << "\"><TT>" << "# total replicates" << "</TT></TH>";
        cout << "</TR>" << endl;

        unsigned int totAll = 0;
        unsigned int totConsensus = 0;
        for (map<string, pair<unsigned int, unsigned int> >::iterator sa = sample.begin(); sa != sample.end(); sa++) {

            cout << "<TR BGCOLOR=\"" << DATACELLCOLOR << "\" ALIGN=\"CENTER\">";
            cout << "<TD BGCOLOR=\"" << DATACELLCOLOR << "\"><TT>" << sa->first << "</TT></TD>";
            cout << "<TD BGCOLOR=\"" << DATACELLCOLOR << "\"><TT>" << sa->second.first << "</TT></TD>";
            cout << "<TD BGCOLOR=\"" << DATACELLCOLOR << "\"><TT>" << sa->second.second << "</TT></TD>";
            cout << "</TR>" << endl;

            totConsensus += sa->second.first;
            totAll += sa->second.second;

        }

        cout << "<TR BGCOLOR=\"" << SUBHEADERCELLCOLOR << "\" ALIGN=\"CENTER\">";
        cout << "<TD BGCOLOR=\"" << SUBHEADERCELLCOLOR << "\"><TT>" << "TOTAL" << "</TT></TD>";
        cout << "<TD BGCOLOR=\"" << SUBHEADERCELLCOLOR << "\"><TT>" << totConsensus << "</TT></TD>";
        cout << "<TD BGCOLOR=\"" << SUBHEADERCELLCOLOR << "\"><TT>" << totAll << "</TT></TD>";
        cout << "</TR>" << endl;

        cout << "</TBODY></TABLE><BR>" << endl;


    }
}

// translateScoreName - translates the two-letter code in the Se= comment field of NIST library to what it means
// (TODO: some are problematic!)
string translateScoreName(string scoreName, string engine) {

    if (engine == "SEQUEST") {
        if (scoreName == "xc") return "XCorr";
        if (scoreName == "sc") return "XCorr";
        if (scoreName == "dc") return "DeltaCn";
        if (scoreName == "ps") return "Sp";
        if (scoreName == "pr") return "SpRank";
        if (scoreName == "fv") return "F-value (PeptideProphet)";
        if (scoreName == "pb") return "Prob (PeptideProphet)";

    } else if (engine == "Mascot") {
        if (scoreName == "sc") return "Ion score";
        if (scoreName == "sd") return "Score(top) - Score(second)";
        if (scoreName == "sr") return "Ion score - Homology score";
        if (scoreName == "td") return "Score(top) - Score(best tryptic)";
        if (scoreName == "is") return "Identity score";
        if (scoreName == "hs") return "Homology score";
        if (scoreName == "fv") return "F-value (PeptideProphet)";
        if (scoreName == "pb") return "Prob (PeptideProphet)";
        if (scoreName == "ex") return "Expectation value";


    } else if (engine == "X!Tandem (K-score)" || engine == "X!Tandem (Native)") {
        if (scoreName == "pr") return "Reported prob";
        if (scoreName == "ex") return "Expectation value";
        if (scoreName == "sd") return "Score(top) - Score(second)";
        if (scoreName == "hs") return "Hyperscore";
        if (scoreName == "ns") return "Nextscore";
        if (scoreName == "td") return "Score(top) - Score(best tryptic)";
        if (scoreName == "fv") return "F-value (PeptideProphet)";
        if (scoreName == "pb") return "Prob (PeptideProphet)";


    } else if (engine == "OMSSA") {
        if (scoreName == "pr") return "Reported prob";
        if (scoreName == "ex") return "Expectation value";
        if (scoreName == "td") return "Score(top) - Score(best tryptic)";
        if (scoreName == "fv") return "F-value (PeptideProphet)";
        if (scoreName == "pb") return "Prob (PeptideProphet)";

    } else if (engine == "Phenyx") {
        if (scoreName == "sc") return "Score";
        if (scoreName == "dc") return "Delta";
        if (scoreName == "fv") return "F-value (PeptideProphet)";
        if (scoreName == "pb") return "Prob (PeptideProphet)";


    } else {
        return (scoreName);
    }
    return (scoreName);
}

// sortLabelsByMzAsc - comparator for sorting peaks
bool sortLabelsByMzAsc(pair<string, Peak> a, pair<string, Peak> b) {
    return (a.second.mz < b.second.mz);
}

#ifdef RUN_AS_CGI

// extractCGI - extracts the CGI query string and sets the environment variables in PlotSpectraSTEnv
void extractCGI(PlotSpectraSTEnv& env) {	
  char* requestType = getenv("REQUEST_METHOD");
  
  if(!requestType || strcmp(requestType, "GET")) {
    cout << " This program needs to be called with CGI GET method." << endl;
    exit(1);
  }
  
  // Decode GET method
  char* queryStr = getenv("QUERY_STRING");
  if (!queryStr) {
    cout << "GET query string empty." << endl;
    exit(1);
  }
  
  
  
  // default parameters
  env.splibFileName = "";
  env.splibFileOffset = 0;
  env.mzXMLFileName = "";
  env.mzXMLScanNum = 0;
  env.labelType = 0;
  env.minMz = 0.0;
  env.maxMz = 0.0;
  env.matchTolerance = DEFAULT_MATCH_TOLERANCE; 
  env.intensityZoom = 1.0;
  for (int c = 0; c < MAX_CHARGE; c++) {
    env.showA[c] = 0;
    env.showB[c] = 1;
    env.showC[c] = 0;
    env.showY[c] = 1;
    env.showZ[c] = 0;
  } 
  env.showPLoss = 1; // set to zero if too messy for publication-quality
  env.showCommonLoss = 0;
  env.showAll = 0;
  env.showMinIntensity = 50.0;
  env.showNumPeaks = 150;
  env.blankNearParentRegion = 0;
  
  env.pngFileName = ""; // will be set later by this script
  env.outputPath = "";
  
  bool firstCall = true; 	
  
  string qs(queryStr);
  
  string::size_type pos = 0;
  string param("");
  while ((param = nextToken(qs, pos, pos, "&\t\r\n")) != "") {
    string::size_type equalPos = 0;
    string name = nextToken(param, equalPos, equalPos, "=\t\r\n");
    name = processCGIStr(name);
    string value = nextToken(param, equalPos + 1, equalPos);
    value = processCGIStr(value);
    
    if (name == "LibFile") {
      env.splibFileName = value;
    } else if (name == "LibFileOffset") {
      env.splibFileOffset = strtoull(value.c_str(), NULL, 10);
    } else if (name == "QueryFile") {
      env.mzXMLFileName = value;
    } else if (name == "QueryScanNum") {
      env.mzXMLScanNum = strtoull(value.c_str(), NULL, 10);
    } else if (name == "LabelType") {
      env.labelType = atoi(value.c_str());
      if (env.labelType < 0 || env.labelType > 2) env.labelType = 0;
    } else if (name == "MinMz") {
      env.minMz = atof(value.c_str());
      if (env.minMz < 0.0 || env.minMz > 5000.0) env.minMz = 0.0; // will set by default based on spectra
    } else if (name == "MaxMz") {
      env.maxMz = atof(value.c_str());
      if (env.maxMz < env.minMz + 10.0 || env.maxMz > 5000.0) env.maxMz = 0.0; // will set by default based on spectra
    } else if (name == "MatchTol") {
      firstCall = false;
      env.matchTolerance = atof(value.c_str());
      if (env.matchTolerance < 0.01 || env.matchTolerance > 100.0) env.matchTolerance = DEFAULT_MATCH_TOLERANCE;
    } else if (name == "IntensityZoom") {
      firstCall = false;
      env.intensityZoom = atof(value.c_str());
      if (env.intensityZoom < 0.01 || env.intensityZoom > 100.0) env.intensityZoom = 1.0; // no zoom
    } else if (name == "ShowPLoss") {
      env.showPLoss = atoi(value.c_str());
      if (env.showPLoss < 0) env.showPLoss = 0;
      if (env.showPLoss > 1) env.showPLoss = 1;
    } else if (name == "ShowCommonLoss") {
      env.showCommonLoss = atoi(value.c_str());
      if (env.showCommonLoss < 0) env.showCommonLoss = 0;
      if (env.showCommonLoss > 1) env.showCommonLoss = 1;
    } else if (name == "ShowNumPeaks") {
      env.showNumPeaks = atoi(value.c_str());
      if (env.showNumPeaks < 0) env.showNumPeaks = 0;
    } else if (name == "ShowMinIntensity") {
      env.showMinIntensity = atof(value.c_str());
      if (env.showMinIntensity < 0.00001) env.showMinIntensity = 0.0;
    } else if (name == "ShowAll") {
      env.showAll = atoi(value.c_str());
      if (env.showAll < 0) env.showAll = 0;
      if (env.showAll > 1) env.showAll = 1;
      
    } else if (name.substr(0, 5) == "ShowB") {
      int ch = atoi(name.substr(5).c_str());
      if (ch > 0 && ch <= MAX_CHARGE) {
        env.showB[ch - 1] = atoi(value.c_str());
        if (env.showB[ch - 1] < 0) env.showB[ch - 1] = 0;
        if (env.showB[ch - 1] > 1) env.showB[ch - 1] = 1;
      }
    } else if (name.substr(0, 5) == "ShowY") { 
      int ch = atoi(name.substr(5).c_str());
      if (ch > 0 && ch <= MAX_CHARGE) {
        env.showY[ch - 1] = atoi(value.c_str());
        if (env.showY[ch - 1] < 0) env.showY[ch - 1] = 0;
        if (env.showY[ch - 1] > 1) env.showY[ch - 1] = 1;
      }
    } else if (name.substr(0, 5) == "ShowZ") { 
      int ch = atoi(name.substr(5).c_str());
      if (ch > 0 && ch <= MAX_CHARGE) {
        env.showZ[ch - 1] = atoi(value.c_str());
        if (env.showZ[ch - 1] < 0) env.showZ[ch - 1] = 0;
        if (env.showZ[ch - 1] > 1) env.showZ[ch - 1] = 1;
      }
    } else if (name.substr(0, 5) == "ShowA") {
      int ch = atoi(name.substr(5).c_str());
      if (ch > 0 && ch <= MAX_CHARGE) {
        env.showA[ch - 1] = atoi(value.c_str());
        if (env.showA[ch - 1] < 0) env.showA[ch - 1] = 0;
        if (env.showA[ch - 1] > 1) env.showA[ch - 1] = 1;
      }
    } else if (name.substr(0, 5) == "ShowC") {
      int ch = atoi(name.substr(5).c_str());
      if (ch > 0 && ch <= MAX_CHARGE) {
        env.showC[ch - 1] = atoi(value.c_str());
        if (env.showC[ch - 1] < 0) env.showC[ch - 1] = 0;
        if (env.showC[ch - 1] > 1) env.showC[ch - 1] = 1;
      }
    } else if (name == "BlankNearParentRegion") {
      env.blankNearParentRegion = atoi(value.c_str());
      if (env.blankNearParentRegion < 0) env.blankNearParentRegion = 0;
      if (env.blankNearParentRegion > 1) env.blankNearParentRegion = 1;
    } else if (name == "OutputPath") {
      env.outputPath = value;
    }
    
    
    pos++;
  }
  
  if (!firstCall) {
    // called when the user click on "Go" in a previous Plotspectrast page

    for (int c = 0; c < MAX_CHARGE; c++) {
      stringstream showass;
      showass << "ShowA" << c + 1;
      if (qs.find(showass.str(), 0) == string::npos) env.showA[c] = 0;
      stringstream showbss;
      showbss << "ShowB" << c + 1;
      if (qs.find(showbss.str(), 0) == string::npos) env.showB[c] = 0;
      stringstream showyss;
      showyss << "ShowY" << c + 1;
      if (qs.find(showyss.str(), 0) == string::npos) env.showY[c] = 0;
      stringstream showzss;
      showyss << "ShowZ" << c + 1;
      if (qs.find(showzss.str(), 0) == string::npos) env.showZ[c] = 0;
      stringstream showcss;
      showyss << "ShowC" << c + 1;
      if (qs.find(showcss.str(), 0) == string::npos) env.showC[c] = 0;
    }
    if (qs.find("ShowPLoss", 0) == string::npos) env.showPLoss = 0;
    if (qs.find("ShowCommonLoss", 0) == string::npos) env.showCommonLoss = 0;
    if (qs.find("ShowAll", 0) == string::npos) env.showAll = 0;
    if (qs.find("ShowMinIntensity", 0) == string::npos) env.showMinIntensity = 50.0;
    if (qs.find("ShowNumPeaks", 0) == string::npos) env.showNumPeaks = 150;
    if (qs.find("BlankNearParentRegion", 0) == string::npos) env.blankNearParentRegion = 0;
    
  }
  
  if (env.splibFileName.empty() || env.mzXMLFileName.empty()) {
    cout << "Missing library or query spectra source file information. Exiting." << endl;
    exit (1);
  }
  
  
}

// processCGIStr - transcode the CGI query string
string processCGIStr(string s) {

  string result("");
  
  for (string::size_type i = 0; i < s.length(); i++) {
    if (s[i] == '%') {
      char c;
      c = s[i+1] >= 'A' ? ((s[i+1] & 0xdf) - 'A') + 10 : (s[i+1] - '0');
      c *= 16;
      c += s[i+2] >= 'A' ? ((s[i+2] & 0xdf) - 'A') + 10 : (s[i+2] - '0');
      result += c;
      i += 2;			
    } else if (s[i] == '+') {
      result += ' ';
    } else {
      result += s[i];
    }
  } 
  return (result);
}
#endif // RUN_AS_CGI



