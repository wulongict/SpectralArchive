#include "SpectraSTLib.hpp"
#include "SpectraSTSearchTask.hpp"
#include "SpectraSTSearchParams.hpp"
#include "SpectraSTCreateParams.hpp"
#include "SpectraSTFileList.hpp"
#include "SpectraSTLog.hpp"
#include "SpectraSTConstants.hpp"
#include "FileUtils.hpp"
#include "Peptide.hpp"

#ifndef STANDALONE_LINUX
#include "common/hooks_tpp.h"
#endif

#include <string>
#include <time.h>
#include <fstream>
#include <sstream>

#include <math.h>

#include "SpectraSTFastaFileHandler.hpp"

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

/* Class: SpectraSTMain
 * 
 * main function for SpectraST.
 * 
 */

using namespace std;

static void doCreate(SpectraSTCreateParams &createParams, vector<string> &fileNames);

static unsigned int doSearch(SpectraSTSearchParams &searchParams, vector<string> &fileNames);

static void printUsage();

static void readUserModFile(string &modFileName);

// Verbose and quiet option flags. Don't want to pass them everywhere, so use global variables
bool g_verbose;
bool g_quiet;

// global pointer to a log object for keeping a log file
SpectraSTLog *g_log;

#ifndef LIB_BUILD

// main
int main(int argc, char *argv[]) {

#ifndef STANDALONE_LINUX
    hooks_tpp handler(argc,argv); // set up install paths etc
#endif

    // remember the time
    time_t startTime = time(NULL);
    char *startTimeStr = ctime(&startTime);
    startTimeStr[strlen(startTimeStr) - 1] = '\0';
    string logFileName("spectrast.log");
    string modFileName("spectrast.usermods");

    Peptide::defaultTables();

    // remember the command line to be put in the log file
    stringstream commandLiness;
    for (int ar = 0; ar < argc; ar++) {
        commandLiness << argv[ar] << ' ';
    }


    char mode = '\0';
    g_verbose = false;
    g_quiet = false;

    SpectraSTCreateParams createParams;
    SpectraSTSearchParams searchParams;

    if (argc < 2) {
        printUsage();
        return (1);
    }

    // get all options
    int i;
    string expectArg = "";

    for (i = 1; i < argc; i++) {

        if (argv[i][0] == '-') {

            if (strlen(argv[i]) <= 1) continue; // a '-' by itself. Ignore.

            // it's an option, the char immediately following the '-' tells us what
            // kind of option it is
            char optionType = argv[i][1];

            if (optionType == 's' || optionType == 'c') {
                // the mode (Search or Create) is set once the program sees an option
                // of that type. e.g. if it sees a -sF option, the program will assume
                // it's in Search mode. Thus, no need to explicitly set the mode in
                // most cases, although an option of -s or -c is also legal and simply tells the program to
                // run in that mode using all defaults.

                if (!mode) {
                    mode = optionType;
                } else if (mode != optionType) {
                    // already chosen another mode!
                    printUsage();
                    return (1);
                }

            } else if (optionType == 'V') {
                g_verbose = true;
                // note: verbose overrides quiet
                if (g_quiet) g_quiet = false;
                expectArg = "";

            } else if (optionType == 'Q' && !g_verbose) {
                g_quiet = true;
                expectArg = "";

            } else if (optionType == 'L') {
                // user-specified log file name
                logFileName = strlen(argv[i]) > 2 ? (argv[i] + 2) : "";
                if (logFileName.empty()) {
                    expectArg = "L";
                } else {
                    expectArg = "";
                    //fixpath(logFileName);
                }
            } else if (optionType == 'M') {
                // user-specified usermods file name
                modFileName = strlen(argv[i]) > 2 ? (argv[i] + 2) : "";
                if (modFileName.empty()) {
                    expectArg = "M";
                } else {
                    expectArg = "";
                    //fixpath(modFileName);
                }
            }

            if (optionType == 's' || optionType == 'c') {
                // sub-options, in the form -s? or -c? where ? is the sub-option of mode s or c
                string subop = strlen(argv[i]) > 2 ? (argv[i] + 2) : "";

                if (subop.empty()) continue; // bare -s or -c. Ignore.

                if (mode == 's') {

                    if (subop[0] == '_' && subop.length() <= 3) {
                        searchParams.printAdvancedOptions(cerr);
                        return (1);
                    }
                    if (!searchParams.addOption(subop)) {
                        expectArg = "s" + subop;
                    } else {
                        expectArg = "";
                    }

                } else if (mode == 'c') {

                    if (subop[0] == '_' && subop.length() <= 3) {
                        createParams.printAdvancedOptions(cerr);
                        return (1);
                    }
                    if (!createParams.addOption(subop)) {
                        expectArg = "c" + subop;
                    } else {
                        expectArg = "";
                    }

                }
            }

        } else {
            // not an option.
            if (expectArg.empty()) {
                // the previous option is not expecting an argument.
                // since all options should preceed files, there should be no more options
                break;
            } else {
                string optionValue(argv[i]);

                if (expectArg == "L") {
                    logFileName = optionValue;
                } else if (expectArg == "M") {
                    modFileName = optionValue;
                } else if (expectArg[0] == 's') {
                    bool dummy = searchParams.addOption(expectArg.substr(1) + optionValue);
                } else if (expectArg[0] == 'c') {
                    bool dummy = createParams.addOption(expectArg.substr(1) + optionValue);
                }
                expectArg = "";

            }
        }
    }

    // create the log object
    g_log = new SpectraSTLog(logFileName);

    // log the command line and the start time
    g_log->log("START", commandLiness.str(), startTimeStr);

    // log the file offset for this machine -- this defines whether the library can be bigger than ~2GB
    stringstream offsetSizess;
    offsetSizess << "File offset size is " << sizeof(fstream::off_type) << " bytes.";
    if (sizeof(fstream::off_type) >= 8) {
        offsetSizess << " Big library supported.";
    }
    g_log->log("GENERAL", offsetSizess.str());

    stringstream ptrSizess;
    ptrSizess << "Pointer size is " << sizeof(int *) << " bytes.";
    g_log->log("GENERAL", ptrSizess.str());

    // read usermod file
    readUserModFile(modFileName);

    if (!g_quiet) {
        cout << "SpectraST started at " << startTimeStr << '.' << endl;
    }
    if (g_verbose) {
        cout << "VERBOSE MODE..." << endl;
    }

    // the index i should now point to the first non-option argument
    int firstNonOption = i;

    // there should be at least one more file following the options
    if (firstNonOption >= argc) {
        printUsage();
        return (1);
    }

    vector<string> fileNames;
    for (int j = firstNonOption; j < argc; j++) {
        string fn(argv[j]);
        //fixpath(fn);
        fileNames.push_back(fn);
    }

    unsigned int searchCount = 0;

    if (mode == 'c') {

        doCreate(createParams, fileNames);

    } else if (mode == 's') {

        searchCount = doSearch(searchParams, fileNames);

    } else {
        // no mode. must be list
        if (SpectraSTFileList::isFileList(fileNames) == 1) {
            // is indeed list. instantiating the SpectraSTFileList object will carry out the required actions (create or search)
            SpectraSTFileList *fileList = new SpectraSTFileList(fileNames);
            // in case it's a search, get the search count
            searchCount = fileList->getSearchCount();
            delete fileList;
        } else {
            g_log->error("GENERAL", "No mode (Create or Search) specified. Exiting.");
        }

    }

    Peptide::deleteTables();

    // Note the time at finish line, calculates how long the program runs.
    time_t endTime = time(NULL);
    char *endTimeStr = ctime(&endTime);
    endTimeStr[strlen(endTimeStr) - 1] = '\0';
    double timeElapsed = difftime(endTime, startTime);


    if (!g_quiet) {

        cout.precision(4);
        if (mode == 's' || searchCount > 0) {
            cout << "Total Number of Searches Performed = " << searchCount << "; Run Time per Search = "
                 << timeElapsed / (double) (searchCount) << " seconds." << endl;
        }
        cout.precision(4);
        cout << "Total Run Time = " << timeElapsed << " seconds." << endl;

        unsigned int numWarning = g_log->getNumWarning();
        if (numWarning != 0) {
            g_log->printWarnings();
        }

        unsigned int numError = g_log->getNumError();
        if (numError == 0) {
            cout << "SpectraST finished at " << endTimeStr << " without error." << endl;
        } else {
            cout << "SpectraST finished at " << endTimeStr << " with " << numError << " error(s):" << endl;
            g_log->printErrors();
            cout << endl;
        }
    }

    stringstream perfss;
    perfss << "Total Run Time = " << timeElapsed << " seconds.";
    if (mode == 's' || searchCount > 0) {
        perfss << " Total Number of Searches Performed = " << searchCount << ". Run Time per Search = "
               << timeElapsed / (double) searchCount << " seconds.";
    }
    // log the end time
    g_log->log("PERFORMANCE", perfss.str());
    g_log->log("END", commandLiness.str(), endTimeStr);
    g_log->log("==========");
    delete (g_log);

    /*
    double retainedAve = g_retained / g_retainedCount;
    cerr << "retained=" << retainedAve << " +/- " << sqrt(g_retainedSq / g_retainedCount - retainedAve * retainedAve) << endl;
    */

    return (0);
}

#endif //LIB_BUILD

// doCreate - performs a Create mode task
void doCreate(SpectraSTCreateParams &createParams, vector<string> &fileNames) {

    // activate all create options
    createParams.finalizeOptions();

    int isList = SpectraSTFileList::isFileList(fileNames);
    if (isList < 0) {
        g_log->error("CREATE", "Input files are of mixed list and non-list types. Exiting.");
        return;
    }

    if (isList == 1) {
        // input files are lists. create a SpectraSTFileList object to do the create task
        SpectraSTFileList *fileList = new SpectraSTFileList(fileNames, createParams);
        delete (fileList);
    } else {
        // input files are individual files to be imported

        // instantiate a library; this will import the file and do everything
        SpectraSTLib *lib = new SpectraSTLib(fileNames, &createParams);
        delete (lib);
    }
}

// doSearch - performs a Search mode task
unsigned int doSearch(SpectraSTSearchParams &searchParams, vector<string> &fileNames) {

    // activate all search options
    searchParams.finalizeOptions();

    int isList = SpectraSTFileList::isFileList(fileNames);

    if (isList < 0) {
        g_log->error("CREATE", "Input files are of mixed list and non-list types. Exiting.");
        return (0);
    }

    if (isList == 1) {
        // input files are lists. create a SpectraSTFileList object to do the search task
        SpectraSTFileList *fileList = new SpectraSTFileList(fileNames, searchParams);
        unsigned int searchCount = fileList->getSearchCount();
        delete (fileList);
        return (searchCount);

    } else {

        // input files are individual files containing spectra to be searched.

        if (searchParams.libraryFile.empty()) {
            // no library file set!!
            g_log->error("SEARCH", "No library file specified. Cannot proceed with search.");
            return (0);
        }

        // instantiate the library
        SpectraSTLib *lib = new SpectraSTLib(searchParams.libraryFile, &searchParams);
        if (!g_quiet) {
            cout << "Library File loaded: \"" << searchParams.libraryFile << "\"." << endl;
        }

        // create the search task and search it
        SpectraSTSearchTask *searchTask = SpectraSTSearchTask::createSpectraSTSearchTask(fileNames, searchParams, lib);
        unsigned int searchCount = 0;
        if (searchTask) {
            searchTask->preSearch();
            searchTask->search();
            searchTask->postSearch();
            searchCount = searchTask->getSearchCount();

            delete (searchTask);
        }

        //if (searchParams.useSWATHScoring) {
        //  lib->printDotTimeProfiles();
        //}

        //if (searchParams.recordDotInLibraryEntry) {
        //  lib->printDotHistograms();
        //}

        delete (lib);

        return (searchCount);
    }
}

// calcTimeElapsed - calculates the time elapsed in seconds
double calcTimeElapsed(clock_t startTime) {
    clock_t endTime = clock();
    return ((double) (endTime - startTime) / CLOCKS_PER_SEC);
}

void readUserModFile(string &modFileName) {

    ifstream fin;
    unsigned int newTypeCount = 0;

    if (!myFileOpen(fin, modFileName)) {
        if (modFileName != "spectrast.usermods") {
            g_log->error("GENERAL", "Cannot open USERMODS file \"" + modFileName +
                                    "\" for reading user-defined modifications. Using all defaults.");
        }
        return;
    }

    g_log->log("GENERAL", "Loading user-defined modifications from \"" + modFileName + "\" .");

    string line("");
    while (nextLine(fin, line, "", "")) {
        if (line == "_EOF_") {
            return;
        }
        string::size_type pos = 0;
        string delimiter(" ,\t\r\n");
        //   string delimiter(" #{}\t\r\n");

        while (pos < line.length()) {

            string mod = nextToken(line, pos, pos, delimiter.c_str(), " \t");

            /*
            if (line[pos] == '{') {
              // {} coming
          delimiter = "#}";
            } else {
          delimiter = " #{}\t\r\n";
            }
            */

            pos++;

            if (mod.empty()) continue;

            string::size_type pipePos = 0;
            string token = nextToken(mod, pipePos, pipePos, " |", " ");
            string modMass = nextToken(mod, pipePos, pipePos, " |", " |");
            string modType = nextToken(mod, pipePos, pipePos, " |", " |");

            string allAA(ALL_AA);
            allAA += "nc";
            if (token.empty() || allAA.find(token[0]) == string::npos) continue;

            if (token.length() > 1 &&
                (token[1] != '[' || token[token.length() - 1] != ']' || (token[2] >= '0' && token[2] <= '9'))) {
                g_log->error("GENERAL", "Illegal user-defined modification token: " + token +
                                        " . It must be of the form A[...], where ... does not start with a number. Ignored.");
                continue;
            }

            if (!modType.empty() && ((modType[0] >= '0' && modType[0] <= '9') ||
                                     (modType.find_first_of(",;|\\[]{}'\"#") != string::npos))) {
                g_log->error("GENERAL", "Illegal user-defined modification type: " + modType +
                                        " . It must not start with a number or contain special characters. Ignored.");
                continue;
            }

            char aa = token[0];
            double deltaMass = 0.0;

            if (modMass.empty()) continue;

            if (modMass[0] == '+' || modMass[0] == '-' || modMass[0] == '.' ||
                (modMass[0] >= '0' && modMass[0] <= '9')) {
                deltaMass = atof(modMass.c_str());
            }

            if (fabs(deltaMass) < 0.0000001) {
                // crappy mass
                g_log->error("GENERAL", "No modification mass. Modification " + token + "|" + modMass + "|" + modType +
                                        " ignored.");
                continue;
            }

            string userToken("");
            if (token.length() > 1) userToken = token;
            double newDeltaMass = deltaMass;
            string newModType(modType);

            if (!(Peptide::processNewMod(aa, newDeltaMass, newModType, userToken))) {
                g_log->error("GENERAL", "Modification " + token + "|" + modMass + "|" + modType +
                                        " collides with existing specifications. Ignored.");
            } else {
                stringstream addModss;
                addModss << "Modification " << token << "|" << modMass << "|" << modType;
                addModss << " successfully added as " << userToken << "|" << newDeltaMass << "|" << newModType << " .";
                g_log->log("GENERAL", addModss.str());
            }
        }
    }

}


// printUsage - prints the usage of the program to console.
static void printUsage() {

    cerr << "SpectraST (version " << SPECTRAST_VERSION << "." << SPECTRAST_SUB_VERSION << ", " << szTPPVersionInfo
         << ") by Henry Lam." << endl;
    cerr << endl;

    SpectraSTCreateParams::printUsage(cerr);
    cerr << endl;

    SpectraSTSearchParams::printUsage(cerr);
    cerr << endl;

    cerr << "Miscellaneous Options:" << endl;
    cerr << "         -V           Verbose mode." << endl;
    cerr << "         -Q           Quiet mode." << endl;
    cerr << "         -L<file>     Specify name of log file. Default is \"spectrast.log\"." << endl;
    cerr << "         -M<file>     Specify name of user-defined modifications file. Default is \"spectrast.usermods\"."
         << endl;
    cerr << endl;

    cerr.flush();
}



