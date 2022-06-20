#include "SpectraSTMgfSearchTask.hpp"
#include "SpectraSTSearch.hpp"
#include "SpectraSTQuery.hpp"
#include "SpectraSTLog.hpp"
#include "SpectraSTConstants.hpp"
#include "FileUtils.hpp"
#include "Peptide.hpp"
#include "ProgressCount.hpp"

extern bool g_quiet;
extern bool g_verbose;
extern SpectraSTLog *g_log;

// constructor
SpectraSTMgfSearchTask::SpectraSTMgfSearchTask(vector<string> &searchFileNames, SpectraSTSearchParams &params,
                                               SpectraSTLib *lib) :
        SpectraSTSearchTask(searchFileNames, params, lib) {
}

// destructor
SpectraSTMgfSearchTask::~SpectraSTMgfSearchTask() {
}

// search - run the searches
void SpectraSTMgfSearchTask::search() {

    for (unsigned int n = 0; n < (unsigned int) m_searchFileNames.size(); n++) {

        searchOneFile(n);
    }

    m_searchTaskStats.logStats();
}

// searchOneFile - search one mgf file
void SpectraSTMgfSearchTask::searchOneFile(unsigned int fileIndex) {

    string searchFileName(m_searchFileNames[fileIndex]);

    ifstream fin;
    if (!myFileOpen(fin, searchFileName)) {
        g_log->error("SEARCH",
                     "Cannot open MGF file \"" + searchFileName + " for reading query spectra. File Skipped.");
        return;
    }

    m_outputs[fileIndex]->openFile();
    m_outputs[fileIndex]->printHeader();

    string line("");
    string::size_type pos = 0;

    if (g_verbose) {
        cout << "Searching query spectra in file \"" << searchFileName << "\" ..." << endl;
    }

    ProgressCount pc(!g_quiet && !g_verbose, 200);
    string msg("Searching query spectra in file \"");
    msg += searchFileName + "\"";
    pc.start(msg);

    while (true) {

        if (line == "_EOF_") {
            // no more record
            pc.done();
            m_outputs[fileIndex]->printFooter();
            m_outputs[fileIndex]->closeFile();
            return;
        }

        // skips over all lines until the line with Name: (will stop when it reaches either "Name:" or the end-of-file)
        if (line.compare(0, 10, "BEGIN IONS") != 0) {
            while (nextLine(fin, line, "BEGIN IONS", ""));
            if (line == "_EOF_") {
                // no more record
                pc.done();
                m_outputs[fileIndex]->printFooter();
                m_outputs[fileIndex]->closeFile();
                return;
            }
        }

        double precursorMz = 0.0;
        string title("");
        string comments("");
        int charge = 0;


        // read the rest of the header fields until TITLE:
        while (nextLine(fin, line, "END IONS", "")) {

            // cerr << line << endl;

            if (line.find('=') == string::npos) {
                // no more headers
                break;
            }

            if (line.compare(0, 6, "TITLE=") == 0) {
                title = nextToken(line, 6, pos, " \t\r\n");
            } else if (line.compare(0, 8, "PEPMASS=") == 0) {
                precursorMz = atof((nextToken(line, 8, pos, " \r\n")).c_str());
                double precursorIntensity = atof((nextToken(line, pos, pos, " \r\n")).c_str());
            } else if (line.compare(0, 7, "CHARGE=") == 0) {
                charge = atoi((nextToken(line, 7, pos, " \t\r\n", "+")).c_str());
            } else if (line == "_EOF_" || line.compare(0, 10, "BEGIN IONS") == 0) {
                cerr << "\nBadly formatted .mgf file! Search task truncated." << endl;
                m_outputs[fileIndex]->printFooter();
                m_outputs[fileIndex]->closeFile();
                return;

            } else {
                // some other header, put the whole thing in comments
                comments += nextToken(line, 0, pos, " \t\r\n") + " ";
            }
        }

        bool hasNoPeaks = false;
        if (line.find("END IONS") != string::npos) {
            // no peaks?
            hasNoPeaks = true;
        }

        bool selected;
        if (!m_searchAll) {
            selected = isInSelectedList(title);
            // not to mess up the parsing, the following fields will still be
            // parsed even though the record is not selected to be searched.
            // the search simply will not be initiated.
        } else {
            selected = true;
        }

        if (precursorMz < 0.0001) {
            // no mass, this can't be searched
            selected = false;
        }

        if (!hasNoPeaks) {


            int numPeaks = 0;

            SpectraSTPeakList *peakList = new SpectraSTPeakList(precursorMz, charge);
            peakList->setNoiseFilterThreshold(m_params.filterRemovePeakIntensityThreshold);

            SpectraSTQuery *query = new SpectraSTQuery(title, precursorMz, charge, comments, peakList);

            // line should contain the first peak

            do { // will stop when it reaches the next "END IONS" or end-of-file

                // here are the peaks
                double mz = atof((nextToken(line, 0, pos)).c_str());
                // pos1 now stores the position of the first space after the mz
                float intensity = atof((nextToken(line, pos, pos)).c_str());
                // pos2 now stores the position of the first space after the intensity
                peakList->insertForSearch(mz, intensity, "");

            } while (nextLine(fin, line, "END IONS", ""));

            if (!selected || !peakList->passFilter(m_params)) {
                // not selected to be searched, or bad spectra, ignore
                delete query;
            } else {

//	double retained = peakList->simplify(m_params.filterMaxPeaksUsed, m_params.filterMaxDynamicRange);

                FileName fn;
                parseFileName(searchFileName, fn);

                // create the search based on what is read, then search
                SpectraSTSearch *s = new SpectraSTSearch(query, m_params, m_outputs[fileIndex]);
                s->search(m_lib);

                m_searchCount++; // counting searches in all mgf files
                m_searchTaskStats.processSearch(s);

                pc.increment();

                // print search result
                s->print();
                delete (s);
            }
        }

    }
    m_outputs[fileIndex]->printFooter();
    m_outputs[fileIndex]->closeFile();

}


